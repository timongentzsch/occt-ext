// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <GeomAlgoExt_Phase.hxx>
#include <GeomAlgoExt_PatchFit.hxx>
#include <GeomAlgoExt_Band.hxx>

namespace BRepExt
{
using GeomAlgoExt::Band;
}

#include <GeomAlgoExt_Failure.hxx>

#include <BSplCLib.hxx>
#include <NCollection_Array1.hxx>
#include <NCollection_Array2.hxx>
#include <math.hxx>
#include <math_Matrix.hxx>
#include <math_Vector.hxx>

#include <algorithm>
#include <cmath>
#include <map>

#if !defined(__wasm_simd128__) || defined(OCCT_GEOMALGOEXT_SCALAR)
namespace BRepExt
{
using GeomAlgoExt::v128_t;
using GeomAlgoExt::wasm_f64x2_add;
using GeomAlgoExt::wasm_f64x2_mul;
using GeomAlgoExt::wasm_f64x2_splat;
using GeomAlgoExt::wasm_v128_load;
using GeomAlgoExt::wasm_v128_store;
} // namespace BRepExt
#endif

namespace BRepExt
{
namespace
{

//! What is added to the diagonal of the scaled matrix, whose entries are one, when it does not
//! factor: from the rounding of a sum of some hundred products of doubles (1e-13) in steps of
//! a hundred up to a millionth, beyond which the pull would show against the conditions.
constexpr double firstHold = 1e-13, holdStep = 100, lastHold = 1e-6;

//! The multipliers are given up when, over this many passes, the misses of the conditions have
//! fallen by less than this share: what is left then is what the net cannot follow.
constexpr int stallPasses = 12;

constexpr double stallShare = 0.02;

//! The most passes of the multipliers: with every one the miss falls or the stall above ends
//! them, and a net that is not there after some hundreds is not brought there by more.
constexpr int mostPasses = 400;

//! The multipliers have nothing left to work on when what drives them has fallen to this share
//! of what it was at first: the square of the rounding of a double.
constexpr double spent = 1e-20;

//! A condition over the unknowns it touches, scaled by the root of its weight; what the fixed
//! poles contribute is taken to its value.
struct Row
{
  std::vector<int>    index;
  std::vector<double> coefficient;
  double              value, within;
  //! Which of the conditions it is.
  int from;
};

//! The energies' terms: the squared second and third derivatives by the parameters, each mixed
//! one as often as it occurs; the coordinates do not meet in them.
const struct
{
  int    du, dv;
  double times;
} parts[] = {{2, 0, 1}, {1, 1, 2}, {0, 2, 1}, {3, 0, 1}, {2, 1, 3}, {1, 2, 3}, {0, 3, 1}};

//! Whether the coordinates are each a system of their own: the energies and the jumps never mix
//! them, and with `dim` 3 no condition does either.
bool separable(const std::vector<FitCondition>& conditions, int dim)
{
  return dim == 3
         && std::all_of(conditions.begin(), conditions.end(), [](const FitCondition& condition) {
              int axes = 0;
              for (int t = 0; t < condition.count; ++t)
                for (int c = 0; c < 3; ++c)
                  if (condition.terms[t].direction.Coord(c + 1) != 0)
                    axes |= 1 << c;
              return (axes & (axes - 1)) == 0;
            });
}

//! The energies are read by Gauss points in every span: as many each way as integrate a product
//! of two derivatives of `degree` exactly. `out` gets each point and its weight.
void gaussSites(const std::vector<double>&              knots,
                int                                     degree,
                std::vector<std::pair<double, double>>& out)
{
  math_Vector points(1, degree + 1), weights(1, degree + 1);
  math::GaussPoints(degree + 1, points);
  math::GaussWeights(degree + 1, weights);
  for (size_t k = degree; k + degree + 1 < knots.size(); ++k)
    for (int g = 1; knots[k + 1] > knots[k] && g <= degree + 1; ++g)
      out.emplace_back((knots[k] + knots[k + 1]) / 2 + (knots[k + 1] - knots[k]) / 2 * points(g),
                       (knots[k + 1] - knots[k]) / 2 * weights(g));
}

//! Per span, the integrals of the products of the basis functions and of their derivatives of
//! each order, from the Gauss sites of the span.
std::vector<std::vector<double>> gramsOf(const std::vector<SplineBasis>&               basis,
                                         const std::vector<std::pair<double, double>>& sites,
                                         int                                           degree)
{
  const int                        q = degree + 1;
  std::vector<std::vector<double>> out(sites.size() / q, std::vector<double>(4 * q * q, 0));
  for (size_t g = 0; g < sites.size(); ++g)
    for (int derivative = 0; derivative <= splineDerivatives; ++derivative)
      for (int a = 0; a < q; ++a)
        for (int b = 0; b < q; ++b)
          out[g / q][(derivative * q + a) * q + b] +=
            sites[g].second * basis[g].value[derivative][a] * basis[g].value[derivative][b];
  return out;
}

double dot(const std::vector<double>& a, const std::vector<double>& b)
{
  double sum = 0;
  for (size_t i = 0; i < a.size(); ++i)
    sum += a[i] * b[i];
  return sum;
}

//! The least squares fit of `fitNet` in its stages: the energies assembled over the spans, the
//! conditions, jumps and prior added, the system scaled and factored, and the multipliers found
//! that bring the conditions to their values.
class NetFit
{
public:
  NetFit(const char*                                  op,
         SplineNet&                                   net,
         const std::vector<char>&                     fixed,
         const std::vector<FitCondition>&             conditions,
         const FitEnergy&                             energy,
         const std::function<double(double, double)>& density,
         const gp_XYZ*                                along,
         const std::vector<FitJump>&                  jumps,
         FitEnergies*                                 again)
      : op(op),
        net(net),
        fixed(fixed),
        conditions(conditions),
        energy(energy),
        density(density),
        along(along),
        jumps(jumps),
        again(again),
        nu(net.countU()),
        nv(net.countV()),
        pu(net.degreeU),
        pv(net.degreeV),
        dim(along ? 1 : 3),
        size(dim * nu * nv),
        separate(separable(conditions, dim)),
        reused(again && again->matrix && again->separate == separate
               && again->matrix->size == size),
        matrix(size, bandHalf(), !jumps.empty()),
        constant(size, 0)
  {
    matrix.threeBlocks = separate;
  }

  //! The Gauss sites and the bases there, and the given poles in the unknowns.
  void prepareEnergy();

  //! The energies over every span; fitted again from the same poles, the energies as they were.
  void assembleEnergy();

  //! The energy copied to the other coordinates, the fixed poles held, and the energies kept for
  //! a fit again.
  void completeEnergy();

  //! Each free pole pulled toward `polePrior` by `poleInfluence`, relative to the mean diagonal.
  void addPrior(const std::vector<gp_XYZ>* polePrior, double poleInfluence);

  //! A row per condition: soft ones go into the energy, watched ones are only read.
  void addConditions();

  //! A wish per jump that the second derivative across a line of the net be the same on both
  //! sides of it.
  void addJumps();

  //! The weights of the conditions and of the energies are many orders apart, and so are the
  //! rows of the matrix. Each unknown is taken in the unit that makes its own diagonal entry
  //! one (the scaling of Jacobi), so that what the factorisation sees is how the unknowns hang
  //! together and not how they were weighed.
  void scaleUnknowns();

  //! Poles that no condition reaches and the energy holds only by rounding (far outside a
  //! loop, at a thousandth of a weight) can still leave the matrix short of definite. They
  //! are then held, as weakly as it takes, to where they were given: a pull that is nothing
  //! against any condition and decides only what nothing else does.
  void factor();

  //! The unknowns: with the multipliers, unless `leastSquares` asks the weighted fit alone.
  void solve(bool leastSquares);

  //! The net moved by the unknowns; returns the watched conditions that are missed.
  std::vector<int> moveNet(bool leastSquares);

private:
  //! A jump is read on both sides of a knot, and so over the poles of two spans.
  int bandHalf() const;

  //! The index of unknown `c` of pole `p`.
  int unknown(int p, int c) const;

  //! What a vector counts for in unknown `c` of a pole, and a pole as it was given in it.
  double part(const gp_XYZ& vector, int c) const;

  //! The energy over one cell into the matrix, by its local poles; with what the given net
  //! bends taken to the other side.
  void addElement(const SplineBasis& bu, const SplineBasis& bv);

  //! The outer product of `values` by `scale` added to the element. Independent SIMD lanes retain
  //! each entry's multiplication and addition order.
  void accumulate(double scale);

  //! The energies read at every pair of Gauss sites.
  void uniformEnergy();

  //! Integrate the tensor factors once per span. Most cells have one density throughout;
  //! only cells crossed by the trim need the full two-dimensional quadrature.
  void tensorEnergy();

  //! Reversing the knots reads the exact left limit at a repeated knot. Reading a
  //! nearby parameter instead introduces tiny coefficients on otherwise absent poles,
  //! widening the Cholesky skyline and every multiplier solve for no useful constraint.
  SplineBasis leftBasis(const std::vector<double>& knots, int degree, double t) const;

  //! The unknowns for given shifts.
  std::vector<double> poles(const std::vector<double>& shift, bool homogeneous);

  //! What the conditions miss by with unknowns `x`.
  std::vector<double> apply(const std::vector<double>& x, bool homogeneous) const;

  //! Whether a watched condition is met by the unknowns as they are.
  bool read(const Row& row) const;

  //! The conditions are weighted, not imposed, so the solution of this system leaves each a
  //! little short by what the energy pulls against it. A shift of each condition's value (a
  //! multiplier) takes that away, down to what the net cannot follow at all. The shifts solve
  //! a symmetric system of their own, whose product with a vector is one solve with the factor
  //! above; it is singular where conditions say the same thing twice, which the method of
  //! conjugate residuals does not mind.
  void conjugateResiduals();

  const char*                                  op;
  SplineNet&                                   net;
  const std::vector<char>&                     fixed;
  const std::vector<FitCondition>&             conditions;
  const FitEnergy&                             energy;
  const std::function<double(double, double)>& density;
  const gp_XYZ* const                          along;
  const std::vector<FitJump>&                  jumps;
  FitEnergies* const                           again;
  const int                                    nu, nv, pu, pv;
  //! Unknowns per pole: its three coordinates, or how far it stands along the one direction.
  const int dim;
  const int size;
  //! The energies and the jumps never mix the coordinates, and where no condition does either
  //! each coordinate is a system of its own: they are then taken one after the other, each in
  //! a band a third as wide, and not side by side pole by pole. Their entries and the order
  //! of every sum are the same; only the zeros between them are not visited.
  const bool                             separate;
  const bool                             reused;
  Band                                   matrix;
  std::vector<double>                    constant;
  std::vector<std::pair<double, double>> us, vs;
  std::vector<SplineBasis>               bus, bvs;
  int                                    local = 0;
  std::vector<double>                    element, values;
  //! The given poles and the local-to-global indices do not depend on the quadrature
  //! point. Read them once, preserving the order of every subsequent energy sum.
  std::vector<double> given;
  std::vector<int>    localPoles;
  std::vector<Row>    rows, watched;
  std::vector<double> unit;
  size_t              count = 0;
  std::vector<double> shifts;
  //! The unknowns are kept in step with the shifts, for the conditions that are only read.
  std::vector<double> solution;
};

int NetFit::bandHalf() const
{
  const int reach = jumps.empty() ? pu * nv + pv : 2 * pu * nv + 2 * pv;
  return separate ? reach : dim * reach + dim - 1;
}

int NetFit::unknown(int p, int c) const
{
  return separate ? c * nu * nv + p : dim * p + c;
}

double NetFit::part(const gp_XYZ& vector, int c) const
{
  return along ? vector.Dot(*along) : vector.Coord(c + 1);
}

void NetFit::prepareEnergy()
{
  if (!reused)
  {
    gaussSites(net.knotsU, pu, us);
    gaussSites(net.knotsV, pv, vs);
  }
  for (const auto& [u, weight] : us)
    bus.push_back(splineBasis(op, net.knotsU, pu, u));
  for (const auto& [v, weight] : vs)
    bvs.push_back(splineBasis(op, net.knotsV, pv, v));
  local = (pu + 1) * (pv + 1);
  element.assign(static_cast<size_t>(local) * local, 0);
  values.assign(local, 0);
  given.assign(size, 0);
  for (int p = 0; p < nu * nv; ++p)
    for (int c = 0; c < dim; ++c)
      given[unknown(p, c)] = part(net.poles[p], c);
  localPoles.assign(local, 0);
}

void NetFit::addElement(const SplineBasis& bu, const SplineBasis& bv)
{
  for (int a = 0; a <= pu; ++a)
    for (int b = 0; b <= pv; ++b)
      localPoles[a * (pv + 1) + b] = (bu.first + a) * nv + bv.first + b;
  for (int i = 0; i < local; ++i)
  {
    const int     p = localPoles[i], row = unknown(p, 0);
    double* const matrixRow = matrix.row(row);
    for (int c = 0; c < dim; ++c)
    {
      const int    pc     = unknown(p, c);
      const double givenP = given[pc];
      for (int k = 0; k <= i; ++k)
      {
        const int    q = localPoles[k], qc = unknown(q, c);
        const double e = element[static_cast<size_t>(i) * local + k];
        // What is looked for is how far each pole moves from where it was
        // given: what the given net bends goes to the other side.
        // Separate coordinates have one energy: it is added to the first.
        // Local poles increase with their index, so q <= p and this is already
        // the lower triangle. Coordinates are independent: visiting one at a
        // time retains the order of each matrix and right-hand-side entry.
        if (!fixed[p] && !fixed[q] && c == 0 && e != 0)
        {
          matrixRow[qc] += e;
          matrix.first[row] = std::min(matrix.first[row], qc);
        }
        const double givenQ = given[qc];
        if (!fixed[p])
          constant[pc] -= e * givenQ;
        if (!fixed[q] && i != k)
          constant[qc] -= e * givenP;
      }
    }
  }
}

void NetFit::accumulate(double scale)
{
  for (int i = 0; i < local; ++i)
  {
    double*      row    = element.data() + static_cast<size_t>(i) * local;
    const v128_t weight = wasm_f64x2_splat(scale * values[i]);
    int          k      = 0;
    for (; k + 1 <= i; k += 2)
      wasm_v128_store(row + k,
                      wasm_f64x2_add(wasm_v128_load(row + k),
                                     wasm_f64x2_mul(weight, wasm_v128_load(values.data() + k))));
    if (k <= i)
      row[k] += scale * values[i] * values[k];
  }
}

void NetFit::assembleEnergy()
{
  const bool derivativeConstraints =
    std::any_of(conditions.begin(), conditions.end(), [](const FitCondition& condition) {
      for (int t = 0; t < condition.count; ++t)
        if (condition.terms[t].du + condition.terms[t].dv > 0)
          return true;
      return false;
    });
  // Free-curve fits can be ill-conditioned enough that quadrature summation order
  // changes the refinement decisions. Preserve their arithmetic; derivative-held
  // quads use the separable integration already used by corner-jump fits.
  // Fitted again from the same poles with other conditions: the energies as they were.
  if (reused)
  {
    matrix   = *again->matrix;
    constant = again->constant;
  }
  else if (jumps.empty() && (!derivativeConstraints || density))
    uniformEnergy();
  else
    tensorEnergy();
}

void NetFit::uniformEnergy()
{
  for (size_t gu = 0; gu < us.size(); ++gu)
    for (size_t gv = 0; gv < vs.size(); ++gv)
    {
      const double here =
        us[gu].second * vs[gv].second * (density ? density(us[gu].first, vs[gv].first) : 1);
      const SplineBasis &bu = bus[gu], &bv = bvs[gv];
      std::fill(element.begin(), element.end(), 0);
      for (const auto& term : parts)
      {
        for (int a = 0; a <= pu; ++a)
          for (int b = 0; b <= pv; ++b)
            values[a * (pv + 1) + b] = bu.value[term.du][a] * bv.value[term.dv][b];
        const double scale =
          here * term.times * (term.du + term.dv == 2 ? energy.bending : energy.variation);
        accumulate(scale);
      }
      addElement(bu, bv);
    }
}

void NetFit::tensorEnergy()
{
  const auto          gramU = gramsOf(bus, us, pu), gramV = gramsOf(bvs, vs, pv);
  const int           qu = pu + 1, qv = pv + 1;
  std::vector<double> densities(qu * qv);
  for (size_t gu = 0; gu < us.size(); gu += qu)
    for (size_t gv = 0; gv < vs.size(); gv += qv)
    {
      const SplineBasis &bu = bus[gu], &bv = bvs[gv];
      bool               uniform = true;
      for (int a = 0; a < qu; ++a)
        for (int b = 0; b < qv; ++b)
        {
          const double value    = density ? density(us[gu + a].first, vs[gv + b].first) : 1;
          densities[a * qv + b] = value;
          uniform               = uniform && value == densities[0];
        }
      std::fill(element.begin(), element.end(), 0);
      if (uniform)
      {
        const auto &u = gramU[gu / qu], &v = gramV[gv / qv];
        for (int i = 0; i < local; ++i)
          for (int k = 0; k <= i; ++k)
          {
            double sum = 0;
            for (const auto& term : parts)
              sum += term.times * (term.du + term.dv == 2 ? energy.bending : energy.variation)
                     * u[(term.du * qu + i / qv) * qu + k / qv]
                     * v[(term.dv * qv + i % qv) * qv + k % qv];
            element[static_cast<size_t>(i) * local + k] = densities[0] * sum;
          }
      }
      else
      {
        for (int a = 0; a < qu; ++a)
          for (int b = 0; b < qv; ++b)
          {
            const double here = us[gu + a].second * vs[gv + b].second * densities[a * qv + b];
            for (const auto& term : parts)
            {
              for (int i = 0; i < qu; ++i)
                for (int j = 0; j < qv; ++j)
                  values[i * qv + j] =
                    bus[gu + a].value[term.du][i] * bvs[gv + b].value[term.dv][j];
              const double scale =
                here * term.times * (term.du + term.dv == 2 ? energy.bending : energy.variation);
              accumulate(scale);
            }
          }
      }
      addElement(bu, bv);
    }
}

void NetFit::completeEnergy()
{
  // Energy is scalar in every spatial coordinate, even when later conditions mix
  // them. Copy the completed sum, retaining its order and the original skyline.
  if (dim == 3 && !reused)
  {
    for (int p = 0; p < nu * nv; ++p)
    {
      const int row = unknown(p, 0), first = matrix.first[row];
      for (int c = 1; c < dim; ++c)
      {
        const int next = unknown(p, c), offset = next - row;
        matrix.first[next] = first + offset;
        for (int q = separate ? first : first / dim; q <= (separate ? row : p); ++q)
        {
          const int column                  = separate ? q : unknown(q, 0);
          matrix.row(next)[column + offset] = matrix.row(row)[column];
        }
      }
    }
  }
  for (int p = 0; !reused && p < nu * nv; ++p)
    if (fixed[p])
      for (int c = 0; c < dim; ++c)
      {
        matrix.add(unknown(p, c), unknown(p, c), 1);
        constant[unknown(p, c)] = 0;
      }
  if (again && !reused)
  {
    again->matrix   = matrix;
    again->constant = constant;
    again->separate = separate;
  }
}

void NetFit::addPrior(const std::vector<gp_XYZ>* polePrior, double poleInfluence)
{
  if (polePrior && poleInfluence > 0)
  {
    if (polePrior->size() != net.poles.size() || !std::isfinite(poleInfluence))
      fail(op, invalid_input, "the pole prior must match the net and have finite influence");
    double diagonal  = 0;
    int    freePoles = 0;
    for (int p = 0; p < nu * nv; ++p)
      if (!fixed[p])
      {
        const int row = unknown(p, 0);
        diagonal += matrix.row(row)[row];
        ++freePoles;
      }
    // Trace normalisation makes influence dimensionless. Unknowns are displacements
    // from this iteration's net, whereas the prior stays at the original net.
    const double weight = freePoles ? poleInfluence * (diagonal / freePoles) : 0;
    for (int p = 0; p < nu * nv; ++p)
      if (!fixed[p])
        for (int c = 0; c < dim; ++c)
        {
          const int row = unknown(p, c);
          matrix.add(row, row, weight);
          constant[row] += weight * part((*polePrior)[p] - net.poles[p], c);
        }
  }
}

void NetFit::addConditions()
{
  rows.reserve(conditions.size());
  for (const auto& condition : conditions)
  {
    const auto bu = splineBasis(op, net.knotsU, pu, condition.u),
               bv = splineBasis(op, net.knotsV, pv, condition.v);
    Row row;
    row.from           = static_cast<int>(&condition - conditions.data());
    const double scale = std::sqrt(condition.weight);
    row.value          = condition.value * scale;
    row.within         = condition.within * scale;
    for (int a = 0; a <= pu; ++a)
      for (int b = 0; b <= pv; ++b)
      {
        const int p = (bu.first + a) * nv + bv.first + b;
        gp_XYZ    direction(0, 0, 0);
        for (int t = 0; t < condition.count; ++t)
          direction += condition.terms[t].direction
                       * (bu.value[condition.terms[t].du][a] * bv.value[condition.terms[t].dv][b]);
        row.value -= direction.Dot(net.poles[p]) * scale;
        for (int c = 0; c < dim && !fixed[p]; ++c)
          if (const double coefficient = part(direction, c) * scale; coefficient != 0)
          {
            row.index.push_back(unknown(p, c));
            row.coefficient.push_back(coefficient);
          }
      }
    // Nothing to be done about it by what is free: a point's place across the direction.
    if (row.index.empty())
      continue;
    if (condition.watch)
    {
      watched.push_back(std::move(row));
      continue;
    }
    for (size_t a = 0; a < row.index.size(); ++a)
      for (size_t b = 0; b <= a; ++b)
        matrix.add(row.index[a], row.index[b], row.coefficient[a] * row.coefficient[b]);
    if (condition.soft)
    {
      for (size_t a = 0; a < row.index.size(); ++a)
        constant[row.index[a]] += row.value * row.coefficient[a];
      continue;
    }
    rows.push_back(std::move(row));
  }
}

SplineBasis NetFit::leftBasis(const std::vector<double>& knots, int degree, double t) const
{
  std::vector<double> reversed(knots.rbegin(), knots.rend());
  for (double& knot : reversed)
    knot = -knot;
  SplineBasis basis = splineBasis(op, reversed, degree, -t);
  basis.first       = static_cast<int>(knots.size()) - 2 * degree - 2 - basis.first;
  for (int d = 0; d <= splineDerivatives; ++d)
  {
    std::reverse(basis.value[d], basis.value[d] + degree + 1);
    if (d % 2)
      for (int k = 0; k <= degree; ++k)
        basis.value[d][k] = -basis.value[d][k];
  }
  return basis;
}

void NetFit::addJumps()
{
  for (const FitJump& jump : jumps)
  {
    std::map<int, double> over;
    const auto&           acrossKnots  = jump.acrossU ? net.knotsU : net.knotsV;
    const int             acrossDegree = jump.acrossU ? pu : pv;
    const auto            alongBasis =
      splineBasis(op, jump.acrossU ? net.knotsV : net.knotsU, jump.acrossU ? pv : pu, jump.along);
    for (const double side : {1.0, -1.0})
    {
      const auto  across = side > 0 ? splineBasis(op, acrossKnots, acrossDegree, jump.at)
                                    : leftBasis(acrossKnots, acrossDegree, jump.at);
      const auto &bu = jump.acrossU ? across : alongBasis, &bv = jump.acrossU ? alongBasis : across;
      for (int a = 0; a <= pu; ++a)
        for (int b = 0; b <= pv; ++b)
          over[(bu.first + a) * nv + bv.first + b] +=
            side * bu.value[jump.acrossU ? 2 : 0][a] * bv.value[jump.acrossU ? 0 : 2][b];
    }
    const double scale = std::sqrt(jump.weight);
    for (int c = 0; c < 3; ++c)
    {
      gp_XYZ axis(0, 0, 0);
      axis.SetCoord(c + 1, 1);
      Row row;
      row.value = 0;
      for (const auto& [p, weight] : over)
      {
        row.value -= weight * axis.Dot(net.poles[p]) * scale;
        for (int k = 0; k < dim && !fixed[p]; ++k)
          if (const double coefficient = part(axis, k) * weight * scale; coefficient != 0)
          {
            // One unknown once: along one direction the three coordinates share it.
            const int  index = unknown(p, k);
            const auto found = std::find(row.index.begin(), row.index.end(), index);
            if (found == row.index.end())
            {
              row.index.push_back(index);
              row.coefficient.push_back(coefficient);
            }
            else
            {
              row.coefficient[found - row.index.begin()] += coefficient;
            }
          }
      }
      for (size_t a = 0; a < row.index.size(); ++a)
      {
        for (size_t b = 0; b <= a; ++b)
          matrix.add(row.index[a], row.index[b], row.coefficient[a] * row.coefficient[b]);
        constant[row.index[a]] += row.value * row.coefficient[a];
      }
    }
  }
}

void NetFit::scaleUnknowns()
{
  unit.assign(size, 1);
  for (int i = 0; i < size; ++i)
    if (const double diagonal = matrix.row(i)[i]; diagonal > 0 && std::isfinite(diagonal))
      unit[i] = 1 / std::sqrt(diagonal);
  for (int i = 0; i < size; ++i)
  {
    double* const ri = matrix.row(i);
    for (int k = std::max(0, i - matrix.half); k <= i; ++k)
      ri[k] *= unit[i] * unit[k];
    constant[i] *= unit[i];
  }
  for (auto* list : {&rows, &watched})
    for (Row& row : *list)
      for (size_t k = 0; k < row.index.size(); ++k)
        row.coefficient[k] *= unit[row.index[k]];
}

void NetFit::factor()
{
  PatchPhase phase("factorization");
  const Band assembled = matrix;
  double     hold      = 0;
  while (!matrix.factor())
  {
    hold = hold == 0 ? firstHold : hold * holdStep;
    if (hold > lastHold)
      fail(op,
           no_solution,
           "the equations of the fit cannot be solved: the net is degenerate, or a "
           "condition is not a number");
    matrix = assembled;
    for (int i = 0; i < size; ++i)
      matrix.row(i)[i] += hold;
  }
}

std::vector<double> NetFit::poles(const std::vector<double>& shift, bool homogeneous)
{
  std::vector<double> x = homogeneous ? std::vector<double>(size, 0) : constant;
  for (size_t r = 0; r < count; ++r)
  {
    const double pull = homogeneous ? shift[r] : rows[r].value - shift[r];
    for (size_t k = 0; k < rows[r].index.size(); ++k)
      x[rows[r].index[k]] += pull * rows[r].coefficient[k];
  }
  {
    PatchPhase phase("linear-solves");
    matrix.solve(x);
  }
  return x;
}

std::vector<double> NetFit::apply(const std::vector<double>& x, bool homogeneous) const
{
  std::vector<double> out(count);
  for (size_t r = 0; r < count; ++r)
  {
    double sum = homogeneous ? 0 : -rows[r].value;
    for (size_t k = 0; k < rows[r].index.size(); ++k)
      sum += rows[r].coefficient[k] * x[rows[r].index[k]];
    out[r] = sum;
  }
  return out;
}

bool NetFit::read(const Row& row) const
{
  double sum = -row.value;
  for (size_t k = 0; k < row.index.size(); ++k)
    sum += row.coefficient[k] * solution[row.index[k]];
  return std::abs(sum) <= row.within;
}

void NetFit::solve(bool leastSquares)
{
  count = rows.size();
  shifts.assign(count, 0);
  solution = poles(shifts, false);
  if (!leastSquares)
    conjugateResiduals();
}

void NetFit::conjugateResiduals()
{
  std::vector<double> residual = apply(solution, false);
  std::vector<double> moved    = poles(residual, true);
  std::vector<double> image = apply(moved, true), direction = residual, step = image;
  double              energy0 = dot(residual, image);
  const double        start   = energy0;
  double              left = dot(residual, residual), leftBefore = left;
  // Every condition as near as it was asked to be: nearer is not better, and costs the
  // surface some of its fairness.
  const auto met = [&] {
    for (size_t r = 0; r < count; ++r)
      if (std::abs(residual[r]) > rows[r].within)
        return false;
    return std::all_of(watched.begin(), watched.end(), [this](const Row& row) {
      return read(row);
    });
  };
  for (int pass = 0; pass < mostPasses && energy0 > spent * start && !met(); ++pass)
  {
    const double alpha = energy0 / dot(step, step);
    for (size_t r = 0; r < count; ++r)
    {
      shifts[r] += alpha * direction[r];
      residual[r] -= alpha * step[r];
    }
    for (int i = 0; i < size; ++i)
      solution[i] -= alpha * moved[i];
    left = dot(residual, residual);
    if (pass % stallPasses == stallPasses - 1)
    {
      if (left > (1 - stallShare) * leftBefore)
        break;
      leftBefore = left;
    }
    const std::vector<double> next = poles(residual, true);
    image                          = apply(next, true);
    const double energy1 = dot(residual, image), beta = energy1 / energy0;
    for (size_t r = 0; r < count; ++r)
    {
      direction[r] = residual[r] + beta * direction[r];
      step[r]      = image[r] + beta * step[r];
    }
    for (int i = 0; i < size; ++i)
      moved[i] = next[i] + beta * moved[i];
    energy0 = energy1;
  }
  solution = poles(shifts, false);
}

std::vector<int> NetFit::moveNet(bool leastSquares)
{
  if (leastSquares && !std::all_of(solution.begin(), solution.end(), [](double value) {
        return std::isfinite(value);
      }))
    fail(op, no_solution, "the initial fair fit is not finite");
  for (int p = 0; p < nu * nv; ++p)
  {
    if (fixed[p])
      continue;
    if (along)
      net.poles[p] += *along * (solution[p] * unit[p]);
    else
      net.poles[p] += gp_XYZ(solution[unknown(p, 0)] * unit[unknown(p, 0)],
                             solution[unknown(p, 1)] * unit[unknown(p, 1)],
                             solution[unknown(p, 2)] * unit[unknown(p, 2)]);
  }
  std::vector<int> missed;
  for (const Row& row : watched)
    if (!read(row))
      missed.push_back(row.from);
  return missed;
}

} // namespace

SplineBasis splineBasis(const char* op, const std::vector<double>& knots, int degree, double t)
{
  const NCollection_Array1<double> flat(knots.front(), 1, static_cast<int>(knots.size()));
  math_Matrix                      values(1, splineDerivatives + 1, 1, degree + 1, 0.0);
  SplineBasis                      basis{};
  if (BSplCLib::EvalBsplineBasis(std::min(splineDerivatives, degree),
                                 degree + 1,
                                 flat,
                                 t,
                                 basis.first,
                                 values))
    fail(op, no_solution, "the patch's knots are degenerate");
  basis.first -= 1;
  for (int d = 0; d <= splineDerivatives; ++d)
    for (int k = 0; k <= degree; ++k)
      basis.value[d][k] = values(d + 1, k + 1);
  return basis;
}

std::vector<int> fitNet(const char*                                  op,
                        SplineNet&                                   net,
                        const std::vector<char>&                     fixed,
                        const std::vector<FitCondition>&             conditions,
                        const FitEnergy&                             energy,
                        const std::function<double(double, double)>& density,
                        const gp_XYZ*                                along,
                        const std::vector<FitJump>&                  jumps,
                        bool                                         leastSquares,
                        FitEnergies*                                 again,
                        const std::vector<gp_XYZ>*                   polePrior,
                        double                                       poleInfluence)
{
  PatchPhase phase("fit");
  if (leastSquares && !along)
    fail(op, no_solution, "the initial fair fit requires a height-graph direction");
  NetFit fit(op, net, fixed, conditions, energy, density, along, jumps, again);
  fit.prepareEnergy();
  PatchPhase assembly("assembly");
  fit.assembleEnergy();
  fit.completeEnergy();
  fit.addPrior(polePrior, poleInfluence);
  fit.addConditions();
  fit.addJumps();
  fit.scaleUnknowns();
  assembly.finish();
  fit.factor();
  fit.solve(leastSquares);
  return fit.moveNet(leastSquares);
}

occ::handle<Geom_BSplineSurface> surfaceOf(const SplineNet& net)
{
  auto distinct =
    [](const std::vector<double>& flat, std::vector<double>& knots, std::vector<int>& counts) {
      for (double knot : flat)
        if (!knots.empty() && knot == knots.back())
        {
          ++counts.back();
        }
        else
        {
          knots.push_back(knot);
          counts.push_back(1);
        }
    };
  std::vector<double> ku, kv;
  std::vector<int>    mu, mv;
  distinct(net.knotsU, ku, mu);
  distinct(net.knotsV, kv, mv);
  const int                  nu = net.countU(), nv = net.countV();
  NCollection_Array2<gp_Pnt> poles(1, nu, 1, nv);
  for (int i = 0; i < nu; ++i)
    for (int j = 0; j < nv; ++j)
      poles.SetValue(i + 1, j + 1, gp_Pnt(net.poles[static_cast<size_t>(i) * nv + j]));
  return new Geom_BSplineSurface(
    poles,
    NCollection_Array1<double>(ku.front(), 1, static_cast<int>(ku.size())),
    NCollection_Array1<double>(kv.front(), 1, static_cast<int>(kv.size())),
    NCollection_Array1<int>(mu.front(), 1, static_cast<int>(mu.size())),
    NCollection_Array1<int>(mv.front(), 1, static_cast<int>(mv.size())),
    net.degreeU,
    net.degreeV);
}

} // namespace BRepExt
