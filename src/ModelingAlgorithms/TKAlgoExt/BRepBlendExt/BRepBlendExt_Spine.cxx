// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepBlendExt_Internal.hxx>
#include <GeomAPI_IntSS.hxx>
#include <Geom_OffsetSurface.hxx>
#include <Precision.hxx>

namespace BRepExt::blending::detail
{

//! The ball at the edge's parameter `t` as the march finds its first one: grown from the edge's
//! points on both faces to its size, so that its contacts are the ones beside the edge.
static bool grownAt(const Rail& rail, double t, Station& s)
{
  s.t = t;
  for (int i = 0; i < 2; ++i)
  {
    double first, last;
    const auto pc = BRep_Tool::CurveOnSurface(rail.edge, rail.s[i].face, first, last);
    if (pc.IsNull())
      return false;
    const gp_Pnt2d uv = pc->Value(t);
    s.uv[i][0] = uv.X();
    s.uv[i][1] = uv.Y();
  }
  for (int i = 1; i <= 8; ++i)
    if (!contactAt(rail, t, s, i / 8.0))
      return false;
  return true;
}

namespace
{
//! A ball on both supports, by its feet's parameters x = (u0, v0, u1, v1): its centre and the unit
//! tangent of the centres' path there, and the parameters' rate along that path (per unit length).
struct Traced
{
  std::array<double, 4> x{};
  gp_XYZ                centre, tangent;
  std::array<double, 4> rate{};
};

//! The ball equations of a rail of constant signed distances r: the centres off the two supports
//! agree, c0(u0, v0) = c1(u1, v1). Their solutions are the offsets' intersection, the centres'
//! path, a curve in the four parameters where the 3 x 4 Jacobian has full rank.
class BallEquations
{
public:
  explicit BallEquations(const Rail& rail)
      : rail_(rail),
        r_{rail.sigma * rail.s[0].distance, rail.sigma * rail.s[1].distance}
  {
  }

  //! c0 - c1 at x and the Jacobian's columns, with the first support's centre.
  bool at(const std::array<double, 4>& x, gp_XYZ& f, gp_XYZ J[4], gp_XYZ& centre) const
  {
    Jet a, b;
    if (!jetAt(rail_.s[0], x[0], x[1], a) || !jetAt(rail_.s[1], x[2], x[3], b))
      return false;
    centre          = a.p.XYZ() + a.n.XYZ() * r_[0];
    f               = centre - (b.p.XYZ() + b.n.XYZ() * r_[1]);
    J[0]            = (a.du + a.ndu * r_[0]).XYZ();
    J[1]            = (a.dv + a.ndv * r_[0]).XYZ();
    J[2]            = (b.du + b.ndu * r_[1]).Reversed().XYZ();
    J[3]            = (b.dv + b.ndv * r_[1]).Reversed().XYZ();
    return std::isfinite(f.Modulus());
  }

  //! The ball at x with its path's tangent: the Jacobian's null vector, whose components are
  //! the signed 3 x 3 minors (a 4 x 4 determinant with a repeated row vanishes), and the
  //! centre's motion along it, J0 t0 + J1 t1. Empty where the path has no tangent (the offsets
  //! touch there) or the ball does not solve to `tolerance`.
  std::optional<Traced> traced(const std::array<double, 4>& x, double tolerance) const
  {
    gp_XYZ f, J[4], centre;
    if (!at(x, f, J, centre) || !(f.Modulus() <= tolerance))
      return {};
    const auto det = [&](int i, int j, int k) { return J[i].Dot(J[j].Crossed(J[k])); };
    const std::array<double, 4> null = {det(1, 2, 3), -det(0, 2, 3), det(0, 1, 3), -det(0, 1, 2)};
    const gp_XYZ                move   = J[0] * null[0] + J[1] * null[1];
    const double                length = move.Modulus();
    if (!(length > 0) || !std::isfinite(length))
      return {};
    Traced out;
    out.x       = x;
    out.centre  = centre;
    out.tangent = move / length;
    for (int k = 0; k < 4; ++k)
      out.rate[k] = null[k] / length;
    return out;
  }

  //! Newton on the ball equations and one more, `extra` (its value and gradient at x).
  template <typename Extra>
  bool solve(std::array<double, 4>& x, double tolerance, Extra extra) const
  {
    return newton(
      x,
      tolerance,
      [&](const std::array<double, 4>& y, std::array<std::array<double, 5>, 4>& m, double& norm) {
        gp_XYZ f, J[4], centre;
        if (!at(y, f, J, centre))
          return false;
        for (int i = 0; i < 3; ++i)
          m[i] = {J[0].Coord(i + 1),
                  J[1].Coord(i + 1),
                  J[2].Coord(i + 1),
                  J[3].Coord(i + 1),
                  -f.Coord(i + 1)};
        double value = 0;
        extra(y, centre, J, m[3], value);
        m[3][4] = -value;
        norm    = f.Modulus() + std::abs(value);
        return std::isfinite(norm);
      });
  }

private:
  const Rail& rail_;
  double      r_[2];
};

//! The centres' path traced through its balls, from the ball at the edge's middle, by
//! pseudo-arclength continuation (Allgower and Georg, Introduction to Numerical Continuation
//! Methods, 1990/2003, ch. 2 and 6): an Euler step along the tangent, corrected by Newton on the
//! ball equations within the hyperplane across the tangent through the step's end.
class SpineTrace
{
public:
  enum class Stop
  {
    Met,     //! through the goal ball
    Closed,  //! back through the ball it started from, so the path closes without the goal
    Left,    //! out of the supports' domains
    Reversed //! through the goal ball against the way the edge runs, or through the other end
             //! ball first: the balls do not follow along the path as along the edge
  };

  SpineTrace(const Rail& rail)
      : rail_(rail),
        equations_(rail),
        tolerance_(rail.resolved(limit::solved)),
        carrier_(rail.resolved(limit::together))
  {
    // The domains are those the supports carry: a parameter that is neither periodic nor
    // unbounded ends at its surface's bound. Where a support runs on without bound, the balls
    // that rest on the faces have their centres within their distance of them.
    for (int i = 0; i < 2; ++i)
    {
      const auto& surface = rail.s[i].surface;
      surface->Bounds(bounds_[i][0], bounds_[i][1], bounds_[i][2], bounds_[i][3]);
      periodic_[i][0] = surface->IsUPeriodic();
      periodic_[i][1] = surface->IsVPeriodic();
      for (int k = 0; k < 4; ++k)
        unbounded_ = unbounded_ || (!periodic_[i][k / 2] && Precision::IsInfinite(bounds_[i][k]));
    }
    if (unbounded_)
    {
      for (const Side& side : rail.s)
        BRepBndLib::Add(side.face, faces_);
      faces_.Enlarge(std::max(rail.s[0].distance, rail.s[1].distance));
    }
  }

  //! The ball `s` on the path, its tangent as the edge runs.
  std::optional<Traced> start(const Station& s) const
  {
    std::array<double, 4> x    = {s.uv[0][0], s.uv[0][1], s.uv[1][0], s.uv[1][1]};
    const gp_XYZ          along = s.along.XYZ(), at = s.centre.XYZ();
    if (!equations_.solve(x, tolerance_, [&](const std::array<double, 4>&, const gp_XYZ& centre,
                                             const gp_XYZ J[4], std::array<double, 5>& row,
                                             double& value) {
          row   = {J[0].Dot(along), J[1].Dot(along), 0, 0, 0};
          value = (centre - at).Dot(along);
        }))
      return {};
    auto out = equations_.traced(x, tolerance_);
    if (out && out->tangent.Dot(along) < 0)
      turn(*out);
    return out;
  }

  static void turn(Traced& p)
  {
    p.tangent.Reverse();
    for (double& r : p.rate)
      r = -r;
  }

  //! From `from`, the way `sense` (+1 as the edge runs) points, until the path passes through
  //! `goal`'s ball (the plane across its path through its centre, at its centre), or through
  //! `origin`'s, or `other`'s first. The balls on the way go into `out`, `from` excluded, the
  //! goal's last, their tangents as the edge runs. Empty where a step cannot be made however
  //! short.
  std::optional<Stop> walk(const Traced&        from,
                           double               sense,
                           const Station&       goal,
                           const Station&       origin,
                           const Station*       other,
                           std::vector<Traced>& out) const
  {
    return march(from, sense, &goal, &origin, other, Precision::Infinite(), out);
  }

  //! On from `from` past `length` along the path, or until it leaves the domains or no step is
  //! made: the carrier then ends there.
  void extend(const Traced& from, double sense, double length, std::vector<Traced>& out) const
  {
    march(from, sense, nullptr, nullptr, nullptr, length, out);
  }

private:
  bool inside(const Traced& p) const
  {
    for (int i = 0; i < 2; ++i)
      for (int k = 0; k < 2; ++k)
      {
        const double value = p.x[2 * i + k], lo = bounds_[i][2 * k], hi = bounds_[i][2 * k + 1];
        if (periodic_[i][k] || (k == 1 && pastPoint(rail_.s[i].surface, value)))
          continue;
        const double rounding =
          tol::rounding(std::max({std::abs(value), std::abs(lo), std::abs(hi)}));
        if ((!Precision::IsInfinite(lo) && value < lo - rounding)
            || (!Precision::IsInfinite(hi) && value > hi + rounding))
          return false;
      }
    return !unbounded_ || !faces_.IsOut(gp_Pnt(p.centre));
  }

  enum class Cross
  {
    None,      //! the path does not cross the plane through the ball across its path
    Elsewhere, //! it crosses it away from the ball
    At,        //! it crosses it at the ball (`met`)
    Failed     //! the crossing could not be solved
  };

  //! Where the path between a and b crosses the plane through `ball`'s centre across its path.
  Cross crossing(const Traced& a, const Traced& b, const Station& ball, Traced* met = nullptr) const
  {
    const gp_XYZ along = ball.along.XYZ(), at = ball.centre.XYZ();
    const double ga = (a.centre - at).Dot(along), gb = (b.centre - at).Dot(along);
    if (!(ga * gb <= 0) || ga == gb)
      return Cross::None;
    std::array<double, 4> x;
    const double          w = ga / (ga - gb);
    for (int k = 0; k < 4; ++k)
      x[k] = a.x[k] + (b.x[k] - a.x[k]) * w;
    if (!equations_.solve(x, tolerance_, [&](const std::array<double, 4>&, const gp_XYZ& centre,
                                             const gp_XYZ J[4], std::array<double, 5>& row,
                                             double& value) {
          row   = {J[0].Dot(along), J[1].Dot(along), 0, 0, 0};
          value = (centre - at).Dot(along);
        }))
      return Cross::Failed;
    auto out = equations_.traced(x, tolerance_);
    if (!out)
      return Cross::Failed;
    if ((out->centre - at).Modulus() > rail_.resolved(limit::together))
      return Cross::Elsewhere;
    if (out->tangent.Dot(a.tangent) < 0)
      turn(*out);
    if (met)
      *met = *out;
    return Cross::At;
  }

  //! The ball on the path across the parameters' tangent `rate` at `through`, from `x`.
  std::optional<Traced> corrected(std::array<double, 4>        x,
                                  const std::array<double, 4>& rate,
                                  const std::array<double, 4>& through,
                                  const Traced&                before) const
  {
    if (!equations_.solve(x, tolerance_, [&](const std::array<double, 4>& y, const gp_XYZ&,
                                             const gp_XYZ*, std::array<double, 5>& row,
                                             double& value) {
          row   = {rate[0], rate[1], rate[2], rate[3], 0};
          value = 0;
          for (int k = 0; k < 4; ++k)
            value += rate[k] * (y[k] - through[k]);
        }))
      return {};
    auto out = equations_.traced(x, tolerance_);
    // The same way along the path.
    if (out && out->tangent.Dot(before.tangent) < 0)
      turn(*out);
    return out;
  }

  std::optional<Stop> march(const Traced&        from,
                            double               sense,
                            const Station*       goal,
                            const Station*       origin,
                            const Station*       other,
                            double               length,
                            std::vector<Traced>& out) const
  {
    // Walked the way it goes, kept as the edge runs.
    const auto kept = [sense](Traced q) {
      if (sense < 0)
        turn(q);
      return q;
    };
    Traced p = kept(from);
    const Traced first = p;
    double                      walked = 0, step = rail_.scale;
    const double                shortest = rail_.resolved(limit::together);
    while (walked <= length)
    {
      if (!(step > shortest))
        return {};
      const double h = step;
      std::array<double, 4> x;
      for (int k = 0; k < 4; ++k)
        x[k] = p.x[k] + p.rate[k] * h;
      const auto q = corrected(x, p.rate, x, p);
      // The step is kept where the cubic through both balls with their tangents, the carrier
      // (spineCurve), stays within `carrier_` of the path: checked at its middle, against the
      // path's ball across the parameters' tangent there.
      std::optional<Traced> middle;
      double                error = Precision::Infinite(), chord = 0;
      if (q && q->tangent.Dot(p.tangent) > 0)
      {
        chord = (q->centre - p.centre).Modulus();
        const gp_XYZ hermite =
          (p.centre + q->centre) / 2 + (p.tangent - q->tangent) * (chord / 8);
        std::array<double, 4> y, rate;
        for (int k = 0; k < 4; ++k)
        {
          y[k]    = (p.x[k] + q->x[k]) / 2 + (p.rate[k] - q->rate[k]) * (chord / 8);
          rate[k] = p.rate[k] + q->rate[k];
        }
        middle = corrected(y, rate, y, p);
        if (middle)
          error = (middle->centre - hermite).Modulus();
      }
      if (!(error <= carrier_))
      {
        step /= 2;
        continue;
      }
      // Fourth order: the cubic's error scales with the step's fourth power.
      step = h * std::min(2.0, error > 0 ? std::pow(carrier_ / error, 0.25) : 2.0);
      // A crossing that cannot be solved leaves the trace unproven either way.
      Traced met;
      if (goal)
        switch (crossing(p, *q, *goal, &met))
        {
          case Cross::Failed:
            return {};
          case Cross::At:
            out.push_back(kept(met));
            return out.back().tangent.Dot(goal->along.XYZ()) > 0 ? Stop::Met : Stop::Reversed;
          default:
            break;
        }
      if (other)
        switch (crossing(p, *q, *other))
        {
          case Cross::Failed:
            return {};
          case Cross::At:
            return Stop::Reversed;
          default:
            break;
        }
      // Back through the first ball, from behind it.
      if (origin && walked > 0 && (p.centre - first.centre).Dot(first.tangent) < 0)
        switch (crossing(p, *q, *origin))
        {
          case Cross::Failed:
            return {};
          case Cross::At:
            return Stop::Closed;
          default:
            break;
        }
      if (!inside(*q))
        return Stop::Left;
      walked += chord;
      out.push_back(kept(*q));
      p = *q;
    }
    return Stop::Met;
  }

  const Rail&   rail_;
  BallEquations equations_;
  double        tolerance_, carrier_;
  double        bounds_[2][4];
  bool          periodic_[2][2];
  bool          unbounded_ = false;
  Bnd_Box       faces_;
};

// Gauss-Legendre on [-1, 1] with 8 points, exact for polynomials up to degree 15 (Abramowitz and
// Stegun, Handbook of Mathematical Functions, 1964, table 25.4): nodes +-x and their weights.
constexpr double gaussNodes[4]   = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
                                    0.9602898564975363};
constexpr double gaussWeights[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745,
                                    0.1012285362903763};

double speedOf(const std::array<double, 5>& q, double tau)
{
  return std::sqrt(std::max(0.0, (((q[4] * tau + q[3]) * tau + q[2]) * tau + q[1]) * tau + q[0]));
}

//! The length of a span from its start to its local parameter tau (negative before its start).
double spanLength(const std::array<double, 5>& q, double tau)
{
  const double half = tau / 2;
  double       sum  = 0;
  for (int k = 0; k < 4; ++k)
    sum += gaussWeights[k] * (speedOf(q, half * (1 - gaussNodes[k])) + speedOf(q, half * (1 + gaussNodes[k])));
  return sum * half;
}

//! The C1 cubic through traced balls with their tangents (Hermite), parameterised by the chords
//! between them, and its arc length.
std::shared_ptr<SpineArc> spineCurve(const std::vector<Traced>& balls)
{
  auto                       out = std::make_shared<SpineArc>();
  const int                  n   = int(balls.size()) - 1;
  NCollection_Array1<gp_Pnt> poles(1, 3 * n + 1);
  NCollection_Array1<double> knots(1, n + 1);
  NCollection_Array1<int>    multiplicities(1, n + 1);
  double                     at = 0;
  knots(1)                      = 0;
  poles(1)                      = gp_Pnt(balls[0].centre);
  out->knots.push_back(0);
  out->lengths.push_back(0);
  for (int i = 0; i < n; ++i)
  {
    const auto&  a = balls[i];
    const auto&  b = balls[i + 1];
    const double chord = (b.centre - a.centre).Modulus();
    const gp_XYZ p[4]  = {a.centre, a.centre + a.tangent * (chord / 3),
                          b.centre - b.tangent * (chord / 3), b.centre};
    for (int k = 1; k < 4; ++k)
      poles(3 * i + 1 + k) = gp_Pnt(p[k]);
    at += chord;
    knots(i + 2) = at;
    // dc/dtau = A + B tau + C tau^2 on the span's Bezier (derived), and |dc/dtau|^2.
    const gp_XYZ A = (p[1] - p[0]) * 3, B = (p[2] - p[1] * 2 + p[0]) * 6,
                 C = (p[3] - p[2] * 3 + p[1] * 3 - p[0]) * 3;
    out->speeds.push_back(
      {A.Dot(A), 2 * A.Dot(B), B.Dot(B) + 2 * A.Dot(C), 2 * B.Dot(C), C.Dot(C)});
    out->knots.push_back(at);
    out->lengths.push_back(out->lengths.back() + spanLength(out->speeds.back(), 1));
  }
  for (int i = 1; i <= n + 1; ++i)
    multiplicities(i) = i == 1 || i == n + 1 ? 4 : 3;
  out->curve = new Geom_BSplineCurve(poles, knots, multiplicities, 3);
  return out;
}
} // namespace

//! An open rail's spine: the centres' path traced from the ball at the edge's middle to the balls
//! the edge's own stations find at its ends (the correspondence of the spine's ends to the
//! edge's), oriented as the edge runs, and walked by its arc length between them
//! (Rail::spineArc). Where the path closes or leaves the supports' domains before it reaches an
//! end ball, that ball is on another branch of the offsets' intersection and the spine does not
//! carry the rail.
static bool openSpine(Rail& rail, const Station& seed, const Station ends[2])
{
  const SpineTrace trace(rail);
  const auto       start = trace.start(seed);
  if (!start)
    return false;
  std::vector<Traced> ahead, behind;
  if (trace.walk(*start, 1, ends[1], seed, &ends[0], ahead) != SpineTrace::Stop::Met
      || trace.walk(*start, -1, ends[0], seed, &ends[1], behind) != SpineTrace::Stop::Met)
    return false;
  // On past the end balls as far as the rail is carried past its ends, at the rate the spine's
  // length between them (by its chords) runs to the edge's range, and a step beyond.
  double length = 0;
  for (const auto* side : {&behind, &ahead})
    for (size_t i = 0; i < side->size(); ++i)
      length += ((*side)[i].centre - (i ? (*side)[i - 1] : *start).centre).Modulus();
  const double rate = length / (rail.last - rail.first);
  std::vector<Traced> before, after;
  trace.extend(behind.back(), -1, std::max(0.0, rail.first - rail.from) * rate, before);
  trace.extend(ahead.back(), 1, std::max(0.0, rail.to - rail.last) * rate, after);
  std::vector<Traced> balls(before.rbegin(), before.rend());
  const size_t        a = balls.size();
  balls.insert(balls.end(), behind.rbegin(), behind.rend());
  balls.push_back(*start);
  balls.insert(balls.end(), ahead.begin(), ahead.end());
  const size_t b = balls.size() - 1;
  balls.insert(balls.end(), after.begin(), after.end());
  const auto arc           = spineCurve(balls);
  rail.spine               = arc->curve;
  rail.spineArc            = arc;
  rail.spineParameterFirst = rail.first;
  rail.spineParameterLast  = rail.last;
  rail.spineFirst          = arc->lengths[a];
  rail.spineLast           = arc->lengths[b];
  return true;
}

std::optional<double> spineStation(const Rail& rail, const gp_Pnt& centre)
{
  const auto& curve = rail.spineArc->curve;
  GeomAPI_ProjectPointOnCurve projection(centre, curve, curve->FirstParameter(), curve->LastParameter());
  if (projection.NbPoints() == 0 || projection.LowerDistance() > rail.resolved(limit::together))
    return {};
  return rail.spineParameterFirst
         + (rail.spineArc->length(projection.LowerDistanceParameter()) - rail.spineFirst)
             * (rail.spineParameterLast - rail.spineParameterFirst)
             / (rail.spineLast - rail.spineFirst);
}

static size_t spanOf(const std::vector<double>& ends, double value)
{
  const auto it = std::upper_bound(ends.begin(), ends.end(), value);
  return it == ends.begin() ? 0 : std::min(size_t(it - ends.begin() - 1), ends.size() - 2);
}

double SpineArc::length(double u) const
{
  const size_t i = spanOf(knots, u);
  return lengths[i] + spanLength(speeds[i], (u - knots[i]) / (knots[i + 1] - knots[i]));
}

double SpineArc::parameter(double distance) const
{
  const size_t i      = spanOf(lengths, distance);
  const double within = distance - lengths[i];
  double       tau    = within / (lengths[i + 1] - lengths[i]);
  // Newton on the span's length, whose derivative is the speed: quadratic from the chord's
  // fraction, which is within the speed's variation of the root.
  for (int iteration = 0; iteration < 16; ++iteration)
  {
    const double speed = speedOf(speeds[i], tau);
    if (!(speed > 0))
      break;
    const double step = (spanLength(speeds[i], tau) - within) / speed;
    tau -= step;
    if (!(std::abs(step) > 4 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(tau))))
      break;
  }
  return knots[i] + tau * (knots[i + 1] - knots[i]);
}

bool intersectSpine(Rail& rail)
{
  // Closed, constant-distance circuits have no stop or size-law correspondence to transfer;
  // an open rail's ends correspond to the balls at the edge's ends (openSpine).
  if ((rail.closed && (rail.from != rail.first || rail.to != rail.last))
      || !rail.sizes.places.empty() || rail.angle != 0
      || (rail.law != Law::Conic && rail.law != Law::G2))
    return false;
  for (const Side& side : rail.s)
    if (side.cliffCurve || !side.trim.IsNull())
      return false;

  const double middle = (rail.first + rail.last) / 2;
  Station seed;
  // Grow a seed on the original branch, without asking its station map to traverse a circuit:
  // a circuit's by the edge's points (its place only sets the carrier's phase), an open rail's
  // by its own stations, whose ends it is matched to.
  Rail seedRail = rail;
  if (rail.closed)
    seedRail.byCentre = false;
  if (!grownAt(seedRail, middle, seed))
    return false;

  // An open rail needs only the branch through its middle ball, traced from there.
  if (!rail.closed)
  {
    Station ends[2];
    return grownAt(rail, rail.first, ends[0]) && grownAt(rail, rail.last, ends[1])
           && openSpine(rail, seed, ends);
  }
  Handle(Geom_Surface) offsets[2];
  for (int i = 0; i < 2; ++i)
  {
    const Handle(Geom_OffsetSurface) offset = new Geom_OffsetSurface(rail.s[i].surface,
      rail.sigma * rail.s[i].distance * rail.s[i].outward);
    // Retain exact canonical offsets where the geometry library supplies them; otherwise
    // the same intersector receives the general offset evaluator.
    offsets[i] = offset->Surface();
    if (offsets[i].IsNull())
      offsets[i] = offset;
  }
  // Use the general intersector, including its analytic and marching cases. No support-family
  // formula is used here. The resulting curve is only a carrier: ballAt corrects on both
  // original offsets at every station, so an approximate SSI curve is not a contact curve.
  GeomAPI_IntSS intersection(offsets[0], offsets[1], rail.resolved(limit::solved));
  if (!intersection.IsDone())
    return false;
  Handle(Geom_Curve) chosen;
  double at = 0;
  for (int i = 1; i <= intersection.NbLines(); ++i)
  {
    const auto curve = intersection.Line(i);
    if (!curve->IsClosed())
      continue;
    GeomAPI_ProjectPointOnCurve projection(seed.centre, curve,
                                           curve->FirstParameter(), curve->LastParameter());
    if (projection.NbPoints() == 0
        || projection.LowerDistance() > rail.resolved(limit::together))
      continue;
    if (!chosen.IsNull())
      return false; // Ambiguous branch identification is not a continuation certificate.
    chosen = curve;
    at = projection.LowerDistanceParameter();
  }
  if (chosen.IsNull())
    return false;
  gp_Pnt point;
  gp_Vec tangent;
  chosen->D1(at, point, tangent);
  const double orientation = tangent.Dot(seed.along);
  if (!(tangent.Magnitude() > 0) || orientation == 0)
    return false;
  gp_Pnt start, end;
  gp_Vec ds, de;
  chosen->D1(chosen->FirstParameter(), start, ds);
  chosen->D1(chosen->LastParameter(), end, de);
  if (start.Distance(end) > rail.resolved(limit::solved)
      || !(ds.Magnitude() > 0 && de.Magnitude() > 0)
      || ds.Normalized().Dot(de.Normalized()) < 1 - tol::angle)
    return false;
  const double span = (chosen->LastParameter() - chosen->FirstParameter())
                      * (orientation < 0 ? -1 : 1);
  rail.spineParameterFirst = rail.first;
  rail.spineParameterLast  = rail.last;
  rail.spine               = chosen;
  rail.spineFirst          = at - span / 2;
  rail.spineLast           = at + span / 2;
  return true;
}

//! Analytic construction for a cylindrical cavity meeting a transverse plane.
//! The circular envelope crosses when r >= (R + r)|n_plane.axis|. The window
//! h=2 asin(r/(R + r)) encloses the complete ellipse double-point domain.
//! This is a construction domain, not an approximation allowance.
bool prepareConicWindow(Rail& rail)
{
  if (rail.spine.IsNull() || rail.law != Law::Conic || rail.rho != 0.5 || rail.sigma != -1)
    return false;
  for (int c = 0; c < 2; ++c)
  {
    GeomAdaptor_Surface cylinder(rail.s[c].surface), plane(rail.s[1 - c].surface);
    if (cylinder.GetType() != GeomAbs_Cylinder || plane.GetType() != GeomAbs_Plane
        || rail.s[c].outward != -1 || rail.s[c].distance != rail.s[1 - c].distance)
      continue;
    const gp_Vec axis(cylinder.Cylinder().Axis().Direction());
    gp_Vec normal(plane.Plane().Axis().Direction());
    normal *= rail.s[1 - c].outward;
    const double cosine = normal.Dot(axis), r = rail.s[c].distance;
    const double R = cylinder.Cylinder().Radius();
    if (!(r > 0 && r < R && std::abs(cosine) > tol::angle
          && r >= (R + r) * std::abs(cosine)))
      return false;
    gp_Vec acute = normal - axis * cosine;
    if (!(acute.Magnitude() > tol::angle))
      return false;
    acute.Normalize();
    rail.conicWindow = Rail::ConicWindow{c,
                                        axis,
                                        acute,
                                        normal,
                                        plane.Plane().Location(),
                                        2 * std::asin(r / (R + r))};
    return true;
  }
  return false;
}

double conicPhase(const Rail& rail, const Station& station)
{
  const auto& window = *rail.conicWindow;
  const gp_Vec radial = station.normal[window.cylinder] * rail.s[window.cylinder].outward;
  return std::atan2(window.axis.Dot(window.acute.Crossed(radial)), window.acute.Dot(radial));
}

} // namespace BRepExt::blending::detail
