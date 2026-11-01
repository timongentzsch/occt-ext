// Copyright (c) 2026 Timon Gentzsch
// SPDX-License-Identifier: LGPL-2.1-only WITH OCCT-exception-1.0

#include <BRepAlgoExt_TrimDomain.hxx>
#include <GeomAlgoExt_Failure.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <BndLib_Add3dCurve.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2dConvert_BSplineCurveToBezierCurve.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <NCollection_Array1.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <algorithm>
#include <array>
#include <cmath>

namespace BRepExt {
namespace {
void unique(std::vector<double>& x) {
    std::sort(x.begin(),x.end());
    x.erase(std::unique(x.begin(),x.end(),[](double a,double b){return std::abs(a-b)<=1e-13*(1+std::max(std::abs(a),std::abs(b)));}),x.end());
}
double choose(int n,int k) {double out=1;for(int i=1;i<=k;++i) out=out*(n-i+1)/i;return out;}
std::vector<double> product(const std::vector<double>& a,const std::vector<double>& b) {
    const int n=int(a.size())-1,m=int(b.size())-1;std::vector<double> out(n+m+1,0);
    for(int i=0;i<=n;++i) for(int j=0;j<=m;++j) out[i+j]+=a[i]*b[j]*choose(n,i)*choose(m,j)/choose(n+m,i+j);
    return out;
}
// The Bernstein coefficients of Y'W-YW', which has the sign of v' along a positive rational
// arc. Y is taken from the first pole, against cancellation far from the origin.
template<class Arc> std::vector<double> slope(const Arc& arc,int axis) {
    const int n=arc->Degree();const double origin=arc->Pole(1).Coord(axis);
    std::vector<double> y,w,dy,dw;
    for(int i=1;i<=n+1;++i) {w.push_back(arc->Weight(i));y.push_back((arc->Pole(i).Coord(axis)-origin)*arc->Weight(i));}
    for(int i=0;i<n;++i) {dy.push_back(n*(y[i+1]-y[i]));dw.push_back(n*(w[i+1]-w[i]));}
    auto out=product(dy,w);const auto other=product(y,dw);
    for(size_t i=0;i<out.size();++i) out[i]-=other[i];
    return out;
}
// Every zero of a Bernstein polynomial, isolated by halving its hull; tangencies included.
// A piece flat to `floor` is one station, not an exponential number of halvings.
void roots(const std::vector<double>& p,double a,double b,double floor,std::vector<double>& out,int depth=0) {
    const auto [lo,hi]=std::minmax_element(p.begin(),p.end());
    if(*lo>0 || *hi<0) return;
    if(depth>=42 || *hi-*lo<=floor) {out.push_back((a+b)/2);return;}
    std::vector<double> q=p,left{p.front()},right{p.back()};
    for(size_t k=1;k<p.size();++k) {
        for(size_t j=0;j<p.size()-k;++j) q[j]=(q[j]+q[j+1])/2;
        left.push_back(q[0]);right.push_back(q[p.size()-k-1]);
    }
    std::reverse(right.begin(),right.end());
    const double m=(a+b)/2;roots(left,a,m,floor,out,depth+1);roots(right,m,b,floor,out,depth+1);
}

using Moments=std::array<double,5>;
void addTo(Moments& a,const Moments& b) {for(size_t i=0;i<a.size();++i) a[i]+=b[i];}
// Gauss-Kronrod 7/15 on area, first moments and optional signed volume together.
// Divide moments and volume by their physical radii to compare errors in area units.
struct Quadrature {
    double scale=1;
    double volumeScale=1;
    size_t calls=0,limit=0;
    bool converged=true;
    static constexpr double x[7]={.9914553711208126,.9491079123427585,.8648644233597691,.7415311855993944,.5860872354676911,.4058451513773972,.2077849550078985};
    static constexpr double wk[8]={.02293532201052922,.06309209262997855,.1047900103222502,.1406532597155259,.1690047266392679,.1903505780647854,.2044329400752989,.2094821410847278};
    static constexpr double wg[4]={.1294849661688697,.2797053914892767,.3818300505051189,.4179591836734694};
    template<class F> std::pair<Moments,double> rule(F& f,double a,double b) {
        const double m=(a+b)/2,h=(b-a)/2;const auto middle=f(m);
        Moments kronrod{},gauss{};
        for(size_t i=0;i<middle.size();++i) {kronrod[i]=h*wk[7]*middle[i];gauss[i]=h*wg[3]*middle[i];}
        for(int j=0;j<7;++j) {
            const auto p=f(m-h*x[j]),q=f(m+h*x[j]);
            for(size_t i=0;i<p.size();++i) {
                const double value=h*(p[i]+q[i]);kronrod[i]+=wk[j]*value;
                if(j%2==1) gauss[i]+=wg[j/2]*value;
            }
        }
        calls+=15;double error=0;
        for(size_t i=0;i<kronrod.size();++i)
            error=std::max(error,std::abs(kronrod[i]-gauss[i])/(i==4?volumeScale:(i?scale:1)));
        return {kronrod,error};
    }
    template<class F> std::pair<Moments,double> adaptive(F& f,double a,double b,double tolerance,int depth=0) {
        auto [sum,error]=rule(f,a,b);
        if(error<=tolerance || depth>=18 || calls>limit) return {sum,error};
        const double m=(a+b)/2;
        auto left=adaptive(f,a,m,tolerance/2,depth+1),right=adaptive(f,m,b,tolerance/2,depth+1);
        addTo(left.first,right.first);
        return {left.first,left.second+right.second};
    }
    template<class F> Moments integrate(F& f,double a,double b,double tolerance) {
        // Subintervals may use unequal shares of the error. Test their summed estimate
        // against this integral's requested budget rather than each leaf's subdivision.
        auto [sum,error]=adaptive(f,a,b,tolerance);
        if(error>tolerance) converged=false;
        return sum;
    }
};
}

std::vector<double> coordinateExtrema(const Handle(Geom_BezierCurve)& curve,int axis) {
    auto derivative=slope(curve,axis);double largest=0;
    for(double value:derivative) largest=std::max(largest,std::abs(value));
    std::vector<double> out{0,1};
    if(largest>0) roots(derivative,0,1,1e-14*largest,out);
    unique(out);return out;
}

TrimDomain::TrimDomain(const TopoDS_Face& input):face(TopoDS::Face(input.Oriented(TopAbs_FORWARD))),classifier(face,1e-12) {
    BRepTools::UVBounds(face,u0,u1,v0,v1);
    for(TopExp_Explorer it(face,TopAbs_EDGE);it.More();it.Next()) {
        // An internal or external edge bounds nothing: a row passes it without entering or leaving.
        const auto orientation=it.Current().Orientation();
        if(orientation!=TopAbs_FORWARD && orientation!=TopAbs_REVERSED) continue;
        double a,b;const auto curve=BRep_Tool::CurveOnSurface(TopoDS::Edge(it.Current()),face,a,b);
        if(curve.IsNull() || !(b>a)) continue;
        Geom2dConvert_BSplineCurveToBezierCurve arcs(Geom2dConvert::CurveToBSplineCurve(new Geom2d_TrimmedCurve(curve,a,b)));
        Piece run{{},0,0};
        const auto flush=[&] {if(!run.segments.empty()) pieces.push_back(std::move(run));run={{},0,0};};
        for(int j=1;j<=arcs.NbArcs();++j) {
            const auto arc=arcs.Arc(j);
            boundaries.push_back(std::make_unique<Boundary>(arc));
            const auto numerator=slope(arc,2);double largest=0;
            for(double c:numerator) largest=std::max(largest,std::abs(c));
            std::vector<double> stations{0,1};roots(numerator,0,1,1e-14*largest,stations);unique(stations);
            for(size_t k=1;k<stations.size();++k) {
                const double y0=arc->Value(stations[k-1]).Y(),y1=arc->Value(stations[k]).Y();
                if(!(std::abs(y1-y0)>1e-13*(1+std::abs(v1-v0)))) {flush();continue;}
                if(!run.segments.empty() && (y1>y0)!=(run.y1>run.y0)) flush();
                if(run.segments.empty()) run.y0=y0;
                run.segments.push_back({arc,stations[k-1],stations[k],y0,y1});run.y1=y1;
            }
        }
        flush();
    }
}
TrimDomain::Boundary::Boundary(const Handle(Geom2d_BezierCurve)& input):curve(input) {
    for(int i=1;i<=curve->NbPoles();++i) box.Add(curve->Pole(i));
    box.SetGap(0);
}
bool TrimDomain::Boundary::intersects(double a,double b,double c,double d,int depth) const {
    double x,y,z,w;box.Get(x,y,z,w);
    if(z<a || x>b || w<c || y>d) return false;
    if((x>=a && z<=b && y>=c && w<=d) || depth>=40) return true;
    if(!left) {
        auto first=Handle(Geom2d_BezierCurve)::DownCast(curve->Copy()),second=Handle(Geom2d_BezierCurve)::DownCast(curve->Copy());
        first->Segment(0,.5);second->Segment(.5,1);
        left=std::make_unique<Boundary>(first);right=std::make_unique<Boundary>(second);
    }
    return left->intersects(a,b,c,d,depth+1) || right->intersects(a,b,c,d,depth+1);
}
bool TrimDomain::boundaryIntersects(double a,double b,double c,double d) const {
    for(const auto& boundary:boundaries) if(boundary->intersects(a,b,c,d)) return true;
    return false;
}
// Where a row crosses a monotone piece: Newton's steps kept inside a bisection bracket.
double TrimDomain::crossing(const Piece& run,double v) const {
    const bool rising=run.y1>run.y0;
    const auto& piece=*std::find_if(run.segments.begin(),run.segments.end()-1,[&](const Segment& s){return rising ? v<=s.y1 : v>=s.y1;});
    double a=piece.t0,b=piece.t1,t=a+(b-a)*std::clamp((v-piece.y0)/(piece.y1-piece.y0),0.0,1.0);
    const double close=1e-15*(1+std::abs(v0)+std::abs(v1));
    gp_Pnt2d p;gp_Vec2d d;
    for(int i=0;i<100;++i) {
        piece.arc->D1(t,p,d);const double f=p.Y()-v;
        if(std::abs(f)<=close) break;
        ((f>0)==rising ? b : a)=t;
        const double next=t-f/d.Y();
        t=next>a && next<b ? next : (a+b)/2;
        if(b-a<=1e-16) {p=piece.arc->Value(t);break;}
    }
    return p.X();
}
std::vector<TrimDomain::Band> TrimDomain::bands(std::vector<double> breaks) const {
    breaks.push_back(v0);breaks.push_back(v1);
    for(const auto& piece:pieces) {breaks.push_back(piece.y0);breaks.push_back(piece.y1);}
    breaks.erase(std::remove_if(breaks.begin(),breaks.end(),[&](double v){return v<v0 || v>v1;}),breaks.end());
    unique(breaks);
    std::vector<Band> out;
    for(size_t i=1;i<breaks.size();++i) {
        Band band{breaks[i-1],breaks[i],{},{}};const double mid=(band.lo+band.hi)/2;
        std::vector<std::pair<double,int>> at;
        for(int k=0;k<int(pieces.size());++k)
            if(std::min(pieces[k].y0,pieces[k].y1)<mid && mid<std::max(pieces[k].y0,pieces[k].y1)) at.emplace_back(crossing(pieces[k],mid),k);
        if(at.size()<2) continue;
        std::sort(at.begin(),at.end());
        for(const auto& [u,k]:at) band.order.push_back(k);
        // A row enters the trim at its first crossing and leaves at its second, and so on; where
        // the classifier disagrees for the first gap, or an odd count tells of pcurves that do
        // not meet, it says for every gap.
        const auto classify=[&](size_t k) {
            const auto state=classifier.Perform(gp_Pnt2d((at[k-1].first+at[k].first)/2,mid),false);
            if(state==TopAbs_UNKNOWN) fail("trimDomain",no_solution,"trim classification is unknown");
            return state==TopAbs_IN || state==TopAbs_ON;
        };
        const bool alternate=at.size()%2==0 && classify(1);
        for(size_t k=1;k<at.size();++k) band.inside.push_back(at[k].first>at[k-1].first && (alternate ? k%2==1 : classify(k)));
        out.push_back(std::move(band));
    }
    return out;
}
bool TrimDomain::inside(const std::vector<Band>& bands,double u,double v) const {
    const auto band=std::upper_bound(bands.begin(),bands.end(),v,[](double at,const Band& b){return at<b.hi;});
    if(band==bands.end() || v<band->lo) return false;
    std::vector<std::pair<double,double>> row;spans(*band,v,row);
    return std::any_of(row.begin(),row.end(),[&](const auto& span){return span.first<=u && u<=span.second;});
}
void TrimDomain::spans(const Band& band,double v,std::vector<std::pair<double,double>>& out) const {
    out.clear();double previous=0;
    for(size_t k=0;k<band.order.size();++k) {
        const double u=crossing(pieces[band.order[k]],v);
        if(k && band.inside[k-1] && u>previous) out.emplace_back(previous,u);
        previous=u;
    }
}

// Rows of constant v across the trim, each integrated between its crossings and the surface's
// knots; the rows in turn across each band. Nothing outside the trim is evaluated, so an
// extrapolated or folded support beyond it costs neither accuracy nor time.
TrimProperties trimmedProperties(const TopoDS_Face& input, const gp_Pnt* volumeOrigin) {
    const auto face=TopoDS::Face(input.Oriented(TopAbs_FORWARD));
    const TrimDomain domain(face);BRepAdaptor_Surface surface(face);
    Bnd_Box physical;
    for(TopExp_Explorer it(face,TopAbs_EDGE);it.More();it.Next()) {
        const auto edge=TopoDS::Edge(it.Current());
        if(BRep_Tool::Degenerated(edge) || !BRep_Tool::IsGeometric(edge)) continue;
        BRepAdaptor_Curve curve(edge);BndLib_Add3dCurve::Add(curve,0,physical);
    }
    if(physical.IsVoid()) fail("surfaceProperties",no_solution,"finite trim boundary required");
    // Moments about the middle of the trim keep large offsets out of the sums.
    double x0,y0,z0,x1,y1,z1;physical.Get(x0,y0,z0,x1,y1,z1);
    const gp_XYZ origin((x0+x1)/2,(y0+y1)/2,(z0+z1)/2);
    Quadrature quad;quad.scale=std::max(1e-300,std::hypot(x1-x0,y1-y0,z1-z0)/2);
    if(volumeOrigin) quad.volumeScale=std::max(1e-300,std::hypot(
        std::max(std::abs(x0-volumeOrigin->X()),std::abs(x1-volumeOrigin->X())),
        std::max(std::abs(y0-volumeOrigin->Y()),std::abs(y1-volumeOrigin->Y())),
        std::max(std::abs(z0-volumeOrigin->Z()),std::abs(z1-volumeOrigin->Z()))));
    std::vector<double> uKnots,vKnots;
    {
        const int nu=surface.NbUIntervals(GeomAbs_CN),nv=surface.NbVIntervals(GeomAbs_CN);
        NCollection_Array1<double> uk(1,nu+1),vk(1,nv+1);surface.UIntervals(uk,GeomAbs_CN);surface.VIntervals(vk,GeomAbs_CN);
        for(int i=1;i<=nu+1;++i) if(uk(i)>domain.u0 && uk(i)<domain.u1) uKnots.push_back(uk(i));
        for(int i=1;i<=nv+1;++i) if(vk(i)>domain.v0 && vk(i)<domain.v1) vKnots.push_back(vk(i));
    }
    const auto bands=domain.bands(vKnots);
    if(bands.empty()) fail("surfaceProperties",no_solution,"trim has no interior");
    auto density=[&](double u,double v) {
        gp_Pnt p;gp_Vec du,dv;surface.D1(u,v,p,du,dv);
        const auto normal=du.Crossed(dv);
        const double j=normal.Magnitude();const auto r=p.XYZ()-origin;
        const double volume=volumeOrigin ? gp_Vec(*volumeOrigin,p).Dot(normal)/3 : 0;
        return Moments{j,j*r.X(),j*r.Y(),j*r.Z(),volume};
    };
    bool coarse=true;double inner=0;
    std::vector<std::pair<double,double>> spans;std::vector<double> cuts;
    auto row=[&](const TrimDomain::Band& band,double v) {
        domain.spans(band,v,spans);cuts.clear();
        for(const auto& [a,b]:spans) {
            cuts.push_back(a);
            for(double k:uKnots) if(k>a && k<b) cuts.push_back(k);
            cuts.push_back(b);cuts.push_back(NAN);
        }
        const double share=inner/std::max<size_t>(1,cuts.size());
        Moments sum{};auto f=[&](double u){return density(u,v);};
        for(size_t i=1;i<cuts.size();++i) if(!std::isnan(cuts[i-1]) && !std::isnan(cuts[i]))
            addTo(sum,coarse ? quad.rule(f,cuts[i-1],cuts[i]).first : quad.integrate(f,cuts[i-1],cuts[i],share));
        return sum;
    };
    // v = lo + width sin^2 t: at a band's end a crossing moves as the square root of v, which
    // this makes smooth instead of halving towards it.
    auto across=[&](const TrimDomain::Band& band,double tolerance) {
        const double width=band.hi-band.lo;
        auto f=[&](double t) {
            const double s=std::sin(t),c=std::cos(t);auto m=row(band,band.lo+width*s*s);
            for(auto& value:m) value*=2*width*s*c;
            return m;
        };
        return quad.integrate(f,0,M_PI/2,tolerance);
    };
    // The scale of the error allowed: each band's middle row times its width.
    double estimate=0,range=0;
    for(const auto& band:bands) {estimate+=row(band,(band.lo+band.hi)/2)[0]*(band.hi-band.lo);range+=band.hi-band.lo;}
    if(!(estimate>0) || !std::isfinite(estimate)) fail("surfaceProperties",no_solution,"trim has no positive area");
    // A tenth of the error to the rows, nine tenths shared equally by the bands.
    // A numerically tiny band must not spend the whole evaluation budget chasing an
    // error proportional to its width and leave the larger bands unrefined. The sum
    // of their assigned errors remains the same requested global budget.
    const double tolerance=(volumeOrigin?1e-11:1e-8)*estimate;
    coarse=false;inner=0.1*tolerance/range;quad.limit=quad.calls+2000000;
    Moments total{};
    for(const auto& band:bands) addTo(total,across(band,0.9*tolerance/bands.size()));
    TrimProperties result;result.area=total[0];
    const bool finite=std::all_of(total.begin(),total.end(),[](double value){return std::isfinite(value);});
    if(!finite || !(total[0]>0)) fail("surfaceProperties",no_solution,"trim integration is not finite");
    result.center=gp_Pnt(origin+gp_XYZ(total[1],total[2],total[3])/total[0]);
    result.volume=total[4];
    result.converged=quad.converged;
    return result;
}
}
