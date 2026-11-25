// 2-D hyperbolic O-grid generator. Algorithm based on Ricky Reusser's
// grid-generation-experiments (2016); independent C++17 implementation.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using std::size_t;
struct Vec { double x=0, y=0; };
Vec operator+(Vec a, Vec b) { return {a.x+b.x,a.y+b.y}; }
Vec operator-(Vec a, Vec b) { return {a.x-b.x,a.y-b.y}; }
Vec operator*(Vec a, double s) { return {a.x*s,a.y*s}; }
Vec operator/(Vec a, double s) { return {a.x/s,a.y/s}; }
double dot(Vec a,Vec b) { return a.x*b.x+a.y*b.y; }
double cross(Vec a,Vec b) { return a.x*b.y-a.y*b.x; }
double norm(Vec a) { return std::hypot(a.x,a.y); }
bool finite(Vec a) { return std::isfinite(a.x)&&std::isfinite(a.y); }
Vec right_normal(Vec t) { return {t.y,-t.x}; }
using Ring=std::vector<Vec>;
using Grid=std::vector<Ring>;
static size_t prev(size_t i,size_t n) {return i==0?n-1:i-1;}
static size_t next(size_t i,size_t n) {return i+1==n?0:i+1;}
static std::string trim(const std::string& s) {
    const auto a=s.find_first_not_of(" \t\r\n"), b=s.find_last_not_of(" \t\r\n");
    return a==std::string::npos?std::string():s.substr(a,b-a+1);
}
static double number(const std::string& s) {
    size_t used=0; double x=std::stod(s,&used);
    if (used!=s.size() || !std::isfinite(x)) throw std::runtime_error("Invalid number: "+s);
    return x;
}
static size_t integer(const std::string& s) {
    if (s.empty() || s.find_first_not_of("0123456789")!=std::string::npos) throw std::runtime_error("Invalid integer: "+s);
    size_t used=0; auto x=std::stoull(s,&used);
    if (used!=s.size() || x>std::numeric_limits<size_t>::max()) throw std::runtime_error("Invalid integer: "+s);
    return static_cast<size_t>(x);
}
struct Options {
    std::string input,output="grid.xyz",mode="height",stencil,integrator="implicit";
    std::string smoothing,first_layer,format="2d",i_order="cw",check="global";
    double h=0,growth=1.15,alpha=-1,diffusion=0.002,cfl=0.25,max_step=0,theta=1,max_skew_angle=90;
    double max_skew_ratio=std::numeric_limits<double>::infinity();
    bool alpha_set=false, multiblock=false, quiet=false;
    size_t layers=0,max_substeps=1000000;
};
static void help() {
    std::cout << R"(Hyperbolic O-grid generator (C++17)
Usage:
  hypergrid --input wall.csv --output grid.xyz --first-height 1e-4 --growth 1.15 --layers 40 [options]
Required:
  --input PATH            CSV x,y; simple closed CCW polyline (closing point optional)
  --first-height H        Positive first layer height, in input length units
  --layers N              Number of cell layers (J=N+1)
Options:
  --output PATH           ASCII Plot3D file (default grid.xyz)
  --growth Q              Geometric growth ratio (>0, default 1.15)
  --mode height|original  Default height; original reproduces JS marching PDE
  --alpha A               >=0; default 0 (height) or 1 (original)
  --diffusion E           Nonnegative dimensionless coefficient (default 0.002)
                          Physical coefficient = E * max(bbox width,height)
  --smoothing tangent|full  Default tangent (height) or full (original)
  --max-skew-angle DEG    Limit tangent/normal RHS speed ratio to tan(DEG)
                          0..90 degrees; 90=disabled (default), 0=pure normal
                          Enabled only for tangent smoothing
  --first-layer exact|march Default exact (height) or march (original)
  --stencil three-point|compact  Default three-point (implicit) or compact (explicit)
  --integrator implicit|euler|rk2|rk4  Default implicit (linearized theta method)
  --theta T               Implicit theta, 0.5<=T<=1 (default 1: backward Euler)
  --cfl C                 Explicit-only substep safety factor (0<C<=1; default 0.25)
  --max-step H            Optional max substep, in input units (0=automatic)
  --max-substeps N        Per-layer limit (default 1000000)
  --format 2d|3d          Plot3D coordinates; 3d uses K=1 and Z=0 (default 2d)
  --multiblock            Prepend block count 1
  --i-order cw|ccw        Output I direction (default cw for positive Jacobian)
  --check global|local    Global contour/nesting + local Jacobian, or local only
  --quiet                 Suppress per-layer log
  --help                  Show this help
Notes: no resampling; input points preserved. Exact first layer uses discrete
       outward normals; subsequent step heights are nominal marching distances.
       No clipping, repair, or partial mesh output on a failed quality check.
)";
}
static Options parse(int argc,char** argv) {
    Options o;
    for (int i=1;i<argc;++i) {
        std::string k=argv[i];
        if(k=="--multiblock") {o.multiblock=true;continue;}
        if(k=="--quiet") {o.quiet=true;continue;}
        if(i+1>=argc) throw std::runtime_error("Missing value for "+k);
        std::string v=argv[++i];
        if(k=="--input") o.input=v; else if(k=="--output") o.output=v;
        else if(k=="--first-height") o.h=number(v); else if(k=="--growth") o.growth=number(v);
        else if(k=="--layers") o.layers=integer(v); else if(k=="--mode") o.mode=v;
        else if(k=="--alpha") {o.alpha=number(v);o.alpha_set=true;} else if(k=="--diffusion") o.diffusion=number(v);
        else if(k=="--max-skew-angle") o.max_skew_angle=number(v);
        else if(k=="--smoothing") o.smoothing=v; else if(k=="--first-layer") o.first_layer=v;
        else if(k=="--stencil") o.stencil=v; else if(k=="--integrator") o.integrator=v;
        else if(k=="--theta") o.theta=number(v); else if(k=="--cfl") o.cfl=number(v); else if(k=="--max-step") o.max_step=number(v);
        else if(k=="--max-substeps") o.max_substeps=integer(v); else if(k=="--format") o.format=v;
        else if(k=="--i-order") o.i_order=v; else if(k=="--check") o.check=v;
        else throw std::runtime_error("Unknown option: "+k);
    }
    if(o.mode!="height"&&o.mode!="original") throw std::runtime_error("Invalid --mode");
    if(!o.alpha_set) o.alpha=o.mode=="height"?0:1;
    if(o.smoothing.empty()) o.smoothing=o.mode=="height"?"tangent":"full";
    if(o.first_layer.empty()) o.first_layer=o.mode=="height"?"exact":"march";
    if(o.input.empty()||o.output.empty()||o.h<=0||o.layers==0) throw std::runtime_error("--input, positive --first-height and --layers are required");
    if(o.growth<=0||o.alpha<0||o.diffusion<0||o.cfl<=0||o.cfl>1||o.max_step<0||o.max_substeps==0) throw std::runtime_error("Invalid numeric control parameter");
    if(o.smoothing!="tangent"&&o.smoothing!="full") throw std::runtime_error("Invalid --smoothing");
    if(o.max_skew_angle<0||o.max_skew_angle>90)throw std::runtime_error("Invalid --max-skew-angle (expected 0 to 90 degrees)");
    if(o.max_skew_angle<90){
        if(o.smoothing!="tangent")throw std::runtime_error("--max-skew-angle below 90 requires --smoothing tangent");
        o.max_skew_ratio=std::tan(o.max_skew_angle*std::acos(-1.0)/180.0);
    }
    if(o.first_layer!="exact"&&o.first_layer!="march") throw std::runtime_error("Invalid --first-layer");
    if(o.stencil.empty())o.stencil=o.integrator=="implicit"?"three-point":"compact";
    if(o.theta<0.5||o.theta>1)throw std::runtime_error("Invalid --theta (expected 0.5 to 1)");
    if(o.integrator=="implicit"&&o.stencil!="three-point")throw std::runtime_error("Implicit integration requires --stencil three-point; use --integrator rk4 for compact");
    if(o.stencil!="compact"&&o.stencil!="three-point") throw std::runtime_error("Invalid --stencil");
    if(o.integrator!="implicit"&&o.integrator!="euler"&&o.integrator!="rk2"&&o.integrator!="rk4") throw std::runtime_error("Invalid --integrator");
    if(o.format!="2d"&&o.format!="3d") throw std::runtime_error("Invalid --format");
    if(o.i_order!="cw"&&o.i_order!="ccw") throw std::runtime_error("Invalid --i-order");
    if(o.check!="local"&&o.check!="global") throw std::runtime_error("Invalid --check");
    return o;
}
struct Wall {Ring raw,p; Vec origin;double scale=0, area=0; std::vector<double> ds;};
static Wall read_wall(const std::string& file) {
    std::ifstream in(file); if(!in) throw std::runtime_error("Cannot read "+file);
    Wall w; std::string s; size_t line=0; bool first=true;
    while(std::getline(in,s)) {
        ++line; if(line==1&&s.compare(0,3,"\xEF\xBB\xBF")==0) s.erase(0,3);
        const auto comment=s.find('#'); if(comment!=std::string::npos) s.erase(comment);
        s=trim(s); if(s.empty()) continue;
        auto comma=s.find(',');
        if(comma==std::string::npos||s.find(',',comma+1)!=std::string::npos) throw std::runtime_error("Expected two CSV columns at line "+std::to_string(line));
        auto x=trim(s.substr(0,comma)), y=trim(s.substr(comma+1));
        if(first&&(x=="x"||x=="X")&&(y=="y"||y=="Y")) {first=false;continue;}
        first=false;
        try {w.raw.push_back({number(x),number(y)});} catch(const std::exception& e){throw std::runtime_error("CSV line "+std::to_string(line)+": "+e.what());}
    }
    if(w.raw.size()<5) throw std::runtime_error("Need at least five distinct wall points");
    Vec lo=w.raw[0],hi=lo;
    for(auto p:w.raw) {lo.x=std::min(lo.x,p.x);lo.y=std::min(lo.y,p.y);hi.x=std::max(hi.x,p.x);hi.y=std::max(hi.y,p.y);}
    w.scale=std::max(hi.x-lo.x,hi.y-lo.y);w.origin=lo;
    if(!(w.scale>0&&std::isfinite(w.scale))) throw std::runtime_error("Degenerate wall bounding box");
    for(auto p:w.raw) w.p.push_back((p-w.origin)/w.scale);
    if(norm(w.p.front()-w.p.back())<=1e-12) {w.p.pop_back();w.raw.pop_back();}
    const size_t n=w.p.size(); if(n<5) throw std::runtime_error("Need at least five distinct wall points after closure removal");
    long double a=0; w.ds.resize(n);
    for(size_t i=0;i<n;++i) {
        auto p=w.p[i],q=w.p[next(i,n)];
        w.ds[i]=norm(q-p);
        if(w.ds[i]<=1e-13) throw std::runtime_error("Duplicate or too-close wall points at index "+std::to_string(i+1));
        a+=static_cast<long double>(p.x)*q.y-static_cast<long double>(p.y)*q.x;
    }
    if(a<=0) throw std::runtime_error("Wall must be counterclockwise with positive signed area");
    w.area=static_cast<double>(a/2); return w;
}
// Global segment checks: sweep over x-intervals, filter by y-intervals.
struct Edge {Vec a,b; double xmin,xmax,ymin,ymax;size_t index;int ring;};
static Edge edge(Vec a,Vec b,size_t i,int r) {return {a,b,std::min(a.x,b.x),std::max(a.x,b.x),std::min(a.y,b.y),std::max(a.y,b.y),i,r};}
static int orient(Vec a,Vec b,Vec c) {
    const long double ux=static_cast<long double>(b.x)-a.x,uy=static_cast<long double>(b.y)-a.y;
    const long double vx=static_cast<long double>(c.x)-a.x,vy=static_cast<long double>(c.y)-a.y;
    const long double d=ux*vy-uy*vx, tol=32*std::numeric_limits<long double>::epsilon()*(std::abs(ux*vy)+std::abs(uy*vx));
    return d>tol?1:(d< -tol?-1:0);
}
static bool on_edge(Vec a,Vec b,Vec c) {
    return c.x>=std::min(a.x,b.x)&&c.x<=std::max(a.x,b.x)&&c.y>=std::min(a.y,b.y)&&c.y<=std::max(a.y,b.y);
}
static bool intersects(const Edge& e,const Edge& f) {
    if(e.xmax<f.xmin||f.xmax<e.xmin||e.ymax<f.ymin||f.ymax<e.ymin) return false;
    const int p=orient(e.a,e.b,f.a),q=orient(e.a,e.b,f.b),r=orient(f.a,f.b,e.a),s=orient(f.a,f.b,e.b);
    return (p*q<0&&r*s<0)||(p==0&&on_edge(e.a,e.b,f.a))||(q==0&&on_edge(e.a,e.b,f.b))||(r==0&&on_edge(f.a,f.b,e.a))||(s==0&&on_edge(f.a,f.b,e.b));
}
static bool contains(const Ring& r,Vec p) {
    bool inside=false;
    for(size_t i=0;i<r.size();++i) {auto a=r[i],b=r[next(i,r.size())]; if((a.y>p.y)!=(b.y>p.y)&&p.x<(b.x-a.x)*(p.y-a.y)/(b.y-a.y)+a.x) inside=!inside;}
    return inside;
}
static void check_contours(const Ring& outer,const Ring* inner=nullptr) {
    const size_t n=outer.size();std::vector<Edge> e; e.reserve(n*(inner?2:1));
    auto add=[&](const Ring& r,int k){for(size_t i=0;i<n;++i){
        const auto a=r[prev(i,n)],b=r[i],c=r[next(i,n)];
        if(orient(a,b,c)==0&&dot(b-a,c-b)<0) throw std::runtime_error("Overlapping adjacent edges at point "+std::to_string(i+1));
        e.push_back(edge(b,c,i,k));
    }};
    add(outer,0);if(inner) add(*inner,1);
    std::sort(e.begin(),e.end(),[](const Edge& a,const Edge& b){return a.xmin<b.xmin;});
    std::vector<size_t> active;active.reserve(n);
    for(size_t i=0;i<e.size();++i) {
        const auto& a=e[i];
        active.erase(std::remove_if(active.begin(),active.end(),[&](size_t j){return e[j].xmax<a.xmin;}),active.end());
        for(auto j:active) {
            const auto& b=e[j];
            if(a.ring==b.ring&&(next(a.index,n)==b.index||prev(a.index,n)==b.index)) continue;
            if(intersects(a,b)) throw std::runtime_error(std::string(a.ring==b.ring?"Self-intersection":"Intersection between adjacent fronts")+" at edges "+std::to_string(a.index+1)+", "+std::to_string(b.index+1));
        }
        active.push_back(i);
    }
    if(inner&&!contains(outer,inner->front())) throw std::runtime_error("Outer front does not contain previous front");
}
// Compact stencil Taylor matching: betaL f^(d)_{i-1}+f^(d)_i+betaR f^(d)_{i+1}
// = sum_{k=-2}^2 w_k f_{i+k}. Seven unknowns, exact for degree 0..6.
static std::array<long double,7> solve7(std::array<std::array<long double,8>,7> a) {
    for(size_t k=0;k<7;++k) {
        size_t p=k;for(size_t j=k+1;j<7;++j) if(std::abs(a[j][k])>std::abs(a[p][k]))p=j;
        if(std::abs(a[p][k])<1e-20L) throw std::runtime_error("Singular compact stencil; try --stencil three-point or redistribute wall points");
        std::swap(a[k],a[p]);
        for(size_t j=k+1;j<7;++j) {const auto f=a[j][k]/a[k][k];for(size_t l=k;l<8;++l)a[j][l]-=f*a[k][l];}
    }
    std::array<long double,7> x{};
    for(int k=6;k>=0;--k){long double s=a[k][7];for(int j=k+1;j<7;++j)s-=a[k][j]*x[j];x[k]=s/a[k][k];}
    return x;
}
static long double monomial(long double x,int k) {if(k<0)return 0;long double v=1;for(int j=1;j<=k;++j)v*=x/j;return v;}
struct StencilRow {double a=0,c=0;std::array<double,5> w{};};
static StencilRow make_stencil(int d,double mm,double m,double p,double pp,bool compact) {
    StencilRow r;
    if(!compact) {
        if(d==1) {r.w[1]=-p/(m*(m+p));r.w[2]=(p-m)/(m*p);r.w[3]=m/(p*(m+p));}
        else {r.w[1]=2/(m*(m+p));r.w[2]=-2/(m*p);r.w[3]=2/(p*(m+p));}
        return r;
    }
    const long double h=std::max(m+mm,p+pp);
    const std::array<long double,5> q{-(static_cast<long double>(m)+mm)/h,-m/h,0,p/h,(static_cast<long double>(p)+pp)/h};
    std::array<std::array<long double,8>,7> a{};
    for(int k=0;k<=6;++k) {
        a[k][0]=-monomial(q[1],k-d);a[k][1]=-monomial(q[3],k-d);
        for(int j=0;j<5;++j)a[k][2+j]=monomial(q[j],k);
        a[k][7]=k==d?1:0;
    }
    const auto v=solve7(a);r.a=static_cast<double>(v[0]);r.c=static_cast<double>(v[1]);
    const long double hd=d==1?h:h*h;
    for(int j=0;j<5;++j)r.w[j]=static_cast<double>(v[2+j]/hd);
    return r;
}
// Factor the cyclic tridiagonal system once. Bordering with last unknown:
// B u + b v = f, c^T u + v = g; u = B^-1 f + w v.
class Derivative {
    size_t n_;std::vector<StencilRow> row_;std::vector<double> diag_,factor_,w_;double schur_=1,bound_=0;bool compact_;
    static void pivot(double x){if(!std::isfinite(x)||std::abs(x)<1e-13)throw std::runtime_error("Ill-conditioned periodic compact system; try --stencil three-point");}
public:
    Derivative(const std::vector<double>& ds,int d,bool compact):n_(ds.size()),row_(n_),diag_(n_-1,1),factor_(n_-1),w_(n_-1),compact_(compact) {
        for(size_t i=0;i<n_;++i){const auto im=prev(i,n_),imm=prev(im,n_);row_[i]=make_stencil(d,ds[imm],ds[im],ds[i],ds[next(i,n_)],compact);}
        std::vector<double> bnd(n_);for(size_t i=0;i<n_;++i){for(auto a:row_[i].w)bnd[i]+=std::abs(a);}
        if(!compact_){bound_=*std::max_element(bnd.begin(),bnd.end());return;}
        w_[0]=-row_[0].a;
        for(size_t i=1;i<n_-1;++i){pivot(diag_[i-1]);factor_[i]=row_[i].a/diag_[i-1];diag_[i]-=factor_[i]*row_[i-1].c;w_[i]=-factor_[i]*w_[i-1];bnd[i]+=std::abs(factor_[i])*bnd[i-1];}
        w_.back()-=row_[n_-2].c;pivot(diag_.back());w_.back()/=diag_.back();bnd[n_-2]/=std::abs(diag_.back());
        for(size_t i=n_-2;i-->0;){pivot(diag_[i]);w_[i]=(w_[i]-row_[i].c*w_[i+1])/diag_[i];bnd[i]=(bnd[i]+std::abs(row_[i].c)*bnd[i+1])/std::abs(diag_[i]);}
        schur_=1+row_.back().c*w_[0]+row_.back().a*w_.back();pivot(schur_);
        bnd.back()=(bnd.back()+std::abs(row_.back().c)*bnd[0]+std::abs(row_.back().a)*bnd[n_-2])/std::abs(schur_);
        for(size_t i=0;i<n_-1;++i)bnd[i]+=std::abs(w_[i])*bnd.back();
        bound_=*std::max_element(bnd.begin(),bnd.end());
        if(!std::isfinite(bound_)||bound_<=0)throw std::runtime_error("Invalid derivative stability bound");
    }
    double bound() const{return bound_;}
    void compute(const Ring& p,Ring& f) const {
        for(size_t i=0;i<n_;++i){const size_t im=prev(i,n_),ip=next(i,n_);const auto& w=row_[i].w;
            // Subtract the center value to avoid cancellation of constant offsets.
            f[i]=(p[prev(im,n_)]-p[i])*w[0]+(p[im]-p[i])*w[1]+(p[ip]-p[i])*w[3]+(p[next(ip,n_)]-p[i])*w[4];}
        if(!compact_)return;
        for(size_t i=1;i<n_-1;++i)f[i]=f[i]-f[i-1]*factor_[i];
        f[n_-2]=f[n_-2]/diag_.back();
        for(size_t i=n_-2;i-->0;)f[i]=(f[i]-f[i+1]*row_[i].c)/diag_[i];
        f.back()=(f.back()-f[0]*row_.back().c-f[n_-2]*row_.back().a)/schur_;
        for(size_t i=0;i<n_-1;++i)f[i]=f[i]+f.back()*w_[i];
    }
};
struct Quality {double min_jac=1,min_h=std::numeric_limits<double>::infinity(),max_h=0,mean_h=0;};
static Quality quality(const Ring& a,const Ring& b,const Ring& tangents) {
    Quality q; const size_t n=a.size();
    auto jac=[&](Vec u,Vec v,size_t i){const double den=norm(u)*norm(v);const double j=-cross(u,v)/den;
        if(!std::isfinite(j)||j<=1e-12) throw std::runtime_error("Nonpositive/degenerate corner Jacobian at cell I="+std::to_string(i+1)+" (input CCW indexing)");
        q.min_jac=std::min(q.min_jac,j);};
    for(size_t i=0;i<n;++i){const size_t ip=next(i,n);const Vec e=a[ip]-a[i],f=b[ip]-b[i],u=b[i]-a[i],v=b[ip]-a[ip];
        jac(e,u,i);jac(e,v,i);jac(f,u,i);jac(f,v,i);
        const double h=dot(u,right_normal(tangents[i])/norm(tangents[i]));
        if(!std::isfinite(h)||h<=0)throw std::runtime_error("Nonpositive normal height at point "+std::to_string(i+1));
        q.min_h=std::min(q.min_h,h);q.max_h=std::max(q.max_h,h);q.mean_h+=h/static_cast<double>(n);
    }
    return q;
}
#include "implicit.hpp"
class Marcher {
    const Options& o_;const Wall& w_;size_t n_;Derivative d1_,d2_;implicit::Stepper implicit_;
    Ring tangent_,curvature_,k1_,k2_,k3_,k4_,tmp_,trial_;
    double rhs(const Ring& r,Ring& k) {
        d1_.compute(r,tangent_);if(o_.diffusion>0)d2_.compute(r,curvature_);
        double smin=std::numeric_limits<double>::infinity();
        for(size_t i=0;i<n_;++i){const Vec t=tangent_[i];const double s=norm(t);if(!(s>1e-12&&std::isfinite(s)))throw std::runtime_error("Degenerate tangent at point "+std::to_string(i+1));smin=std::min(smin,s);
            const Vec unit=t/s; k[i]=right_normal(unit)/std::pow(s,o_.alpha);
            if(o_.diffusion>0){
                if(o_.smoothing=="tangent"&&o_.max_skew_angle<90){
                    const double limit=o_.max_skew_ratio*std::pow(s,-o_.alpha);
                    const double vt=std::clamp(o_.diffusion*dot(curvature_[i],unit),-limit,limit);
                    k[i]=k[i]+unit*vt;
                }else{Vec cur=curvature_[i];if(o_.smoothing=="tangent")cur=unit*dot(cur,unit);k[i]=k[i]+cur*o_.diffusion;}
            }
            if(!finite(k[i]))throw std::runtime_error("Nonfinite marching RHS");
        }return smin;
    }
    void trial(const Ring& r,double h) {
        if(o_.integrator=="euler") {for(size_t i=0;i<n_;++i)trial_[i]=r[i]+k1_[i]*h;}
        else {
            for(size_t i=0;i<n_;++i)tmp_[i]=r[i]+k1_[i]*(h/2);
            rhs(tmp_,k2_);
            if(o_.integrator=="rk2")for(size_t i=0;i<n_;++i)trial_[i]=r[i]+k2_[i]*h;
            else {
                for(size_t i=0;i<n_;++i)tmp_[i]=r[i]+k2_[i]*(h/2);
                rhs(tmp_,k3_);
                for(size_t i=0;i<n_;++i)tmp_[i]=r[i]+k3_[i]*h;
                rhs(tmp_,k4_);
                for(size_t i=0;i<n_;++i)trial_[i]=r[i]+(k1_[i]+k2_[i]*2+k3_[i]*2+k4_[i])*(h/6);
            }
        }
        for(auto p:trial_)if(!finite(p))throw std::runtime_error("Nonfinite marching coordinate");
    }
public:
    Marcher(const Options& o,const Wall& w):o_(o),w_(w),n_(w.p.size()),d1_(w.ds,1,o.stencil=="compact"),d2_(w.ds,2,o.stencil=="compact"),implicit_(o,w,d1_,d2_),tangent_(n_),curvature_(n_),k1_(n_),k2_(n_),k3_(n_),k4_(n_),tmp_(n_),trial_(n_){}
    Ring tangents(const Ring& r) {Ring t(n_);d1_.compute(r,t);for(size_t i=0;i<n_;++i)if(!(norm(t[i])>1e-12&&finite(t[i])))throw std::runtime_error("Degenerate tangent at point "+std::to_string(i+1));return t;}
    size_t advance(Ring& r,double h) {
        double done=0;size_t steps=0;
        while(done<h) {
            if(++steps>o_.max_substeps)throw std::runtime_error("Substep limit exceeded; reduce diffusion or clustering, or increase --max-substeps");
            const double smin=rhs(r,k1_);
            // Conservative operator-norm CFL estimates. First derivative factor
            // bounds the Jacobian of the nonlinear normal marching term.
            const double adv=std::max(1.0,o_.alpha)/std::pow(smin,1+o_.alpha);
            double dt=std::min(h-done,o_.cfl/(d1_.bound()*adv));
            if(o_.diffusion>0) dt=std::min(dt,2*o_.cfl/(o_.diffusion*d2_.bound()));
            if(o_.max_step>0)dt=std::min(dt,o_.max_step/w_.scale);
            if(!(dt>0&&std::isfinite(dt)&&done+dt>done))throw std::runtime_error("Marching substep underflow");
            bool accepted=false;
            for(int attempt=0;attempt<16;++attempt){
                try {trial(r,dt);accepted=true;break;} catch(const std::runtime_error&) {dt*=0.5;}
                if(!(dt>0&&done+dt>done))break;
            }
            if(!accepted)throw std::runtime_error("Runge-Kutta stages failed after reducing substep");
            r.swap(trial_);done+=dt;
            if(h-done<=8*std::numeric_limits<double>::epsilon()*h)done=h;
        }return steps;
    }
    Grid generate() {
        Grid grid;grid.push_back(w_.p);double h=o_.h/w_.scale,xi=0;
        if(o_.layers==std::numeric_limits<size_t>::max()||n_>std::numeric_limits<size_t>::max()/(o_.layers+1))throw std::runtime_error("Grid dimensions overflow");
        for(size_t j=1;j<=o_.layers;++j){
            try{
                if(!(h>0&&std::isfinite(h)))throw std::runtime_error("Invalid layer increment");
                const Ring& a=grid.back();Ring b=a;const Ring t=tangents(a);size_t steps=0;
                if(j==1&&o_.first_layer=="exact") {for(size_t i=0;i<n_;++i)b[i]=a[i]+right_normal(t[i])/norm(t[i])*h;steps=1;}
                else if(o_.integrator=="implicit")steps=implicit_.advance(b,h,t);
                else steps=advance(b,h);
                const auto q=quality(a,b,t);if(o_.check=="global")check_contours(b,&a);
                xi+=h;
                if(!o_.quiet)std::cout<<"layer "<<std::setw(4)<<j<<" target="<<h*w_.scale<<" normal[min/mean/max]="<<q.min_h*w_.scale<<"/"<<q.mean_h*w_.scale<<"/"<<q.max_h*w_.scale<<" minJ="<<q.min_jac<<" steps="<<steps<<'\n';
                grid.push_back(std::move(b));h*=o_.growth;
            }catch(const std::exception& e){throw std::runtime_error("Layer "+std::to_string(j)+": "+e.what());}
        }
        std::cout<<"Nominal total marching distance: "<<xi*w_.scale<<'\n';return grid;
    }
};
static void write_plot3d(const Grid& g,const Wall& wall,const Options& o) {
    const size_t n=wall.p.size();std::ofstream out(o.output);if(!out)throw std::runtime_error("Cannot write "+o.output);
    out.imbue(std::locale::classic());out<<std::scientific<<std::setprecision(16);
    if(o.multiblock)out<<"1\n";
    out<<n+1<<' '<<g.size();if(o.format=="3d")out<<" 1";out<<'\n';
    size_t col=0;const int nc=o.format=="3d"?3:2;
    for(int c=0;c<nc;++c){
        for(size_t j=0;j<g.size();++j)for(size_t i=0;i<=n;++i){
            const size_t k=o.i_order=="cw"?(n-i)%n:i%n;
            // Preserve input double coordinates on the wall, including seam.
            const Vec p=j==0?wall.raw[k]:g[j][k]*wall.scale+wall.origin;
            if(!finite(p))throw std::runtime_error("Output coordinate overflow");
            out<<(c==0?p.x:c==1?p.y:0.0);if(++col==5){out<<'\n';col=0;}else out<<' ';
        }if(col){out<<'\n';col=0;}
    }
    out.close();if(!out)throw std::runtime_error("Failed writing Plot3D file");
}
int main(int argc,char** argv){
    try{
        if(argc==1){help();return 1;}if(argc==2&&std::string(argv[1])=="--help"){help();return 0;}
        const auto o=parse(argc,argv);namespace fs=std::filesystem;
        if(fs::absolute(o.input).lexically_normal()==fs::absolute(o.output).lexically_normal()||(fs::exists(o.output)&&fs::equivalent(o.input,o.output)))throw std::runtime_error("Input and output must be different files");
        auto w=read_wall(o.input);check_contours(w.p);
        const auto minmax=std::minmax_element(w.ds.begin(),w.ds.end());
        std::cout<<std::setprecision(7)<<"Wall points: "<<w.p.size()<<"; reference length: "<<w.scale<<"; min/max wall segment: "<<*minmax.first*w.scale<<"/"<<*minmax.second*w.scale<<'\n';
        std::cout<<"mode="<<o.mode<<" alpha="<<o.alpha<<" diffusion="<<o.diffusion<<" smoothing="<<o.smoothing<<" first-layer="<<o.first_layer<<" stencil="<<o.stencil<<" integrator="<<o.integrator;
        if(o.integrator=="implicit")std::cout<<" theta="<<o.theta;
        std::cout<<" max-skew-angle="<<o.max_skew_angle;
        std::cout<<'\n';
        Marcher m(o,w);auto g=m.generate();write_plot3d(g,w,o);
        std::cout<<"Wrote "<<o.output<<" (I="<<w.p.size()+1<<", J="<<g.size()<<", I-order="<<o.i_order<<")\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<'\n';return 2;}
}
