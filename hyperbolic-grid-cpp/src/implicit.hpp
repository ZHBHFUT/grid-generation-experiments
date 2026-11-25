#pragma once

// Linearized theta method with a periodic 2x2 block-tridiagonal Jacobian.
// Included after Vec, Derivative and the mesh quality routines in main.cpp.
namespace implicit {
struct Mat {
    double a=0,b=0,c=0,d=0;
};
inline Mat eye(){return {1,0,0,1};}
inline Mat operator+(Mat x,Mat y){return {x.a+y.a,x.b+y.b,x.c+y.c,x.d+y.d};}
inline Mat operator-(Mat x,Mat y){return {x.a-y.a,x.b-y.b,x.c-y.c,x.d-y.d};}
inline Mat operator*(Mat x,double k){return {x.a*k,x.b*k,x.c*k,x.d*k};}
inline Mat operator*(Mat x,Mat y){return {x.a*y.a+x.b*y.c,x.a*y.b+x.b*y.d,x.c*y.a+x.d*y.c,x.c*y.b+x.d*y.d};}
inline Vec operator*(Mat x,Vec y){return {x.a*y.x+x.b*y.y,x.c*y.x+x.d*y.y};}
inline Mat outer(Vec x,Vec y){return {x.x*y.x,x.x*y.y,x.y*y.x,x.y*y.y};}
inline Mat inverse(Mat x){
    const double s=std::max({std::abs(x.a),std::abs(x.b),std::abs(x.c),std::abs(x.d)});
    if(!(s>0&&std::isfinite(s)))throw std::runtime_error("Invalid implicit block pivot");
    x=x*(1/s);const double det=x.a*x.d-x.b*x.c;
    if(!std::isfinite(det)||std::abs(det)<1e-13)throw std::runtime_error("Ill-conditioned implicit block pivot");
    return Mat{x.d,-x.b,-x.c,x.a}*(1/det/s);
}
// Block bordering: eliminate the first M-1 unknowns while carrying the
// last unknown as two additional RHS columns, then solve its 2x2 Schur block.
// Coefficient diagonal and RHS are overwritten. Work arrays are reused.
inline void solve_cyclic(const std::vector<Mat>& l,std::vector<Mat>& d,
                         const std::vector<Mat>& u,Ring& f,std::vector<Mat>& w){
    const size_t n=f.size(),m=n-1;
    w[0]=l[0]*(-1);
    for(size_t i=1;i<m;++i){
        const Mat g=l[i]*inverse(d[i-1]);
        d[i]=d[i]-g*u[i-1];f[i]=f[i]-g*f[i-1];w[i]=g*w[i-1]*(-1);
    }
    w[m-1]=w[m-1]-u[m-1];
    for(size_t i=m;i-->0;){
        if(i+1<m){f[i]=f[i]-u[i]*f[i+1];w[i]=w[i]-u[i]*w[i+1];}
        const Mat inv=inverse(d[i]);f[i]=inv*f[i];w[i]=inv*w[i];
    }
    const Mat schur=d[m]+u[m]*w[0]+l[m]*w[m-1];
    f[m]=inverse(schur)*(f[m]-u[m]*f[0]-l[m]*f[m-1]);
    for(size_t i=0;i<m;++i)f[i]=f[i]+w[i]*f[m];
}
// F(v,z) = R t / s^alpha + eps * P(t) z, v=r_eta, z=r_etaeta.
// A=dF/dv, B=dF/dz; include the derivative of the tangent projector.
inline Vec linearize(Vec v,Vec z,double alpha,double eps,bool tangent,Mat& A,Mat& B,
                     double max_skew_ratio=std::numeric_limits<double>::infinity()){
    const double s=norm(v);
    if(!(s>1e-12&&std::isfinite(s)))throw std::runtime_error("Degenerate tangent in implicit step");
    const Vec t=v/s,n=right_normal(t);
    const Mat P=outer(t,t),Q=eye()-P,R{0,1,-1,0};
    const double speed=std::pow(s,-alpha);
    A=(R*(eye()-P*(1+alpha)))*(speed/s);
    // Active cap: vt = sign(vt_raw) * max_skew_ratio * speed.
    // Differentiate both the speed and the moving tangent; the clipped
    // branch has zero derivative with respect to curvature z.
    if(tangent&&std::isfinite(max_skew_ratio)){
        const double raw=eps*dot(t,z),limit=max_skew_ratio*speed;
        if(max_skew_ratio==0||std::abs(raw)>limit){
            const double signed_ratio=std::copysign(max_skew_ratio,raw);
            A=A+(eye()-P*(1+alpha))*(signed_ratio*speed/s);
            B=Mat{};
            return (n+t*signed_ratio)*speed;
        }
    }
    B=(tangent?P:eye())*eps;
    if(tangent&&eps>0)A=A+(Q*dot(t,z)+outer(t,Q*z))*(eps/s);
    return n*speed+B*z;
}
// Local length-based tangential redistribution, with exact gradients of
// Vt_raw = eps * (1/l_minus - 1/l_plus). Positive means move forward.
struct Spacing {double raw;Vec gm,gp;};
inline Spacing spacing(Vec rm,Vec r,Vec rp,double eps){
    const Vec em=r-rm,ep=rp-r;
    const double lm=norm(em),lp=norm(ep);
    if(!(lm>1e-13&&lp>1e-13&&std::isfinite(lm)&&std::isfinite(lp)))
        throw std::runtime_error("Degenerate edge in current-spacing smoothing");
    const double im=1/lm,ip=1/lp;
    return {eps*(im-ip),em*(eps*im*im*im),ep*(eps*ip*ip*ip)};
}
// gm and gp are derivatives with respect to previous and next coordinates;
// the center derivative follows by translation invariance. A is dF/dv.
inline Vec linearize_current(Vec v,const Spacing& e,double alpha,double ratio,
                             Mat& A,Mat& Gm,Mat& Gp){
    const double s=norm(v);
    if(!(s>1e-12&&std::isfinite(s)))throw std::runtime_error("Degenerate tangent in implicit step");
    const Vec t=v/s,n=right_normal(t);
    const Mat P=outer(t,t),Q=eye()-P,R{0,1,-1,0},K=eye()-P*(1+alpha);
    const double speed=std::pow(s,-alpha);
    A=R*K*(speed/s);
    Gm=Gp=Mat{};
    if(std::isfinite(ratio)&&(ratio==0||std::abs(e.raw)>ratio*speed)){
        const double c=std::copysign(ratio,e.raw);
        A=A+K*(c*speed/s);
        return (n+t*c)*speed;
    }
    A=A+Q*(e.raw/s);
    Gm=outer(t,e.gm);Gp=outer(t,e.gp);
    return n*speed+t*e.raw;
}
class Stepper {
    const Options& o_;const Wall& w_;const Derivative& d1_;const Derivative& d2_;
    size_t n_;std::vector<StencilRow> s1_,s2_;std::vector<Mat> l_,d_,u_,work_;
    Ring v_,z_,delta_,trial_,current_;
    void step(const Ring& r,double h,double eps){
        const bool current=o_.smoothing_metric=="current"&&eps>0;
        d1_.compute(r,v_);if(eps>0&&!current)d2_.compute(r,z_);
        for(size_t i=0;i<n_;++i){
            const double k=-o_.theta*h;
            Mat A,B;
            Vec F;
            if(current){
                const auto e=spacing(r[prev(i,n_)],r[i],r[next(i,n_)],eps);
                Mat Gm,Gp;
                F=linearize_current(v_[i],e,o_.alpha,o_.max_skew_ratio,A,Gm,Gp);
                l_[i]=(A*s1_[i].w[1]+Gm)*k;
                d_[i]=eye()+(A*s1_[i].w[2]-Gm-Gp)*k;
                u_[i]=(A*s1_[i].w[3]+Gp)*k;
            }else{
                F=linearize(v_[i],z_[i],o_.alpha,eps,o_.smoothing=="tangent",A,B,o_.max_skew_ratio);
                l_[i]=(A*s1_[i].w[1]+B*s2_[i].w[1])*k;
                d_[i]=eye()+(A*s1_[i].w[2]+B*s2_[i].w[2])*k;
                u_[i]=(A*s1_[i].w[3]+B*s2_[i].w[3])*k;
            }
            delta_[i]=F*h;
        }
        solve_cyclic(l_,d_,u_,delta_,work_);
        for(size_t i=0;i<n_;++i){trial_[i]=r[i]+delta_[i];if(!finite(trial_[i]))throw std::runtime_error("Nonfinite implicit coordinate");}
        // Reject invalid intermediate fronts before they can contaminate later steps.
        quality(r,trial_,v_);
        if(o_.check=="global")check_contours(trial_,&r);
    }
public:
    Stepper(const Options& o,const Wall& w,const Derivative& d1,const Derivative& d2):
        o_(o),w_(w),d1_(d1),d2_(d2),n_(w.p.size()),s1_(n_),s2_(n_),l_(n_),d_(n_),u_(n_),work_(n_-1),v_(n_),z_(n_),delta_(n_),trial_(n_),current_(n_){
        for(size_t i=0;i<n_;++i){const auto im=prev(i,n_);s1_[i]=make_stencil(1,0,w.ds[im],w.ds[i],0,false);s2_[i]=make_stencil(2,0,w.ds[im],w.ds[i],0,false);}
    }
    size_t advance(Ring& r,double h,const Ring& start_tangent,double eps){
        size_t subdivisions=1,attempted=0;
        if(o_.max_step>0){const double count=std::ceil(h/(o_.max_step/w_.scale));
            if(!std::isfinite(count)||count>static_cast<double>(o_.max_substeps)||count>=static_cast<double>(std::numeric_limits<size_t>::max()))throw std::runtime_error("Implicit substep limit exceeded by --max-step");
            subdivisions=std::max<size_t>(1,static_cast<size_t>(count));}
        std::string reason;
        // Start with one solve per output layer. On failure retry the entire
        // layer with doubled subdivisions, keeping the saved front unchanged.
        for(int retry=0;retry<=10;++retry){
            current_=r;
            try{
                for(size_t j=0;j<subdivisions;++j){
                    if(attempted==o_.max_substeps)throw std::runtime_error("Implicit substep limit exceeded");
                    ++attempted;step(current_,h/static_cast<double>(subdivisions),eps);current_.swap(trial_);
                }
                // Even if all intermediate fronts are valid, the final coarse
                // quadrilaterals must also pass the same checks.
                if(subdivisions>1){quality(r,current_,start_tangent);if(o_.check=="global")check_contours(current_,&r);}
                r.swap(current_);return subdivisions;
            }catch(const std::runtime_error& e){reason=e.what();}
            if(subdivisions>o_.max_substeps/2||attempted>=o_.max_substeps)break;
            subdivisions*=2;
        }
        throw std::runtime_error("Implicit layer retry failed: "+reason);
    }
};
} // namespace implicit
