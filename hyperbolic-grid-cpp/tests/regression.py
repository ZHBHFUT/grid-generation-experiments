"""Dependency-free regression checks for the hypergrid command-line program."""
import csv
import math
import pathlib
import re
import subprocess
import sys
import tempfile

EXE = str(pathlib.Path(sys.argv[1]).resolve())
EXAMPLES = pathlib.Path(sys.argv[2]).resolve()

def read_grid(path, multiblock=False, dim=2):
    tok = path.read_text().split()
    pos = 0
    if multiblock:
        assert tok[pos] == '1'
        pos += 1
    ni, nj = map(int, tok[pos:pos+2]); pos += 2
    if dim == 3:
        assert tok[pos] == '1'; pos += 1
    n = ni*nj
    a = list(map(float, tok[pos:]))
    assert len(a) == dim*n and all(math.isfinite(v) for v in a)
    if dim == 3:
        assert all(v == 0.0 for v in a[2*n:])
    g = [[(a[j*ni+i], a[n+j*ni+i]) for i in range(ni)] for j in range(nj)]
    assert all(r[0] == r[-1] for r in g)
    return g

def cross(a,b):
    return a[0]*b[1]-a[1]*b[0]

def sub(a,b):
    return a[0]-b[0], a[1]-b[1]

def min_scaled_j(g, sign=1):
    z=1.0
    for r,s in zip(g,g[1:]):
        for i in range(len(r)-1):
            e,f=sub(r[i+1],r[i]),sub(s[i+1],s[i]); u,v=sub(s[i],r[i]),sub(s[i+1],r[i+1])
            for a,b in ((e,u),(e,v),(f,u),(f,v)):
                z=min(z,sign*cross(a,b)/(math.hypot(*a)*math.hypot(*b)))
    assert z>0.0, z
    return z

def run(wall,out,*args,success=True):
    p=subprocess.run([EXE,'--input',str(wall),'--output',str(out),*args],capture_output=True,text=True)
    if success:
        assert p.returncode==0,p.stdout+p.stderr
    else:
        assert p.returncode!=0 and not out.exists(),p.stdout+p.stderr
    return p

def write_ring(path,r):
    with path.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['x','y']);w.writerows(r)

with tempfile.TemporaryDirectory() as d:
    d=pathlib.Path(d)
    # Constant-normal-speed analytic solution, with periodic seam and geometric growth.
    wall=EXAMPLES/'circle.csv'; out=d/'circle.xyz'; h=0.001; q=1.15; n=20
    run(wall,out,'--first-height',str(h),'--growth',str(q),'--layers',str(n),'--quiet')
    g=read_grid(out);assert (len(g[0]),len(g))==(121,21)
    xi=0;err=0
    for j,r in enumerate(g):
        if j: xi+=h*q**(j-1)
        err=max(err,max(abs(math.hypot(*p)-(0.5+xi)) for p in r))
    assert err<1e-10,err
    print(f'circle / physical: max radius error {err:.3e}; min scaled Jacobian {min_scaled_j(g):.6f}')
    # Original alpha=1 equation, no diffusion: r^2=R^2+2R_eff xi.
    # R_eff accounts for initial eta defined by polygon edge length instead of true arc.
    out=d/'original.xyz';h=0.002;q=1.12;n=20
    run(wall,out,'--first-height',str(h),'--growth',str(q),'--layers',str(n),'--mode','original','--diffusion','0','--integrator','rk4','--quiet')
    g=read_grid(out);reff=0.5*math.sin(math.pi/120)/(math.pi/120);xi=0;err=0
    for j,r in enumerate(g):
        if j:xi+=h*q**(j-1)
        rad=math.sqrt(0.25+2*reff*xi)
        err=max(err,max(abs(math.hypot(*p)-rad) for p in r))
    assert err<2e-9,err
    print(f'circle / original: max radius error {err:.3e}; min scaled Jacobian {min_scaled_j(g):.6f}')
    # Sharp trailing edge; preserve coordinates and enforce first height.
    out=d/'naca.xyz';wall=EXAMPLES/'naca0012.csv';h=1e-4
    run(wall,out,'--first-height',str(h),'--growth','1.15','--layers','40','--quiet')
    g=read_grid(out);assert (len(g[0]),len(g))==(161,41)
    with wall.open() as f: raw=[tuple(map(float,r)) for r in list(csv.reader(f))[1:]]
    assert g[0]==[raw[0],*raw[-2:0:-1],raw[0]]
    err=max(abs(math.dist(a,b)-h) for a,b in zip(g[0],g[1]))
    assert err<1e-12,err
    print(f'NACA0012: max first-height error {err:.3e}; min scaled Jacobian {min_scaled_j(g):.6f}')
    # Nonuniform periodic grid, alternate integrators/stencils and Plot3D 3D multiblock.
    wall=d/'nonuniform.csv';t=[2*math.pi*i/100 for i in range(100)]
    r=[(2+0.8*math.cos(x+0.35*math.sin(x)),3+0.8*math.sin(x+0.35*math.sin(x))) for x in t]
    write_ring(wall,r)
    for stencil,integ in [('compact','rk2'),('three-point','euler')]:
        out=d/f'{stencil}.xyz'
        run(wall,out,'--first-height','0.0001','--growth','1.05','--layers','6','--format','3d','--multiblock','--stencil',stencil,'--integrator',integ,'--i-order','ccw','--quiet')
        g=read_grid(out,True,3);assert g[0]==r+[r[0]]
        assert len(g)==7 and len(g[0])==101
        print(f'nonuniform / {stencil} / {integ}: min scaled Jacobian {min_scaled_j(g,-1):.6f}')
    # Linearized implicit theta scheme: test temporal convergence against the
    # exact radius for the *three-point spatially discretized* circle equation.
    for theta,order in [(1.0,1),(0.5,2)]:
        errors=[]
        for n in (10,20,40):
            out=d/f'implicit-circle-{theta}-{n}.xyz'
            run(EXAMPLES/'circle.csv',out,'--mode','original','--diffusion','0','--theta',str(theta),'--first-height',str(0.2/n),'--growth','1','--layers',str(n),'--quiet')
            g=read_grid(out);rad=math.sqrt(0.25+0.2/math.cos(math.pi/120))
            errors.append(max(abs(math.hypot(*p)-rad) for p in g[-1]))
            min_scaled_j(g)
        assert all(0.8*2**order<a/b<1.3*2**order for a,b in zip(errors,errors[1:])),errors
        print(f'implicit circle theta={theta}: errors {errors}')
    # The clustered 60-layer wing should need only one solve per marched layer.
    wall=EXAMPLES/'naca0012_clustered.csv'
    for theta in ('1','0.5'):
        out=d/f'clustered-{theta}.xyz'
        p=run(wall,out,'--first-height','0.0001','--growth','1.15','--layers','60','--theta',theta,'--smoothing-start','0','--smoothing-end','0')
        g=read_grid(out);assert (len(g[0]),len(g))==(161,61)
        steps=list(map(int,re.findall(r'steps=(\d+)',p.stdout)))
        assert len(steps)==60 and sum(steps)==60,steps
        err=max(abs(math.dist(a,b)-1e-4) for a,b in zip(g[0],g[1]))
        assert err<1e-12 and min_scaled_j(g)>0.1
        print(f'clustered NACA0012 theta={theta}: 60 advances, min scaled Jacobian {min_scaled_j(g):.6f}')
    # Tangential speed cap: optional, independent of theta, and not a hard
    # bound on finite-cell corner angles. Compare actual outer grid quality.
    def outer_mean_skew(g):
        angles=[]
        for r,s in zip(g[-11:-1],g[-10:]):
            for i in range(len(r)-1):
                e,f=sub(r[i+1],r[i]),sub(s[i+1],s[i]);u,v=sub(s[i],r[i]),sub(s[i+1],r[i+1])
                for a,b in ((e,u),(e,v),(f,u),(f,v)):
                    c=abs(a[0]*b[0]+a[1]*b[1])/(math.hypot(*a)*math.hypot(*b))
                    angles.append(math.degrees(math.asin(min(1.0,c))))
        return sum(angles)/len(angles)
    args=['--first-height','0.0001','--growth','1.15','--layers','60','--theta','0.5','--smoothing-start','0','--smoothing-end','0']
    baseline=read_grid(d/'clustered-0.5.xyz')
    for angle in ('90','10','5','0'):
        out=d/f'skew-{angle}.xyz'
        p=run(wall,out,*args,'--max-skew-angle',angle)
        g=read_grid(out);assert (len(g[0]),len(g))==(161,61)
        assert sum(map(int,re.findall(r'steps=(\d+)',p.stdout)))==60
        min_scaled_j(g)
        assert g[0]==baseline[0] and g[1]==baseline[1]
        avg=outer_mean_skew(g)
        if angle=='90': assert g==baseline
        if angle=='10': assert avg<0.8*outer_mean_skew(baseline)
        if angle=='5': assert avg<0.5*outer_mean_skew(baseline)
        print(f'clustered NACA0012 angle cap {angle}: outer mean skew {avg:.6f} degrees')
    out=d/'zero-diffusion.xyz';run(wall,out,*args,'--diffusion','0')
    g=read_grid(out);g0=read_grid(d/'skew-0.xyz')
    assert max(math.dist(a,b) for r,s in zip(g,g0) for a,b in zip(r,s))<1e-11
    # Explicit RK4 accepts the same cap. With cap zero, the ODE is pure normal.
    args=['--first-height','0.0001','--layers','5','--integrator','rk4','--stencil','three-point']
    run(wall,d/'explicit-cap0.xyz',*args,'--max-skew-angle','0')
    run(wall,d/'explicit-no-diff.xyz',*args,'--diffusion','0')
    g=read_grid(d/'explicit-cap0.xyz');g0=read_grid(d/'explicit-no-diff.xyz')
    assert max(math.dist(a,b) for r,s in zip(g,g0) for a,b in zip(r,s))<1e-9
    for angle in ('-1','90.1'):
        run(wall,d/'bad-angle.xyz',*args,'--max-skew-angle',angle,success=False)
    run(wall,d/'bad-smoothing.xyz',*args,'--smoothing','full','--max-skew-angle','5',success=False)
    # Current physical spacing smoothing: independent of wall arc-length metric.
    args=['--first-height','0.0001','--growth','1.15','--layers','60','--theta','0.5','--smoothing-start','0','--smoothing-end','0']
    for eps in ('0.002','0.03'):
        out=d/f'current-{eps}.xyz'
        p=run(wall,out,*args,'--smoothing-metric','current','--diffusion',eps)
        g=read_grid(out);assert (len(g[0]),len(g))==(161,61)
        assert g[:2]==baseline[:2]
        min_scaled_j(g)
        assert sum(map(int,re.findall(r'steps=(\d+)',p.stdout)))==60
        assert outer_mean_skew(g)<0.5*outer_mean_skew(baseline)
        spacing=[math.dist(a,b) for a,b in zip(g[-1],g[-1][1:])]
        if eps=='0.03': assert max(spacing)/min(spacing)<25
        print(f'current spacing eps={eps}: outer mean skew {outer_mean_skew(g):.6f}, spacing ratio {max(spacing)/min(spacing):.6f}')
    for angle in ('10','0'):
        out=d/f'current-cap-{angle}.xyz'
        run(wall,out,*args,'--smoothing-metric','current','--max-skew-angle',angle)
        g=read_grid(out);min_scaled_j(g)
        assert g[:2]==baseline[:2]
        if angle=='0':
            g0=read_grid(d/'zero-diffusion.xyz')
            assert max(math.dist(a,b) for r,s in zip(g,g0) for a,b in zip(r,s))<1e-11
    # Explicit mode uses the same RHS and current metric CFL estimate.
    args=['--first-height','0.0001','--layers','5','--integrator','rk4','--stencil','three-point']
    run(wall,d/'explicit-current.xyz',*args,'--smoothing-metric','current','--max-skew-angle','5')
    min_scaled_j(read_grid(d/'explicit-current.xyz'))
    run(wall,d/'bad-metric.xyz',*args,'--smoothing-metric','invalid',success=False)
    run(wall,d/'bad-metric-full.xyz',*args,'--smoothing-metric','current','--smoothing','full',success=False)
    # Near-wall protection: ramp the coefficient in nominal marching distance.
    args=['--first-height','0.0001','--growth','1.15','--layers','60','--theta','0.5']
    out=d/'ramp.xyz'
    p=run(wall,out,*args,'--smoothing-metric','current','--diffusion','0.03','--smoothing-start','0.005','--smoothing-end','0.05')
    g=read_grid(out);normal=read_grid(d/'zero-diffusion.xyz')
    assert g[:16]==normal[:16], 'First 15 cell layers must match pure normal marching'
    assert sum(map(int,re.findall(r'steps=(\d+)',p.stdout)))==60
    assert min_scaled_j(g)>0.1 and outer_mean_skew(g)<2.5
    spacing=[math.dist(a,b) for a,b in zip(g[-1],g[-1][1:])]
    assert max(spacing)/min(spacing)<22
    eps=list(map(float,re.findall(r'effective-diffusion=([0-9.eE+-]+)',p.stdout)))
    assert len(eps)==60 and eps[:15]==[0]*15 and eps[-1]==0.03
    assert all(a<=b for a,b in zip(eps,eps[1:]))
    print(f'near-wall ramp: 15 pure-normal layers, outer mean skew {outer_mean_skew(g):.6f}, spacing ratio {max(spacing)/min(spacing):.6f}')
    out=d/'ramp-off.xyz'
    run(wall,out,*args,'--smoothing-metric','current','--diffusion','0.03','--smoothing-start','0','--smoothing-end','0')
    assert read_grid(out)==read_grid(d/'current-0.03.xyz')
    # Auto near-wall protection is now the default in height mode, including
    # current and wall metrics. Explicit zero/zero recovers the old behavior.
    autoargs=['--first-height','0.0001','--growth','1.15','--layers','60','--theta','0.5']
    for metric,eps0 in [('wall','0.002'),('current','0.03')]:
        out=d/f'auto-{metric}.xyz'
        p=run(wall,out,*autoargs,'--smoothing-metric',metric,'--diffusion',eps0)
        g=read_grid(out);assert g[:16]==normal[:16]
        min_scaled_j(g)
        assert 'smoothing-ramp=auto smoothing-start=0.005 smoothing-end=0.05' in p.stdout
        schedule=list(map(float,re.findall(r'effective-diffusion=([0-9.eE+-]+)',p.stdout)))
        assert schedule[:15]==[0]*15 and schedule[-1]==float(eps0)
        nodal=float(re.search(r'Wall orthogonality \(deg\): nodal max=([0-9.eE+-]+)',p.stdout)[1])
        assert nodal<1e-7,nodal
        manual=d/f'manual-{metric}.xyz'
        run(wall,manual,*autoargs,'--smoothing-metric',metric,'--diffusion',eps0,'--smoothing-start','0.005','--smoothing-end','0.05')
        assert g==read_grid(manual)
        print(f'auto {metric}: first 15 layers pure normal, wall nodal deviation {nodal:.3e} deg')
    # Scale the wall and the first layer by the same factor. Automatic distances
    # and normalized mesh should be scale-invariant (up to round-off).
    with wall.open() as f: pts=[tuple(map(float,r)) for r in list(csv.reader(f))[1:]]
    scaled_wall=d/'scaled-wall.csv';write_ring(scaled_wall,[(4*x,4*y) for x,y in pts])
    out=d/'scaled-auto.xyz'
    p=run(scaled_wall,out,'--first-height','0.0004','--growth','1.15','--layers','60','--theta','0.5','--smoothing-metric','current','--diffusion','0.03')
    assert 'smoothing-ramp=auto smoothing-start=0.02 smoothing-end=0.2' in p.stdout
    g=read_grid(out);ref=read_grid(d/'auto-current.xyz')
    assert max(math.dist((a[0]/4,a[1]/4),b) for r,s in zip(g,ref) for a,b in zip(r,s))<1e-10
    # Original mode continues to opt out unless explicit ramp distances are given.
    circle=EXAMPLES/'circle.csv';origargs=['--mode','original','--first-height','0.001','--growth','1','--layers','4']
    p=run(circle,d/'original-auto.xyz',*origargs)
    assert 'smoothing-ramp=off' in p.stdout
    run(circle,d/'original-off.xyz',*origargs,'--smoothing-start','0','--smoothing-end','0')
    assert read_grid(d/'original-auto.xyz')==read_grid(d/'original-off.xyz')
    # Explicit integrator also honors the same coefficient schedule.
    args=['--first-height','0.0001','--layers','5','--integrator','rk4','--stencil','three-point']
    run(wall,d/'explicit-ramp.xyz',*args,'--smoothing-metric','current','--diffusion','0.03','--smoothing-start','0.005','--smoothing-end','0.05')
    assert read_grid(d/'explicit-ramp.xyz')==read_grid(d/'explicit-no-diff.xyz')
    for a,b in [('0.05','0.005'),('0.01','0.01'),('-0.01','0.05'),('0.01','0')]:
        run(wall,d/'bad-ramp.xyz',*args,'--smoothing-start',a,'--smoothing-end',b,success=False)
    # Explicit max-step still enforces substeps in implicit mode (except exact first layer).
    out=d/'subdivided.xyz'
    p=run(EXAMPLES/'circle.csv',out,'--first-height','0.001','--growth','1','--layers','2','--max-step','0.0003')
    assert list(map(int,re.findall(r'steps=(\d+)',p.stdout)))==[1,4]
    min_scaled_j(read_grid(out))
    run(wall,d/'bad-implicit.xyz','--first-height','0.005','--growth','1.2','--layers','2','--max-substeps','16','--smoothing-start','0','--smoothing-end','0',success=False)
    run(wall,d/'bad-stencil.xyz','--first-height','0.0001','--layers','1','--stencil','compact',success=False)
    # Reject malformed input, clockwise boundaries and oversized invalid first layer.
    bad=d/'clockwise.csv';write_ring(bad,r[::-1]);run(bad,d/'bad.xyz','--first-height','0.001','--layers','1',success=False)
    bad=d/'duplicate.csv';write_ring(bad,r[:2]+r[1:]);run(bad,d/'duplicate.xyz','--first-height','0.001','--layers','1',success=False)
    bad=d/'crossing.csv';write_ring(bad,[(0,0),(1,0),(0,1),(1,1),(0.2,0.5),(0,2)])
    run(bad,d/'crossing.xyz','--first-height','0.001','--layers','1',success=False)
    print('invalid CSV/clockwise/self-intersection rejection: passed')
print('All regression checks passed.')
