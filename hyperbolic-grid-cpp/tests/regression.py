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
        p=run(wall,out,'--first-height','0.0001','--growth','1.15','--layers','60','--theta',theta)
        g=read_grid(out);assert (len(g[0]),len(g))==(161,61)
        steps=list(map(int,re.findall(r'steps=(\d+)',p.stdout)))
        assert len(steps)==60 and sum(steps)==60,steps
        err=max(abs(math.dist(a,b)-1e-4) for a,b in zip(g[0],g[1]))
        assert err<1e-12 and min_scaled_j(g)>0.1
        print(f'clustered NACA0012 theta={theta}: 60 advances, min scaled Jacobian {min_scaled_j(g):.6f}')
    # Explicit max-step still enforces substeps in implicit mode (except exact first layer).
    out=d/'subdivided.xyz'
    p=run(EXAMPLES/'circle.csv',out,'--first-height','0.001','--growth','1','--layers','2','--max-step','0.0003')
    assert list(map(int,re.findall(r'steps=(\d+)',p.stdout)))==[1,4]
    min_scaled_j(read_grid(out))
    run(wall,d/'bad-implicit.xyz','--first-height','0.005','--growth','1.2','--layers','2','--max-substeps','16',success=False)
    run(wall,d/'bad-stencil.xyz','--first-height','0.0001','--layers','1','--stencil','compact',success=False)
    # Reject malformed input, clockwise boundaries and oversized invalid first layer.
    bad=d/'clockwise.csv';write_ring(bad,r[::-1]);run(bad,d/'bad.xyz','--first-height','0.001','--layers','1',success=False)
    bad=d/'duplicate.csv';write_ring(bad,r[:2]+r[1:]);run(bad,d/'duplicate.xyz','--first-height','0.001','--layers','1',success=False)
    bad=d/'crossing.csv';write_ring(bad,[(0,0),(1,0),(0,1),(1,1),(0.2,0.5),(0,2)])
    run(bad,d/'crossing.xyz','--first-height','0.001','--layers','1',success=False)
    print('invalid CSV/clockwise/self-intersection rejection: passed')
print('All regression checks passed.')
