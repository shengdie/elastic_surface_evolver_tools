#!/usr/bin/env python3

import numpy as np
from scipy.spatial import Delaunay

ls = []
a = (1,0)
b = (1/2, np.sqrt(3)/2)
use_fix_len = False

def create_VEF(R,r):
    v=[]
    e=[]
    f=[]
    a=1/2
    b=np.sqrt(3)/2
    R_=2*R
    r_=2*r
    numheight=int(2/np.sqrt(3)*R)+1
    # add a hexagon
    for i in range(R+1):
        if i %2 ==0:
            for j in range(-R+i//2,R-i//2+1):
                v.append((2*j*a, i*b))
        else:
            for j in range(-R+i//2,R-i//2):
                v.append(((2*j+1)*a, i*b))
    for i in range(1,R+1):
        if i %2 ==0:
            for j in range(-R+i//2,R-i//2+1):
                v.append((2*j*a, -i*b))
        else:
            for j in range(-R+i//2,R-i//2):
                v.append(((2*j+1)*a, -i*b))
    # add extra triangles to make it a circle

# square of vector
def norm2(vec):
    return vec[0]**2 + vec[1]**2 + vec[0]*vec[1]
def create_VFE_part(R,r):
    a = (1,0)
    b = (1/2, np.sqrt(3)/2)
    n_ext_out=int(2/np.sqrt(3)*R)+1
    n_ext_in = int(2/np.sqrt(3)*r)+1
    R2 = R**2
    r2 = r**2
    v=[]
    e=[]
    f=[]
    # 1/6 of circle
    for i in range(r+1, R+1):
        for j in range(0, i):
            v.append((i-j, j))
    v.remove((R,0))
    #if n_ext_out - R > 1:
    for i in range(R+1, n_ext_out):
        for j in range(i):
            vec = (i-j,j)
            if norm2(vec) < R2:
                v.append(vec)
    #if n_ext_in - r > 1:
    for i in range(r+1, n_ext_in):
        for j in range(i):
            vec = (i-j, j)
            if norm2(vec) < r2:
                v.remove(vec)
    vv=[(v[i][0]+b[0]*v[i][1], b[1]*v[i][1]) for i in range(len(v))]
    num_except_outter= len(vv)
    # add circle outter boundary
    angle=np.pi/3/(R+1)
    for a in range(0, R+1):
        aa = a*angle
        vv.append((R*np.cos(aa), R*np.sin(aa)))
    
    # add circle inner boundary
    num_except_inner=len(vv)
    num_all_part = num_except_inner +r
    angle=np.pi/3/r
    for a in range(0, r):
        aa = a*angle
        vv.append((r*np.cos(aa), r*np.sin(aa)))
    vv = np.array(vv)
    #"""
    vo = [vv]
    vvt=np.transpose(vv)
    for a in range(1, 6):
        aa = np.pi/3*a
        T = np.array([
            [np.cos(aa), -np.sin(aa)],
            [np.sin(aa), np.cos(aa)]])
        vo.append(np.transpose(np.dot(T, vvt)))
    vv = np.concatenate(vo)
    # id of inner boundary verts
    inner_verts_id = [num_except_inner+i+j*num_all_part for i in range(0,r) for j in range(6)]
    outter_verts_id = [num_except_outter+i+j*num_all_part for i in range(0,R+1) for j in range(6)]
    
    # construct mesh
    tri = Delaunay(vv)
    simp = []
    # remove the inner simplices
    for s in tri.simplices:
        if any(i not in inner_verts_id for i in s):
            simp.append(s)
    
    #e = []
    ee = {}
    for s in simp:
        for ed in [(s[0]+1, s[1]+1), (s[1]+1, s[2]+1), (s[2]+1, s[0]+1)]: 
            if not any(em in e for em in [ed, (ed[1], ed[0])]):
                e.append(ed)
    for i in range(len(e)):
        ee[e[i]] = i+1
        ee[(e[i][1], e[i][0])] = -i-1 
    for s in simp:
        f.append((ee[(s[0]+1, s[1]+1)], ee[(s[1]+1, s[2]+1)], ee[(s[2]+1, s[0]+1)]))
    
    ls.append("vertices\n")
    for i in range(len(vv)):
        if i in inner_verts_id:
            if use_fix_len:
                l = '{0}\t{1}*c_ratio {2}*c_ratio {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, vv[i][0], vv[i][1], 0)
                if vv[i][0] == r and vv[i][1] == 0:
                    l+=' fixed'
            else:
                if vv[i][0] != r or vv[i][1] != 0:
                    l = '{0}\t{1}*c_ratio {2}*c_ratio {3} ref_coord {{ {1} {2} {3} }} constraint 1'.format(i+1, vv[i][0], vv[i][1], 0)
                else:
                    l = '{0}\t{3} boundary 1 ref_coord {{ {1} {2} {3} }}'.format(i+1, vv[i][0], vv[i][1], 0)
                #l+=' '
        elif i in outter_verts_id:
            l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, vv[i][0], vv[i][1], 0)
        else:
            l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }} bend gbend'.format(i+1, vv[i][0], vv[i][1], 0)
        l+='\n'
        ls.append(l)
    ls.append('\nedges\n')
    for i in range(len(e)):
        v1 = e[i][0]-1
        v2 = e[i][1]-1
        l = '{0}\t{1} {2}'.format(i+1, e[i][0], e[i][1])
        if v1 in inner_verts_id and v2 in inner_verts_id:
            l += ' border_in 1'
            if use_fix_len:
                l+= ' fixlen'
        elif v1 in outter_verts_id and v2 in outter_verts_id:
            l += ' border_out 1'
        l+='\n'
        ls.append(l)
    ls.append('\nfaces\n')
    for i in range(len(f)):
        l = '{0}\t {1} {2} {3}\n'.format(i+1, f[i][0], f[i][1], f[i][2])
        ls.append(l)

    
    """
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.triplot(vv[:,0], vv[:,1], simp)
    ax.plot(vv[:,0], vv[:,1], 'o')
    ax.set_aspect(1)
    plt.show()
    """

def create_header(R, r):
    #l.append("PARAMETER lngth = ", n)
    ls.append("PARAMETER refine_times = 0\nPARAMETER run_time = 240\nPARAMETER submit_times = 0\nPARAMETER thseek = 0\n\n")
    #l.append("PARAMETER tstep = 0.0005 // this the step of changing the tension\n\n")
    ls.append("PARAMETER radout = {}\n".format(R))
    ls.append("PARAMETER radin = {}\n".format(r))
    ls.append("PARAMETER gridsize = {}\n".format(b[1]))

    header ="""PARAMETER c_ratio = 0.98
PARAMETER c_r = c_ratio * radin
PARAMETER thicknezz = (radin+radout)/2/1000
//OPTIMIZING_PARAMETER leftEndX = 0

// For resetting form_factors in refining
define vertex attribute ref_coord real[3]
//define vertex attribute on_border integer
//define vertex attribute corner integer
define vertex attribute old_vid  integer
define edge   attribute old_eid  integer

define facet attribute poisson_ratio real
define facet attribute form_factors real[3]
quantity stretch energy method linear_elastic global

quantity bend energy method star_perp_sq_mean_curvature
quantity gbend energy method star_gauss_curvature

//quantity fixlen conserved method edge_length modulus 1

define edge attribute border_in integer //this is real border
define edge attribute border_out integer
//define edge attribute righttension_mark integer

//define edge attribute divide   integer
//define edge attribute rheight  integer


boundary 1 parameters 1
x1: c_r
x2: 0
x3: p1

constraint 1
formula: x^2 + y^2 = c_r^2

//y_mirror
//constraint 2
//    formula: y = 0
//boundary 2 parameters 1
//x1: rightEndX
//x2: p1*cos(delta)
//x3: p1*sin(delta)

"""
    ls.append(header)

def creatTail():
    tail = """
read

set facet tension 0
set facet poisson_ratio 0.4


read "refinement_large.cmd"

set_thickness(thicknezz)

//transform_expr "a"
read "diagnostics.cmd"

U
hessian_normal off
check_increase off
suppress_warning 1825  // suppress the warning about hessian not positively defined
printf "Outter radius := %g\\n", radout
printf "Inner radius := %g\\n", radin
printf "radout/radin := %g\\n", radout/radin
printf "constricted r (ratio) := %g\\n", c_ratio

    """
    ls.append(tail)

def main(R, r, oname, use_fix_len_=False):
    #use_fix_len = use_fix_len_
    use_fix_len = False
    ls.append("PARAMETER origin_name = \"{0}\"\n".format(oname))
    create_header(R, r)
    create_VFE_part(R, r)
    creatTail()
    s = ''.join(ls)
    with open(oname + '.fe', 'w') as f:
        f.write(s) 

#main(2, 3, 'test')
if __name__ == "__main__":
    import sys
    usage="""usage: SCRIPT R r "filename_without_extension" use_fix_len
R: the outter radius (int)
r: the inner radius (int)
          """
    if len(sys.argv) == 1:
        print(usage)
        sys.exit(1) 
    try:
        R = int(sys.argv[1])
        r = int(sys.argv[2])
    except (TypeError, IndexError):
        sys.stderr.write("must has 3 variables\n")
        print(usage)
        sys.exit(2)
    if len(sys.argv) == 5:
        main(R, r, str(sys.argv[3]), use_fix_len_=True)
    else:
        main(R, r, str(sys.argv[3]))


#create_VFE_part(10,3)

        