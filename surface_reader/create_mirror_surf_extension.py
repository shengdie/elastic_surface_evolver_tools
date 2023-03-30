#!/usr/bin/env python3
#mport numpy as np
import math

ls = []
a = 0.5 # shorter leg
b = math.sqrt(3)/2 # longer leg
c = 1 # hypotenuse
def create_VEF(nx,ny):
    """Create vertex, edges, faces
    
    Arguments:
        nx {int} -- number of equilateral triangle one row
        ny {int} -- number of equilateral triangle one column, (height of e-triangle)
    """

    #for y in range(1,ny+2):
    #
    n_whole_edge = (ny+2)//2
    n_whole_row_vers = nx+1
    n_half_row_vers = nx + 2 
    ls.append('vertices\n')
    for y in range(n_whole_edge):
        for x in range(nx+1):
            l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(x+1 + n_half_row_vers*y+ n_whole_row_vers*y,
                                                                     x*c, 2*y*b, 0)
            if x == 0:
                l += ' constraint 1'
            elif x == n_whole_row_vers -1:
                l += ' boundary 1 fixed'
            if y == 0 and x < n_whole_row_vers-1:
                l += ' constraint 2'
            if ny %2 ==0: # the boundary has whole edges
                if y < n_whole_edge - 1 and x < n_whole_row_vers -1: # the boundary no bend and gbend
                    l += ' bend gbend'
            elif x < n_whole_row_vers:
                l += ' bend gbend'
            ls.append(l)
    n_half_edge = (ny+1)//2
    
    for y in range(n_half_edge):
        for x in range(n_half_row_vers):
            if x == 0:
                l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(x+1 + n_whole_edge*(y+1) + n_half_edge*y,
                                                                            x, (2*y+1)*b, 0)
            elif x == n_half_row_vers - 1:
                l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(x+1 + n_whole_edge*(y+1) + n_half_edge*y, n_half_row_vers-2, (2*y+1)*b, 0)
            else:
                l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(x+1 + n_whole_edge*(y+1) + n_half_edge*y, x*c-a, (2*y+1)*b, 0)
            if x == 0:
                l += ' constraint 1'
            elif x == n_half_row_vers - 1:
                l += ' boundary 1 fixed'
            if ny % 2 == 1: # it is boundary
                if y < n_half_edge -1 and x < n_half_row_vers-1:
                    l+=' bend gbend'
            elif x < n_half_row_vers -1:
                l+= ' bend gbend'
            ls.append(l)
    ls.append('\nedges\n')

def create_VEF_m(nx, ny):
    nwl = (ny+2)//2
    nhl = (ny+1)//2
    nwv = nx+1
    nhv = nx+2
    n_v = nwv * nwl + nhl * nhv
    #v = np.zeros((n_v, 3)) # vertices
    #v = [0 for i in range(n_v)] # v has the coordinates of the vertices
    v = []
    e = [] # edges
    f = [] # faces
    #a = 0.5 # shorter leg
    #b = math.sqrt(3)/2 # longer leg
    #c = 1 # hypotenuse
    # vortices coor
    for y in range(nwl):
        for x in range(nwv):
            #v[x + nhv*y+ nwv*y] = (x*2, 2*y)
            v.append((x*2, 2*y))
    for y in range(nhl):
        for x in range(nhv):
            #v[x + nwl*(y+1) + nhl*y] = (2*x-1, 2*y+1)
            v.append((2*x-1, 2*y+1))
        #    if x == 0:
        #        v[x + nwl*(y+1) + nhl*y] = (x, (2*y+1), 0)
        #    elif x == nhv - 1:
        #        v[x + nwl*(y+1) + nhl*y] = ((nhv-2)*2, (2*y+1), 0)
        #    else:
        #        v[x + nwl*(y+1) + nhl*y] = (2*x-1, (2*y+1), 0)
    vv = {v[i]:i+1 for i in range(n_v)} # dict for id of vertices

    # edges
    for y in range(nwl):
        for x in range(nwv):
            if 2*y+1 <= ny:
                e.append((vv[(2*x, 2*y)], vv[(2*x+1, 2*y+1)]))
            if 2*y-1 >= 0:
                e.append((vv[(2*x, 2*y)], vv[(2*x+1, 2*y-1)]))
            if x < nx:
                e.append((vv[(2*x, 2*y)], vv[(2*(x+1), 2*y)]))
    for y in range(nhl):
        for x in range(nwv):
            if 2*y+2 <= ny:
                e.append((vv[(2*x-1, 2*y+1)], vv[(2*x, 2*y+2)]))
            e.append((vv[(2*x-1, 2*y+1)], vv[(2*x, 2*y)]))
            e.append((vv[(2*x-1, 2*y+1)], vv[(2*x+1, 2*y+1)]))
    ee = {e[i]:i+1 for i in range(len(e))}

    # faces
    for y in range(nwl):
        #if ny % 2 == 0: # boundary if whole line
        for x in range(nwv-1):
            pa = (2*x, 2*y)
            pb = (2*x+2, 2*y)
            if 2*y+1 <= ny:
                #pa = (2*x, 2*y)
                #b = (2*x+2, 2*y)
                pc = (2*x+1, 2*y+1)
                f.append((ee[(vv[pa], vv[pb])],
                            -ee[(vv[pc], vv[pb])],
                            -ee[(vv[pa], vv[pc])]))
            if 2*y-1 >= 0:
                pc = (2*x+1, 2*y-1)
                f.append((-ee[(vv[pa], vv[pb])],
                        ee[(vv[pa], vv[pc])],
                        ee[(vv[pc], vv[pb])]))
    for y in range(nhl):
        for x in range(nhv-1):
            pa = (2*x-1, 2*y+1)
            pc = (2*x, 2*y)
            pb = (2*x+1, 2*y+1)
            f.append((-ee[(vv[pa], vv[pb])],
                    ee[(vv[pa], vv[pc])],
                    ee[(vv[pc], vv[pb])]))
            if 2*y+2 <= ny:
                pc = (2*x, 2*y+2)
                f.append((ee[(vv[pa], vv[pb])],
                        -ee[(vv[pc], vv[pb])],
                        -ee[(vv[pa], vv[pc])]))
    ls.append("vertices\n")
    for i in range(len(v)):
        if v[i][0] <= 0:
            l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, 0, v[i][1]*b, 0)
        elif v[i][0] >= 2*nx:
            l = '{0}\t{2} boundary 1 fixed ref_coord {{ {1} {2} {3} }}'.format(i+1, 2*nx*a, v[i][1]*b, 0)
        else:
            l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, v[i][0]*a, v[i][1]*b, 0)
        if v[i][0] <= 0:
            l+= ' constraint 1'
        if v[i][1] == 0 and v[i][0] < 2* nx:
            l+= ' constraint 2'
        #if v[i][0] >= 2*nx:
        #    l += ' fixed boundary 1'
        if v[i][1] < ny and v[i][0] < 2*nx:
            l += ' bend gbend'
        #if v[i][0] == 2*nx and v[i][1] == 0:
        #    l += ' righttension'
        l+='\n'
        ls.append(l)
    ls.append("\nedges\n")
    for i in range(len(e)):
        v1 = e[i][0] -1
        v2 = e[i][1] -1
        l = '{0}\t{1} {2}'.format(i+1, v1+1, v2+1)
        if v[v1][0] <= 0 and v[v2][0] <= 0:
            l += ' constraint 1 mborder 1'
        if v[v1][1] == 0 and v[v2][1] == 0:
            l += ' constraint 2 mborder 1'
        if v[v1][0] >= 2*nx and v[v2][0] >= 2*nx:
            l += ' border 1 mborder 1 boundary 1'
        elif v[v1][1] ==ny and v[v2][1] == ny:
            l += ' border 1 mborder 1'
        l += '\n'
        ls.append(l)
    ls.append('\nfaces\n')
    for i in range(len(f)):
        l = '{0}\t {1} {2} {3}\n'.format(i+1, f[i][0], f[i][1], f[i][2])
        ls.append(l)

def create_header(nx, ny):
    #l.append("PARAMETER lngth = ", n)
    ls.append("PARAMETER refine_times = 0\nPARAMETER run_time = 240\nPARAMETER submit_times = 0\nPARAMETER thseek = 0\n\n")
    #l.append("PARAMETER tstep = 0.0005 // this the step of changing the tension\n\n")

    ls.append("PARAMETER wd =  {}\n".format(b* ny))
    ls.append("PARAMETER gridsize = {}\n".format(b))
    ls.append("PARAMETER lngth = {}\n".format(nx * c))

    header ="""PARAMETER hypertenuse = sqrt(wd**2+(wd)**2/3)
PARAMETER thicknezz = 2*wd/10000
PARAMETER oth = thicknezz
PARAMETER ow = wd
PARAMETER pratio = 0.4

PARAMETER tensionG = 0.01
PARAMETER eratio = 1/sqrt(1-2*pratio*tensionG)
OPTIMIZING_PARAMETER rightEndX = lngth
//OPTIMIZING_PARAMETER leftEndX = 0

// For resetting form_factors in refining
define vertex attribute ref_coord real[3]
define vertex attribute corner integer
define vertex attribute old_vid  integer
define edge   attribute old_eid  integer

define facet attribute poisson_ratio real
define facet attribute form_factors real[3]
quantity stretch energy method linear_elastic global

quantity bend energy method star_perp_sq_mean_curvature
quantity gbend energy method star_gauss_curvature

define edge attribute border integer //this is real border
define edge attribute righttension_mark integer
define edge attribute mborder integer //mirror border
//define edge attribute lefttension_mark integer
define edge attribute sqcurve_string_mark integer
quantity freeedgebend energy method sqcurve_string_marked
    parameter_1: 2

//quantity righttension energy method vertex_scalar_integral
//    scalar_integrand: -(x - ref_coord[1]) * ow


define edge attribute divide   integer
define edge attribute rheight  integer

boundary 1 parameters 1
x1: rightEndX
x2: p1 * eratio
x3: 0

//x_mirror
constraint 1
    formula: x = 0

//y_mirror
constraint 2
    formula: y = 0
//boundary 2 parameters 1
//x1: rightEndX
//x2: p1*cos(delta)
//x3: p1*sin(delta)

view_transform_generators 1
 1 0 0 0
 0 -1 0 0
 0 0 1 0
 0 0 0 1

"""
    ls.append(header)

def creatTail():
    tail = """
read

set facet tension 0
set facet poisson_ratio pratio

//lefttension.modulus := tensionG
//righttension.modulus := tensionG
freeedgebend.modulus := 0

read "refinement_large.cmd"

set_thickness(thicknezz)

transform_expr "a"
read "diagnostics.cmd"

U
hessian_normal off
check_increase off
suppress_warning 1825  // suppress the warning about hessian not positively defined
printf "Current thickness: %g\\n",thicknezz/2/wd
printf "Current tension: %g\\n", tensionG
printf "Aspect ratio: %g\\n", lngth/wd

    """
    ls.append(tail)

def main(nx, ny, oname):
    ls.append("PARAMETER origin_name = \"{0}\"\n".format(oname))
    create_header(nx,ny)
    create_VEF_m(nx, ny)
    creatTail()
    s = ''.join(ls)
    with open(oname + '.fe', 'w') as f:
        f.write(s)
            
            
        #x = nwv-1
        #e.append((v[2*x, 2*y, 0], v[2*x, 2*y+1, 0]))

#main(2, 3, 'test')
if __name__ == "__main__":
    import sys
    usage="""usage: SCRIPT nx ny "filename_without_extension"
nx: the total number of equilateral-triangle edges in one row
ny: number of heights
          """
    if len(sys.argv) == 1:
        print(usage)
        sys.exit(1) 
    try:
        nx = int(sys.argv[1])
        ny = int(sys.argv[2])
    except (TypeError, IndexError):
        sys.stderr.write("must has 3 variables\n")
        print(usage)
        sys.exit(2)
    main(nx, ny, str(sys.argv[3]))

            

