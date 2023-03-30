#!/usr/bin/env python3
import sys
import math
import argparse
from helper import deal_cmd

ls = []
a = 0.5 # shorter leg
b = math.sqrt(3)/2 # longer leg
c = 1 # hypotenuse


# for torus mode, periodic in y-direction, ny must be even, since the minimum period is 2 heights

def create_VEF_m(nx, ny, gap=False, edge_gap=False, energy_gap=False, 
                 edge_egap=False, edge_pr=False, fix_bwavlen=False, edge_slope=False, stretch=False, stretch_every=False, **kw):
    # the real vertices we actually need is ny - 1
    ny = ny - 1
    nwl = (ny+2)//2  # number of whole lines
    nhl = (ny+1)//2  # number of half lines
    nwv = nx+1       # number of vertex of whole lines
    nhv = nx+2       # number of vertex of half lines
    n_v = nwv * nwl + nhl * nhv   # number of all vertices

    v = []
    e = [] # edges
    f = [] # faces

    # the idea is whole edge lines and half edge lines, here shows the x coords, nx = number of whole edges, ny = number of heights
    # thus the nwl = ()
    # whole lines         *(0)    *(2)    *(4)
    #                   /    \   /   \   /   \
    # half lines      *(-1)   *(1)    *(3)    *(5)
    #                   \    /   \   /   \   /
    # whole lines         *(0)    *(2)    *(4)
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

    period_edges_id_p = [] # store edges that is in periodic boundary * + *
    #period_edges_id_n = [] # * - *
    # edges
    for y in range(nwl):
        for x in range(nwv):
            if 2*y+1 <= ny:
                e.append((vv[(2*x, 2*y)], vv[(2*x+1, 2*y+1)]))
            if 2*y-1 >= 0:
                e.append((vv[(2*x, 2*y)], vv[(2*x+1, 2*y-1)]))
            #else: # perioidc (* - *)
            #    e.append((vv[(2*x, 2*y)], vv[(2*x-1, ny)]))
            #    period_edges_id_n.append(len(e)-1)
            if x < nx:
                e.append((vv[(2*x, 2*y)], vv[(2*(x+1), 2*y)]))
    for y in range(nhl):
        for x in range(nwv):
            if 2*y+2 <= ny:
                e.append((vv[(2*x-1, 2*y+1)], vv[(2*x, 2*y+2)]))
            else: # periodic (* + *)
                e.append((vv[(2*x-1, 2*y+1)], vv[(2*x, 0)]))
                period_edges_id_p.append(len(e)-1)
                e.append((vv[(2*x+1, 2*y+1)], vv[(2*x, 0)])) # left up edges
                period_edges_id_p.append(len(e)-1)
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
            #else:
            #    pc = (2*x+1, ny)
                f.append((-ee[(vv[pa], vv[pb])],
                        ee[(vv[pa], vv[pc])],
                        ee[(vv[pc], vv[pb])]))
            else:
                pc = (2*x+1, ny)
                f.append((-ee[(vv[pa], vv[pb])],
                        -ee[(vv[pc], vv[pa])],
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
            else: # period, next is 
                pc = (2*x, 0)
                f.append((ee[(vv[pa], vv[pb])],
                    ee[(vv[pb], vv[pc])],
                    -ee[(vv[pa], vv[pc])]))
            
    
    # extra for period
    

    #ff = {f[i]:i+1 for i in range(len(f))}
    ls.append("vertices\n")
    for i in range(len(v)):
        # if v[i][1] == 0 and v[i][0] < 2*nx:
        #     l = '{0}\t{1} boundary 1 ref_coord {{ {1} {2} {3} }}'.format(i+1, v[i][0] * a, 0, 0)
        # elif v[i][1] == 0 and v[i][0] >= 2*nx:
        #     if not fix_bwavlen:
        #         l = '{0}\t{1} boundary 1 ref_coord {{ {1} {2} {3} }}'.format(i+1, v[i][0] * a, 0, 0)
        #     else:
        #         l = '{0}\t{2} boundary 4 ref_coord {{ {1} {2} {3} }} fixed'.format(i+1, v[i][0] * a, 0, 0)
        # elif v[i][1] == ny and v[i][0] < 2*nx:
        #     l = '{0}\t{1} boundary 2 ref_coord {{ {1} {2} {3} }}'.format(i+1, v[i][0] * a if v[i][0] > 0 else 0, ny * b, 0)
        # elif v[i][1] == ny and v[i][0] >= 2*nx:
        #     if not fix_bwavlen:
        #         l = '{0}\t{1} boundary 2 ref_coord {{ {1} {2} {3} }}'.format(i+1, 2*nx*a, ny * b, 0)
        #     else:
        #         l = '{0}\t upEndY boundary 4 ref_coord {{ {1} {2} {3} }} fixed'.format(i+1, 2*nx*a, ny *b, 0)
        # elif v[i][0] <=0:
        #     if not fix_bwavlen:
        #         l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, 0, v[i][1]*b, 0)
        #     else:
        #         l = '{0}\t{2} boundary 3 ref_coord {{ {1} {2} {3} }}'.format(i+1, 0, v[i][1]*b, 0)
            
        #     if edge_gap:
        #         l += ' constraint 1'
        #     if energy_gap:
        #         l += ' egap'
        #     if edge_egap:
        #         l += ' egap'
        #     if gap:
        #         l += ' constraint 1'
        # elif v[i][0] >= 2 * nx:
        #     if not fix_bwavlen:
        #         l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, 2 * nx * a, v[i][1]*b, 0)
        #     else:
        #         l = '{0}\t{2} boundary 4 ref_coord {{ {1} {2} {3} }}'.format(i+1, 2 * nx * a, v[i][1]*b, 0)
        #     if edge_gap:
        #         l += ' constraint 1'
        #     if edge_egap:
        #         l += ' egap'
        #     if energy_gap:
        #         l += ' egap'
        #     if gap:
        #         l += ' constraint 1'
        if v[i][0] <= 0:
            if fix_bwavlen:
                l = '{0}\t{2} boundary 1 ref_coord {{ {1} {2} {3} }}'.format(i+1, 0, v[i][1]*b, 0)
            elif stretch or edge_slope:
                l = '{0}\t{2} {3} boundary 1 ref_coord {{ {1} {2} {3} }}'.format(i+1, 0, v[i][1]*b, 0)
            else:
                l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, 0, v[i][1]*b, 0)
            # if not fix_bwavlen:
            #     l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, 0, v[i][1]*b, 0)
            # else:
            #     l = '{0}\t{2} boundary 1 ref_coord {{ {1} {2} {3} }}'.format(i+1, 0, v[i][1]*b, 0)
                
        elif v[i][0] >= 2 * nx:
            if fix_bwavlen:
                l = '{0}\t{2} boundary 2 ref_coord {{ {1} {2} {3} }}'.format(i+1, 2 * nx * a, v[i][1]*b, 0)
            elif stretch or edge_slope:
                l = '{0}\t{2} {3} boundary 2 ref_coord {{ {1} {2} {3} }}'.format(i+1, 2 * nx * a, v[i][1]*b, 0)
            else:
                l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, 2 * nx * a, v[i][1]*b, 0)
            # if not fix_bwavlen:
            #     l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, 2 * nx * a, v[i][1]*b, 0)
            # else:
            #     l = '{0}\t{2} boundary 2 ref_coord {{ {1} {2} {3} }}'.format(i+1, 2 * nx * a, v[i][1]*b, 0)
            # if stretch_every:
            #     l += ' righttension'
        else:
            l = '{0}\t{1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, v[i][0]*a, v[i][1]*b, 0)
        if energy_gap:
            l += ' egap'
        if gap:
            l += ' constraint 1'
        #if v[i][0] <= 0:
        #    l+= ' constraint 1'
        #if v[i][1] == 0 and v[i][0] < 2* nx:
        #    l+= ' constraint 2'


        if v[i][0]  <= 0 and v[i][1] == 0:
            l+= ' fixed'
        #elif v[i][0] <= 0  and v[i][1] == ny:
        #    l+= ' fixed'


        #if fix_bwavlen
        

        #if v[i][1] < ny and v[i][0] < 2*nx and v[i][1] > 0 and v[i][0] > 0:
        #    l += ' bend' # gbend'
        if stretch:
            if v[i][0] == 2*nx and v[i][1] == 0:
                l += ' righttension'
        l+='\n'
        ls.append(l)
    ls.append("\nedges\n")
    for i in range(len(e)):
        v1 = e[i][0] -1
        v2 = e[i][1] -1
        l = '{0}\t{1} {2}'.format(i+1, v1+1, v2+1)
        #if v[v1][1] == ny and v[v2][1] == 0:
        #    l += ' * + *'
        if i in period_edges_id_p:
            l += ' * + *'
        #elif i in period_edges_id_n:
        #    l += ' * - *'
        else:
            l += ' * * *'
        if v[v1][0] <= 0 and v[v2][0] <= 0:
            l += ' border 1'
            if fix_bwavlen or stretch or edge_slope:
                l += ' boundary 1'
            if edge_slope:
                l += ' ese'
        if v[v1][0] >= 2*nx and v[v2][0] >= 2*nx:
            l += ' border 1'
            if fix_bwavlen or stretch or edge_slope:
                l += ' boundary 2'
            if edge_slope:
                l += ' ese'
        # if v[v1][1] == 0 and v[v2][1] == 0:
        #     l += ' border 1 boundary 1'
        # elif v[v1][1] ==ny and v[v2][1] == ny:
        #     l += ' border 1 boundary 2'
        l += '\n'
        ls.append(l)
    ls.append('\nfaces\n')
    for i in range(len(f)):
        l = '{0}\t {1} {2} {3}\n'.format(i+1, f[i][0], f[i][1], f[i][2])
        ls.append(l)
    #ls.append('\nbodies\n')
    #ls.append('1 ' + ' '.join(str(i) for i in range(1, len(f)+1)) + '\n')

def create_header(nx, ny, substrate='linear', cmd=None, gap=False, edge_gap=False, 
                  prevent_reverse=False, edge_pr=False, energy_gap=False, 
                  edge_egap=False, fix_bwavlen=False, edge_slope=False, 
                  stretch=False, stretch_every=False, gbend=False, nonuni=True, pos=0.5):
    #l.append("PARAMETER lngth = ", n)
    #if edge_slope:
    #    "function real cryz(real y1, real z1, real y2, real z2);"
    if cmd is not None:
        # (rfc, rpc, rpcn), (dfc, dpc, dpcn) = cmd
        rfc, rpc, rpcn = cmd
        ls.extend([k + ';\n' for k in rfc.keys()])
        ls.extend([k + ';\n' for k in rpc.keys()])
        # ls.extend([k + ';\n' for k in dfc.keys()])
        # ls.extend([k + ';\n' for k in dpc.keys()])
        ls.append('\n')
    if nonuni:
        ls.append('PARAMETER max_refine_times = 0\n')
    ls.append('PARAMETER last_stime = 0\n')
    ls.append("PARAMETER refine_times = 0\nPARAMETER save_period = 10800\nPARAMETER run_time = 120\nPARAMETER submit_times = 1\nPARAMETER ost = submit_times - 1\nPARAMETER thseek = 0\n\n")
    #l.append("PARAMETER tstep = 0.0005 // this the step of changing the tension\n\n")

    ls.append("PARAMETER wd =  {}\n".format(b* ny))
    ls.append("PARAMETER gridsize = {}\n".format(b))
    ls.append("PARAMETER lngth = {}\n".format(nx * c))
    if substrate.casefold() == 'jump'.casefold():
        ls.append("PARAMETER pos = {}\n".format(pos))
    if edge_egap:
        ls.append("PARAMETER edgegap = 1 // if only edge egap\n")
    else:
        ls.append("PARAMETER edgegap = 0 // if only edge egap\n")
    if edge_pr:
        ls.append('PARAMETER edgepr = 1 // if prevent edge facet revert only\n')
    else:
        ls.append('PARAMETER edgepr = 0 // prevent edge facet revert only\n')
    if edge_gap:
        ls.append('PARAMETER edrgap = 1 // if real edge gap with constraint method\n')
    else:
        ls.append('PARAMETER edrgap = 0 // if real edge gap with constraint method\n')
    #if fix_bwavlen:
    #    ls.append('OPTIMIZING_PARAMETER a0 = 0\nOPTIMIZING_PARAMETER al\n')

    header ="""PARAMETER hypertenuse = sqrt(wd**2+(wd)**2/3)
PARAMETER thicknezz = wd/500
PARAMETER oth = thicknezz
PARAMETER ow = wd

PARAMETER tensionG = 0.01
PARAMETER fLen = sqrt(1+2*tensionG) * lngth
PARAMETER delta = -0.001
PARAMETER ww = wd * (1 + delta)
//PARAMETER rightEndX = lngth * (1 + delta)
PARAMETER upEndY = wd * (1 + delta)
PARAMETER num0 = 1
PARAMETER numl = 2
PARAMETER wl0 = ww/num0
PARAMETER wll = ww/numl
PARAMETER lf = (wll - wl0) /lngth
PARAMETER ls = wl0
PARAMETER k0 = 0
PARAMETER kl = 0
PARAMETER kf = (kl - k0) / lngth
"""
    if fix_bwavlen:
        header += 'OPTIMIZING_PARAMETER a0 = 0\nOPTIMIZING_PARAMETER al = 0\nOPTIMIZING_PARAMETER rEndX=lngth\n'
    elif stretch or edge_slope:
        header += 'OPTIMIZING_PARAMETER rEndX=lngth\n'
    else:
        header += 'PARAMETER a0 = wd/10\nPARAMETER al = wd/10\n'
    header +="""
//PARAMETER tti = a0 * num0 * pi/ww
//PARAMETER tta = al * numl * pi/ww
PARAMETER ttt = 2 * sqrt( - delta / (1+delta))
PARAMETER af = 1
PARAMETER as = 1
PARAMETER fa = gridsize * 0.5
PARAMETER shift = pi 

// periodic definition
torus
// periodic in three direction, only y-direction is useful, but we have to define the other two directions even we don't use it
// the period in y-direction is the confined width wd * (1+ delta)
"""
    if fix_bwavlen or stretch or edge_slope:
        header += """
periods
rEndX 0 0
0 upEndY 0
0 0 1
"""
    else:
        header += """
periods
lngth 0 0
0 upEndY 0
0 0 1
"""
    header += """
//OPTIMIZING_PARAMETER rightEndX = lngth
//OPTIMIZING_PARAMETER leftEndX = 0

// For resetting form_factors in refining
define vertex attribute ref_coord real[3]
define vertex attribute corner integer
define vertex attribute old_vid  integer
define edge attribute old_eid  integer
define edge attribute divide   integer
define edge attribute rheight  integer
define edge attribute border integer //this is real border

define facet attribute poisson_ratio real
define facet attribute form_factors real[3]
//define facet attribute ksub real
quantity stretch energy method linear_elastic global

quantity bend energy method star_perp_sq_mean_curvature global
"""
    if gbend:
        header += 'quantity gbend energy method star_gauss_curvature //global\n'
    if substrate.casefold() == 'linear'.casefold():
        # if stretch or fix_bwavlen:
        #     header += 'quantity subenergy energy method facet_scalar_integral global\n scalar_integrand: 1/(wl0 + (wll - wl0)/rEndX *x)^4 * z^2\n'
        # else:
        header += 'quantity subenergy energy method facet_scalar_integral global\n scalar_integrand: 1/(ls + lf*x)^4 * z^2\n'
        #if edge_egap:
        #    header += 'quantity egap energy modulus 10 method vertex_scalar_integral\n  scalar_integrand: floor(abs(z)/a0) * (abs(z) - a0)^2\n'
        if gap or edge_gap:
            header += 'constraint 1 nonpositive\n  formula: abs(z) - af * (x+as) = 0\n'
        if energy_gap or edge_egap:
            header += 'quantity egap energy modulus 10 method vertex_scalar_integral\n  scalar_integrand: floor(abs(z)/(af * (x+as))) * (abs(z) - (af * (x+as)))^2\n'
    elif substrate.casefold() == 'nlinear'.casefold():
        header += 'quantity subenergy energy method facet_scalar_integral global\n scalar_integrand: (lf * x + ls)^4  * z^2\n'
    elif substrate.casefold() == 'klinear'.casefold():
        header += 'quantity subenergy energy method facet_scalar_integral global\n scalar_integrand: (kf * x + k0) * z^2\n'
    elif substrate.casefold() == 'uniform'.casefold():
        header += 'quantity subenergy energy method facet_scalar_integral global\n scalar_integrand: z^2\n'
        if gap or edge_gap:
            header += 'constraint 1 nonpositive\n  formula: abs(z) = a0\n'
        if energy_gap or edge_egap:
            header += 'quantity egap energy modulus 10 method vertex_scalar_integral\n  scalar_integrand: floor(abs(z)/a0) * (abs(z) - a0)^2\n'
    elif substrate.casefold() == 'jump'.casefold():
        header += 'quantity subenergy energy method facet_scalar_integral global\n scalar_integrand: (k0 + kf * ceil((x - pos*{0})/{0})) * z^2\n'.format('rEndX' if fix_bwavlen or stretch else 'lngth')
    elif substrate.casefold() == 'sqrt':
        header += 'quantity subenergy energy method facet_scalar_integral global\n scalar_integrand: (ls + lf*x)^(-2) * z^2\n'
    else:
        raise Exception('No such K variance')
    
    if prevent_reverse:
        header += 'quantity cenergy energy modulus 1e9 method facet_scalar_integral global\n scalar_integrand: floor(1 - facet_normal[3]) * (facet_normal[3] / fa)^2\n'
    if edge_pr:
        header += 'quantity cenergy energy modulus 1e9 method facet_scalar_integral\n scalar_integrand: floor(1 - facet_normal[3]) * (facet_normal[3] / fa)^2\n'
    if edge_slope:
        header += 'quantity ese energy modulus 5e-3 method edge_scalar_integral\n scalar_integrand: abs(edge_vector[2]) * ttt > abs(edge_vector[3]) ? 0 : ((abs(edge_vector[2]) * ttt - abs(edge_vector[3]))/gridsize)^2 \n'
        #header += 'quantity es0 energy modulus 5e-3 method edge_scalar_integral\n scalar_integrand: abs(edge_vector[2]) * tti > abs(edge_vector[3]) ? 0 : ((abs(edge_vector[2]) * tti - abs(edge_vector[3]))/gridsize)^2 \n'
        #header += 'quantity esl energy modulus 53-3 method edge_scalar_integral\n scalar_integrand: abs(edge_vector[2]) * tta > abs(edge_vector[3]) ? 0 : ((abs(edge_vector[2]) * tta - abs(edge_vector[3]))/gridsize)^2 \n'
    if stretch:
        header += 'quantity righttension energy method vertex_scalar_integral\n scalar_integrand: -(x - fLen) * ow\n'
    elif stretch_every:
        header += 'quantity righttension energy method vertex_scalar_integral\n scalar_integrand: -(x - fLen) * gridsize\n'
    #if edge_egap:
    #    header += 'quantity egap energy modulus 10 method vertex_scalar_integral\n  scalar_integrand: floor(abs(z)/al) * (abs(z) - al)^2\n'
    #if gap:
    #    header += 'constraint 1 nonpositive\n  formula: abs(z) = a0\n'

    # we don't need boundary in periodic model
#     header += """boundary 1 parameters 1
# x1: p1
# x2: 0
# x3: 0

# boundary 2 parameters 1
# x1: p1
# x2: upEndY
# x3: 0
# """
    if fix_bwavlen:
        header += """boundary 1 parameters 1
x1: 0
x2: p1
x3: abs(a0) * sin(num0 * pi *2/ww * p1)

boundary 2 parameters 1
x1: rEndX
x2: p1
x3: abs(al) * sin(numl * pi *2/ww * p1 + shift)
"""
    elif stretch or edge_slope:
        header += """boundary 1 parameters 2
x1: 0
x2: p1
x3: p2

boundary 2 parameters 2
x1: rEndX
x2: p1
x3: p2
"""
    header += """

//constraint 1 nonpositive
//    formula: x = rightEndX
//constraint 2 nonnegative
//    formula: x = 0

//x_mirror
//constraint 1
//    formula: x = 0

//y_mirror
//constraint 2
//    formula: y = 0

/*
view_transform_generators 1
 -1 0 0 0
 0 1 0 0
 0 0 1 0
 0 0 0 1
 */

"""
    ls.append(header)

def creatTail(edge_egap=False, stretch=False, substrate='linear', fix_bwavlen=False, cmd=None, nonuni=False, edge_slope=False, **kw):
    tail = """
read

set facet tension 0
set facet poisson_ratio 0.4

ost := submit_times - 1

//lefttension.modulus := tensionG
//righttension.modulus := tensionG
//freeedgebend.modulus := 0

"""
    if cmd is not None:
        # (rfc, rpc, rpcn), (dfc, dpc, dpcn) = cmd
        rfc, rpc, rpcn = cmd
        
        # pcn.update(dpcn)
        if substrate.casefold() == 'jump'.casefold():
            for k, v in rpc.items():
                if 'set_ksub' in k:
                    rpc[k] = r"""procedure set_ksub(real lnum0, real lnuml) {
    num0 := lnum0;
    numl := lnuml;
    wl0 := ww / lnum0;
    wll := ww / lnuml;
    kl := get_ksub(numl); 
    k0 := get_ksub(num0);
    kf := kl - k0;
    recalc;
}
"""
        elif substrate.casefold() == 'sqrt'.casefold():
            for k, v in rpc.items():
                # print(k)
                if 'set_ksub' in k:
                    rpc[k] = r"""procedure set_ksub(real lnum0, real lnuml) {
    num0 := lnum0;
    numl := lnuml;
    wl0 := ww / lnum0;
    wll := ww / lnuml;
    kl := get_ksub(numl); 
    k0 := get_ksub(num0);
    ls := wl0^2 / 2 / sqrt(bend.modulus) / pi^2;
    lf := (num0^2 - numl^2) * ww^2 / 2 / sqrt(bend.modulus) / lngth / num0^2 / numl^2 / pi^2;
    recalc;
}
"""
            for k, v in rpcn.items():
                if 'set_awave' in k:
                    #print('hs')
                    rpcn[k] = r"""set_awave := {
    //local ww;
    local aa;
    local lnm;
    //ww := wd * (1 + delta);
    foreach vertex vv do {
        lnm := ((ls + lf*x)^(-2) * 4 / bend.modulus)^(1/4);
        //lnm := ww / (lf * (vv.x+ls));
        aa := 2/lnm * sqrt(- delta );
        vv.z := aa / 2 * (sin(lnm * vv.y) + sin(lnm * (ww - vv.y)));
    };
}
"""
        elif substrate.casefold() == 'nlinear'.casefold():
            for k, v in rpc.items():
                # print(k)
                if 'set_ksub' in k:
                    rpc[k] = r"""procedure set_ksub(real lnum0, real lnuml) {
    num0 := lnum0;
    numl := lnuml;
    wl0 := ww / lnum0;
    wll := ww / lnuml;
    lf := (numl - num0) / lngth;
    ls := num0;
    subenergy.modulus := (2 * pi / ww)^4 * bend.modulus/4; // remember K/2
    recalc;
}
"""
        elif substrate.casefold() == 'klinear'.casefold():
            for k, v in rpc.items():
                # print(k)
                if 'set_ksub' in k:
                    rpc[k] = r"""procedure set_ksub(real lnum0, real lnuml) {
    num0 := lnum0;
    numl := lnuml;
    wl0 := ww / lnum0;
    wll := ww / lnuml;
    kl := get_ksub(numl); 
    k0 := get_ksub(num0);
    kf := (kl - k0) / lngth;
    recalc;
}
"""
#         elif (fix_bwavlen or stretch) and substrate.casefold() == 'linear'.casefold():
#             for k, v in rpc.items():
#                 if 'set_ksub' in k:
#                     rpc[k] = r"""procedure set_ksub(real lnum0, real lnuml) {
#     num0 := lnum0;
#     numl := lnuml;

#     if lnum0 == lnuml then { 
#         print "uniform";
#         wl0 := ww / lnum0;
#         wll := ww / lnuml;
#         subenergy.modulus := get_ksub(num0);
#     } else {
#         wl0 := ww / lnum0;
#         wll := ww / lnuml;
#         kl := get_ksub(numl); 
#         k0 := get_ksub(num0);
#         kf := (kl - k0) / rEndX;
#         ls := wl0;
#         lf := (wll - wl0) / rEndX;
#         // E_s = 1/2 * K z^2, K = 1/2 B / (wave_len)^4 * (2 *pi)^4 = B / (l0 + lf * x)^4 * 16 * pi^4/2 -> sube.modulus := B *8 *pi^4
#         subenergy.modulus := bend.modulus * 4 * pi^4; // B=bend.modulus/2
#     };
#     //set_a(lnum0, lnuml, delta);
#     recalc;
# }
# """

        fcpc, pcn = {}, {}
        fcpc.update(rfc)
        fcpc.update(rpc)
        # fcpc.update(dfc)
        # fcpc.update(dpc)
        pcn.update(rpcn)
        for k in pcn.keys() : tail += k + '\n'
        # for k in rpcn.keys(): tail += k + '\n'
        # for k in dpcn.keys(): tail += k + '\n'
        tail += '// procedure without parameters\n'
        for v in pcn.values() : tail += v + '\n'
        # for v in rpcn.values(): tail += v + '\n'
        # for v in dpcn.values(): tail += v + '\n'
        tail += '// function and procedures\n'
        for v in fcpc.values():
            tail += v + '\n'
    else:
        refine_f = 'refinement_nom_uniform.cmd' if not nonuni else 'refine_nuni.cmd'
        tail += 'read "{}"\n'.format(refine_f)

    if stretch:
        tail += """procedure set_tension(real ten) {
    tensionG := ten;
    fLen := sqrt(1+2*tensionG) * lngth;
    foreach vertex vv do {
        vv.y := vv.y *  sqrt(1 - facet[1].poisson_ratio * ten *2) / sqrt(1 - facet[1].poisson_ratio * tensionG*2);
        //vv.z := vv.z * (1+dtt) /(1+delta) * sqrt(dtt/delta);
    };
    foreach vertex vv where on_boundary 1 || on_boundary 2 do vv.p1 := vv.p1 *  sqrt(1 - facet[1].poisson_ratio * ten*2) /sqrt(1 - facet[1].poisson_ratio * tensionG*2);
    
    righttension.modulus := tensionG;
    rEndX := sqrt(1+2*tensionG) * lngth;
    upEndY := wd * (1+ delta) * sqrt(1 - facet[1].poisson_ratio * tensionG * 2);
    ww := upEndY;
    set_ksub(num0,numl);
    recalc;    
    printf "Setted to new tension: %g\\n", tensionG;
}
set_tension(tensionG)

set_t := {
    foreach vertex vv do vv.x := (tensionG * lngth + lngth) / rEndX * vv.x;
    rEndX := tensionG * lngth + lngth;
}

set_smallt := {
    local lratio;
    k0 := get_ksub(num0);
    lratio := (0.5901888630769108 * (numl - num0) * wd * abs(delta)**(1 / 4) / (bend.modulus * k0)**(1 / 8) / num0 / numl / lngth)**4;
    printf "lratio: %g, sqrt ratio: %g\n", 1/lratio, 1/sqrt(lratio);
    set_tension(6.236976698022494 * sqrt(lratio * bend.modulus * k0) / (wd/lngth * (1/num0 - 1/numl))**2);
    recalc;
}

fix_end := {
    rEndX := sqrt(1+2*tensionG) * lngth;
    righttension.modulus := 0;
    fix rEndX;
}

unfix_end := {
    righttension.modulus := tensionG;
    unfix rEndX;
}

"""
    if substrate.casefold() == 'jump'.casefold():
        tail += """
procedure set_wave0(real lnum0, real lrat) {{
    local aa0;
    local aal;
    local scaa;
    
    //ww := wd * (1 + delta);
    scaa := 8/(wd/num0*2);
    aa0 := 1/pi/lnum0 * sqrt(- delta * wd^2 * (1 + delta) ) * lrat;
    aal := 1/pi/numl * sqrt(- delta * wd^2 * (1 + delta) );
    foreach vertex vv do {{
        vv.z := aal * 1/(1+exp(-(vv.x-{0} * pos) * scaa)) * sin(numl * pi * 2 / ww * vv.y) + aa0 * 1/(1+exp((vv.x-{0} * pos) * scaa)) * sin(lnum0 * pi * 2 / ww * vv.y);
    }};
    
    //foreach vertex vv where x <= {0} * pos do {{
    //    vv.z := aa * 1/(1+exp((vv.x-{0} * pos) * scaa)) * sin(num0 * pi * 2 / ww * vv.y);
    //}};
}}
        """.format('rEndX' if fix_bwavlen or stretch or edge_slope else 'lngth')


    tail += """
set_form_factors
set_thickness(thicknezz)
set_ksub(num0,numl)

read "diagnostics.cmd"

U
hessian_normal off
check_increase off
suppress_warning 1825  // suppress the warning about hessian not positively defined
printf "Current thickness: %g\\n",thicknezz/2/wd
//printf "Current tension: %g\\n", tensionG
printf "Aspect ratio: %g\\n", lngth/wd

//cc0
procedure auto(real lnum0, real lnuml, real ltimes, real rtimes) {
    if lnum0 > lnuml then {
        print "Error: num0 should smaller than numl";
    } else {
        set_thd(lnum0, lnuml, ltimes);
        exec "metis_factor";
        cc0;
        r rtimes;
        set_wave;
        print "Saddling: \\n";
        saddle; hgb(0);
    }
}

    """
    ls.append(tail)

def main(nx, ny, oname, substrate='linear', load_cmd=True, **kw):
    if ny < 4:
        print('Current only work on ny >= 4')
        sys.exit(1)
    if ny % 2 != 0:
        print('Error: ny must be even for periodic model')
        sys.exit(1)
    ls.append("PARAMETER origin_name = \"{0}\"\n".format(oname))
    if load_cmd:
        # rfc, rpc, rpcn = deal_cmd('refine_nuni.cmd')
        # dfc, dpc, dpcn = deal_cmd('diagnostics.cmd')
        cmd = deal_cmd('refine_nuni.cmd')  # , deal_cmd('diagnostics.cmd')
    else:
        cmd = None
    create_header(nx,ny, substrate, cmd=cmd, **kw)
    create_VEF_m(nx, ny, **kw)
    creatTail(substrate=substrate, cmd=cmd, **kw)
    s = ''.join(ls)
    with open(oname + '.fe', 'w') as f:
        f.write(s)


def get_aspect(nx,ny):
    if ny == 0:
        return 1e12
    return nx/(ny * math.sqrt(3)/2)
def generate_aspect(aspect, range_nx, periodic=False):
    err = 1e9
    cy = -1
    cx = -1
    for x in range_nx:
        y = int(x/(math.sqrt(3)/2)/aspect)
        apr = abs(get_aspect(x, y) - aspect)
        apr1 = abs(get_aspect(x, y+1) - aspect)
        if periodic:
            if y % 2 == 0:
                if apr < err:
                    cy = y
                    cx = x
                    err = apr
            else:
                if apr1 < err:
                    cy = y + 1
                    cx = x
                    err = apr1
        else:
            if apr > apr1:
                if apr1 < err:
                    cy = y + 1
                    cx = x
                    err = apr1
            else:
                if apr < err:
                    cy = y
                    cx = x
                    err = apr
    print('relative err:', err/aspect, ' ,(nx, ny) = ({}, {})'.format(cx, cy))
    return(cx,cy)



if __name__ == "__main__":
#     import sys
#     usage="""usage: SCRIPT nx ny "filename_without_extension"
# nx: the total number of equilateral-triangle edges in one row
# ny: number of heights
#           """
#     if len(sys.argv) == 1:
#         print(usage)
#         sys.exit(1) 
#     try:
#         nx = int(sys.argv[1])
#         ny = int(sys.argv[2])
#     except (TypeError, IndexError):
#         sys.stderr.write("must has 3 variables\n")
#         print(usage)
#         sys.exit(2)
#     main(nx, ny, str(sys.argv[3]))
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--substrate', '-ss', dest='substrate', type=str, default='linear' #, nargs='+'
                        , help='if include stubstrate in the surface, "uniform" for a constant K_sub, "linear" for a changing K_sub linearly with wave_len = lf * (ls+x)')
    parser.add_argument('--gap', '-g', dest='gap', action='store_true', help='if include a gap constraint, make amplitude not larger amax')
    parser.add_argument('--energygap', '-eg', dest='energygap', action='store_true', help='if include a gap as energy, make amplitude not larger amax')
    parser.add_argument('--edgeegap', '-ee', dest='edgeegap', action='store_true', help='if set the edge only gap, make amplitude not larger amax')
    parser.add_argument('--edgergap', dest='edge_gap', action='store_true', help='if set the edge only gap, (constraint method), make amplitude not larger amax')
    parser.add_argument('--edgeslope', '-es', dest='edge_slope', action='store_true', help='if set edge slope cannot exceed the max theory')
    parser.add_argument('--preventreverse', '-pr', dest='prevent_reverse', action='store_true', help='if prevent the facet totally reverseded')
    parser.add_argument('--fixbwavelen', '-fb', dest='fix_bwavlen', action='store_true', help='if fix the wavelenght of the two free boundaries')
    parser.add_argument('--stretch', dest='stretch', action='store_true', help='add a stretch tension on a boundary')
    parser.add_argument('--uni', dest='uni', action='store_false', help='use uniform mesh')
    parser.add_argument('--severy', dest='stretch_every', action='store_true', help='add a stretch tension on every boundary vertex, works for periodic')
    
    parser.add_argument('--gbend', dest='gbend', action='store_true', help='if consider gaussian curvature')
    parser.add_argument('--preventedgereverse', dest='edge_pr', action='store_true', help='if prevent the facet totally reverseded on the edge only')
    parser.add_argument('--numedgex', '-nx', dest='nx', type=int, help='number of equilateral-triangle edges in x direction')
    parser.add_argument('--numheighty', '-ny', dest='ny', type=int, help='number of equilateral-triangle hegits in y direction')
    parser.add_argument('--aspectratio', '-ar', dest='ar', type=float, help='generate the given aspect ratio, AR = length_x / length_y')
    parser.add_argument('--nxrangeaspect', '-xr', dest='xr', default='200', help='generate the given aspect ratio, pick NX in XR')
    parser.add_argument('--pos', dest='pos', type=float, help='position at which the k jump, deault: 0.5 (middle)', default=0.5)
    parser.add_argument('filename', nargs='?', metavar='filename', default=None, help='"File_Name_Without_ext"')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    pargs = parser.parse_args()
    
    if pargs.ar is not None:
        xr = eval(pargs.xr)
        if not hasattr(xr, '__len__'):
            xr = range(1, xr+1)
        nx, ny = generate_aspect(pargs.ar, xr, periodic=True)
        if (ny & 1):
            nx *= 2
            ny *= 2
        if ny == 2:
            nx *= 2
            ny *= 2
    else:
        if pargs.nx is None or pargs.ny is None:
            raise ValueError('You have to specify (nx, ny) or aspect ratio, see the help!')
        nx = pargs.nx
        ny = pargs.ny
    if pargs.filename is None:
        filename = '{}_{}'.format(nx, ny)
    else:
        filename = pargs.filename
    main(nx, ny, filename, pargs.substrate, gap=pargs.gap, prevent_reverse=pargs.prevent_reverse, 
         energy_gap=pargs.energygap, edge_egap=pargs.edgeegap, edge_pr=pargs.edge_pr, 
         edge_gap=pargs.edge_gap, fix_bwavlen=pargs.fix_bwavlen, edge_slope=pargs.edge_slope,
         stretch=pargs.stretch, stretch_every=pargs.stretch_every, gbend=pargs.gbend, nonuni=pargs.uni, pos=pargs.pos)


            

