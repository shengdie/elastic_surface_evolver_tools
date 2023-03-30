import re
import math
import numpy as np
from scipy.spatial.distance import pdist

def edge_chain(edges, bound_es):
    """Get one boundary chains of edges in order
    
    Arguments:
        edges {mesh.edges, (Num_edges, 2) array} -- edges construct with vid list [(vid0, vid1), ....]
        bound_es {mesh.boudn_es, 1d array} -- boundary edges id 
    
    Returns:
        circ -- the chain
    """
    bound = edges[bound_es]  # [(startvid, endvid), ....]
    circ = []
    added_edge_bid = []
    start = 0  # first bound edge id
    startvid = bound[start][0]  # start vid
    endvid = bound[start][1]  # end vid
    edge_bound_id = start
    added_edge_bid.append(start)
    circ.append(bound_es[start])  # edge id
    while endvid != startvid:
        pos = np.where(bound == endvid)[0]
        pos = pos[pos != edge_bound_id][0]
        circ.append(bound_es[pos])
        ed = bound[pos]
        endvid = ed[ed != endvid][0]
        edge_bound_id = pos
        added_edge_bid.append(pos)
    return circ, added_edge_bid

def get_edge_chains(edges, bound_es):  # get all chain, return all chains, every chian contains the edge ids in that chain in order
    chains = []
    bound_es = bound_es
    while len(bound_es) != 0:
        chain, idx = edge_chain(edges, bound_es)
        bound_es = np.delete(bound_es, idx)
        chains.append(np.array(chain))
    return chains

def get_max_dist(bound_es, edges, verts):  # get max dist of boudanry verts
    vids = np.unique(edges[bound_es])
    D = pdist(verts[vids])
    return np.nanmax(D)

def read_withoutcomment(file, rmspace=True):
    """Read file without comments of c style
    
    Arguments:
        file {str} -- file path
    
    Returns:
        str -- readed text
    """
    single = '//'
    multi_s, multi_e = '/\*', '*/'
    txt = open(file).read()
    # txt = file
    delrange = [0]
    for m in re.finditer(multi_s, txt):
        start = m.start()
        delrange.extend([start, txt.index(multi_e, start) + 2])
    delrange.append(len(txt))
    txt = ''.join(txt[delrange[i]:delrange[i + 1]] for i in range(0, len(delrange), 2))
    delrange = [0]
    for m in re.finditer(single, txt):
        start = m.start()
        delrange.extend([start, txt.index('\n', start)])
    delrange.append(len(txt))
    txt = ''.join(txt[delrange[i]:delrange[i + 1]] for i in range(0, len(delrange), 2))   
    if rmspace: txt = re.sub(r'^\s*$', '', txt, flags=re.MULTILINE)
    return txt


def find_mbracket(txt, s, r='}', l='{'):
    """find match right bracket give current postion of left bracket
    
    Arguments:
        txt {str} -- string target
        s {int} -- position of current left bracket
    
    Keyword Arguments:
        r {str} -- right bracket (default: {'}'})
        l {str} -- left bracket (default: {'{'})
    
    Returns:
        int -- the postion of right bracket which matches current left bracket,  actually + 1
    """
    stack = [s]
    i = s + 1
    e = -1
    while len(stack) > 0:
        c = txt[i]
        if c == r: 
            stack.pop()
            e = i
        elif c == l:
            stack.append(i)
        i += 1
    
    return e + 1


def deal_cmd(file):
    cmd = read_withoutcomment(file)
    match = {'pc': r'procedure .+{', 'func': 'function .+{', 'pcn': r'^.+:=\s*{'}
    # function
    fc_dec = []
    fc_dict = {}
    for m in re.finditer(match['func'], cmd):
        start, end = m.start(), m.end() - 1
        fc_dec.append(cmd[start:end].strip())
        end_p = find_mbracket(cmd, end)
        fc_dict[fc_dec[-1]] = cmd[start:end_p]
        # fc_body.append(cmd[start:end_p])
    # procedure
    pc_dec = []  # declare of procedure with parameters
    pc_dict = {}  # full body
    for m in re.finditer(match['pc'], cmd):
        start, end = m.start(), m.end() - 1
        pc_dec.append(cmd[start:end].strip())
        end_p = find_mbracket(cmd, end)
        pc_dict[pc_dec[-1]] = cmd[start:end_p]
        # pc_body.append(cmd[start:end_p])
    # procedure without parameters
    pcn_dec = []
    pcn_dict = {}
    for m in re.finditer(match['pcn'], cmd, re.MULTILINE):
        start, end = m.start(), m.end()
        pcn_dec.append(cmd[start:end].strip() + '}')
        end_p = find_mbracket(cmd, end - 1)
        pcn_dict[pcn_dec[-1]] = cmd[start:end_p]
        # pcn_body.append(cmd[start:end_p])
    # return (fc_dec, fc_dict), (pc_dec, pc_dict), (pcn_dec, pcn_dict)
    return fc_dict, pc_dict, pcn_dict


def get_proper_length(num0, numl, w, ltimes, pratio=0.4, lratio=0.004):
    # dl = (w / num0 - w / numl) / l
    delta = 0.9 / numl**2
    # delta / (2 * np.sqrt(B K)) == ltimes
    # B = (delta / ltimes / 2)^2 / K
    # h = np.sqrt(3 / 2 * delta * (1 - pratio**2) / ltimes) / numl / np.pi * w
    B = delta * w**2 / 8 / ltimes / numl**2 / math.pi**2
    # Kl = 2 * delta * numl**2 * np.pi**2 / ltimes / w**2
    K0 = (2 * math.pi * num0 / w)**4 * B 
    l = 0.5901888630769108 * (numl - num0) * w * delta**(1 / 4) / (B * K0)**(1 / 8) / num0 / numl / lratio**(1 / 4)
    dl = (w / num0 - w / numl) / l
    T = 6.236976698022494 * math.sqrt(lratio * B * K0) / dl**2
    return l, T 


def create_header(wd, b, lngth, oname):
    ls = []
    ls.append('PARAMETER origin_name = "{}"\nPARAMETER edrgap = 0\nPARAMETER edgepr = 0\n\n'.format(oname))
    ls.append("PARAMETER refine_times = 0\nPARAMETER run_time = 120\nPARAMETER submit_times = 1\nPARAMETER ost = submit_times - 1\n\n")
    #l.append("PARAMETER tstep = 0.0005 // this the step of changing the tension\n\n")

    ls.append("PARAMETER wd =  {}\n".format(wd))
    ls.append("PARAMETER gridsize = {}\n".format(b))
    ls.append("PARAMETER lngth = {}\n".format(lngth))

    header ="""PARAMETER hypertenuse = sqrt(wd**2+(wd)**2/3)
PARAMETER thicknezz = wd/100
PARAMETER oth = thicknezz
PARAMETER ow = wd

PARAMETER delta = -0.05
PARAMETER ww = wd * (1 + delta)
//PARAMETER rightEndX = lngth * (1 + delta)
PARAMETER upEndY = wd * (1 + delta)
PARAMETER kf = 0
PARAMETER ks = 0
PARAMETER minnum = 1
PARAMETER maxnum = 2
//PARAMETER tti = amax * minnum * pi/ww
//PARAMETER tta = amin * maxnum * pi/ww
PARAMETER af = 1
PARAMETER as = 1
PARAMETER amax = 0
PARAMETER amin = 0
PARAMETER fa = gridsize * 0.5

//OPTIMIZING_PARAMETER rightEndX = lngth
//OPTIMIZING_PARAMETER leftEndX = 0

// For resetting form_factors in refining
define vertex attribute ref_coord real[3]
define vertex attribute corner integer
define vertex attribute old_vid  integer
define edge attribute old_eid  integer
define edge attribute divide   integer
define edge attribute rheight  integer

define facet attribute poisson_ratio real
define facet attribute form_factors real[3]
//define facet attribute ksub real
quantity stretch energy method linear_elastic global

quantity bend energy method star_perp_sq_mean_curvature //global
quantity gbend energy method star_gauss_curvature //global

quantity subenergy energy method facet_scalar_integral global\n scalar_integrand: z^2\n

define edge attribute border integer //this is real border

"""
    header += """

//constraint 1 nonpositive
//    formula: x = rightEndX
//constraint 2 nonnegative
//    formula: x = 0


"""
    ls.append(header)
    return ls

def creatTail(edge_egap=False, **kw):
    tail = """
read

set facet tension 0
set facet poisson_ratio 0.4

ost := submit_times - 1

//lefttension.modulus := tensionG
//righttension.modulus := tensionG
//freeedgebend.modulus := 0
read "refinement_nom.cmd"
set_thickness(thicknezz)
set_ksub(minnum,maxnum)

read "diagnostics.cmd"

U
hessian_normal off
check_increase off
suppress_warning 1825  // suppress the warning about hessian not positively defined
printf "Current thickness: %g\\n",thicknezz/2/wd
//printf "Current tension: %g\\n", tensionG
printf "Aspect ratio: %g\\n", lngth/wd

cc0
evolve := {
    set_delta(-0.1);
    set_thickness(lngth/85);
    set_ksub(5,15);
    cc0;
    r 3;
    cc0;
    saddle;
    cc0;
    r;
    cc0; saddle;
    cc0;
    r; //saddle;
    cc0
}
evotest := {//metis_factor; 
    set_delta(-0.05);
    set_thickness(lngth/1000);
    set_ksub(10,30);
    cc0; r 6; cc0; saddle; cc0
}
evotest1 := {
    //metis_factor;
    set_delta(-0.1);
    set_thickness(0.00666666666666667);
    set_ksub(5,15);
    cc0; r 6; cc0; saddle; cc0
}

    """
    return [tail]

def write_fe(ls, file):
    with open(file, 'w') as f:
        f.write(''.join(ls))
