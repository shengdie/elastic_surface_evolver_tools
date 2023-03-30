import numpy as np
import dmsh
from scipy.spatial import Delaunay
from helpers import get_edge_chains, get_max_dist, write_fe, deal_cmd

class Mesh(object):
    def __init__(self):
        self.multibound = False
        self.faces = None
        self.edges = None
    # @property
    # def verts(self):
    #     return self._verts
    # @property
    # def cells(self):
    #     return self._cells

    def generate_delaunay(self):
        raise NotImplementedError

    def generate_VEF(self): # generate Vertices, edges, facets, facet has edge id from 1, which is  idx+1
        edges = {}
        faces = []
        i = 1 # edge id
        for c in self.cells:
            fe=[]
            for c0, c1 in zip(c, np.roll(c, -1)):
                if ((c0, c1) not in edges.keys()): # (c0, c1) not in edges
                    if ((c1, c0) not in edges.keys()): # check (c1, c0) is in edges
                        edges[(c0, c1)] = i # not in, then add
                        fe.append(i) # add one edge of fe
                        i +=1
                    else:
                        fe.append(-edges[(c1, c0)]) # (c1,c0) already in, then add - id
                else: # (c0, c1) already in
                    fe.append(edges[(c0,c1)])
            faces.append(fe)
        self.edges = np.array(list(edges.keys()))
        self.faces = np.array(faces)
        uni, count = np.unique(np.abs(self.faces) - 1, return_counts=True)
        #self.bound_es = self.edges[uni[count == 1]]
        self.bound_es = uni[count == 1] # edge id on boundary
        self.bound_vs = np.unique(self.edges[self.bound_es]) # verts id on boundary
        if self.multibound:
            chains = get_edge_chains(self.edges, self.bound_es)
            if len(chains) == 1:
                print("Error: No multi boundary  chain\n")
                # return 0;
                return self.edges, self.faces
            elif len(chains) >2:
                print("Error: You have more than 2 boounday, cannot handel, now\n")
                return self.edges, self.faces
            d1 = get_max_dist(chains[0], self.edges, self.verts)
            d2 = get_max_dist(chains[1], self.edges, self.verts)
            self.bounds_es = {}
            if d1 > d2: 
                self.bounds_es['inner'] = chains[1]
                self.bounds_es['outer'] = chains[0]
            else:
                self.bounds_es['inner'] = chains[0]
                self.bounds_es['outer'] = chains[1]
            self.bounds_vs = {}
            self.bounds_vs['inner'] = np.unique(self.edges[self.bounds_es['inner']])
            self.bounds_vs['outer'] = np.unique(self.edges[self.bounds_es['outer']])
        return self.edges, self.faces

    def write_VFE(self, non_gloabl_energy=[], boundary=None, constraint=None, fixed=None, 
                 defined_bound=None, defined_bound_fixed=True, fixonebound=False, fixoneinner=False, 
                 fixatcenter=False, twopoint=False, gbend=None):
        """write VFE, 
        
        Arguments:
            energy {list} -- {energy names append on verts not on border}
            boundary {dict} -- {name: [boundary free parameter, boundary equation]}
            constraint {dict} -- {name: constraint_eq_func}
            defined_bound {tuple} -- [(('*num1', '*num2', '*num3'), boundary equation), ...]
            fixed {list} -- [fixed eq]
        """
        if self.faces is None: self.generate_VEF()
        ls = []
        ls.append('\nvertices\n')
        #if boundary is None:
        bisnotNone = (boundary is not None)
        dbisnotNone = (defined_bound is not None)
        cisnotNone = (constraint is not None)
        fisnotNone = (fixed is not None)
        if fixatcenter:
            centeri = np.argmin(np.linalg.norm(self.verts, axis=1))
        pid = []
        if twopoint:
            inr = np.min(np.linalg.norm(self.verts, axis=1))
            inrid = np.argwhere(np.linalg.norm(self.verts, axis=1) / inr - 1 <= 0.01).flatten()
            ps = self.verts[inrid, ]
            a1 = np.argwhere(ps[:, 0] < 0).flatten()
            id1 = np.argmin(np.abs(ps[a1, ][:, 1]))
            
            pid.append(inrid[a1[id1]])
            a2 = np.argwhere(ps[:, 0] > 0).flatten()
            id2 = np.argmin(np.abs(ps[a2, ][:, 1]))
            pid.append(inrid[a2[id2]])

        for i in range(len(self.verts)):
            if bisnotNone:
                inb = False
                for bname, (bp, beq) in boundary.items():
                    if beq(self.verts[i]):
                        l = '{0} {5} boundary {4} ref_coord {{ {1} {2} {3} }}'.format(i+1, self.verts[i,0], self.verts[i,1], self.verts[i,2], bname, self.verts[i, bp])
                        inb = True
                        break
                if not inb:
                    l = '{0} {1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, self.verts[i,0], self.verts[i,1], self.verts[i,2])
            elif dbisnotNone:
                inb = False
                for bp, beq in defined_bound:
                    if beq(self.verts[i]):
                        l = '{0} {1}{4} {2}{5} {3}{6} ref_coord {{ {1} {2} {3} }}'.format(i+1, self.verts[i,0], self.verts[i,1], self.verts[i,2], *bp)
                        if defined_bound_fixed:
                            l += ' fixed'
                        inb = True
                        break
                if not inb:
                    l = '{0} {1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, self.verts[i,0], self.verts[i,1], self.verts[i,2])
            else:
                l = '{0} {1} {2} {3} ref_coord {{ {1} {2} {3} }}'.format(i+1, self.verts[i,0], self.verts[i,1], self.verts[i,2])
            if twopoint:
                if i in pid:
                    l += ' constraint 1'
            if cisnotNone:
                for cname, ceq in constraint.items():
                    if ceq(self.verts[i]):
                        l += ' constraint ' + str(cname)
            if fisnotNone:
                for f in fixed:
                    if f(self.verts[i]): l += ' fixed'
            if not i in self.bound_vs:
                l += ' ' + ' '.join(non_gloabl_energy) 
            if fixoneinner:
                if i in self.bounds_vs['inner']:
                    l += ' fixed'
                    fixoneinner = False
            elif (fixonebound):
                if i in self.bound_vs:
                    l += ' fixed'
                    fixonebound = False
            if fixatcenter and i == centeri:
                l += ' fixed'
                fixatcenter = False
            ls.append(l + '\n')
        ls.append('\nedges\n')
        for i, ed in enumerate(self.edges):
            l = '{0} {1} {2}'.format(i+1, ed[0] + 1, ed[1] + 1)
            if boundary is not None: # add edges to boundary
                for bname, (bp, beq) in boundary.items():
                    if beq(self.verts[ed[0]]) and beq(self.verts[ed[1]]): l+= ' boundary ' + str(bname)
            if constraint is not None: # add edges to constraint
                for cname, ceq in constraint.items():
                    if ceq(self.verts[ed[0]]) and ceq(self.verts[ed[1]]): l+= ' constraint ' + str(cname)
            if i in self.bound_es: # set border to 1 if on border
                l += ' border 1'
                if self.multibound:
                    if i in self.bounds_es['inner']:
                        l += ' inborder 1'
                    else:
                        l += ' outborder 1' 
                else:
                    l += ' outborder 1' 
            
            ls.append(l + '\n')
        ls.append('\nfaces\n')
        for i, fa in enumerate(self.faces):
            ls.append('{} {} {} {}\n'.format(i+1, fa[0], fa[1], fa[2]))
        return ls
                
    # def write_fe(self, verts=None, edges=None, faces=None):
    #     if verts is None: verts = self.verts
    #     if edges is None: edges = self.edges
    #     if faces is None: faces = self.faces
        
class Rectangle(Mesh):
    def __init__(self, l, w, edgesize, dynamicsize=None, dim=3):
        super().__init__()
        self.l = l
        self.w = w
        self.edgesize = edgesize
        self.dynamicsize = dynamicsize
        self.dim = dim

    def generate_delaunay(self, dynamicsize=None, **kw):
        geo = dmsh.Rectangle(0, self.l, 0, self.w)
        esize = eval(dynamicsize, locals()) if dynamicsize is not None else self.edgesize
        X, self.cells = dmsh.generate(geo, esize, **kw)
        if self.dim == 3: 
            self.verts = np.concatenate((X, np.zeros(len(X))[:, None]), axis=1)
        else:
            self.verts = X
        return self.verts, self.cells
    
    def create_header(self, extra=None, cmd=None):
        ls = []
        if cmd is not None:
            # (rfc, rpc, rpcn), (dfc, dpc, dpcn) = cmd
            rfc, rpc, rpcn = cmd
            ls.extend([k + ';\n' for k in rfc.keys()])
            ls.extend([k + ';\n' for k in rpc.keys()])
            # ls.extend([k + ';\n' for k in dfc.keys()])
            # ls.extend([k + ';\n' for k in dpc.keys()])
            ls.append('\n')
        ls.append('PARAMETER submit_times = 1\n')
        ls.append('PARAMETER origin_name = "{}"\n'.format(self.oname))
        ls.append("PARAMETER refine_times = 0\n")
        ls.append('PARAMETER last_stime = 0\nPARAMETER save_period = 10800\nPARAMETER run_time = 120\n\n')
        # l.append("PARAMETER tstep = 0.0005 // this the step of changing the tension\n\n")
        ls.append('PARAMETER in_r = {}\nPARAMETER out_r = {}\n'.format(self.in_r, self.out_r))
        
        ls.append("PARAMETER gridsize = {}\n".format(self.edgesize))

        header = """
PARAMETER thicknezz = in_r/100
PARAMETER oth = thicknezz

PARAMETER delta = -0.05
PARAMETER dratio = 1 + delta

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

quantity bend energy method star_perp_sq_mean_curvature global
//quantity gbend energy method star_gauss_curvature //global

//quantity subenergy energy method facet_scalar_integral global\n //scalar_integrand: z^2\n

define edge attribute border integer //this is real border
define edge attribute inborder integer
define edge attribute outborder integer

"""
        if self.gbend:
            header += "quantity gbend energy method star_gauss_curvature global\n"
        elif self.egbend:
            header += "quantity gbend energy method star_gauss_curvature \n"
        
        if extra is not None:
            header += extra
        ls.append(header)
        if self.in_r == 0:
            ls.append("""PARAMETER c_r = out_r*dratio
constraint 1
formula: x^2 + y^2 = c_r^2
    """)
        elif self.in_free:
            ls.append("""PARAMETER c_r = in_r*dratio
constraint 1
formula: x^2 + y^2 = c_r^2
    """)
        return ls
    
    def creatTail(self, edge_egap=False, cmd=None, **kw):
        tail = """
read

set facet tension 0
set facet poisson_ratio 0.33

ost := submit_times - 1
U
hessian_normal off
check_increase off
suppress_warning 1825  // suppress the warning about hessian not positively defined
"""
        # if self.gbend:
        #     tail += """
        #     foreach vertex vv where old_vid == 0 do {if max(vv.edge, inborder) == 0 && max(vv.edge, outborder) == 0 then set vv gbend}
        #     """
        if self.egbend:
            tail += "\n foreach vertex vv do {if max(vv.edge, inborder) == 0 && max(vv.edge, outborder) == 0 then set vv gbend}\n"
        if cmd is not None:
            rfc, rpc, rpcn = cmd
            if self.in_free:
                for k, v in rpc.items():
                    if 'set_delta' in k:
                        rpc[k] = v.replace("set vv fixed;", "")
                        if self.twopoint:
                            v = rpc[k]
                            v = v.replace("vv.x := vv.ref_coord[1] * dratio;", "")
                            v = v.replace("vv.y := vv.ref_coord[2] * dratio;", "")
                            rpc[k] = v
                    if 'refinemarked' in k:
                        rpc[k] = v.replace("set vv fixed;", "")
                        if self.twopoint:
                            v = rpc[k]
                            v = v.replace("vv.x := vv.ref_coord[1] * dratio;", "")
                            v = v.replace("vv.y := vv.ref_coord[2] * dratio;", "")
                            rpc[k] = v
            
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
            tail += 'set_form_factors'
        else:
            tail += """
read "refinement_uni.cmd"
"""
        tail += """
set_thickness(thicknezz)
set_delta(delta)
read "diagnostics.cmd"
"""
        return [tail]

    def write_input(self, load_cmd=False, twopoint=False):
        cmd = None
        if load_cmd:
            cmd = deal_cmd('refinement_uni.cmd')

        ls = self.create_header(cmd=cmd)
        if self.in_r == 0:
            ls += self.write_VFE(constraint={1: lambda x: abs((x**2).sum() / self.out_r**2 - 1) <= 0.01}, 
                                 fixonebound=True, defined_bound_fixed=False)
        elif twopoint:
            ls += self.write_VFE(fixonebound=True, defined_bound_fixed=False, twopoint=twopoint)
        elif self.in_free: 
            ls += self.write_VFE(constraint={1: lambda x: abs((x**2).sum() / self.in_r**2 - 1) <= 0.01},
                                 fixoneinner=True, defined_bound_fixed=False)
        else:
            ls += self.write_VFE(defined_bound=[(('*dratio', '*dratio', ''), 
                                 lambda x: abs((x**2).sum()/ self.in_r**2 - 1) <= 0.01)])
        ls += self.creatTail(cmd=cmd)
        write_fe(ls, self.oname + '.fe')


class Hexagon(Mesh):
    def __init__(self, elen, snum, dim=3):
        super().__init__()
        self.elen = elen
        self.snum = snum
        self.dim = dim
        self.verts, self.cells = self.generate_delaunay()

    def generate_delaunay(self):
        edge_len = self.elen / self.snum
        self.edge_len = edge_len
        a = edge_len / 2
        b = a * np.sqrt(3)
        self.grid_size = b

        all_p = []
        numr = self.snum * 2
        for j in range(0, self.snum+1):
            for i in range(-numr + j, numr - j + 1, 2):
                all_p.append([i, j])
                if j > 0: all_p.append([i, -j])
        all_p = np.array(all_p, dtype=float)
        all_p[:,0] *= a
        all_p[:,1] *= b
        tri = Delaunay(all_p)
        #print(all_p)
        if self.dim == 3: return np.concatenate((tri.points, np.zeros(len(all_p))[:,None]), axis=1), tri.simplices
        return tri.points, tri.simplices

class Disk(Mesh):
    def __init__(self, radis, edgesize, **kw):
        super().__init__()
        self.radis = radis
        self.edgesize = edgesize
        self.kw = kw

    def generate_delaunay(self):
        geo = dmsh.Circle([0.0, 0.0], self.radis)
        X, cells = dmsh.generate(geo, self.edgesize, **self.kw)
        return X, cells

class Ring(Mesh):
    def __init__(self, in_r, out_r, edgesize, dynamicsize=None, dim=3, oname=None, in_free=False, twopoint=False, gbend=False, egbend=False, **kw):
        """init, create a ring mesh
        
        Arguments:
            in_r {float} -- iner radius
            out_r {float} -- outer radius
            edgesize {float} -- edge size
        """
        super().__init__()
        self.in_r = in_r
        self.out_r = out_r
        self.in_free = in_free
        self.twopoint = twopoint
        self.gbend = gbend
        self.egbend = egbend
        # if edgesize is str:
        #     self.edgesize = eval(edgesize)
        # else:
        self.edgesize = edgesize
        self.dim = dim
        if self.in_r > 0:
            self.multibound = True

        self.generate_delaunay(dynamicsize, **kw)
        if oname is None:
            oname = 'ring{}_{}'.format(in_r, out_r)
            if dynamicsize is not None:
                oname += 'ds'
        self.oname = oname
        
    def generate_delaunay(self, dynamicsize=None, **kw):
        if (self.in_r == 0):
            geo = dmsh.Circle([0, 0.0], self.out_r)
        else:
            ic = dmsh.Circle([0, 0.0], self.in_r)
            geo = dmsh.Difference(dmsh.Circle([0, 0.0], self.out_r), ic)
        esize = eval(dynamicsize, locals()) if dynamicsize is not None else self.edgesize
        X, self.cells = dmsh.generate(geo, esize, **kw)
        if self.dim == 3: 
            self.verts = np.concatenate((X, np.zeros(len(X))[:, None]), axis=1)
        else:
            self.verts = X
        return self.verts, self.cells
    
    def create_header(self, extra=None, cmd=None):
        ls = []
        if cmd is not None:
            # (rfc, rpc, rpcn), (dfc, dpc, dpcn) = cmd
            rfc, rpc, rpcn = cmd
            ls.extend([k + ';\n' for k in rfc.keys()])
            ls.extend([k + ';\n' for k in rpc.keys()])
            # ls.extend([k + ';\n' for k in dfc.keys()])
            # ls.extend([k + ';\n' for k in dpc.keys()])
            ls.append('\n')
        ls.append('PARAMETER submit_times = 1\n')
        ls.append('PARAMETER origin_name = "{}"\n'.format(self.oname))
        ls.append("PARAMETER refine_times = 0\n")
        ls.append('PARAMETER last_stime = 0\nPARAMETER save_period = 10800\nPARAMETER run_time = 120\n\n')
        # l.append("PARAMETER tstep = 0.0005 // this the step of changing the tension\n\n")
        ls.append('PARAMETER in_r = {}\nPARAMETER out_r = {}\n'.format(self.in_r, self.out_r))
        
        ls.append("PARAMETER gridsize = {}\n".format(self.edgesize))

        header = """
PARAMETER thicknezz = in_r/100
PARAMETER oth = thicknezz

PARAMETER delta = -0.05
PARAMETER dratio = 1 + delta

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

quantity bend energy method star_perp_sq_mean_curvature global
//quantity gbend energy method star_gauss_curvature //global

//quantity subenergy energy method facet_scalar_integral global\n //scalar_integrand: z^2\n

define edge attribute border integer //this is real border
define edge attribute inborder integer
define edge attribute outborder integer

"""
        if self.gbend:
            header += "quantity gbend energy method star_gauss_curvature global\n"
        elif self.egbend:
            header += "quantity gbend energy method star_gauss_curvature \n"
        
        if extra is not None:
            header += extra
        ls.append(header)
        if self.in_r == 0:
            ls.append("""PARAMETER c_r = out_r*dratio
constraint 1
formula: x^2 + y^2 = c_r^2
    """)
        elif self.in_free:
            ls.append("""PARAMETER c_r = in_r*dratio
constraint 1
formula: x^2 + y^2 = c_r^2
    """)
        return ls
    
    def creatTail(self, edge_egap=False, cmd=None, **kw):
        tail = """
read

set facet tension 0
set facet poisson_ratio 0.33

ost := submit_times - 1
U
hessian_normal off
check_increase off
suppress_warning 1825  // suppress the warning about hessian not positively defined
"""
        # if self.gbend:
        #     tail += """
        #     foreach vertex vv where old_vid == 0 do {if max(vv.edge, inborder) == 0 && max(vv.edge, outborder) == 0 then set vv gbend}
        #     """
        if self.egbend:
            tail += "\n foreach vertex vv do {if max(vv.edge, inborder) == 0 && max(vv.edge, outborder) == 0 then set vv gbend}\n"
        if cmd is not None:
            rfc, rpc, rpcn = cmd
            if self.in_free:
                for k, v in rpc.items():
                    if 'set_delta' in k:
                        rpc[k] = v.replace("set vv fixed;", "")
                        if self.twopoint:
                            v = rpc[k]
                            v = v.replace("vv.x := vv.ref_coord[1] * dratio;", "")
                            v = v.replace("vv.y := vv.ref_coord[2] * dratio;", "")
                            rpc[k] = v
                    if 'refinemarked' in k:
                        rpc[k] = v.replace("set vv fixed;", "")
                        if self.twopoint:
                            v = rpc[k]
                            v = v.replace("vv.x := vv.ref_coord[1] * dratio;", "")
                            v = v.replace("vv.y := vv.ref_coord[2] * dratio;", "")
                            rpc[k] = v
            
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
            tail += 'set_form_factors'
        else:
            tail += """
read "refinement_uni.cmd"
"""
        tail += """
set_thickness(thicknezz)
set_delta(delta)
read "diagnostics.cmd"
"""
        return [tail]

    def write_input(self, load_cmd=False, twopoint=False):
        cmd = None
        if load_cmd:
            cmd = deal_cmd('refinement_uni.cmd')

        ls = self.create_header(cmd=cmd)
        if self.in_r == 0:
            ls += self.write_VFE(constraint={1: lambda x: abs((x**2).sum() / self.out_r**2 - 1) <= 0.01}, 
                                 fixonebound=True, defined_bound_fixed=False)
        elif twopoint:
            ls += self.write_VFE(fixonebound=True, defined_bound_fixed=False, twopoint=twopoint)
        elif self.in_free: 
            ls += self.write_VFE(constraint={1: lambda x: abs((x**2).sum() / self.in_r**2 - 1) <= 0.01},
                                 fixoneinner=True, defined_bound_fixed=False)
        else:
            ls += self.write_VFE(defined_bound=[(('*dratio', '*dratio', ''), 
                                 lambda x: abs((x**2).sum()/ self.in_r**2 - 1) <= 0.01)])
        ls += self.creatTail(cmd=cmd)
        write_fe(ls, self.oname + '.fe')
        
                

            


