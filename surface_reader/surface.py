#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 01:17:38 2017

@author: xinmm
"""

import os
import sys, warnings
import json
import numpy as np
from functools import cached_property
from matplotlib import transforms

from numpy.linalg import norm, inv, eigh

from scipy.interpolate import griddata, Rbf, interp1d
from scipy import integrate, optimize
from scipy.signal import argrelmax, argrelmin
import matplotlib.pyplot as plt
import plotly.offline as ply
import plotly.graph_objs as go
from other_funcs import contour_plot_ply, contour_plot, plt_curve, plt_figure_mpl, sample_function, plt_func1d, plt_any_mpl, plt_mesh, car2pol, pol2car

warnings.filterwarnings("ignore")

class EvolveSurf(object):
    """Abstracting a surface from Surface Evolver."""
    def __init__(self,filename=None, mirror=False, newdata=True, periodic=False, period=None, polar=False):
        """Args:
                filename: file name to load
                mirror: mirror the surface by x=0 and y=0
        """
        self.periodic = periodic
        self.period = period
        self.polar = polar
        np.set_printoptions(precision=15)
        self._work = None  # work done by tension, positive
        self._fake_gbende = None  # fake gbend
        
        self.plabel = ((((os.path.basename(filename)).split('_'))[0]).split('.'))[0]
        self.mirror = mirror
        self.mverts = None
        self.mrefverts = None
        self.fstrain = None  # strain tensor {{strain_11,strain_12},{strain_21,strain_22}}
        self.fstress = None  # stress, not tensor {stress_11,stress_12,stress_22}
        self.farea = None    # area of every facet
        self.frefarea = None  # ref area of every facet
        self.vstress_y = None  # stress_y value of every vertex
        self.vstress_x = None  # stress_x value of every vertex
        self.vstress_xy = None  # shear stress
        self.vstrainE = None  # vertex strain energy
        self.vbendE = None  # vtexbend energy
        self.btangel = None
        self.cstress_y = None  # stress y value of centroids of facets
        self.fcentroid = None  # centroids of facets
        self.reffcentroid = None  # centroids of ref facets
        self.interp_stress_y = None  # interpolation of stress y
        self.gstressy = None  # interpolated meshgrids of stress y
        self.sverts = None  # sorted verts
        self.tensionfiled = None
        self.ybound = None
        #self.strainE = None # strain energy per facet
        self._min_rx = None  # min max ref coord
        self._max_rx = None
        self._min_ry = None
        self._max_ry = None
        self._min_x = None
        self._min_y = None
        self._max_x = None
        self._max_y = None
        self._w = None  # width
        self._deltaLN = None  # (Lf-L)/L
        self._delta = None
        self._h_n = None  # normalized thickness, h/w
        self._L_n = None  # normalized length by width
        self._e_sheet = None  # sheet energy = total_energy + work
        self._total_strainE = None  # total strain energy
        self._largest_compression = None
        self._maxyf = None  # interp of ybound
        self._bmod = None  # bending modulus
        if filename:
            self.load(filename, newdata)
            if self.polar and hasattr(self, 'in_r'):
                self.rin_r = self.in_r * (1 + self.delta)
                self.tension = 0
                self.vstress_r = None
                self.vstress_t = None  # \sigma_theta
                self.vstress_rt = None  # \sigma_r\theta
                self.verts_rt = None  # coordinates in polar 

    #@jit
    def load(self, filename, newdata):
        type_, out = load_output(filename)
        print(type_)
        if type_ == 'txt':
            d = np.loadtxt(filename)
        
        # the input file should have 8 columns
        # nv: number of verts
        # nf: number of facets
        # pratio: poisson ratio
        # verts: the verties
        # refverts: the reference coords of verts
        # fverts: verties of every facet
        # tension: tension excerted
            self.nv, self.nf = map(int, d[0,:2])
            self.thickness, self.gridsize, self.energy, self.refine_times, self.pratio, self.period = d[0,2:]
            self.tension = self.period
            if newdata:
                self._fake_gbende, self.stretche, self.bende, self.gbende, self._work = d[1,0:5]
                self.sube = self._fake_gbende
                self.min_wavnum, self.max_wavnum = d[1, -2:]
                #self.total_bende = self.bende + self.gbende - self._fake_gbende
                self.verts = d[2:self.nv+2,:3]
                self.bendE = d[2:self.nv+2, 3]
                self.refverts = d[2:self.nv+2,-3:-1] # ref verts only has two coords(x,y), z=0
                self.fverts = (d[-self.nf:,:3] - 1).astype(int)
                self.strainE = d[-self.nf:, 3]
                if self.periodic:
                    self.in_period = d[-self.nf:, -1].astype(bool)
                
            else:
                self.verts = d[1:self.nv+1,:3]
                self.refverts = d[1:self.nv+1,-3:-1] # ref verts only has two coords(x,y), z=0
                self.fverts = (d[-self.nf:,:3] - 1).astype(int)
                self.strainE = d[-self.nf:, 3]
        elif type_ == 'dict':
            if out.find("'") >= 0:
                out = out.replace("'", '"')
            out = out.replace('], ]', ']]')
            d = json.loads(out)
            for k, v in d.items():
                if k not in ['vlist', 'flist']:
                    self.__dict__[k] = v
            vlist = np.array(d['vlist'])
            flist = np.array(d['flist'])
            self.verts = vlist[:, :3]
            self.bendE = vlist[:, 3]
            self.refverts = vlist[:, -3:-1]
            self.fverts = (flist[:, :3] - 1).astype(int)
            self.strainE = flist[:, 3]
            if self.periodic:
                self.in_period = flist[:, -1].astype(bool)
        if self.periodic:
            y0_idx = np.argwhere(self.refverts[:, 1] == 0).flatten()
            #print(self.refverts[:,1].shape)
            #if self.period is None:
            if self._w is None:
                self._w = self.max_ry - self.min_ry + self.gridsize
            self.refverts = np.concatenate((self.refverts, np.array([self.refverts[y0_idx, 0].flatten(), np.repeat(self._w, len(y0_idx))]).T))
            self.nv = len(self.refverts)
            self.verts = np.concatenate((self.verts, np.array([self.verts[y0_idx, 0].flatten(), self.verts[y0_idx, 1].flatten() + self.period,
                                                               self.verts[y0_idx, 2].flatten()]).T))
            self.bendE = np.concatenate((self.bendE, self.bendE[y0_idx]))                                                    
            # print(y0_idx)
            dic_y0_yb = {y0_idx[i]: (self.nv - len(y0_idx) + i) for i in range(len(y0_idx))}
            for i in np.argwhere(self.in_period).flatten():
                for j in range(3):
                    if self.fverts[i][j] in y0_idx:
                        self.fverts[i][j] = dic_y0_yb[self.fverts[i][j]]
    
    
    @cached_property
    def all_energy(self):
        return {'stretche': self.stretche, 'bende': self.bende, 'gbende': self.gbende}
    
    @property
    def min_rx(self):
        if self._min_rx is None:
            self._min_rx = np.min(self.refverts[:,0])
        return self._min_rx

    @property
    def min_ry(self):
        if self._min_ry is None:
            self._min_ry = np.min(self.refverts[:,1])
        return self._min_ry

    @property
    def max_rx(self):
        if self._max_rx is None:
            self._max_rx = np.max(self.refverts[:,0])
        return self._max_rx

    @property
    def max_ry(self):
        if self._max_ry is None:
            self._max_ry = np.max(self.refverts[:,1])
        return self._max_ry

    @property
    def min_x(self):
        if self._min_x is None:
            self._min_x = np.min(self.verts[:,0])
        return self._min_x

    @property
    def max_x(self):
        if self._max_x is None:
            self._max_x = np.max(self.verts[:,0])
            if self._max_x > 10* self.max_rx: self._max_x = np.sort(self.verts[:,0])[-2]
        return self._max_x

    @property
    def min_y(self):
        if self._min_y is None:
            self._min_y = np.min(self.verts[:,1])
        return self._min_y

    @property
    def max_y(self):
        if self._max_y is None:
            self._max_y = np.max(self.verts[:,1])
        return self._max_y

    # width w,
    @property
    def w(self):
        if self._w is None:
            if self.mirror:
                self._w = 2 * self.max_ry
            elif self.polar:
                self._w = self.in_r
            else:
                self._w = self.max_ry - self.min_ry
        return self._w
    
    @cached_property
    def l(self):
        return self.max_rx - self.min_rx

    @property
    def delta(self):
        if self._delta is None:
            self._delta = (self.w - self.period) / self.w if self.periodic else (self.w - self.max_y) / self.w
        return self._delta
        
    # normalized delta L, (L-L0)/W
    @property
    def deltaLN(self):
        if self._deltaLN is None:
            if self.mirror:
                self._deltaLN = (self.max_x  - self.max_rx)/self.max_rx
            else:
                self._deltaLN = (self.max_x - self.min_x)/(self.max_rx - self.min_rx) - 1
        return self._deltaLN

    # normalized thcikness h, h/width
    @property
    def h_n(self):
        if self._h_n is None:
            self._h_n = self.thickness / self.w
            # if self.mirror:
            #     self._h_n = self.thickness / 2 / self.max_ry
            # else:
            #     self._h_n = self.thickness / (self.max_ry-self.min_ry)
        return self._h_n

    # normalized L, length/width
    @property
    def L_n(self):
        if self._L_n is None:
            if self.mirror:
                self._L_n = self.max_rx/self.w
            else:
                self._L_n = (self.max_rx - self.min_rx)/self.w
        return self._L_n

    @property
    def work(self):
        if self._work is None:
            self._work = self.tension * (self.max_ry - self.min_ry) * ((self.max_x - self.min_x) - (self.max_rx - self.min_rx))
        return self._work

    @property
    def e_sheet(self):
        if self._e_sheet is None:
            self._e_sheet = self.energy + self.work
        return self._e_sheet

    @property
    def fake_gbende(self):
        if self._fake_gbende is None:
            if not self.mirror:
                return 0
            numx0 = 0
            numy0 = 0
            for i in range(self.nv):
                if self.refverts[i,0] == self.min_rx:
                    numx0 += 1
                if self.refverts[i,1] == self.min_ry:
                    numy0 += 1
            gmod = - self.thickness**2/12/(1 + self.pratio)
            self._fake_gbende = ((numx0+numy0-3) * np.pi + np.pi/2)*gmod
        return self._fake_gbende
    
    @property
    def total_strainE(self):
        if self._total_strainE is None:
            self._total_strainE = np.sum(self.strainE)
        return self._total_strainE

    @property
    def total_bende(self):
        return self.bende + self.gbende - self.fake_gbende
    
    @property
    def total_e(self):
        return self.total_bende + self.stretche

    @property
    def bmod(self):
        if self._bmod is None:
            self._bmod = self.thickness**2/ 12 /(1 - self.pratio**2)
        return self._bmod
    
    @property
    def epsilon(self):
        return self.bmod / self.tension / self.w**2

    def kmod(self, n):
        return (2 * np.pi * n / self.period)**4 * self.bmod
        
    # don't work, don't know the reason
    #def get_strainE(self):
    #    if self.strainE is None:
    #        self.get_strain_area()
    #        self.strainE = 0.5 /(1+self.pratio) * ( np.trace(self.fstrain @ self.fstrain, axis1=1, axis2=2) 
    #                        + self.pratio/(1-self.pratio) * np.trace(self.fstrain, axis1=1, axis2=2)**2 )
        #return self.strainE
    
    def get_ybound(self):  # works for rectangular sheet only, works for mirror only
        """Get the vertexes on the real boundary (sorted)
        """
        if self.ybound is None:
            self.ybound = (self.refverts[:,1]==self.max_ry).nonzero()[0]
            #yb = np.asarray([i for i in range(self.nv) if self.refverts[i,1] == self.max_ry])
            #for i in range(self.nv):
            #    if self.refverts[i,1] == self.max_ry:
            #        yb.append(i)
            #so = np.argsort((self.refverts[yb])[:,0])
            so = np.argsort(self.refverts[self.ybound, 0])
            #self.ybound = np.array(yb)
            self.ybound = self.ybound[so]
        #return self._ybound

    #@jit
    def maxy(self, x, interp=False, kind='cubic'):
        """Get the maxy at x = x

        Arguments:
            x {float} -- the position to get max y
        """
        self.get_ybound()
        if interp:
            if self._maxyf is None:
                vs = self.verts[self.ybound]
                x = vs[:, 0]
                y = vs[:, 1]
                self._maxyf = interp1d(x, y, kind=kind)
            return self._maxyf(x)

        if x ==  self.verts[self.ybound[0],0]:
            return self.verts[self.ybound[0],1]
        for i in range(1, len(self.ybound)):
            if x <= self.verts[self.ybound[i],0]:
                return (x - self.verts[self.ybound[i-1],0])*(self.verts[self.ybound[i],1]-self.verts[self.ybound[i-1],1])/(self.verts[self.ybound[i],0]-self.verts[self.ybound[i-1],0]) + self.verts[self.ybound[i-1],1]

    #vc = np.zeros(self.nv).astype(int)
    #for i in range(self.nf):
    #    for j in range(3):
    #        vc[self.fverts[i,j]] += 1
    #for i in range(self.nv):
    #    if vc[i] == 3 and self.verts[i]

    #@jit
    def get_mirror_coords(self):
        """get the mirror coords of the surface, remenber that we cannot duplicate the verts at the boundary (x=0 or y=0)
        """
        if self.mirror and (self.mverts is None):
            # self.mverts = np.concatenate((self.verts, \
            #                              [[self.verts[x,0],-self.verts[x,1],self.verts[x,2]] for x in range(len(self.verts)) if self.verts[x,1] != 0],\
            #                              [[-self.verts[x,0],self.verts[x,1],self.verts[x,2]] for x in range(len(self.verts)) if self.verts[x,0] !=0],\
            #                              [[-self.verts[x,0],-self.verts[x,1],self.verts[x,2]] for x in range(len(self.verts)) if self.verts[x,0] != 0 and self.verts[x,1] != 0]\
            #                             ))
            a = self.verts[:,1] != 0
            b = self.verts[:, 0] != 0
            c = np.logical_and(a, b)
            self.mverts = np.concatenate((self.verts, self.verts[a] * np.array([[1, -1, 1]]),
                                          self.verts[b] * np.array([[-1, 1, 1]]), self.verts[c]  * np.array([[-1, -1, 1]])))
        #return self.mverts
    def get_mirror_refcoords(self):
        if self.mirror and (self.mrefverts is None):
            # self.mrefverts = np.concatenate((self.refverts, \
            #                              [[self.refverts[x,0],-self.refverts[x,1]] for x in range(len(self.refverts)) if self.refverts[x,1] != 0],\
            #                              [[-self.refverts[x,0],self.refverts[x,1]] for x in range(len(self.refverts)) if self.refverts[x,0] !=0],\
            #                              [[-self.refverts[x,0],-self.refverts[x,1]] for x in range(len(self.refverts)) if self.refverts[x,0] != 0 and self.refverts[x,1] != 0]\
            #                             ))
            a = self.refverts[:,1] != 0
            b = self.refverts[:, 0] != 0
            c = np.logical_and(a, b)
            self.mrefverts = np.concatenate((self.refverts, self.refverts[a] * np.array([[1, -1]]), 
                                        self.refverts[b] * np.array([[-1, 1]]), -self.refverts[c]))


    def get_strain_area(self):  # remember strain and stress are both defined in ref coords
        if self.fstrain is None:
            self.fstrain = np.zeros((self.nf, 2, 2))
            self.frefarea = np.zeros(self.nf)
            self.farea = np.zeros(self.nf)

            # e1 = self.verts[self.fverts[:,1]]-self.verts[self.fverts[:,0]]
            # e2 = self.verts[self.fverts[:,2]]-self.verts[self.fverts[:,0]]
            # e12 = (e1 * e2).sum(axis=1)
            # F = 
            for i in range(self.nf):
                e1 = self.verts[self.fverts[i, 1]] - self.verts[self.fverts[i, 0]]  # edge 1 of facet i, which is vector
                e2 = self.verts[self.fverts[i, 2]] - self.verts[self.fverts[i, 0]]  # edge 2 of facet i, which is vector
                e12 = np.dot(e1,e2)
                F = np.matrix([[np.dot(e1,e1),e12],[e12,np.dot(e2,e2)]]) # F matrix
                re1 = self.refverts[self.fverts[i,1]] - self.refverts[self.fverts[i,0]] # edge 1 before strain
                re2 = self.refverts[self.fverts[i,2]] - self.refverts[self.fverts[i,0]] # edge 2 before strain
                self.frefarea[i] = norm(np.cross(re1,re2))*0.5
                # self.farea[i] = norm(np.cross(e1, e2)) * 0.5
                try:
                    QT = np.matrix([re1,re2]).I # transpose of Q matrix
                except:
                    print(np.matrix([re1,re2]))
                self.fstrain[i] = (QT * F * QT.T - np.identity(2))/2
                #ST= np.matrix([re1, re2]) # unstrained edges , transpose of S
                #QT = ST.I # inverse of ST
                #WT = np.matrix([e1, e2]) # strained edges, transpose of W
                #self.fstrain.append((QT * F * QT.T - np.identity(2))/2)
    
    def get_rt_coordinates(self):
        if self.verts_rt is None:
            self.verts_rt = np.concatenate((np.hypot(self.verts[:, 0], self.verts[:,1])[:, None], np.angle(self.verts[:,0] + self.verts[:,1]*1j)[:,None]), axis=-1)

    def get_stress(self): # on ref coords, not work for large strain or take it as PK2 stress
        if self.fstress is None:
            self.get_strain_area()
            fs = []  # remember that fstress is (stress_11, stress_12, stress_22)
            for i in range(self.nf):
                fa1 = 1/(1-self.pratio**2)
                s11 = fa1*self.fstrain[i][0,0]
                s22 = fa1*self.fstrain[i][1,1]
                fs.append([s11 + self.pratio*s22,\
                            1/(1+self.pratio)*self.fstrain[i][0,1],\
                            s22 + self.pratio*s11])
            self.fstress = np.array(fs)

    def get_vstrainE(self, density=True):
        if self.vstrainE is None:
            self.get_stress()
            t1 = np.zeros(self.nv)
            t2 = np.zeros(self.nv)
            if density:
                self.get_strain_area()
                t3 = np.zeros(self.nv)
                for i in range(self.nf):
                    for j in range(3):
                        #print(self.fverts[i,j])
                        t1[self.fverts[i,j]] += self.strainE[i]
                        #t2[self.fverts[i,j]] += 1
                        t3[self.fverts[i,j]] += self.frefarea[i]
                self.vstrainE = np.divide(t1, t3)
            else:
                for i in range(self.nf):
                    for j in range(3):
                        t1[self.fverts[i,j]] += self.strainE[i]
                        t2[self.fverts[i,j]] += 1
                self.vstrainE = np.divide(t1, t2)

    def get_vbendE(self, density=True):
        if self.vbendE is None:
            if not density:
                self.vbendE = self.bendE
                return
            # t1 = np.zeros(self.nv)
            self.get_strain_area()
            t2 = np.zeros(self.nv)
            
            for i in range(self.nf):
                for j in range(3):
                    #t2[self.fverts[i,j]] += 1
                    t2[self.fverts[i,j]] += self.frefarea[i]
            self.vbendE = self.bendE / t2 / 3
        
    #@jit
    def get_vstress_y(self): # stress y, average of surfaces value
        if self.vstress_y is None:
            self.get_stress()
            t1 = np.zeros(self.nv)
            t2 = np.zeros(self.nv)
            for i in range(self.nf):
                for j in range(3):
                    #print(self.fverts[i,j])
                    t1[self.fverts[i,j]] += self.fstress[i][2]
                    t2[self.fverts[i,j]] += 1
            self.vstress_y = np.divide(t1,t2)
    #@jit
    def get_vstress_x(self): # stress x, average of surfaces value
        if self.vstress_x is None:
            self.get_stress()
            t1 = np.zeros(self.nv)
            t2 = np.zeros(self.nv)
            for i in range(self.nf):
                for j in range(3):
                    #print(self.fverts[i,j])
                    t1[self.fverts[i,j]] += self.fstress[i][0]
                    t2[self.fverts[i,j]] += 1
            self.vstress_x = np.divide(t1,t2)
    
    def get_vstress_xy(self): # shear stress
        if self.vstress_xy is None:
            self.get_stress()
            t1 = np.zeros(self.nv)
            t2 = np.zeros(self.nv)
            for i in range(self.nf):
                for j in range(3):
                    #print(self.fverts[i,j])
                    t1[self.fverts[i,j]] += self.fstress[i][1]
                    t2[self.fverts[i,j]] += 1
            self.vstress_xy = np.divide(t1,t2)
    
    def get_vstress_r(self):  # \singma_rr in polar coordinates
        if self.vstress_r is None:
            self.get_vstress_y()
            self.get_vstress_x()
            self.get_vstress_xy()
            self.get_rt_coordinates()
            cos = np.cos(self.verts_rt[:, 1])**2
            self.vstress_r = self.vstress_x * cos + self.vstress_y * (1 - cos) + np.sin(2 * self.verts_rt[:, 1]) * self.vstress_xy
    
    def get_vstress_t(self):
        if self.vstress_t is None:
            self.get_vstress_y()
            self.get_vstress_x()
            self.get_vstress_xy()
            self.get_rt_coordinates()
            cos = np.cos(self.verts_rt[:, 1])**2
            self.vstress_t = self.vstress_x * (1 - cos) + self.vstress_y * cos - np.sin(2 * self.verts_rt[:, 1]) * self.vstress_xy
    
    def get_vstress_rt(self):
        if self.vstress_rt is None:
            self.get_vstress_y()
            self.get_vstress_x()
            self.get_vstress_xy()
            self.get_rt_coordinates()
            cos = np.cos(self.verts_rt[:, 1])
            sin = np.sin(self.verts_rt[:, 1])
            sincos = sin * cos
            c2s2 = cos**2 - sin**2
            self.vstress_rt = -self.vstress_x * sincos + self.vstress_y * sincos + c2s2 * self.vstress_xy

    def get_vbound_angle(self):
        if self.btangel is None:
            self.get_vstress_x()
            self.get_vstress_xy()
            self.get_vstress_y()
            self.get_ybound()
            vboundsy = self.vstress_xy[self.ybound]
            vboundsx = self.vstress_x[self.ybound]
            self.btangel = (vboundsy + np.nan_to_num(np.sqrt(vboundsy**2 - vboundsx * self.vstress_y[self.ybound]))) / vboundsx

    #def get_bound_x
    #@jit
    def get_crossdata_rbf(self, data, c, ci, cf, num=500, a=0, k='cubic', wd=5, polar=False, retfunc=False, **kw): # data should be shape of Nx3
        """Get cross data of surface using interpolation of Rbf method.

        Get cross data of surface defined by 'data', return like (x[], z[])
        returns a=0: x or r=c, a=1: y or theta = c , only origin at 0 for polar 

        Arguments:
            data {array like with Nx3 shape or (base[Nx2],value[N])} -- contains the list of verts which define the surface.
            c {float} -- cross position, like x=1.2 or y=5
            ci {float} -- start position of the cross line
            cf {float} -- final position of the cross line
            **kw {args} -- extra kwargs for Rbf method

        Keyword Arguments:
            num {int} -- how many points you want get in the cross line (default: {500})
            a {int} -- (0,1), axis to cut, 0 is x or r, 1 is y or theta(default: {0})
            wd {int} -- area width for interpolation (wd * gridsize) (default: {5})
            polar {bool} -- if the coordinates are polar
        """
        w = wd * self.gridsize
        if not polar:
            if len(data) == 2:
                d = np.array([[data[0][i][0],data[0][i][1],data[1][i]] for i in range(len(data[1])) \
                                                                        if np.abs(data[0][i][a]-c) <= w])
            else:
                #print('True')
                d = np.array([data[i] for i in range(len(data)) if np.abs(data[i][a]-c) <= w])
            cd = np.linspace(ci, cf, num)
            rbf = Rbf(d[:,0],d[:,1],d[:,2], function=k, **kw)
            if retfunc:
                return rbf
        
            if a==0:
                return (cd, rbf([c]*num, cd))
            elif a == 1:
                return (cd, rbf(cd, [c]*num))
        else:
            if len(data) == 2:
                data = np.concatenate((data[0], data[1][:,None]), axis=1)
            if a == 0:  # return r=c
                d = data[np.abs(np.hypot(data[:,0], data[:,1]) - c) <= w]
                # print(d)
                ang = np.linspace(ci, cf, num)[:, None]
                cd = np.concatenate((c* np.cos(ang), c * np.sin(ang)), axis=1)
                rbf = Rbf(d[:,0],d[:,1],d[:,2], function=k, **kw)
                if retfunc:
                    return rbf
                
                return (ang[:, 0], rbf(cd[:,0], cd[:, 1]))
                #cd = np.concatenate((np.repeat(c, num), np.linspace(ci, cf, num)) 
            elif a==1:
                p1 = np.array((0,0))[None, :]
                p2 = np.array((1, np.tan(c)))[None, :]
                dist = np.abs(np.cross(p2, data[:, :2])) / norm(p2)
                d = data[dist <= w]
                po = np.linspace(ci, cf, num)[:, None]
                cd = np.concatenate((po * np.cos(c), po * np.sin(c)), axis=1)
        
                rbf = Rbf(d[:,0],d[:,1],d[:,2], function=k, **kw)
                if retfunc:
                    return rbf
                return (po[:, 0], rbf(cd[:,0], cd[:, 1]))
            
            

    def get_grid_stress_y(self, numx =500, numy=500, k='cubic', **kw):  #using griddata method to interpolate
        if self.gstressy is None:
            self.get_vstress_y()
            #so = np.lexsort((self.refverts[:,1],self.refverts[:,0]))
            X = np.linspace(self.min_rx, self.max_rx, numx)
            Y = np.linspace(self.min_ry, self.max_ry, numy)
            grid_x, grid_y = np.meshgrid(X,Y)
            self.gstressy = griddata(self.refverts, self.vstress_y, (grid_x, grid_y), method=k, **kw)
            return (X,Y)
    #@jit
    def get_xcross_stressy(self, x, ynum=500, m='rbf', k='cubic', yi=None, yf=None, **kw): # get a x cross view of stress y, ynum is the num of values
        if yi is None: yi = self.min_ry
        if yf is None: yf = self.max_ry
        if x >= self.min_rx and x <= self.max_rx:
            self.get_vstress_y()
            if m == 'griddata':
                Y = np.linspace(self.min_ry, self.max_ry, ynum)
                gx, gy = np.meshgrid(x,Y)
                gz = griddata(self.refverts, self.vstress_y, (gx,gy), method=k, **kw)
                return (gy.flatten(), gz.flatten())
            else:
                return self.get_crossdata_rbf((self.refverts,self.vstress_y), x, yi, \
                                                yf, num=ynum, k=k, **kw)
    
    def get_ycross_stress(self, xy, y, xnum=500, k='cubic', xi=None, xf=None, **kw):
        if xi is None:
            xi=self.min_rx
        if xf is None:
            xf=self.max_rx

        if xy == 'x':
            self.get_vstress_x()
            return self.get_crossdata_rbf((self.refverts, self.vstress_x), y, xi, xf,
                                            num=xnum, k=k, a=1, **kw)
        elif xy == 'y':
            self.get_vstress_y()
            return self.get_crossdata_rbf((self.refverts, self.vstress_y), y, xi, xf,
                                            num=xnum, k=k, a=1, **kw)
    
    def get_xcross_stressxy(self, x, yi=None, yf=None, ynum=500, k='cubic', **kw):
        if yi is None: yi = self.min_ry
        if yf is None: yf = self.max_ry
        if x >= self.min_rx and x <= self.max_rx:
            self.get_vstress_xy()
            return self.get_crossdata_rbf((self.refverts,self.vstress_xy), x, yi, \
                                                yf, num=ynum, k=k, **kw)

    def get_xcross_stressx(self, x, ynum=500, k='cubic', **kw): # get a x cross view of stress y, ynum is the num of values
        if x >= self.min_rx and x <= self.max_rx:
            self.get_vstress_x()
            return self.get_crossdata_rbf((self.refverts,self.vstress_x), x, self.min_ry, \
                                            self.max_ry, num=ynum, k=k, **kw)

    def get_fcentroid(self, fid=-1):  # fid <0: get all centroids
        if fid<0:
            if self.fcentroid is None:
                self.fcentroid = np.array([(self.verts[self.fverts[i,0]] + self.verts[self.fverts[i,1]]\
                                   + self.verts[self.fverts[i,2]])/3 for i in range(self.nf)])
        else:
            fid = int(fid)
            return (self.verts[self.fverts[fid,0]] + self.verts[self.fverts[fid,1]]\
                                   + self.verts[self.fverts[fid,2]])/3

    def get_reffcentroid(self, fid=-1):
        if fid<0:
            if self.reffcentroid is None:
                self.reffcentroid = np.array([(self.refverts[self.fverts[i,0]] + self.refverts[self.fverts[i,1]]\
                                   + self.refverts[self.fverts[i,2]])/3 for i in range(self.nf)])
        else:
            fid = int(fid)
            return (self.refverts[self.fverts[fid,0]] + self.refverts[self.fverts[fid,1]]\
                                   + self.refverts[self.fverts[fid,2]])/3

    def get_cstress_y(self): # stress y, centorids of facets
        """Get the the stress y on every centroid of ref facet
        """
        if self.cstress_y is None:
            self.get_reffcentroid()
            self.get_stress()
            self.cstress_y = []
            #print(self.cstress_y)
            #print(self.reffcentroid[0])
            #self.cstress_y = np.zeros((self.nf,3))
            minrx=self.min_rx
            minry=self.min_ry
            maxrx= self.max_rx
            maxry = self.max_ry
            bverts = {} # boundary verts
            for i in range(self.nf):
                self.cstress_y.append([self.reffcentroid[i,0], self.reffcentroid[i,1], self.fstress[i,2]])
                for j in range(3):
                    if self.refverts[self.fverts[i,j],0] in (minrx, minry) or self.refverts[self.fverts[i,j],1] in (minry,maxry):
                        try:
                            bverts[self.fverts[i,j]][0] += 1
                            bverts[self.fverts[i,j]][1] += self.fstress[i,2]
                        except KeyError:
                            bverts[self.fverts[i,j]] = [0,0]
            for k,v in bverts.items():
                self.cstress_y.append([self.refverts[k,0],self.refverts[k,1], v[1]/v[0]])
            self.cstress_y = np.array(self.cstress_y)

    def get_xcross_stressy_centroid(self, x, ynum=500, m='', k='cubic', **kw):
        self.get_cstress_y()
        if m == 'griddata':
            Y = np.linspace(self.min_ry, self.max_ry, ynum)
            gx, gy = np.meshgrid(x,Y)
            gz = griddata(self.cstress_y[:,:2], self.cstress_y[:,2], (gx,gy), method=k, **kw)
            return (gy.flatten(), gz.flatten())
        else:
            return self.get_crossdata_rbf(self.cstress_y, x, self.min_ry, \
                                                self.max_ry, num=ynum, k=k, **kw)

    @property
    def largest_compression(self):
        if self._largest_compression is None:
            self.get_stress()
            ind = np.argmin(self.fstress[:,2])
            self._largest_compression = (self.fstress[ind,2],ind, self.get_reffcentroid(ind))
        return self._largest_compression

    def plt_contour_stress(self, xy='y', numx=500, numy=500, k='cubic', usempl=False, ret=False, **kw):
        X = np.linspace(self.min_rx, self.max_rx, numx)
        Y = np.linspace(self.min_ry, self.max_ry, numy)
        grid_x, grid_y = np.meshgrid(X,Y)
        if xy == 'y':
            self.get_vstress_y()
            gstress = griddata(self.refverts, self.vstress_y, (grid_x, grid_y), method=k, **kw)
        elif xy == 'x':
            self.get_vstress_x()
            gstress = griddata(self.refverts, self.vstress_x, (grid_x, grid_y), method=k, **kw)
        elif xy == 'xy':
            self.get_vstress_xy()
            gstress = griddata(self.refverts, self.vstress_xy, (grid_x, grid_y), method=k, **kw)
        else:
            print('must be x, y, or xy')
            exit(1)
        if not usempl:
            contour_plot_ply(X,Y,gstress, title=r'$\text{{Contour plot of }}\sigma_{3},\ \tilde{{h}}={0:.3g},\ \tilde{{T}}={1},\ \tilde{{L}}={2:.4f} $'.format(self.h_n, self.tension, np.around(self.L_n,4), xy),
                            xlabel='$x$', ylabel='$y$')
        else:
            contour_plot(X,Y, gstress, title=r'Contour plot of $\sigma_{3},\ \tilde{{h}}={0:.3g},\ \tilde{{T}}={1},\ \tilde{{L}}={2:.4f} $'.format(self.h_n, self.tension, np.around(self.L_n,4), xy),
                            xlabel='$x$', ylabel='$y$')
        if ret:
            return (X,Y,gstress)

    def plt_xcross_stressy(self, x, ynum=500, m='rbf', k='cubic', vs=True, elabel='', **kw):
        if vs:
            gy, gz = self.get_xcross_stressy(x, ynum, m, k, **kw)
        else:
            gy, gz = self.get_xcross_stressy_centroid(x, ynum, m, k, **kw)
        plt.plot(gy, gz)
        plt.ylabel(r'$\sigma_y$')
        plt.xlabel('y')
        plt.title(r'$\sigma_y(x={0},y),\ h={1},\ T={2}$, {3}{4}'.format(x,self.thickness, self.tension, self.plabel, elabel))
        plt.show()

    def plt_xcross_stressy_ply(self, x, ynum=500, m='rbf', k='cubic', vs=True, elabel='', **kw):
        """Plot x cross line of stress_y

        Arguments:
            x {float} -- the cross line position
            **kw {extra kargs} -- extra arg for rbf or griddata

        Keyword Arguments:
            ynum {int} -- how many sample points in the cross line (default: {500})
            m {string} -- 'griddata' using griddata to intepolate, otherwise using rbf
            k {string} -- method for griddata or rbf
            vs {bool} -- True using verts stress, False using centroids stress (default: {True})
            elabel {string} -- extra string in figure title
        """

        if vs:
            gy, gz = self.get_xcross_stressy(x, ynum, m, k, **kw)
        else:
            gy, gz = self.get_xcross_stressy_centroid(x, ynum, m, k, **kw)
        t = go.Scatter(
            x=gy,
            y=gz
        )
        layout= dict(
            title=r'$\sigma_y(x={0},y),\ \tilde{{h}}={1:.3g},\ \tilde{{T}}={2},\ \tilde{{L}}={3:.4f},\ \text{{{4}}} {5}$'.format(x, self.h_n, self.tension, np.around(self.max_rx/self.max_ry,4), self.plabel, elabel),
            xaxis=dict(
                title='y'
            ),
            yaxis=dict(
                title='$\\sigma_y$'
            )
        )
        ply.init_notebook_mode(connected=False)
        f=go.Figure(data=[t],layout=layout)
        ply.iplot(f)
    
    def plt_xcross_stressxy_ply(self, x, ynum=500, m='rbf', k='cubic', elabel='', **kw):
        """Plot x cross line of stress_y

        Arguments:
            x {float} -- the cross line position
            **kw {extra kargs} -- extra arg for rbf or griddata

        Keyword Arguments:
            ynum {int} -- how many sample points in the cross line (default: {500})
            m {string} -- 'griddata' using griddata to intepolate, otherwise using rbf
            k {string} -- method for griddata or rbf
            vs {bool} -- True using verts stress, False using centroids stress (default: {True})
            elabel {string} -- extra string in figure title
        """

        gy, gz = self.get_xcross_stressxy(x, ynum, k, **kw)
        #else:
        #    gy, gz = self.get_xcross_stressxy_centroid(x, ynum, m, k, **kw)
        t = go.Scatter(
            x=gy,
            y=gz
        )
        title=r'$\sigma_y(x={0},y),\ \tilde{{h}}={1:.3g},\ \tilde{{T}}={2},\ \tilde{{L}}={3:.4f},\ \text{{{4}}} {5}$'.format(x, \
            self.h_n, self.tension, np.around(self.max_rx/self.max_ry,4), self.plabel, elabel)
        plt_curve(gy, gz, title=title, xaxis='y', yaxis='$\\sigma_{xy}$')

    def get_itg_data(self, x, ynum=513, data='straine', k='cubic', **kw):
        zdata = np.zeros(self.verts[:,2].size) if not self.mirror else np.zeros(self.mverts[:,2].size)
        data = data.casefold().split()
        if len(data) == 0: raise Exception('have to spcify a data to plot')
        if 'bende' in data:
            zdata += self.bendE
        if 'straine' in data:
            self.get_vstrainE()
            zdata += self.vstrainE
        if 'sube' in data:
            l0 = self.period / self.min_wavnum 
            ll = self.period / self.max_wavnum
            def kkl(x): return (2 * np.pi / (l0 + (ll - l0)/self.max_x * x))**4 * self.bmod / 2  # K/2
            zdata += kkl(self.verts[:,0]) * self.verts[:,2]**2
       
        
        gy, gz = self.get_crossdata_rbf((self.verts[:,:2],zdata), x, self.min_ry, \
                                            self.max_ry, num=ynum, k=k)
        sp = (self.max_ry-self.min_ry)/(ynum-1)
        intz = (2 * integrate.romb(gz,sp)) if self.mirror else integrate.romb(gz,sp) # we need times 2 if this is mirror
        return intz
    
    def plt_itg_data(self, data='straine', xi=None, xf=None, tol=0.001, samplekw=None, normalize = True, 
                        normalize_x =False, flip_x=False, elabel='', plt=True, ret=False, **kw):
        """Plot the \int stress dy over x.
        Arguments:
            whichstress {str} -- 'xx' for stress_xx, 'yy' for stress_yy, 'xy' for stress_xy.
        Keyword Arguments:
            xi {float} -- start point. (default: None)
            xf {float} -- end point. (default: None)
            tol {float} -- tolerance for adaptive sampling. (default: 0.001)
            samplekw {dict} -- sample function keywords. (default: None)
            normalize {bool} -- if True, norm stress by width and tension. (default: True)
            normalize_x {bool} -- if True, norm x axis by width. (default: True)
            flip_x {bool} -- flip the plot, let edge start from x=0. (default: False)
            ret {bool} -- if True, return x,y coords. (default: False)
            **kw -- kargs for plt_figure_mpl.
        """
        if xi is None:
            xi = self.min_rx
        if xf is None:
            xf = self.max_rx
        if samplekw is None:
            samplekw = {}
        X, Y = sample_function(lambda x: self.get_itg_data(x, data=data), [xi, xf], tol=tol, **samplekw)
        
        #itgs =  [self.get_itg_stressy(x) for x in X]
        if flip_x:
            X = self.max_rx - X
        if normalize_x:
            X = X/self.w
        if normalize:
            Y = Y/self.w
        if plt:
            plt_figure_mpl(X,Y, xlabel='$x$', ylabel=data, **kw)
        if ret:
            return X,Y
        
    #@jit
    def get_itg_stressy(self, x, ynum=513, **kw): # ynum should be 2**k + 1, using romb method to get high precision
        gy, gz = self.get_xcross_stressy(x, ynum=ynum, **kw)
        #intz = integrate.cumtrapz(gz.flatten(), gy.flatten(), initial=0)
        sp = (self.max_ry-self.min_ry)/(ynum-1)
        intz = (2 * integrate.romb(gz,sp)) if self.mirror else integrate.romb(gz,sp) # we need times 2 if this is mirror
        return intz
    
    def get_itg_stressxy(self, x, ynum=513, **kw): # ynum should be 2**k + 1, using romb method to get high precision
        gy, gz = self.get_xcross_stressxy(x, ynum=ynum, **kw)
        #intz = integrate.cumtrapz(gz.flatten(), gy.flatten(), initial=0)
        sp = (self.max_ry-self.min_ry)/(ynum-1)
        intz = (2 * integrate.romb(gz,sp)) if self.mirror else integrate.romb(gz,sp) # we need times 2 if this is mirror
        return intz

    def get_itg_stressx(self, x, ynum=513, **kw): # ynum should be 2**k + 1, using romb method to get high precision
        gy, gz = self.get_xcross_stressx(x, ynum=ynum, **kw)
        #intz = integrate.cumtrapz(gz.flatten(), gy.flatten(), initial=0)
        sp = (self.max_ry-self.min_ry)/(ynum-1)
        intz = (2 * integrate.romb(gz,sp)) if self.mirror else integrate.romb(gz,sp) # we need times 2 if this is mirror
        return intz
    
    def get_itg_all_stressy(self, xi, xf=None, xnum=513, ynum=513, norm=False, **kw):
        if xf is None:
            xf = optimize.brentq(self.get_itg_stressy, self.max_rx - 0.4 * self.w, self.max_rx-0.65*self.w, args=257)
        print(xf)
        xid = np.linspace(xi, xf, xnum)
        itgy = [self.get_itg_stressy(x, ynum=ynum) for x in xid]
        sp = (xf-xi)/(xnum-1)
        itg = integrate.romb(itgy, sp)
        
        return itg/self.tension/self.w/(xf-xi) if norm else itg

    def plt_itg_stress(self, whichstress, xi=None, xf=None, yi=None, yf=None, tol=0.001, samplekw=None, normalize = True, 
                        normalize_x =True, flip_x=False, elabel='', plt=True, ret=False, **kw):
        """Plot the \int stress dy over x.
        Arguments:
            whichstress {str} -- 'xx' for stress_xx, 'yy' for stress_yy, 'xy' for stress_xy.
        Keyword Arguments:
            xi {float} -- start point. (default: None)
            xf {float} -- end point. (default: None)
            tol {float} -- tolerance for adaptive sampling. (default: 0.001)
            samplekw {dict} -- sample function keywords. (default: None)
            normalize {bool} -- if True, norm stress by width and tension. (default: True)
            normalize_x {bool} -- if True, norm x axis by width. (default: True)
            flip_x {bool} -- flip the plot, let edge start from x=0. (default: False)
            ret {bool} -- if True, return x,y coords. (default: False)
            **kw -- kargs for plt_figure_mpl.
        """
        if xi is None:
            xi = self.min_rx
        if xf is None:
            xf = self.max_rx
        if samplekw is None:
            samplekw = {}
        
        if whichstress == 'yy':
            if yi is not None or yf is not None:
                X, Y = sample_function(lambda x : self.get_itg_stressy(x, yi=yi, yf=yf), [xi, xf], tol=tol, **samplekw)
            else:
                X, Y = sample_function(self.get_itg_stressy, [xi, xf], tol=tol, **samplekw)
        elif whichstress == 'xy':
            X, Y = sample_function(self.get_itg_stressxy, [xi, xf], tol=tol, **samplekw)
        elif whichstress == 'xx':
            X, Y = sample_function(self.get_itg_stressx, [xi, xf], tol=tol, **samplekw)
        else:
            raise ValueError('should be xx, yy or xy.')
        #itgs =  [self.get_itg_stressy(x) for x in X]
        if flip_x:
            X = self.max_rx - X
        if normalize_x:
            X = X/self.w
        if normalize:
            Y = Y/self.w/self.tension
        if plt:
            plt_figure_mpl(X,Y, xlabel='$x$', ylabel=r'$\frac{1}{TW} \int \sigma_{' + whichstress + '} dy$',
                            title=r'$\bar{{\sigma}}_{{{5}}},\ \tilde{{h}}={0:.3g},\ T={1},\ \tilde{{L}}={2:.4f},$ {3} {4}'.format(self.h_n, self.tension, np.around(self.L_n,4), self.plabel, elabel, whichstress), **kw)
        if ret:
            return X,Y

    def plt_itg_stressy_ply(self, xi=None, xf=None, pn=50, normalize = True, 
                            normalize_x =False, flip_x=False, elabel='', ret=False, **kw):
        """Plot the \int stress_y dy over x

        Keyword Arguments:
            xi {float} -- start point of plot, None to start from most begining (default: {None})
            xf {float} -- end point of plot, None to end at most ending. (default: {None})
            pn {int or array} -- number of samples of the plotted curve, 
                                        can be array like [[[x0,x1], pn1],[[x1,x2], pn2],...],
                                        if array, xi and xf will not be effect  (default: {50})
            normalize {bool} -- if normalize the integrated stress by width and tension (default: {True})
            normalize_x {bool} -- if normalize x axis by width [description] (default: {False})
            flip_x {bool} -- if flip x axis so that the clamped edge locates at x=0 [description] (default: {False})
            elabel {str} -- extra plot title [description] (default: {''})
        """

        if xi is None:
            xi = self.min_rx
        if xf is None:
            xf = self.max_rx
        if hasattr(pn, '__len__'):
            X = []
            for i in range(len(pn)):
                r = pn[i][0]
                pni = pn[i][1]
                X.extend(np.linspace(r[0],r[1],pni)[:-1])
            X.append(r[1])
            X = np.array(X)
        else:
            X = np.linspace(xi,xf, pn)
        itgs =  [self.get_itg_stressy(x) for x in X]
        if flip_x:
            X = self.max_rx - X
        if normalize_x:
            X = X/self.w
        if normalize:
            itgs = itgs/self.w/self.tension

        t = go.Scatter(
            x=X,
            y=itgs
        )
        layout = dict(
            title=r'$\int \sigma_y,\ \tilde{{h}}={0:.3g},\ T={1},\ \tilde{{L}}={2:.4f},\ \text{{{3}}} {4}$'.format(self.h_n, self.tension, np.around(self.L_n,4), self.plabel, elabel),
            xaxis=dict(
                title='x'
            ),
            yaxis=dict(
                title='$\\int \\sigma_y$'
            ),
            **kw
        )
        f = go.Figure(data=[t], layout=layout)
        ply.init_notebook_mode(connected=True)
        ply.iplot(f)
        if ret:
            return f

    def plt_itg_stressx(self, xi=None, xf=None, pn=50, normalize=True, elabel='', **kw):
        if xi is None:
            xi = self.min_rx
        if xf is None:
            xf = self.max_rx
        X = np.linspace(xi,xf, pn)
        itgs =  [self.get_itg_stressx(x) for x in X]
        if normalize:
            itgs = itgs/self.w
        t = go.Scatter(
            x=X,
            y=itgs
        )
        layout = dict(
            title='$\\int \\sigma_x,\\ h={0},\\ T={1},\\ \\text{{{2}{3}}}$'.format(self.thickness, self.tension, self.plabel, elabel),
            xaxis=dict(
                title='x'
            ),
            yaxis=dict(
                title='$\\int \\sigma_x$'
            ),
            **kw
        )
        f = go.Figure(data=[t], layout=layout)
        ply.init_notebook_mode(connected=True)
        ply.iplot(f)

    def plt_delta_x(self, xi=None, xf=None, tol=0.001, samplekw=None, normalize = 'arclength', mpl=True,
                        normalize_x =True, flip_x=False, elabel='', plt=True, ret=False, **kw):
        if xi is None:
            xi = self.min_x
        if xf is None:
            xf = self.max_x
        if samplekw is None:
            samplekw = {}

        X, Y = sample_function(lambda x: self.get_delta(x, normalize=normalize), [xi, xf], tol=tol, **samplekw)

        #itgs =  [self.get_itg_stressy(x) for x in X]
        if flip_x:
            X = self.max_x - X
        if normalize_x:
            X = X/self.w
            xlabel='$x/W$'
        else:
            xlabel='$x$'
        if plt:
            if mpl:
                plt_figure_mpl(X, Y, xlabel=xlabel, ylabel='$\\tilde{\\Delta}$'
                                , title=r'$\tilde{{\Delta}}, \tilde{{T}}={0}, \tilde{{h}}={1:.2E}$'.format(self.tension, self.h_n), **kw)
            else:
                plt_curve(X,Y, xaxis=xlabel, yaxis='$\\tilde{\\Delta}$'
                                , title=r'$\tilde{{\Delta}}, \tilde{{T}}={0}, \tilde{{h}}={1:.2E}$'.format(self.tension, self.h_n), **kw)
        if ret:
            return X, Y

    def get_crosssurf(self, xy, c, ci=None, cf=None, cnum=500,
                         k = 'cubic', inref=False, scaled = False,
                        usempl=False, ret=False, **kw):
        if scaled:
            c = c * self.w
        if xy == 'x':
            #X = c
            if ci is None:
                if not inref:
                    ci = self.min_y if not self.mirror else (-self.maxy(c))
                else:
                    ci = self.min_ry if not self.mirror else (-self.max_ry)
            if cf is None:
                if not inref:
                    cf = self.maxy(c)
                else:
                    cf = self.max_ry
            #Y = np.linspace(ci,cf, cnum)
            xt = 'y'
            a = 0
        elif xy == 'y':
            if ci is None:
                if not inref:
                    ci = self.min_x if not self.mirror else (-self.max_x)
                else:
                    ci = self.min_rx if not self.mirror else (-self.max_rx)
            if cf is None:
                if not inref:
                    cf = self.max_x
                else:
                    cf = self.max_rx
            #X = np.linspace(ci,cf,cnum)
            #Y = c
            xt = 'x'
            a = 1
        else:
            raise ValueError("xy must be 'x' or 'y'")
        #grid_x, grid_y = np.meshgrid(X,Y)
        #print(ci, cf)
        if self.mirror:
            self.get_mirror_coords()
            #Z = griddata(self.mverts[:,:2], self.mverts[:,2], (grid_x, grid_y), method=k, **kw)
            if not inref:
                xx, yy = self.get_crossdata_rbf(self.mverts, c, ci, cf, num=cnum, a=a, wd=3, k=k, **kw)
            else:
                self.get_mirror_refcoords()
                xx, yy = self.get_crossdata_rbf((self.mrefverts, self.mverts[:,2]), c, ci, cf, num=cnum, a=a, wd=3, k=k, **kw)
        else:
            #Z = griddata(self.verts[:,:2], self.verts[:,2], (grid_x, grid_y), method=k, **kw)
            xx, yy = self.get_crossdata_rbf(self.verts, c, ci, cf, num=cnum, a=a, wd = 3, k=k, **kw)
        return xx, yy, xt

    def get_curvature(self, xy, c, num=500, edge_order=(2,2), **kw):
        xx, yy, xt = self.get_crosssurf(xy, c, cnum=num)
        g1 = np.gradient(yy, xx, edge_order=edge_order[0])
        g2 = np.gradient(g1, xx, edge_order=edge_order[1])
        cvat = g2/(1+g1**2)**(3/2)
        return xx, cvat, xt
        
    def plt_curvature(self, xy, c, num=500, edge_order=(2,2), **kw):
        xx, cv, xt = self.get_curvature(xy, c, num=num, edge_order=edge_order)
        plt_figure_mpl(xx, cv, xlabel='${}$'.format(xt), ylabel='Curvature', title='Curvature')
        


    def plt_xcross_stressx(self, x, ynum=500, k='cubic', elabel='', **kw):
        """Plot x cross line of stress_x

        Arguments:
            x {float} -- the cross line position
            **kw {extra kargs} -- extra arg for rbf or griddata

        Keyword Arguments:
            ynum {int} -- how many sample points in the cross line (default: {500})
            k {string} -- method for rbf
            elabel {string} -- extra string in figure title
        """
        gy, gz = self.get_xcross_stressx(x, ynum, k, **kw)
        t = go.Scatter(
            x=gy,
            y=gz
        )
        layout= dict(
            title='$\\sigma_x(x={0},y),\\ \\tilde{{h}}={1:.3g},\\ T={2} \\text{{{3}}} {4}$'.format(x,self.thickness, self.tension,self.plabel, elabel),
            xaxis=dict(
                title='y'
            ),
            yaxis=dict(
                title='$\\sigma_x$'
            )
        )
        ply.init_notebook_mode(connected=False)
        f=go.Figure(data=[t],layout=layout)
        ply.iplot(f)

    def get_wavenumber(self, x='m', thr= 0.001, numsample=500, method='peak', 
                        retcoord=False):
        """get the wave number of the wrinkle
        Arguments:
            x {float} -- Position to get wavenumber.
        Karguments:
            thr {float} -- peak smaller than this threshold, will not be counted.
            numsample {int} -- num of sample in y direction
            method {str} -- 'peak' count peak, 'npeak' count n-peak, 'all' count both
            revert {bool} -- if True, central is always peak 
            retcoord {bool} -- Return coords of peaks if True
        Returns:
            method == 'peak': 
            num_peak, <peak_coords>, 
            'npeak'
            num_npeak, <npeak_coords>
            'all'
            num_peak+num_npeak, num_peak, <peak_coords>, num_npeak, <npeak_coords>
        """
        if x == 'm':
            x = self.get_wrinkle_peak()[0]
        yy, zz = self.get_crossdata_rbf(self.verts, x, self.min_ry, self.maxy(x), num=numsample, a=0, wd = 3)
        if zz[0] < 0:
            zz = -zz
        if self.mirror:
            yy = np.concatenate((-yy[::-1][:-1], yy))
            zz = np.concatenate((zz[::-1][:-1], zz))
        if method == 'peak':
            up = True
            unp = False
        elif method == 'npeak':
            unp = True
            up = False
        elif method == 'all':
            up = True
            unp = True
        #npeaks = 0
        ret= []
        if up:
            ma = max(zz)
            mask=zz>=thr*ma
            cz = zz[mask]
            cy = yy[mask]
            arm = argrelmax(cz)[0]
            npeaks = len(arm)
            ret.append(npeaks)
            if retcoord:
                peakc = (cy[arm],cz[arm])
                ret.append(peakc)
            #return npeaks, peakc
        if unp:
            mi = min(zz)
            maski = zz<=thr*mi
            cz = zz[maski]
            cy = yy[maski]
            arm = argrelmax(cz)[0]
            nnpeaks = len(arm)
            ret.append(nnpeaks)
            if retcoord:
                npeakc = (cy[arm],cz[arm])
                ret.append(npeakc)
            #return npeaks, peakc
        if up and unp:
            ret.insert(0, npeaks+nnpeaks)
        return ret

    def plt_surface(self, xi=None, xf=None, yi=None, yf=None, numx=500, numy=500, normbywidth=False, data='surf',
                     usemesh=False, mirror=None, revert=False, straindensity=True, bendedensity=True,
                    k='cubic', ar=None, ret=False, onlysurf=False, colorscale='Rainbow', camera=None,
                    retXYZ=False, contour=False, reverse=False,
                    deletebound=False,  removebase=None, **kw):
        """plot the surface

        plot the surface from in {xi, xf}

        Arguments:
            **kw  -- extra args used for griddata

        Keyword Arguments:
            xi {float} -- initial x (default: min_x)
            xf {float} -- final x (default: max_x)
            normbywidth {bool} -- if normlize the scale by width. (default: False)
            ar {array} -- aspect ratio (x,y,z) (default: None)
            ret {bool} -- return data if True (default: False)
            onlysurf {bool} -- don't show the gird and axis if True. (default: False)
        """
        if data.casefold() != 'surf':
            zdata = np.zeros(self.verts[:,2].size) if not self.mirror else np.zeros(self.mverts[:,2].size)
            data = data.casefold().split()
            if len(data) == 0: raise Exception('have to spcify a data to plot')
            if 'bende' in data:
                self.get_vbendE(density=bendedensity)
                zdata += self.vbendE
            if 'straine' in data:
                self.get_vstrainE(density=straindensity)
                zdata += self.vstrainE
            if 'sube' in data:
                l0 = self.period / self.min_wavnum 
                ll = self.period / self.max_wavnum
                def kkl(x): return (2 * np.pi / (l0 + (ll - l0)/self.max_x * x))**4 * self.bmod / 2  # K/2
                zdata += kkl(self.verts[:,0]) * self.verts[:,2]**2
                # assert zdata.shape[0] == self.verts[:,2].size
            if 'stressxy' in data:
                self.get_vstress_xy()
                zdata = self.vstress_xy
            elif 'stressy' in data:
                self.get_vstress_y()
                zdata = self.vstress_y
            elif 'stressx' in data:
                self.get_vstress_x()
                zdata = self.vstress_x
            elif 'stressr' in data:
                self.get_vstress_r()
                zdata = self.vstress_r
            elif 'stresst' in data:
                self.get_vstress_t()
                zdata = self.vstress_t
            elif 'stressrt' in data:
                self.get_vstress_rt()
                zdata = self.vstress_rt

        elif data.casefold() == 'surf':
            if self.mirror:
                self.get_mirror_coords()
                zdata = self.mverts[:,2]
            else:
                zdata = self.verts[:,2]
            #print(zdata.shape)
        else:
            raise Exception('No such data')

        if removebase is not None:
            zdata -= removebase

        if not usemesh:
            if xi is None:
                xi = self.min_x if not self.mirror else - self.max_x
            if xf is None:
                xf = self.max_x
            X = np.linspace(xi, xf, numx)
            if yi is None:
                yi = self.min_y if not self.mirror else - self.max_y
            if yf is None:
                yf = self.max_y
            Y = np.linspace(yi, yf, numy)
            grid_x, grid_y = np.meshgrid(X,Y)
            if self.mirror:
                self.get_mirror_coords()
                Z = griddata(self.mverts[:,:2], self.mverts[:,2], (grid_x, grid_y), method=k, **kw)
            else:
                Z = griddata(self.verts[:,:2], zdata, (grid_x, grid_y), method=k, **kw)
            
            if self.polar and hasattr(self, 'in_r'):
                sel = np.hypot(grid_x, grid_y) < self.rin_r
                Z[sel] = np.nan
                    
            #h = self.thickness
            if normbywidth:
                X = X / self.w
                Y = Y / self.w
                Z = Z / self.w
                h = '\\tilde{h}=' + str(self.h_n)
            else:
                h = 'h=' + str(self.thickness)
            if revert:
                Z = -Z 

            if deletebound:
                pass
            
            if retXYZ:
                return X, Y, Z
            
            data = [go.Surface(x=X,y=Y,z=Z, colorscale=colorscale)]
            
            if ar is None:
                ar = (5,2,1) if not self.mirror else (10,2,1)   
            layout = dict(
                title='$\\text{{Surface}},\\ {0},\\ \\text{{{1}}}$'.format(h, self.plabel),
                scene=dict(
                    aspectratio=dict(
                        x=ar[0],
                        y=ar[1],
                        z=ar[2]
                    )#,
                    #xaxis=dict(title='$x$'),
                    #axis=dict(title='$y$'),
                    #zaxis=dict(title='$z$')
                )
            )
            if onlysurf:
                layout['scene'].update(dict(
                    xaxis=dict(
                        visible=False,
                        title='',
                        showgrid=False,
                        showline=False
                    ),
                    yaxis=dict(
                        visible=False,
                        title='',
                        showgrid=False,
                        showline=False
                    ),
                    zaxis=dict(
                        visible=False,
                        title='',
                        showgrid=False,
                        showline=False
                    )))
                
            fig = go.Figure(data=data,layout=layout)
            if camera is not None:
                camera = dict(eye=dict(x=camera[0], y=camera[1], z=camera[2]))
                fig.update_layout(scene_camera=camera)
            ply.init_notebook_mode(connected=False)
            ply.iplot(fig)
        else:
            plt_mesh(np.concatenate((self.verts[:,:2], zdata[:, None]), axis=-1), self.fverts, mirror=mirror, ar=ar, onlysurf=onlysurf, camera=camera, **kw)
        if ret:
            return data
    
    def plt_contour(self, title=None, xlabel='$x$', ylabel='$y$', xi=None, xf=None, yi=None, yf=None, numx=500, numy=500, normbywidth=False,revert=False, cmap='binary', rotate=False,
                    vmin=None, vmax=None, k='cubic', ar=1, interpolation='quadric', bad_color='k', extent=None, straindensity=True, bendedensity=True,
                     dpi=150, figname=None, ret=False, colorbar=True, data='surf', removebase=None, onlysurf=False, rotatedeg=0, **kw):
        """plot contour of surface
        
        Keyword Arguments:
            cmap {str} -- camp of plot (default: {'binary'})
            vmin {float} -- min scale value of color (default: {None})
            vmax {float} -- max scale value of color (default: {None})
            ar {int or str} -- aspect ratio (default: {1})
            interpolation {str} -- interpolation method  (default: {'quadric'})
            bad_color {str} -- nan value color (default: {'k'})
        """
        plt.rcParams['mathtext.fontset'] = 'cm'
        X, Y, Z = self.plt_surface(retXYZ=True, xi=xi, xf=xf, yi=yi, yf=yf, numx=numx, numy=numy, 
                                   normbywidth=normbywidth, revert=revert, k=k, data=data,
                                   straindensity=straindensity, bendedensity=bendedensity)
        fig, ax = plt.subplots()
        if type(cmap) is str:
            cmap = plt.cm.get_cmap(cmap)
            cmap.set_bad(color=bad_color, alpha=0)
        if removebase is not None:
            Z -= removebase
        if vmin is not None: vmin=np.nanmin(Z) * vmin
        if vmax is not None: vmax = np.nanmax(Z) * vmax
        #if rotatedeg is not None:
        tr = transforms.Affine2D().rotate_deg(rotatedeg)
        if rotate:
            Z = np.rot90(Z, axes=(0,1))
            if extent is None: extent = [Y.min(), Y.max(), X.min(), X.max()]
        else:
            if extent is None: extent = [X.min(), X.max(), Y.min(), Y.max()]
        # print(Z)
        im = ax.imshow(Z, cmap=cmap, extent=extent,
                       aspect=ar, vmax=vmax, vmin=vmin, interpolation=interpolation, origin='lower', transform=tr + ax.transData, **kw)
        if onlysurf:
            colorbar = False
            xlabel=None
            ylabel=None
            ax.axis('off')
        if colorbar:
            fig.colorbar(im)
        if title is not None: ax.set_title(title)
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        # print(Z)
        fig.set_dpi(dpi)
        # fig
        if figname is not None:
            fig.savefig(figname, bbox_inches='tight')
        if ret:
            return X,Y,Z

            
    def plt_cross_surf(self, xy, c, ci=None, cf=None, cnum=500, mirror=None, reverse=False,
                        normbywidth=True, k='cubic', inref=False, scaled = False, wd=3,
                        usempl=False, ret=False, plt=True, data='surf', **kw):
        """Plot cross view of surface, where x or y = c

        Arguments:
            xy {str} -- 'x' or 'y', specify x cross or y cross
            c {float} -- cross position like x=5 or y=3.2
            ci {float} -- begining of cross view
            cf {float} -- end of cross view
            inref {bool} -- if inref, plot in reference coords.
            scaled {bool} -- if use scaled length, like values = c * self.w (default: {False})
            **kw {kargs} -- other args for rbf method

        Keyword Arguments:
            cnum {int} -- the number of sample points in this plot (default: {500})
        """
        if mirror is None:
            mirror = self.mirror
        if scaled:
            c = c * self.w
        if xy == 'x':
            #X = c
            if ci is None:
                if not inref:
                    ci = self.min_y if not mirror else (-self.maxy(c))
                else:
                    ci = self.min_ry if not mirror else (- self.max_ry)
            if cf is None:
                if not inref:
                    cf = self.maxy(c)
                else:
                    cf = self.max_ry
            #Y = np.linspace(ci,cf, cnum)
            xt = 'y'
            a = 0
        elif xy == 'y':
            if ci is None:
                if not inref:
                    ci = self.min_x if not mirror else (-self.max_x)
                else:
                    ci = self.min_rx if not mirror else (-self.max_rx)
            if cf is None:
                if not inref:
                    cf = self.max_x
                else:
                    cf = self.max_rx
            #X = np.linspace(ci,cf,cnum)
            #Y = c
            xt = 'x'
            a = 1
        elif xy == 'r':
            if ci is None:
                ci = 0
            if cf is None:
                cf = 2 * np.pi
            a = 0
            xt=r'\theta'
            normbywidth = False
        elif xy == 't':
            if ci is None:
                ci = self.in_r
            if cf is None:
                cf = self.out_r
            a = 1
            xt='r'
        else:
            raise ValueError("xy must be 'x' or 'y', 'r', 't'")
        #grid_x, grid_y = np.meshgrid(X,Y)
        #print(ci, cf)
        if 'surf' not in data:
            if 'stressxy' in data:
                self.get_vstress_xy()
                zdata = self.vstress_xy
            elif 'stressy' in data:
                self.get_vstress_y()
                zdata = self.vstress_y
            elif 'stressx' in data:
                self.get_vstress_x()
                zdata = self.vstress_x
            elif 'stressr' in data:
                self.get_vstress_r()
                zdata = self.vstress_r
            elif 'stresst' in data:
                self.get_vstress_t()
                zdata = self.vstress_t
            else:
                print('No such data option, surf, stressr, stresst')
                return
            fdata = (self.verts[:, :2], zdata)
        else:
            if mirror:
                self.get_mirror_coords()
                # zdata = self.mverts[:, 2]
                if not inref:
                    fdata = self.mverts
                else:
                    fdata = (self.mrefverts, self.mverts[:,2])
            else:
                fdata = self.verts
            
        xx, yy = self.get_crossdata_rbf(fdata, c, ci, cf, num=cnum, a=a, wd = wd, k=k, polar=self.polar, **kw)
        # if mirror:
        #     self.get_mirror_coords()
        #     #Z = griddata(self.mverts[:,:2], self.mverts[:,2], (grid_x, grid_y), method=k, **kw)
        #     if not inref:
        #         xx, yy = self.get_crossdata_rbf(self.mverts, c, ci, cf, num=cnum, a=a, wd=wd, k=k, **kw)
        #     else:
        #         self.get_mirror_refcoords()
        #         xx, yy = self.get_crossdata_rbf((self.mrefverts, self.mverts[:,2]), c, ci, cf, num=cnum, a=a, wd=wd, k=k, **kw)
        # else:
        #     #Z = griddata(self.verts[:,:2], self.verts[:,2], (grid_x, grid_y), method=k, **kw)
        #     xx, yy = self.get_crossdata_rbf(self.verts, c, ci, cf, num=cnum, a=a, wd = wd, k=k, polar=self.polar, **kw)
        #     print(xx.shape, yy.shape)
        #print(xx, Z)
        #xx, yy = self.get_crossdata_rbf(self.mverts, c, ci, cf, num=cnum, a=a)
        if normbywidth and 'surf' in data:
            xx= xx/self.w
            yy = yy/self.w
        if reverse:
            yy = -yy
        if plt:
            if not usempl:
                data = [go.Scatter(
                    x = xx,
                    y = yy,
                    mode = 'lines'
                )]
                layout = dict(
                    title='$\\text{{Surface cross}},\\ h={0},\\ T={1},\\ \\text{{{2}}}$'.format(self.thickness, self.tension, self.plabel),
                    scene=dict(
                        aspectratio=dict(
                            x=5 if not mirror else 10,
                            y=1,
                        )
                    ),
                    xaxis=dict(title=xt),
                    yaxis=dict(title='z')
                )
                fig = go.Figure(data=data,layout=layout)
                ply.init_notebook_mode(connected=False)
                ply.iplot(fig)
            else:
                plt_figure_mpl(xx, yy, title='Surface cross, $h={0},\\ T={1}$, {2}'.format(self.thickness, self.tension, self.plabel),
                                xlabel='${}$'.format(xt),
                                ylabel='$z$')
        if ret:
            return (xx, yy)

    def plt_cross_fft(self, xy, c, ci=None, cf=None, cnum=1000, mirror=None,
                        normbywidth=False, k = 'cubic', wavlen=False, ret=False, plt=True, **kw):
        """plot fouier transformation of cross section at axis, x=c or y=c, i.e. amplitude of each wave number or wave len
        
        Arguments:
            xy {str} -- 
            c {[type]} -- [description]
        
        Keyword Arguments:
            ci {[type]} -- [description] (default: {None})
            cf {[type]} -- [description] (default: {None})
            cnum {int} -- [description] (default: {500})
            wavlen {bool} -- xaxis is wavlen, default is wave number (default: {False})
        
        Returns:
            [type] -- [description]
        """
        x, y = self.plt_cross_surf(xy, c, ci=ci, cf=cf, cnum=cnum, mirror=mirror,
                              normbywidth=normbywidth, k=k, ret=True, plt=False)
        sp = np.fft.rfft(y)
        amp = np.abs(sp/y.size) * 2 # get amplitude
        fr = np.fft.rfftfreq(y.size, d=x[1]-x[0])
        if not wavlen: # x axis default is wav num
            xa = fr * (x[-1] - x[0])
        else: # xaxis is wave length
            xa = 1/fr
        if plt:
            plt_figure_mpl(xa, amp, **kw)
        if ret:
            return xa, amp
    
    def get_fmode_a_k(self, c, k, ci=None, cf=None, cnum=1000, xy='x', mirror=None):
        """get a_k at x=c, k is a array
        """
        x, y = self.plt_cross_surf(xy, c, ci=ci, cf=cf, cnum=cnum, mirror=mirror,
                              normbywidth=False, k='cubic', ret=True, plt=False)
        # if simple:
        #     
        #     fx = y*np.exp(-2j*np.pi*k * x/w) 
        #     return np.abs((fx[1:-1].sum() + (fx[0] + fx[-1])/2) * (x[1] - x[0])) / w * 2
        # else:
        w = x[-1] - x[0]
        return np.abs(integrate.trapz(y*np.exp(-2j*np.pi* k *x/w), dx=x[1]-x[0])/ w)*2
    
    def get_fmode_phase(self, c, k, ci=None, cf=None, cnum=1000, xy='x', mirror=None, th=0):
        x, y = self.plt_cross_surf(xy, c, ci=ci, cf=cf, cnum=cnum, mirror=mirror,
                                   normbywidth=False, k='cubic', ret=True, plt=False)

        N = y.size - 1
        if type(k) is int: 
            return np.angle((y[:-1] * np.exp(-2j * np.pi * np.arange(N)/N * k)).sum())
        else:
            k = np.array(k)
            if np.max(k) > (cnum-1)/2: raise Exception('max k should be smaller than (cnum-1)/2')
            ak = (y[:-1] * np.exp(-2j * np.pi * np.arange(N)/N * k[:,None])).sum(axis=1)
            amp = np.abs(ak)**2 * k**2
            return np.array([np.angle(a) if p > th else None for a, p in zip(ak, amp)])
    
    def get_fmode_complex(self, c, k, ci=None, cf=None, cnum=1000, xy='x', mirror=None, th=0):
        x, y = self.plt_cross_surf(xy, c, ci=ci, cf=cf, cnum=cnum, mirror=mirror,
                                   normbywidth=False, k='cubic', ret=True, plt=False)
        N = y.size - 1
        if type(k) is int: 
            return (y[:-1] * np.exp(-2j * np.pi * np.arange(N)/N * k)).sum()
        else:
            k = np.array(k)
            if np.max(k) > (cnum-1)/2: raise Exception('max k should be smaller than (cnum-1)/2')
            return (y[:-1] * np.exp(-2j * np.pi * np.arange(N)/N * k[:,None])).sum(axis=1)
    
    def plt_fft_mode_complex(self, xy, k, ci=None, cf=None, num=50, ret=False, plt=True, **kw):
        """Plot fft mode a[k] amplitude changes over c
        
        Arguments:
            xy {[type]} -- [description]
            k {int or list of int} -- a[k]
        
        Keyword Arguments:
            ci {[type]} -- [description] (default: {None})
            cf {[type]} -- [description] (default: {None})
            tol {float} -- [description] (default: {0.001})
            samplekw {[type]} -- [description] (default: {None})
            ret {bool} -- [description] (default: {False})
            plt {bool} -- [description] (default: {True})
        
        Returns:
            [type] -- [description]
        """
        if xy == 'x':
            if ci is None: ci = self.min_x
            if cf is None: cf = self.max_x
            
        else:
            if ci is None:  ci = self.min_y #if not self.mirror else -self.maxy(c)
            if cf is None: cf = self.max_y
        if not hasattr(k, '__len__'): k = [k]
        X=[]
        Yr=[]
        Yi=[]
        names= []
        # if not hasattr(k, '__len__'):
        #     k = [k]
        # k = np.array(k)
        xx= np.linspace(ci, cf, num)
        for kk in k:
            yy = np.array([self.get_fmode_complex(x, kk, xy=xy) for x in xx])
            X.append(xx)
            Yr.append(yy.real)
            Yi.append(yy.imag)
            names.append('$f_{{{}}}$'.format(kk))
        if plt:
            plt_figure_mpl(X+X,Yr+Yi, name=names + names, xlabel='${}$'.format(xy), ylabel='$A$')
        if ret:
            return (X,Yr), (X,Yi), names

    def plt_fft_phase(self, xy, k, ci=None, cf=None, tol=0.001, samplekw=None, ret=False, plt=True, th=0, num=100, **kw):
        if xy == 'x':
            if ci is None: ci = self.min_x
            if cf is None: cf = self.max_x
            
        else:
            if ci is None:  ci = self.min_y #if not self.mirror else -self.maxy(c)
            if cf is None: cf = self.max_y
        if samplekw is None:
            samplekw = {}
        if not hasattr(k, '__len__'): k = [k]
        X= np.linspace(ci, cf, num)
        names= ['$f_{{{}}}$'.format(kk) for kk in k]

        Y= np.transpose([self.get_fmode_phase(x, k, th=th, **kw) for x in X])
        xx = [X] * len(k)
        # for kk in k:
        #     xx, yy = sample_function(lambda x: self.get_fmode_phase(x, kk, xy=xy, **kw), [ci, cf], tol=tol, **samplekw)
        #     X.append(xx)
        #     Y.append(yy)
        #     names.append('$f_{{{}}}$'.format(kk))
        if plt:
            plt_figure_mpl(xx,Y, name=names, xlabel='${}$'.format(xy), ylabel='$A$')
        if ret:
            return xx,Y, names
        

    def plt_fft_mode(self, xy, k, ci=None, cf=None, tol=0.001, samplekw=None, ret=False, plt=True, **kw):
        """Plot fft mode a[k] amplitude changes over c
        
        Arguments:
            xy {[type]} -- [description]
            k {int or list of int} -- a[k]
        
        Keyword Arguments:
            ci {[type]} -- [description] (default: {None})
            cf {[type]} -- [description] (default: {None})
            tol {float} -- [description] (default: {0.001})
            samplekw {[type]} -- [description] (default: {None})
            ret {bool} -- [description] (default: {False})
            plt {bool} -- [description] (default: {True})
        
        Returns:
            [type] -- [description]
        """
        if xy == 'x':
            if ci is None: ci = self.min_x
            if cf is None: cf = self.max_x
            
        else:
            if ci is None:  ci = self.min_y #if not self.mirror else -self.maxy(c)
            if cf is None: cf = self.max_y
        if samplekw is None:
            samplekw = {}
        if not hasattr(k, '__len__'): k = [k]
        X=[]
        Y=[]
        names= []
        # if not hasattr(k, '__len__'):
        #     k = [k]
        # k = np.array(k)

        for kk in k:
            xx, yy = sample_function(lambda x: self.get_fmode_a_k(x, kk, xy=xy), [ci, cf], tol=tol, **samplekw)
            X.append(xx)
            Y.append(yy)
            names.append('$f_{{{}}}$'.format(kk))
        if plt:
            plt_figure_mpl(X,Y, name=names, xlabel='${}$'.format(xy), ylabel='$A$')
        if ret:
            return X,Y, names
        

    def get_arclength(self, x, cnum=1000, k = 'cubic', **kw):
        """Get the 1/2 arc length at x = x

        Arguments:
            x {float} -- the position to get arc length
        """
        yy, zz = self.get_crossdata_rbf(self.verts, x, self.min_ry, self.maxy(x), num=cnum, a=0, wd = 3, k=k, **kw)
        # al = 0
        return np.hypot(yy[1:] - yy[:-1], zz[1:] - zz[:-1]).sum()
        #for i in range(len(yy)-1):
        #    al += norm((yy[i+1] - yy[i], zz[i+1] - zz[i]))
        # return al

    def get_delta(self, x, cnum=1000, k='cubic', normalize=None, **kw):
        if normalize == 'arclength':
            return 1 - self.maxy(x)/self.get_arclength(x, cnum=cnum, k=k, **kw)
        elif normalize == 'projectedlength':
            return self.get_arclength(x, cnum=cnum, k=k, **kw)/self.maxy(x) - 1
        elif normalize == 'width':
            return 2*(self.get_arclength(x, cnum=cnum, k=k, **kw) - self.maxy(x))/self.w
        return 2*(self.get_arclength(x, cnum=cnum, k=k, **kw) - self.maxy(x))

    def plt_ybound(self,tol=0.001, samplekw={}, normalize_xy=True, flip_x=False, plt=True, mpl=False,  displacement=False):
        """Plot the y boundary

        Keyword Arguments:
            xnum {int} -- num of points (default: {500})
        """
        X, Y = sample_function(self.maxy, [self.min_x, self.max_x], tol=tol, **samplekw)
        if flip_x:
            X = (self.max_x - X)[::-1]
            Y = Y[::-1]
        if displacement:
            Y = Y - self.max_ry

        if normalize_xy:
            X = X/self.w
            Y = Y/self.w
        if plt:
            if not mpl:
                plt_curve(X,Y, ar=(5,1), title='$\\text{{Boundary}},\\ h={0},\\ T={1},\\ \\text{{{2}}}$'.format(self.thickness, self.tension, self.plabel),
                        xaxis='x', yaxis='y')
            else:
                plt_figure_mpl(X,Y, xlabel='$x$', ylabel='$y$', 
                            title='Boundary, $h={0},\\ T={1},$ {2}'.format(self.thickness, self.tension, self.plabel))
        else:
            return X,Y
    
    def plt_defined_delta(self, tol=0.001, samplekw={}, normalize_xy=True, 
                          flip_x=True, plt=True, mpl=False, shift=False, svk=False):
        xb, yb = self.plt_ybound(tol=tol, samplekw=samplekw, normalize_xy=normalize_xy, flip_x=flip_x, plt=False) 
        if shift:
            if flip_x:
                delta = (yb[-1] + yb[-2]) - yb * 2
            else:
                delta = yb[0]+ yb[1] - yb * 2
        else:
            exx = np.sqrt(1-2 * self.pratio * self.tension) if not svk else np.sqrt(1- 2 * self.pratio * self.tension/np.sqrt(1+2 * self.tension))
            delta = exx - yb * 2 
        if plt:
            plt_figure_mpl(xb, delta, title=r'Defined $\tilde{\Delta}$', xlabel='$x$', ylabel=r'$\tilde{\Delta}$')
        else:
            return xb, delta
            
    
    def get_ybound_derivative(self, x, dx=None, originx=False, normalize_xy=True):
        if not originx:
            x *= self.w
        if dx is None:
            dx = (self.max_x - self.min_x) / 100000

        X = [x-2*dx, x-dx, x+dx, x+2*dx]
        #print(dx*self.w)
        Y = [self.maxy(i) for i in X]
        return (-Y[-1] + 8* (Y[-2] - Y[1]) + Y[0])/12/dx

    def plt_wavlength_x(self, tol=0.001, threshold=0.01, m='cubic', normlize=True, verbose=True, plt=True, ret=False):
        #self.get_crosssurf('y', 0)
        X,Y = self.get_crossdata_rbf(self.verts, 0, self.min_x, self.max_x, num=500, a=1, wd = 3)
        #argm = np.argmax(np.abs(Y))
        # X[argm], Y[argm]
        if Y[0] <0:
            Y = -Y
        Y = Y - threshold * self.thickness
        argneg = np.argmax(Y < 0)
        
        if X[argneg-1] == 0:
            max_x = X[argneg-1]
        else:
            f = interp1d(X,Y, kind=m)
            max_x = optimize.brentq(f, X[argneg-1], X[argneg+1])
            
        if verbose:
            print(max_x/self.w)
        xx, yy = sample_function(lambda x: self.get_wavelength(x)[-1], [self.min_x, max_x], tol=tol)
        if normlize:
            xx = xx/self.w
            yy = yy/self.w
        if plt:
            plt_curve(xx,yy, name='wavelength vs x', xaxis='$x$', yaxis='wavelength')
        if ret:
            return xx,yy
        

    def plt_arclength(self,xnum=500, theory=False):
        """Plot the 1/2 arclength changing with x

        Keyword Arguments:
            xnum {int} -- num of points (default: {500})
            theory {bool} -- plot out theory W/2 * (1- pratio * T)
        """
        X = np.linspace(self.min_x, self.max_x, xnum)
        Y = list(map(self.get_arclength, X))
        data = [go.Scatter(
            x = X,
            y = Y,
            mode = 'lines',
            name='1/2 Arclength'
        )]
        if theory:
            wd = 2 * self.max_ry if self.mirror else self.max_ry - self.min_ry
            Yt = [self.max_ry * np.sqrt(1 - 2 * self.pratio * self.get_itg_stressx(0)/wd)] * xnum
            data.append(go.Scatter(
                x=X,
                y=Yt,
                mode='lines',
                name='\\ (1-\\nu T)*W/2'
            ))
        layout = dict(
            title='$1/2\\ \\text{{Arclength}},\\ h={0},\\ T={1},\\ \\text{{{2}}}$'.format(self.thickness, self.tension, self.plabel),
            scene=dict(
                aspectratio=dict(
                    x=5,
                    y=1,
                )
            ),
            xaxis=dict(title='x'),
            yaxis=dict(title='Arclength')
        )
        fig = go.Figure(data=data,layout=layout)
        ply.init_notebook_mode(connected=False)
        ply.iplot(fig)



    def get_tensionfiled(self, fid=None):
        if fid is None:
            if self.tensionfiled is None:
                self.get_stress()
                #self.tensionfiled = [eigh([self.fstress[i,:2],self.fstress[i,1:]])
                #                                    for i in range(self.nf)]
                fsm = np.array([self.fstress[:, :2], self.fstress[:, 1:]])
                self.tensionfiled = eigh(np.transpose(fsm, axes=(1,2,0))) # it is [eigenvales, eigenvectors]
        else:
            self.get_stress()
            return eigh([self.fstress[fid,:2],self.fstress[fid,1:]])

    def get_plt_data_tensionfiled(self, xy='x',xi=None, xf=None, yi=None, yf=None, num=5000):
        self.get_tensionfiled()
        self.get_reffcentroid()
        if xi is None:
            xi = self.min_rx
        if xf is None:
            xf = self.max_rx
        if yi is None:
            yi = self.min_ry
        if yf is None:
            yf = self.max_ry
        ind= [i for i in range(num) if self.reffcentroid[i,0] <= xf and self.reffcentroid[i,0] >= xi\
                                        and self.reffcentroid[i,1] <= yf and self.reffcentroid[i,1] >= yi]
        tf = np.array([self.tensionfiled[0][i] * self.tensionfiled[1][i] for i in ind])
        
        refc = np.array([self.reffcentroid[i] for i in ind])
        X = refc[:, 0]
        Y = refc[:, 1]
        U2 = tf[:, :, 0][:,0]
        V2 = tf[:, :, 0][:,1]
        U1 = tf[:, :, 1][:,0]
        V1 = tf[:, :, 1][:,1]
        #return np.append(X,X), np.append(Y,Y), np.append(U1, U2), np.append(V1, V2)
        return X, Y, U1, U2, V1, V2
    
    def plt_tensionfield_mpl(self, xy='xy', xi=None, xf=None, yi=None, yf=None, num=5000,
                             scale=1, units='xy', angles='xy',
                             quiverprop={}, **kw):
        X, Y, U1, U2, V1, V2 = self.get_plt_data_tensionfiled(xi=xi, xf=xf, yi=yi, yf=yf, num=num)
        l = X.shape[0]
        #color = (['r']*l).extend(['b']*l)
        fig, ax = plt.subplots()
        if xy == 'x':
            #X = X[l:]
            #Y = Y[l:]
            U = -U1
            V = -V1
            color=['b']*l
        elif xy == 'y':
            #X = X[:l]
            #Y = Y[:l]
            U = -U2
            V = -V2
            color = ['r']*l
        elif xy == 'xy':
            X = np.append(X, X)
            Y = np.append(Y, Y)
            U = -np.append(U1, U2)
            V = -np.append(V1, V2)
            color = (['b']*l).extend(['l']*l)
        ax.quiver(X, Y, U, V, color=color, scale=scale, units=units, angles=angles, **quiverprop)
        ax.set_aspect(1)
        plt_any_mpl(ax, fig, **kw)
        plt.show()
    

    def plot_theta(self, xy='xy', xi=None, xf=None, yi=None, yf=None, num=5000):
        X, Y, U1, U2, V1, V2 = self.get_plt_data_tensionfiled(xi=xi, xf=xf, yi=yi, yf=yf, num=num)
        #tan = V1/U1
        theta = np.arctan2(V1, U1)
        


    def plt_tensionfield(self,xy='x',xi=None, xf=None, yi=None, yf=None, num=5000, h=300, data=[], **kw):
        self.get_tensionfiled()
        self.get_reffcentroid()
        if xi is None:
            xi = self.min_rx
        if xf is None:
            xf = self.max_rx
        if yi is None:
            yi = self.min_ry
        if yf is None:
            yf = self.max_ry
        ind= [i for i in range(num) if self.reffcentroid[i,0] <= xf and self.reffcentroid[i,0] >= xi\
                                        and self.reffcentroid[i,1] <= yf and self.reffcentroid[i,1] >= yi]
        tf = np.array([self.tensionfiled[i][0]* self.tensionfiled[i][1] for i in ind])
        refc = np.array([self.reffcentroid[i] for i in ind])
        #data=[]
        if 'x' in xy:
            bx = []
            by = []
            for i in range(len(ind)):
                bx.extend((refc[i,0],refc[i,0]+tf[i,0,1], None))
                by.extend((refc[i,1],refc[i,1]+tf[i,1,1], None))
            data.append(go.Scatter(x=bx, y=by, mode='lines', name='$\\sigma_x\'$'))
        if 'y' in xy:
            bx = []
            by = []
            for i in range(len(ind)):
                bx.extend((refc[i,0],refc[i,0]+tf[i,0,0], None))
                by.extend((refc[i,1],refc[i,1]+tf[i,1,0], None))
            data.append(go.Scatter(x=bx, y=by, mode='lines', name='$\\sigma_y\'$'))
        #data = [q]
        #xri =
        #h = h
        w = (xf-xi)/(yf-yi) * h
        #print(h,w)
        #print([xi,xf])
        layout = dict(
            title='$\\text{{Tension field}},\\ h={0},\\ T={1},\\ \\text{{{2}}}$'.format(self.thickness, self.tension, self.plabel),
            #hovermode='closest',
            autosize=False,
            width= w,
            height= h,
            xaxis=dict(
                title='x',
                range=[xi,xf]
            ),
            yaxis=dict(
                title='y',
                range=[yi,yf]
            ),
            **kw
        )
        fig = go.Figure(data=data,layout=layout)
        ply.init_notebook_mode(connected=False)
        ply.iplot(fig)

    def get_wrinkle_peak(self):
        X,Y = self.get_crossdata_rbf(self.verts, 0, self.min_x, self.max_x, num=500, a=1, wd = 3)
        argm = np.argmax(np.abs(Y))
        return X[argm], Y[argm]    

    def get_wavelength(self, x='m',  usetwowav=True, onlycentralwav=2, num=500, usefmode=False,):
        """Get the wave length at x=x

        Arguments:
            onlycentraltwo {bool} -- use only central two waves
            x {float} -- position to calculate the wave length
            usefmode {bool} -- use fourier modes
        """
        if x == 'm':
            x = self.get_wrinkle_peak()[0]
        if usefmode:
            xa, amp = self.plt_cross_fft('x', x, ret=True, plt=False)
            so = np.argmax(amp)
            return self.w / xa[so]
        xx, yy = self.get_crossdata_rbf(self.verts, x, self.min_ry, self.max_ry, num=num, a=0, wd = 3)
        maxima = (argrelmax(yy))[0]
        minima = (argrelmin(yy))[0]
        if yy[1] < yy[0]: # peak in the middle
            maxima=np.insert(maxima,0,0,axis=0)
        elif yy[1] > yy[0]:
            minima=np.insert(minima,0,0,axis=0)
        #print(maxima, minima)
        wlmaxX = xx[maxima]
        wlminX = xx[minima]
        #print(maxima, minima, wlmaxX,wlminX)
        l = len(maxima)
        if l >= 3 and usetwowav:
            wlmax = (wlmaxX[2]-wlmaxX[0])/2
        elif l >= 2:
            wlmax = (wlmaxX[1]-wlmaxX[0])
        else:
            wlmax = None
        l = len(minima)
        if l >= 3 and usetwowav:
            wlmin = (wlminX[2]-wlminX[0])/2
        elif l >= 2:
            wlmin = (wlminX[1]-wlminX[0])
        else:
            wlmin = None
        if not (wlmin is None or wlmax is None):
            wlave = (wlmin+wlmax)/2
        else:
            wlave = None
        
        #if onlycentralwav == 1:
        if yy[1] < yy[0]:
            wlcen1 = np.abs(wlminX[0]*2)
        else:
            wlcen1 = np.abs(wlmaxX[0]*2)
        #elif onlycentralwav ==2:
        if yy[1] < yy[0]:
            wlcen2 = np.abs(wlmaxX[1])
        else:
            wlcen2 = np.abs(wlminX[1])

        return (wlmax, wlmin, wlave, wlcen1, wlcen2)

    def get_k(self, x, ynum=513, **kw):
        """sigma_yy / lambda**2 + (2*pi)**2 * B/lambda**4
        """
        gy, gz = self.get_xcross_stressy(x, ynum=ynum, **kw)
        #intz = integrate.cumtrapz(gz.flatten(), gy.flatten(), initial=0)
        gz = np.clip(gz, a_max=0, a_min=None)
        plt_figure_mpl(gy,gz)
        sp = (self.max_ry-self.min_ry)/(ynum-1)
        intz = (2 * integrate.romb(gz,sp)) if self.mirror else integrate.romb(gz,sp) # we need times 2 if this is mirror
        #sy_ave = intz / self.
        #print(sy_ave)
        wavlen = self.get_wavelength(x)[-1]
        k = gz[0] / wavlen**2 + (2*np.pi)**2 * self.thickness**2/12/(1-self.pratio**2) / wavlen**4
        print((gz[0] / wavlen**2) / ((2*np.pi)**2 * self.thickness**2/12/(1-self.pratio**2) / wavlen**4))
        return k

    def get_amplitude(self,x = None, inte=True, o=2):
        """Get the amplitude of wrinkle at x=x
        use the maximum of the amplitude
        Keyword Arguments:
            x {float} -- position to get amplitude, if None, get the largest (default: {None})
            inte {bool} -- using integral as Amplitude
            o {int} -- order in A = (\int A(y)^o dy)^(1/o)
        """
        if x is None:
            xx, zz = self.get_crossdata_rbf(self.verts, 0, self.min_x, self.max_x, num=1000, a=1, wd = 3)
            #maxima = (argrelmax(zz))[0]
            absz=np.abs(zz)
            ind = np.argmax(absz)
            x=xx[ind]
            if not inte:
                return (x, absz[ind])
        ynum=513
        maxy= self.maxy(x)
        sp = (maxy-self.min_ry)/(ynum-1)
        gx, gz = self.get_crossdata_rbf(self.verts, x, self.min_ry, maxy, num=ynum, a=0, wd = 3)
        gz = np.abs(gz)
        if inte:
            gzo = gz**o
            inta = (2 * integrate.romb(gzo,sp))
            return (x, inta**(1/o))
        else:
            return (x, np.max(gz))


        #minima = (argrelmin(yy))[0]



def fix_plot_latex(): # fix the issue of latex does not work
    from IPython.core.display import display, HTML
    # The polling here is to ensure that plotly.js has already been loaded before
    # setting display alignment in order to avoid a race condition.
    display(HTML(
        '<script>'
            'var waitForPlotly = setInterval( function() {'
                'if( typeof(window.Plotly) !== "undefined" ){'
                    'MathJax.Hub.Config({ SVG: { font: "STIX-Web" }, displayAlign: "center" });'
                    'MathJax.Hub.Queue(["setRenderer", MathJax.Hub, "SVG"]);'
                    'clearInterval(waitForPlotly);'
                '}}, 250 );'
        '</script>'
    ))

def load_output(file):
    with open(file, 'r') as f:
        d = f.readlines()
    if 'type' in d[0]:
        type = d[0].split()[-1]
    else:
        type = 'txt'
    if type == 'txt':
        return type, None
    #out = []
    # for l in d[1:]:
    #     if not l.startswith('#') and l.strip():
    #         out.append(l)
    return type, ' '.join(d[1:])

