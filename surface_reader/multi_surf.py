# This script is for dealing with multi surface data in a folder
# created by Meng Xin
import sys
import os
import glob
import numpy as np
import plotly.offline as ply
import plotly.graph_objs as go
from scipy.optimize import curve_fit
from surface import EvolveSurf, fix_plot_latex
from other_funcs import plt_curve, plt_figure_mpl, plt_func1d, sample_function

class SurfSet(object):
    """Get the set of surface data from a folder
    """
    def __init__(self, filename, newdata=True, sort=None, plabel='', mirror=True):
        """Get the set of surface data from path=folder
        sort {str} -- 'h' sort by thickness, 't' for tension, None don't sort. 
        """
        #if not folder: # folder is empty string, use current path
        self.fl = []
        if hasattr(filename, '__len__') and not isinstance(filename, str):
            for f in filename:
                self.fl.extend(glob.glob(f)) # file list
        else:
            #print(glob.glob(filename))
            self.fl.extend(glob.glob(filename))
        if len(self.fl) == 0:
            print("Error: there is no specified file in this folder")
            sys.exit(2)
        self.fl.sort()
        #print(self.fl)
        self.s = [] # keep all the data of surfaces
        self.hs = None # sorted thickness
        self.ts = None # sorted tension
        for f in self.fl:
            self.s.append(EvolveSurf(f, newdata=newdata, mirror=mirror))
        self.numsurf = len(self.s)
        self.s = np.array(self.s, dtype=object)
        if sort is not None:
            self.sort(sort)
        self.plabel = plabel
        fix_plot_latex()
        ply.init_notebook_mode(connected=False)
    
    def sort(self, sort):
        if sort == 'h':
            if self.hs is None:
                hs = np.array([ss.h_n for ss in self.s])
                so = np.argsort(hs)
                self.s = self.s[so]
                self.hs = hs[so]
                self.epsilons = [ss.epsilon for ss in self.s]
        elif sort == 't':
            if self.ts is None:
                ts = np.array([ss.tension for ss in self.s])
                so = np.argsort(ts)
                self.s = self.s[so]
                self.ts = ts[so]
    
    def get_torh(self, torh):
        if torh is None:
            if self.hs is not None:
                xx = self.hs
                label='h'
            elif self.ts is not None:
                xx = self.ts
                label='T'
            else:
                raise Exception("Surface is not sorted, specify 'torh'.")
        elif torh == 'h':
            if self.hs is None:
                self.sort(torh)
            xx = self.hs
            label = 'h'
        elif torh == 't' or torh == 'T':
            if self.ts is None:
                self.sort(torh)
            xx = self.ts
            label='T'
        else:
            raise Exception('sort by T or h')
        return xx, label

    def get_amplitude_t_dL(self, x=None,inte=True, o=2, tl=None, tu=None, usedL=False, retargsort=False):
        """return (tension, amplitude) of the set, if tl or tu non-zero, only return 
        the amplitude when tl <= tension <= tu
        
        Keyword Arguments:
            inte {bool} -- if use the integral to get a average amplitude, false to use largest amplitude (default: {True})
            o {int} -- the order for integral A=(\int A**o dy)**(1/o) (default: {1})
            tl {int} -- lower tension bound to get the amp (default: {None})
            tu {int} -- upper tension bound to get the amp (default: {None})
            usedL {bool} -- return delta_L/L instead of tension (default: {False})
            retargsort {bool} -- return argsort (default: {False})
        """
        ten = []
        amp = []
        for i in range(len(self.s)):
            if ((tl is None or self.s[i].tension >= tl) and (tu is None or self.s[i].tension <=tu)):
                if not usedL:
                    ten.append(self.s[i].h_n)
                else:
                    ten.append(self.s[i].deltaLN) # deltaLN = delta_L/L
                amp.append((self.s[i].get_amplitude(x=x,inte=inte,o=o))[1])
        ten = np.array(ten)
        amp = np.array(amp)
        arg = np.argsort(ten)
        #ten = ten[arg]
        #amp = amp[arg]
        if retargsort:
            return (ten[arg], amp[arg], arg)
        else:
            return (ten[arg], amp[arg])
    

    def plt_wavenum(self, x='m', thr= 0.001, numsample=500, method='peak', 
                    torh=None, ret=False, **kw):
        xx, xlabel = self.get_torh(torh)
        wns = np.zeros(self.numsurf)
        for i in range(self.numsurf):
            wns[i] = self.s[i].get_wavenumber(x, thr, numsample, method, retcoord=False)[0]
        
        plt_figure_mpl(xx, wns, title=r'Wave num v.s. $\tilde{{{}}}$'.format(xlabel),
                        xlabel=r'$\tilde{{{}}}$'.format(xlabel), 
                        ylabel='Wave number', marker='.', **kw)
        if ret:
            return xx, wns

    def plt_cross_surf(self, xy, c, ci=None, cf=None, cnum=500, 
                    k = 'cubic',label='t', reverse=False, ftr=None, scaled = False,
                    normbywith=True, ret = False,
                    **kw):
        """Plot cross view of surface, where x or y = c

        Arguments:
            xy {str} -- 'x' or 'y', specify x cross or y cross
            c {float} -- cross position like x=5 or y=3.2
            ci {float} -- begining of cross view
            cf {float} -- end of cross view
            **kw {kargs} -- other args for rbf method

        Keyword Arguments:
            cnum {int} -- the number of sample points in this plot (default: {500})
            k {str} -- method to use for rbf interpolation (default: {'cubic'})
            label {str} -- plot label, 't' means tension is different, 'h' means thickness is different
            ftr {list like} -- just plot the indices specified in ftr (sorted indices)
        """

        if xy == 'x':
            #X = c
            if ci is None:
                # ci = self.s[0].min_ry if not self.s[0].mirror else (-self.s[0].max_ry)
                ci = [-ss.maxy(c) if self.s[0].mirror else ss.min_y for ss in self.s]
            if cf is None:
                # cf = self.s[0].max_ry
                cf = [ss.maxy(c) for ss in self.s]
            #Y = np.linspace(ci,cf, cnum)
            xa = 'y'
            a = 0
        elif xy == 'y':
            if ci is None:
                ci = [self.s[0].min_x if not self.s[0].mirror else (-self.s[0].max_x)] * len(self.s)
            if cf is None:
                cf = [self.s[0].max_x] * len(self.s)
            #X = np.linspace(ci,cf,cnum)
            #Y = c
            xa = 'x'
            a = 1
        else:
            raise ValueError("xy must be 'x' or 'y'")
        #grid_x, grid_y = np.meshgrid(X,Y)
        #print(ci, cf)
        xx=[]
        yy=[]
        if label == 't':
            tens = [self.s[i].tension for i in range(len(self.s))]
            sorto = np.argsort(tens)
        else:
            ths = [self.s[i].thickness for i in range(len(self.s))]
            sorto = np.argsort(ths)
        if ftr is not None and (not isinstance(ftr,str)):
            sorto = sorto[ftr]
        for i in sorto:
            if scaled:
                c = c* self.s[i].w
            if self.s[i].mirror:
                self.s[i].get_mirror_coords()
                xt, yt = self.s[i].get_crossdata_rbf(self.s[i].mverts, c, ci[i], cf[i], num=cnum, a=a, wd=3, k=k, **kw)
                if normbywith:
                    xt = xt/self.s[i].w
                    yt = yt/self.s[i].w
                xx.append(xt)
                yy.append(yt)
            else:
                #Z = griddata(self.verts[:,:2], self.verts[:,2], (grid_x, grid_y), method=k, **kw)
                xt, yt = self.s[i].get_crossdata_rbf(self.s[i].verts, c, ci[i], cf[i], num=cnum, a=a, wd = 3, k=k, **kw)
                if normbywith:
                    xt = xt/self.s[i].w
                    yt = yt/self.s[i].w
                xx.append(xt)
                yy.append(yt)

        if reverse:
            nid = []
            mid = int(cnum/2)
            #print len(yy[i])
            for i in range(len(yy)):
                if yy[i][mid]<0:
                    nid.append(i)
            yy = np.array(yy)
            yy[nid] = -yy[nid]

        data = [go.Scatter(
            x = xx[i],
            y = yy[i],
            mode = 'lines',
            name = 'T={0}'.format(self.s[sorto[i]].tension) if label == 't' else '$\\tilde{{h}}={0}$'.format(self.s[sorto[i]].h_n)
        ) for i in range(len(xx))]
        layout = dict(
            title='$\\text{{Surface cross}},\\ h={0},\\ \\text{{{1}}}$'.format(self.s[0].thickness, self.plabel) \
                if label == 't' else '$\\text{{Surface cross}},\\ T={0},\\ \\text{{{1}}}$'.format(self.s[0].tension, self.plabel),
            scene=dict(
                aspectratio=dict(
                    x=5 if not self.s[0].mirror else 10,
                    y=1,
                )
            ),
            xaxis=dict(title=xa),
            yaxis=dict(title='z')
        )
        fig = go.Figure(data=data,layout=layout)
        ply.init_notebook_mode(connected=False)
        ply.iplot(fig)
        if ret:
            return (xx, yy)

    def plt_curvature(self, xy, c, num=500, torh=None, edge_order=(2,2), ret=False, **kw):
        names, xlabel =self.get_torh(torh)
        xx = []
        ct = []
        for ss in self.s:
            xxx, cct, _ = ss.get_curvature(xy, c, num=num, edge_order=edge_order)
            xx.append(xxx)
            ct.append(cct)
        name=[r'$\tilde{{{}}} = {}$'.format(xlabel, hh) for hh in names]
        plt_figure_mpl(xx, ct, xlabel='${}$'.format('x' if xy=='y' else 'y')
                        , ylabel='curvature', title='curvature', name=name, **kw)
        if ret:
            return xx, ct
    
    def plt_wavlength_x(self, torh=None, tol=0.001, threshold=0.01, m='cubic', normlize=True, verbose=True, ret=False, **kw):
        names, xlabel =self.get_torh(torh)
        xx = []
        wl = []
        for ss in self.s:
            xxx, wls = ss.plt_wavlength_x(tol=tol, threshold=threshold, m=m, normlize=normlize, verbose=verbose, plt=False, ret=True)
            xx.append(xxx)
            wl.append(wls)
        name=[r'$\tilde{{{}}} = {}$'.format(xlabel, hh) for hh in names]
        plt_figure_mpl(xx, wl, xlabel='$x$'
                        , ylabel='$\lambda$', title='wavelength vs $x$', name=name, **kw)
        if ret:
            return xx, wl

    
    def plt_itg_stressy(self, xi=None, xf=None, pn=50, normalize = True, 
                        normalize_x =True, flip_x=False, legend = 'h', 
                        filter = None, elabel='', usempl=True, ret=False, **kw):
        """Plot the \int stress_y dy over x
        
        Keyword Arguments:
            xi {float} -- start point of plot, None to start from most begining (default: {None})
            xf {float} -- end point of plot, None to end at most ending. (default: {None})
            pn {int or array of int} -- number of samples of the plotted curve, 
                                        can be array to specify every curve  (default: {50})
            normalize {bool} -- if normalize the integrated stress by width and tension (default: {True})
            normalize_x {bool} -- if normalize x axis by width [description] (default: {False})
            flip_x {bool} -- if flip x axis so that the clamped edge locates at x=0 [description] (default: {False})
            lengend {char} -- 'L' for l as lengend, 'T' for tension, 'h' for thickness (default: {'L'})
            filter {list} -- just plot the surfaces listed in the filter, indices are gave in the filter (defalt: {None})
            elabel {str} -- extra plot title [description] (default: {''})
        """
        xx = []
        yy = []
        name = []
        if filter is None:
            ss = self.s
        else:
            ss = [self.s[i] for i in filter]
            
        if (not hasattr(pn, '__len__')) or hasattr(pn[0], '__len__'):
                #(hasattr(pn[0],'__len__') and len(pn[0]) == 2 and not hasattr(pn[0][1], '__len__')):
            pn = [pn]*len(ss)
        #print(pn)
        for i in range(len(ss)):
            if xi is None:
                xi = ss[i].min_rx
            if xf is None:
                xf = ss[i].max_rx
            if hasattr(pn[i], '__len__'):
                X = []
                for j in range(len(pn[i])):
                    r = pn[i][j][0]
                    pni = pn[i][j][1]
                    #print(r, pni)
                    X.extend(np.linspace(r[0],r[1],pni)[:-1])
                X.append(r[1])
                X = np.array(X)
            else:
                X = np.linspace(xi, xf, pn[i])
            itgs = [ss[i].get_itg_stressy(x) for x in X]
            if flip_x:
                X = ss[i].max_rx - X
            if normalize_x:
                X = X/ss[i].w
            if normalize:
                itgs = itgs/ss[i].w/ss[i].tension
            if legend=='L' or legend == 'l':
                name.append(r'$\tilde{{L}}={:3g}$'.format(np.around(ss[i].L_n, 3)))
            elif legend == 't' or legend == 'T':
                name.append(r'$\tilde{{T}} = {:3g}$'.format(ss[i].tension))
            elif legend == 'h' or legend == 'H':
                name.append(r'$\tilde{{h}} = {:3g}$'.format(ss[i].h_n))
            xx.append(X)
            yy.append(itgs)
        if usempl:
            plt_figure_mpl(xx, yy, name=name, title=r'$\frac{1}{TW} \int \sigma_{yy} dy$',
                            xlabel='$x$', ylabel=r'$\frac{1}{TW} \int \sigma_{yy} dy$', **kw)
        else:
            plt_curve(xx, yy, title=r'$\int \sigma_y dy$', name=name, xaxis='$x$', yaxis=r'$\int \sigma_y dy$')
        if ret:
            return xx,yy
    
    def plt_itg_stress(self, whichstress, xi=None, xf=None, tol=0.001, samplekw={}, normalize = True, 
                        normalize_x =True, flip_x=False, torh=None, around=6, elabel='', usempl=True, ret=False, **kw):
        names, label = self.get_torh(torh)
        xx = []
        yy = []
        for ss in self.s:
            if xi is None:
                xi = ss.min_rx
            if xf is None:
                xf = ss.max_rx
            if whichstress == 'y' or whichstress == 'yy':
                X, Y = sample_function(ss.get_itg_stressy, [xi, xf], tol=tol, **samplekw)
            elif whichstress == 'x' or whichstress == 'xx':
                X, Y = sample_function(ss.get_itg_stressx, [xi, xf], tol=tol, **samplekw)
            if flip_x:
                X = ss.max_rx - X
            if normalize_x:
                X = X/ss.w
            if normalize:
                Y = Y/ss.w/ss.tension
            xx.append(X)
            yy.append(Y)
        name = [r'$\tilde{{{}}}={}$'.format(label, np.around(na,around)) for na in names]
        if usempl:
            plt_figure_mpl(xx, yy, name=name, title=r'$\frac{1}{TW} \int \sigma_{yy} dy$',
                            xlabel='$x$', ylabel=r'$\frac{1}{TW} \int \sigma_{yy} dy$', **kw)
        else:
            plt_curve(xx, yy, title=r'$\int \sigma_y dy$', name=name, xaxis='$x$', yaxis=r'$\int \sigma_y dy$')
        if ret:
            return xx, yy, name


    def plt_amplitude_t_dL(self,x=None, inte=True, o=2, tl=None, tu=None, 
                normalize=True, normbyh=False, usedL=False, ret=False):
        """plot amplitude vs tension, if tl or tu non-zero, only plot 
        the amplitude when tl <= tension <= tu
        
        Keyword Arguments:
            inte {bool} -- if using integral to get the amplitude (default: {True})
            o {int} -- the order for integral A=(\int A**o dy)**(1/o) (default: {1})
            tl {int} -- lower tension bound to plot the amp (default: {None})
            tu {int} -- upper tension bound to plot the amp (default: {None})
            normalize {bool} -- if normalize by width (default: {True})
            usedL {bool} -- use delta_L/L as xaxis (default: {False})
            ret {bool} -- return data if True (default: {False})
        """
        #tension = [s[i].tension for i in range(len(s))]
        #amplitude = [(s[i].get_amplitude())[1] for i in range(len(s))]
        ten, amp, arg = self.get_amplitude_t_dL(x=x,inte=inte, o=o, tl=tl, tu=tu, usedL=usedL, retargsort=True)
        if not usedL:
            xlabel = r'\tilde{T}'
        else:
            xlabel = r'\Delta \tilde{L}'
        if normalize:
            if normbyh: # norm by thickness
                norm = np.array([self.s[i].thickness for i in range(len(self.s))])
                norm = norm[arg]
            else: # norm by width
                norm = np.array([self.s[i].max_ry * 2 for i in range(len(self.s))])
                norm = norm[arg] 
            if inte:
                amp = np.divide(amp, norm**((o+1)/o))
                title='$ \\tilde{{A}}\\text{{ vs }} {2},\\ h={0},\\ \\tilde{{A}}=(\\int dy (A/W)^{1}/)^{{1/{1}}}/W$'.format(self.s[0].thickness, o, xlabel)
            else:
                amp = np.divide(amp, norm)
                title='$ \\tilde{{A}}\\text{{ vs }} {1},\\ h={0},\\ \\tilde{{A}}=A_p/W$'.format(self.s[0].thickness, xlabel)
            yaxis = '$\\tilde{A}$'
        else:
            if inte:
                title='$\\text{{A vs }} {2},\\ h={0},\\ A=(\\int dy A^{1}/)^{{1/{1}}}$'.format(self.s[0].thickness, o, xlabel)
            else:
                title='$\\text{{A vs }} {1},\\ h={0},\\ A=A_p$'.format(self.s[0].thickness, xlabel)
            yaxis='$A$'

        data = [go.Scatter(
            x=ten,
            y=amp,
            mode='lines+markers'
        )]
        layout = dict(
            title=title,
            scene=dict(
                aspectratio=dict(
                    x=10,
                    y=1,
                )
            ),
            xaxis=dict(title='${}$'.format(xlabel)),
            yaxis=dict(title=yaxis)
        )
        fig = go.Figure(data=data, layout=layout)
        ply.iplot(fig)
        if ret:
            return (ten,amp)
    
    def plt_delta_h_t(self, x, torh=None, delta_norm=None, plt=True, ret=False, **kw):
        xx, label = self.get_torh(torh)
        delta = [ss.get_delta(x, normalize=delta_norm) for ss in self.s]
        if plt:
            plt_figure_mpl(xx, delta, xlabel=r'$\tilde{{{}}}$'.format(label), ylabel=r'$\tilde{\Delta}$',
                        title=r'$\tilde{{\Delta}}$ vs $\tilde{{{}}}$'.format(label), **kw)
        if ret:
            return xx, delta
            

    def curv_fit_amp_ten(self, initialguess=None, delta=None, inte=True, o=1, tl=None,tu=None):
        """We expect A = b*(T/T_c-)^d
        We do logA = b + \\delta*log(T/T_c-1)
        use integral A=(\int A**o dy)**(1/o) to get the amplitude
        initialguess {array} -- The initial guess of (T_c, \\delta)
        """
        ten, amp = self.get_amplitude_t_dL(inte=inte, o=o, tl=tl, tu=tu)
        #lten = np.log(ten)
        lamp = np.log(amp)
        def func(x, t,d,b): #t <-> T_c, x <-> T
            return b+d*np.log(x/t-1)
        def fd(x, t,b):
            return b+delta*np.log(x/t-1)
        if delta is None:
            popt, pcov = curve_fit(func, ten, lamp, p0=initialguess)
            fitname = '$\\text{{fitting: }}\\log A = b + d \\log (T/T_c-1)$'
        else:
            popt, pcov = curve_fit(fd, ten, lamp, p0=initialguess)
            fitname = '$\\text{{fitting: }}\\log A = b + {0} \\log (T/T_c-1)$'.format(delta)
        data = [go.Scatter(
            x=ten,
            y=lamp,
            mode='markers',
            name='data'
        ),
        go.Scatter(
            x=ten,
            y=func(ten, *popt) if delta is None else fd(ten,*popt),
            mode='lines',
            name=fitname
        )]
        layout = dict(
            title='log(A) vs T and fitting',
            scene=dict(
                aspectratio=dict(
                    x=10,
                    y=1,
                )
            ),
            xaxis=dict(title='T'),
            yaxis=dict(title='log(A)')
        )
        fig = go.Figure(data=data,layout=layout)
        ply.iplot(fig)
        return popt
    
    def curv_fit_amp_ten_nonlog(self, initialguess=None, delta=None, inte=True, o=1, tl=None,tu=None, 
                                usedL=False):
        """We expect A = b*(T/T_c-1)^d
        use integral A=(\int A**o dy)**(1/o) to get the amplitude

        initialguess {array} -- The initial guess of (T_c, \\delta)
        delta {folat} -- d==delta, if delta is None, will fit it.
        """
        ten, amp = self.get_amplitude_t_dL(inte=inte, o=o, tl=tl, tu=tu, usedL=usedL)
        #lten = np.log(ten)
        #lamp = np.log(amp)
        def func(x, t,d,b): #t <-> T_c, x <-> T
            return b*(x/t-1)**d
        def fd(x, t,b):
            return b*(x/t-1)**0.5
        if delta is None:
            popt, pcov = curve_fit(func, ten, amp, p0=initialguess)
            fitname='$\\text{{fitting: }}A = b (T/T_c-1)^d$'
        else:
            popt, pcov = curve_fit(fd, ten, amp, p0=initialguess)
            fitname = '$\\text{{fitting: }}A = b (T/T_c-1)^{0}$'.format(delta)
        data = [go.Scatter(
            x=ten,
            y=amp,
            mode='markers',
            name='data'
        ),
        go.Scatter(
            x=ten,
            y=func(ten, *popt) if delta is None else fd(ten,*popt),
            mode='lines',
            name=fitname
        )]
        layout = dict(
            title='A vs T and fitting',
            scene=dict(
                aspectratio=dict(
                    x=10,
                    y=1,
                )
            ),
            xaxis=dict(title='T'),
            yaxis=dict(title='A')
        )
        fig = go.Figure(data=data,layout=layout)
        ply.iplot(fig)
        return popt

    def plt_energy_tension(self, dfakegbend=True, sheetonly=False, log=False,f=False, correct_t=False):
        """Plot energy vs tension, with fitting
        
        Keyword Arguments:
            eshift {float} -- energy shift for some reason (default: {0})
            dfakegbend {bool} -- delete the fake gbend (default: {True})
            sheetonly {bool} -- delete the work from the system, count the energy of sheet only (default: {False})
            log {bool} -- use log plot (default: {False})
            f {bool} -- use fitting (default: {False})
            correct_t {bool} -- use corrected tension
        """
        energy=[]
        ts = []
        for i in range(len(self.s)):
            if dfakegbend:
                if sheetonly:
                    energy.append(self.s[i].energy-self.s[i].fake_gbende + self.s[i].work)
                else:
                    energy.append(self.s[i].energy-self.s[i].fake_gbende)
            else:
                if sheetonly:
                    energy.append(self.s[i].energy + self.s[i].work)
                else:
                    energy.append(self.s[i].energy)
            ts.append(self.s[i].tension)
        soa = np.argsort(ts)
        energy = np.array(energy)
        ts = np.array(ts)
        if log: 
            lenergy =np.log(np.abs(energy[soa]))
            lts=np.log(ts[soa])
            title= 'log(E) vs log(T) (delete fake gbend)' if dfakegbend else 'log(E) vs log(T)'
            xtitle='log(T)'
            ytitle='log(E)'
        else:
            lenergy =energy[soa]
            lts=ts[soa]
            title='E vs T (delete fake gbend)' if dfakegbend else 'E vs T'
            xtitle='T'
            ytitle='E'
        if sheetonly:
            title += ' (sheet only)'
        parg=np.polyfit(lts,lenergy,1)
        fit=np.poly1d(parg)
        data = [go.Scatter(
            x=lts,
            y=lenergy,
            name='data',
            mode='markers'
        )]
        if f:
            parg=np.polyfit(lts,lenergy,1)
            fit=np.poly1d(parg)
            data.append(go.Scatter(
                x=lts,
                y=fit(lts),
                name='{:.2f}*log(T)+{:.2f}'.format(parg[0],parg[1]),
                mode='lines'
            ))
            
        layout = dict(
            title=title,
            scene=dict(
                aspectratio=dict(
                    x=10,
                    y=1,
                )
            ),
            xaxis=dict(title=xtitle),
            yaxis=dict(title=ytitle)
        )
        fig = go.Figure(data=data,layout=layout)
        ply.iplot(fig)
        return (lts,lenergy)

    def plt_delta_x(self, torh=None,  mask=None,
                 xi=None, xf=None, tol=0.001, samplekw=None, normalize = 'arclength', 
                        normalize_x =True, flip_x=False, elabel='', plt=True, ret=False, **kw):
        """plot delta vs x
        """
        names, label = self.get_torh(torh)
        xxs = []
        deltas = []
        if mask is None:
            mask = range(self.numsurf)
        for i in mask:
            x, y = self.s[i].plt_delta_x(plt=False, ret=True, xi=xi, xf=xf, tol=tol, samplekw=samplekw,
                                    normalize=normalize, normalize_x=normalize_x, flip_x=flip_x)
            xxs.append(x)
            deltas.append(y)
        lables = [r'$\tilde{{{0}}}={1:.2E}$'.format(label, names[i]) for i in mask]
        if plt:
            plt_figure_mpl(xxs, deltas, name=lables,
                            xlabel='$x/W$' if normalize_x else '$x$', ylabel=r'$\tilde{\Delta}$'
                            , title=r'$\tilde{\Delta}$', **kw)
        if ret:
            return xxs, deltas, lables
    
    def plt_ybound(self, torh=None,  mask=None, displacement=False,
                  tol=0.001, samplekw={},
                        normalize_xy =True, plt=True, flip_x=True, ret=False, **kw):
        """plot delta vs x
        """
        names, label = self.get_torh(torh)
        xxs = []
        deltas = []
        if mask is None:
            mask = range(self.numsurf)
        for i in mask:
            x, y = self.s[i].plt_ybound(plt=False, tol=tol, samplekw=samplekw, normalize_xy=normalize_xy, flip_x=flip_x, displacement=displacement)
            xxs.append(x)
            # if displacement:
            #     y = y - self.s[i].max_ry
            deltas.append(y)
        lables = [r'$\tilde{{{0}}}={1:.2E}$'.format(label, names[i]) for i in mask]
        if plt:
            plt_figure_mpl(xxs, deltas, name=lables,
                            xlabel='$x/W$' if normalize_xy else '$x$', ylabel='$y/W$' if normalize_xy else '$y$'
                            , title=r'Projected boundary', **kw)
        if ret:
            return xxs, deltas, lables
    
    def plt_defined_delta(self, torh=None, mask=None, tol=0.001, samplekw={}, flip_x=True, shift=False, svk=False,
                        normalize_xy =True, plt=True, ret=False, **kw):
        names, label = self.get_torh(torh)
        xxs = []
        deltas = []
        if mask is None:
            mask = range(self.numsurf)
        for ss in self.s[mask]:
            x, y = ss.plt_defined_delta(plt=False, tol=tol, samplekw=samplekw, normalize_xy=normalize_xy, flip_x=flip_x, shift=shift, svk=svk)
            xxs.append(x)
            deltas.append(y)
        lables = [r'$\tilde{{{0}}}={1:.2E}$'.format(label, n) for n in names[mask]]
        if plt:
            plt_figure_mpl(xxs, deltas, name=lables,
                            xlabel='$x/W$' if normalize_xy else '$x$', ylabel=r'$\tilde{\Delta}$' if normalize_xy else r'$\Delta$'
                            , title=r'Defined $\Delta$', **kw)
        if ret:
            return xxs, deltas, lables



    def plt_delta_tension_thick(self, x, normalize='arclength', torh=None, ret=False):
        """plot delta vs tension, delta = arclength - projected length
        
        Arguments:
            x {float} -- position to get the delta
        
        Keyword Arguments:
            normalize {str} -- normalize the delta, None-> no normalize, arclength-> n by arclength, projectedlength-> n by pl, width-> n by width (default: {None})
        """

        delta = []
        xx, label = self.get_torh(torh)
        for i in range(len(self.s)):
            delta.append(self.s[i].get_delta(x, normalize=normalize))
            #ts.append(self.s[i].tension)
        if normalize == 'width':
            width = [ss.w for ss in self.s]
            #width = np.array([self.s[i].max_ry * 2 for i in range(len(self.s))])
            delta = np.divide(delta, width)
            title = '$\\Delta/W\\ \\text{{vs}}\\ \\tilde{{{1}}}\\ \\text{{at}}\\ x={0}$'.format(x, label)
            yaxis= '$\\Delta/W$'
        elif normalize is None:
            title = '$\\Delta\\ \\text{{vs}}\\ \\tilde{{{1}}}\\ \\text{{at}}\\ x={0}$'.format(x, label)
            yaxis= '$\\Delta$'
        else:
            title = '$\\tilde{{\\Delta}}\\text{{ (normalized by {1}) }} \\text{{vs}}\\ \\tilde{{{2}}}\\ \\text{{at}}\\ x={0}$'.format(x, normalize, label)
            yaxis= '$\\tilde{\\Delta}$'
        #soa = np.argsort(ts)
        #ts = np.array(ts)
        #delta = np.array(delta)
        #ts = ts[soa]
        #delta = delta[soa]
        data = [go.Scatter(
            x=xx,
            y=delta
        )]
        #if normalize is None:
        #    title = '$\\Delta\\ \\text{{vs}}\\ T\\ \\text{{at}}\\ x={0}$'.format(x)
        #else:
        #    title = '$\\tilde{{\\Delta}}\\text{{ (normalized by {1}) }} \\text{{vs}}\\ T\\ \\text{{at}}\\ x={0}$'.format(x, normalize)
        layout = dict(
            title=title,
            scene=dict(
                aspectratio=dict(
                    x=10,
                    y=1,
                )
            ),
            xaxis=dict(title='$\\tilde{{{}}}$'.format(label)),
            yaxis=dict(title=yaxis)
        )
        fig = go.Figure(data=data, layout=layout)
        ply.iplot(fig)
        if ret:
            return xx, delta
    
    def plt_projectedlength_tension(self, x, ret=False, normalize=True):
        pl = np.zeros(len(self.s))
        ts = np.zeros(len(self.s))
        for i in range(len(self.s)):
            pl[i] = self.s[i].maxy(x) * 2
            ts[i] = self.s[i].tension
        if normalize:
            width = np.array([self.s[i].max_ry * 2 for i in range(len(self.s))])
            pl = np.divide(pl, width)
            title = 'Projected length/W vs T'
            yaxis = 'Projected length/W'
        else:
            title = 'Projected length vs T'
            yaxis = 'Projected length'
        soa = np.argsort(ts)
        ts = ts[soa]
        pl = pl[soa]
        return plt_curve(ts, pl, title=title, xaxis='T', yaxis=yaxis, ret=ret)
        
    def plt_arclength_tension(self, x, theory=True, ret=False, normalize=True):
        al = np.zeros(len(self.s))
        ts = np.zeros(len(self.s))
        for i in range(len(self.s)):
            al[i] = self.s[i].get_arclength(x) * 2
            ts[i] = self.s[i].tension
        if normalize:
            width = np.array([self.s[i].max_ry * 2 for i in range(len(self.s))])
            al = np.divide(al, width)
            title = 'Arclength/W vs T at x = {0}'.format(x)
            yaxis= 'Arclength/W'
            namet = '$ \\sqrt{1-2 \\nu S_{xx}} $'
        else:
            title = 'Arclength vs T at x = {0}'.format(x)
            yaxis = 'Arclength' 
            namet = '$ W \\sqrt{1-2 \\nu S_{xx}} $'
        soa = np.argsort(ts)
        ts = ts[soa]
        al = al[soa]
        xx=[ts]
        yy=[al]
        if theory:
            def arcl(t, n=normalize):
                if n:
                    return np.sqrt(1 - 2 * self.s[0].pratio * t/ \
                    (1/3*np.sqrt(3)*np.sqrt((1+(-1+54*t**2+6*np.sqrt(-3*t**2+81*t**4))**(1/3))**2/(-1+54*t**2+6*np.sqrt(-3*t**2+81*t**4))**(1/3))))
                else:
                    return (self.s[0].max_ry - self.s[0].min_ry) * 2 * np.sqrt(1 - 2 * self.s[0].pratio * t/ \
                    (1/3*np.sqrt(3)*np.sqrt((1+(-1+54*t**2+6*np.sqrt(-3*t**2+81*t**4))**(1/3))**2/(-1+54*t**2+6*np.sqrt(-3*t**2+81*t**4))**(1/3))))
            alt = np.real(arcl(ts + 0j))
            xx.append(ts)
            yy.append(alt)
        return plt_curve(xx,yy, title=title, name=['Data', namet] ,xaxis='T', yaxis=yaxis, ret=ret)
    
    def plt_wavelength_torh(self, x, torh = 'h', usetwowav=True, num=500, methodwav='ave', yrange=None, normalize=True,
                            ret=False):
        """Plot wavelength (at x=x) vs tension

        Arguments:
            x {float} -- the position to calculate the wavelength
        
        Keyword Arguments:
            usetwowav {bool} -- use two waves to calculate the wavelength (default: {True})
            methodwav {str} -- the method to get wavelength, 'peak'-> distanse between peaks, 'npeak'-> d between negtive peaks, 'ave'-> (peak+npeak)/2, 'fmode'-> major fmodes  (default: {'ave'})
            yrange {list} -- the range of y, [y_i, y_f], if None will auto (default: {None})
            normalize {bool} -- if normalize the wavelength by width (default: {True})
        """
        wl = np.zeros(len(self.s))
        ts = np.zeros(len(self.s))
        for i in range(len(self.s)):
            if methodwav == 'ave':
                wl[i] = (self.s[i].get_wavelength(x, usetwowav=usetwowav, num=num))[2]
            elif methodwav == 'peak':
                wl[i] = (self.s[i].get_wavelength(x, usetwowav=usetwowav, num=num))[0]
            elif methodwav == 'central1':
                wl[i] = (self.s[i].get_wavelength(x, usetwowav=usetwowav, num=num))[3]
            elif methodwav == 'central2':
                wl[i] = (self.s[i].get_wavelength(x, usetwowav=usetwowav, num=num))[4]
            elif methodwav == 'fmode':
                wl[i] = self.s[i].get_wavelength(x, usefmode=True)
            if torh == 't':
                ts[i] = self.s[i].tension  
            elif torh == 'h':
                ts[i] = self.s[i].h_n
        if normalize:
            width = np.array([ss.w for ss in self.s])  # np.array([self.s[i].max_ry * 2 for i in range(len(self.s))])
            wl = np.divide(wl, width)
            title = '$\\tilde{{\\lambda}}\\text{{ ({1}) }} \\text{{vs}}\\ T\\ \\text{{at}}\\ x={0}$'.format(x, methodwav)
            yaxis = '$\\tilde{{\\lambda}}$'
        else:
            title = '$\\lambda\\text{{ ({1}) }} \\text{{vs}}\\ T\\ \\text{{at}}\\ x={0}$'.format(x, methodwav)
            yaxis = '$\\lambda$'
        soa = np.argsort(ts)
        ts = ts[soa]
        wl = wl[soa]
        data = [go.Scatter(
            x=ts,
            y=wl
        )]
        layout = dict(
            title=title,
            scene=dict(
                aspectratio=dict(
                    x=10,
                    y=1,
                )
            ),
            xaxis=dict(title='$T$'),
            yaxis=dict(title=yaxis) if yrange is None else dict(title=yaxis,range=yrange)
        )
        fig = go.Figure(data=data, layout=layout)
        ply.iplot(fig)
        if ret:
            return ts, wl

    def plt_deltaLN_tension(self, ret=False):
        """Plot delta_LN vs tension, delta_LN = (strained L (length) - original L)/original L
        """
        dL = np.zeros(len(self.s))
        ts = np.zeros(len(self.s))
        for i in range(len(self.s)):
            dL[i] = self.s[i].deltaLN
            ts[i] = self.s[i].tension
        soa = np.argsort(ts)
        ts = ts[soa]
        dL = dL[soa]
        data = [go.Scatter(
            x=ts,
            y=dL
        )]
        layout = dict(
            title='$\\Delta \\tilde{L}\\ \\text{vs}\\ T$',
            scene=dict(
                aspectratio=dict(
                    x=10,
                    y=1,
                )
            ),
            xaxis=dict(title='$T$'),
            yaxis=dict(title='$\\Delta \\tilde{L}$')
        )
        fig = go.Figure(data=data, layout=layout)
        ply.iplot(fig)
        if ret:
            return (ts,dL)
    
    def plt_A_Delta_wavlength(self, x=0, inte=True, o=2, usetwowav=True, wavlengthind=2,
                            onlycentralwav=2, normbywidth=True,
                            log = False, fit= False, controledfit=False, 
                            initguess=None, normalize='arclength', ret=False):
        """Plot A ~ Delta^(1/2) wavlength, compare with E. Cerda's approach
        
        Arguments:
            x {float} -- position to get Delta, and wavlength
        
        Keyword Arguments:
            inte {bool} -- If use integral as amplitude (default: {True})
            o {int} -- order to get A (default: {2})
            log {bool} -- If use Log (default: {True})
            normalize {str} -- normalize Delta/normalize, 'arclength' or 'projectedlength' (default: {None})
        """
        #if normalize is not None:
        #    estr = '(normalized)'
        #else:
        #    estr = ''
        l = len(self.s)
        A = np.zeros(l)
        Delta = np.zeros(l)
        wavlength = np.zeros(l)
        xx = np.zeros(l)
        for i in range(l):
            A[i] = (self.s[i].get_amplitude(x, inte=inte, o=o))[1] 
            wavlength[i] = (self.s[i].get_wavelength(x, usetwowav=usetwowav, onlycentralwav=onlycentralwav))[wavlengthind]
            if normbywidth:
                A[i] = A[i] / self.s[i].w
                wavlength[i] = wavlength[i] / self.s[i].w
            Delta[i] = self.s[i].get_delta(x, normalize=normalize)
            xx[i] = Delta[i]**(1/2) * wavlength[i]
        if log:
            A = np.log(A)
            #Delta = np.log(Delta)
            #wavlength = np.log(wavlength)
            xx = np.log(xx)
            if normalize is None:
                title= '$ \\log(A) \\sim 1/2 \\log(\\Delta) + \\log(\\lambda) $'
                xaxis = '$1/2 \\log(\\Delta) + \\log(\\lambda) $'
            else:
                title= '$ \\log(A) \\sim 1/2 \\log(\\tilde{\\Delta}) + \\log(\\lambda) $'
                xaxis = '$1/2 \\log(\\tilde{\\Delta}) + \\log(\\lambda) $'
            yaxis = '$ \\log(A) $'
        else:
            if normalize is None:
                title = '$ A \\sim \\Delta^{1/2} \\lambda $'
                xaxis = '$\\Delta^{1/2} \\lambda $'
            else:
                title = '$ A \\sim \\tilde{\\Delta}^{1/2} \\lambda $'
                xaxis = '$\\tilde{\\Delta}^{1/2} \\lambda $'
            yaxis = '$ A $'
        soa = np.argsort(xx)
        A = A[soa]
        xx = xx[soa]
        data = [go.Scatter(
            x=xx,
            y=A,
            name = 'data'
        )]
        if fit:
            if not controledfit:
                p1 = np.polyfit(xx,A,1)
                f = np.poly1d(p1)
                yf = f(xx)
                namef = '{0}*x + {1}'.format(np.around(p1[0],4),np.around(p1[1],4))
            else:
                if log:
                    def fc(x,b):
                        return x+b
                else:
                    def fc(x,k):
                        return k*x
                #f = lambda x, k: k*x if not log else lambda x, b: x + b
                popt, pcrv = curve_fit(fc, xx, A, p0=initguess)
                yf = fc(xx, *popt)
                namef = '{0} x'.format(popt[0]) if not log else 'x + {0}'.format(np.around(popt[0],4))
            data.append(go.Scatter(
                x=xx,
                y=yf,
                name=namef
            ))
        layout = dict(
            title=title,
            scene=dict(
                aspectratio=dict(
                    x=10,
                    y=1,
                )
            ),
            xaxis=dict(title=xaxis),
            yaxis=dict(title=yaxis)
        )
        fig = go.Figure(data=data, layout=layout)
        ply.iplot(fig)
        if ret:
            return data







