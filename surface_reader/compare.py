from multi_surf import SurfSet
from surface import EvolveSurf
from other_funcs import plt_figure_mpl, plt_any_mpl, plt_func1d
import numpy as np

class Compare(object):
    def __init__(self, fname1, fname2):
        self.s1 = EvolveSurf(fname1)
        self.s2 = EvolveSurf(fname2)
        self.aspect1 = '{:.2g}'.format(self.s1.L_n)
        self.aspect2 = '{:.2g}'.format(self.s2.L_n)
        self.legend_name = (r'$\tilde{{L}}={{{0}}}$'.format(self.aspect1), 
                             r'$\tilde{{L}}={{{0}}}$'.format(self.aspect2))
    def amplitude(self):
        cx1, cy1 = self.s1.plt_cross_surf('y', 0, ret=True, mirror=False)
        cx2, cy2 = self.s2.plt_cross_surf('y', 0, ret=True, mirror=False)
        if cy1[0] < 0:
            cy1 = -cy1
        if cy2[0] < 0:
            cy2 = -cy2
        plt_figure_mpl([self.s1.max_x/self.s1.w - cx1, self.s2.max_x/self.s2.w - cx2], [cy1, cy2], 
                       name=self.legend_name, xlabel='$x$', ylabel='$z$'
                      , title='Compare Cross section at $y=0$')
        self.center_amp1 = np.array((cx1, cy1))
        self.center_amp2 = np.array((cx2, cy2))
        
    def delta(self):
        x1, dt1 = self.s1.plt_delta_x(flip_x=True, plt=False, ret=True)
        x2, dt2 = self.s2.plt_delta_x(flip_x=True, plt=False, ret=True)
        plt_figure_mpl([x1, x2], [dt1, dt2], name=self.legend_name
                      , xlabel='$x$', ylabel=r'$\tilde{\Delta}$'
                      , title=r'Compare $\tilde{\Delta}$')
        self.delta1 = np.array((x1, dt1))
        self.delta2 = np.array((x2, dt2))
        
    def boundary(self):
        x1, y1 = self.s1.plt_ybound(flip_x=True, plt=False)
        x2, y2 = self.s2.plt_ybound(flip_x=True, plt=False)
        plt_figure_mpl([x1, x2], [y1, y2], name=self.legend_name
                      , xlabel='$x$', ylabel=r'$y$'
                      , title=r'Compare boundary')
        self.bound1 = np.array((x1, y1))
        self.bound2 = np.array((x2, y2))
        
    def stress(self, whichstress='yy'):
        x1, sy1 = self.s1.plt_itg_stress(whichstress, plt=False, ret=True, flip_x=True)
        x2, sy2 = self.s2.plt_itg_stress(whichstress, plt=False, ret=True, flip_x=True)
        plt_figure_mpl([x1, x2], [sy1, sy2]
                      , xlabel='$x$', ylabel=r'$\sigma_{{{0}}}$'.format(whichstress)
                      , title=r'Compare stress')
        self.stress1 = np.array((x1, sy1))
        self.stress2 = np.array((x2, sy2))