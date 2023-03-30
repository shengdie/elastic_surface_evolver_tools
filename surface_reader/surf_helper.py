import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from styles import rev_slim


# class RevPlots(object):
#     # fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
#     # inches_per_pt = 1.0 / 72.27               # Convert pt to inches
#     # golden_mean = (np.sqrt(5) - 1.0) / 2.0         # Aesthetic ratio
#     # fig_width = fig_width_pt * inches_per_pt  # width in inches
#     # fig_height = fig_width * golden_mean       # height in inches
#     # aesthetic_size = np.array([fig_width, fig_height])
#     @staticmethod

def subplots(*args, **kw):
    mp.rcParams.update(rev_slim)
    return plt.subplots(*args, **kw)

def plt_sharex():
    pass

def set_figlabel(axes, labels, pos=[-0.1, 1.1], **kw):
    # prop = {'fontsize': size} if size is not None else {}
    for ax, l in zip(axes, labels):
        ax.text(pos[0], pos[1], l, transform=ax.transAxes, va='top', ha='right', **kw)

def savefig(fig, file, **kw):
    fig.tight_layout()
    fig.savefig(file, bbox_inches='tight', **kw)
# class RevPlots(object):
#     @staticmethod
#     def plt_rev_share_x(data1, data2, alllw=0.3, lw=0.3, borderlw=0.3, xticklw=0.3, yticklw=0.3, names=None,
#                         mathfont='cm', fontfamily='serif', usetex=True, axislabelsize=10, ticklabelsize=8, figlabelsize=8,
#                         cycler=['r', 'b', 'c'], figlabel=['(a)', '(b)'], figlabelpos=[-0.1, 1.15], figname=None):
#         fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
#         inches_per_pt = 1.0 / 72.27               # Convert pt to inches
#         golden_mean = (np.sqrt(5) - 1.0) / 2.0         # Aesthetic ratio
#         fig_width = fig_width_pt * inches_per_pt  # width in inches
#         fig_height =fig_width * golden_mean       # height in inches
#         aesthetic_size = np.array([fig_width, fig_height])
        
#         fig, (ax1, ax2) = plt.subplots(2,1, figsize=[fig_size[0], fig_size[1]*2.2], sharex=True)

#         if hasattr(data1[0][0], '__len__'): 
#             for x, y in zip(data1[0], data1[1]):
#                 ax1.plot(x, y)
#         if len(data1) == 2:
#             x1 = data1[0]
#             y1 = data1[1]
#         else:
#             x1 = data1[:, 0]
#             y1 = data1[:, 1]
#         ax1.plot(x1, y1)


def _rect_inter_inner(x1,x2):
    n1=x1.shape[0]-1
    n2=x2.shape[0]-1
    X1=np.c_[x1[:-1],x1[1:]]
    X2=np.c_[x2[:-1],x2[1:]]
    S1=np.tile(X1.min(axis=1),(n2,1)).T
    S2=np.tile(X2.max(axis=1),(n1,1))
    S3=np.tile(X1.max(axis=1),(n2,1)).T
    S4=np.tile(X2.min(axis=1),(n1,1))
    return S1,S2,S3,S4

def _rectangle_intersection_(x1,y1,x2,y2):
    S1,S2,S3,S4=_rect_inter_inner(x1,x2)
    S5,S6,S7,S8=_rect_inter_inner(y1,y2)

    C1=np.less_equal(S1,S2)
    C2=np.greater_equal(S3,S4)
    C3=np.less_equal(S5,S6)
    C4=np.greater_equal(S7,S8)

    ii,jj=np.nonzero(C1 & C2 & C3 & C4)
    return ii,jj

def intersection(x1,y1,x2,y2):
    """
INTERSECTIONS Intersections of curves.
   Computes the (x,y) locations where two curves intersect.  The curves
   can be broken with NaNs or have vertical segments.
usage:
x,y=intersection(x1,y1,x2,y2)
    Example:
    a, b = 1, 2
    phi = np.linspace(3, 10, 100)
    x1 = a*phi - b*np.sin(phi)
    y1 = a - b*np.cos(phi)
    x2=phi
    y2=np.sin(phi)+2
    x,y=intersection(x1,y1,x2,y2)
    plt.plot(x1,y1,c='r')
    plt.plot(x2,y2,c='g')
    plt.plot(x,y,'*k')
    plt.show()
    """
    ii,jj=_rectangle_intersection_(x1,y1,x2,y2)
    n=len(ii)

    dxy1=np.diff(np.c_[x1,y1],axis=0)
    dxy2=np.diff(np.c_[x2,y2],axis=0)

    T=np.zeros((4,n))
    AA=np.zeros((4,4,n))
    AA[0:2,2,:]=-1
    AA[2:4,3,:]=-1
    AA[0::2,0,:]=dxy1[ii,:].T
    AA[1::2,1,:]=dxy2[jj,:].T

    BB=np.zeros((4,n))
    BB[0,:]=-x1[ii].ravel()
    BB[1,:]=-x2[jj].ravel()
    BB[2,:]=-y1[ii].ravel()
    BB[3,:]=-y2[jj].ravel()

    for i in range(n):
        try:
            T[:,i]=np.linalg.solve(AA[:,:,i],BB[:,i])
        except:
            T[:,i]=np.NaN


    in_range= (T[0,:] >=0) & (T[1,:] >=0) & (T[0,:] <=1) & (T[1,:] <=1)

    xy0=T[2:,in_range]
    xy0=xy0.T
    return xy0[:,0],xy0[:,1]


# if __name__ == '__main__':

#     # a piece of a prolate cycloid, and am going to find
#     a, b = 1, 2
#     phi = np.linspace(3, 10, 100)
#     x1 = a*phi - b*np.sin(phi)
#     y1 = a - b*np.cos(phi)

#     x2=phi
#     y2=np.sin(phi)+2
#     x,y=intersection(x1,y1,x2,y2)
#     plt.plot(x1,y1,c='r')
#     plt.plot(x2,y2,c='g')
#     plt.plot(x,y,'*k')
#     plt.show()