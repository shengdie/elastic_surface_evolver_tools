import sys
import os
import glob
import numpy as np
import h5py
import plotly.offline as ply
import plotly.graph_objs as go
import plotly.figure_factory as FF
from collections import Sequence
from scipy.interpolate import interp1d, griddata
from scipy.optimize import fmin, root, brentq, curve_fit
from scipy.integrate import trapz
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from functools import reduce

ply.init_notebook_mode(connected=False)

def import_hdf5(filename):
    """import hdf5 formate file which is exported by mathematica.
    """
    f = h5py.File(filename, 'r')
    a_group_key = list(f.keys())[0]
    data = list(f[a_group_key])
    return np.array(data)

def get_griddata(x,y,z, numx =500, numy=500, k='cubic', retgridxy=False, **kw):  #using griddata method to interpolate
    X = np.linspace(np.min(x),np.max(x), numx)
    Y = np.linspace(np.min(y), np.max(y), numy)
    grid_x, grid_y = np.meshgrid(X,Y)
    verts = np.transpose([x,y])
    Z = griddata(verts, z, (grid_x, grid_y), method=k, **kw)
    if not retgridxy:
        return (X,Y, Z) 
    else:
        return (grid_x, grid_y, Z)

def get_mathematica_data_for_contour(data, mirror=True, a=1, unstructured=False):
    """give the data from mathematica, return x,y,z for contour plot with matplotlib,
    for example, in Mathematica, for a function f[x,y], export the data with
    Table[f[x,y], {x,xi,xf}, {y,yi,yf}]
    
    Arguments:
        data {array} -- 3d array shape of (n,m,3), exported from mathematica, with
        Table[f[x,y], {x,xi,xf}, {y,yi,yf}], if unstructured, then it should be [n,3]
    
    Keyword Arguments:
        mirror {bool} -- if mirror the data, with axis of x=0 or y=0 (default: {True})
        a {int} -- a=0, mirror on x=0 (y axis), a=1 mirror on y=0 (y axis) (default: {1})
    
    Returns:
        (x,y,z) -- coordinates for contour plot with matplotlib
    """
    if unstructured:
        data = np.array(data)
        x, y, z = get_griddata(data[:,0], data[:,1], data[:,2])
        return x,y,z

    if mirror:
        sx, sy, _ = np.shape(data)
        if a == 1:
            data = np.array([[[data[x,y,0],-data[x,y,1],data[x,y,2]] for y in range(sy-1,0,-1)] +
                            [[data[x,y,0],data[x,y,1],data[x,y,2]] for y in range(sy)]
                            for x in range(sx)])
        elif a == 0:
            data = np.array([[[-data[x,y,0],data[x,y,1],data[x,y,2]] for x in range(sx-1,0,-1)] +
                            [[data[x,y,0],data[x,y,1],data[x,y,2]] for x in range(sx)]
                            for y in range(sy)])
    x = data[:,:,0][:,0]
    y = data[:,:,1][0]
    z = np.transpose(data[:,:,2])
    return (x,y,z)

def map_z2color(zval, colormap, vmin, vmax):
    #map the normalized value zval to a corresponding color in the colormap
    colormap = cm.get_cmap(colormap)
    if vmin>vmax:
        raise ValueError('incorrect relation between vmin and vmax')
    t=(zval-vmin)/float((vmax-vmin))#normalize val
    R, G, B, alpha=colormap(t)
    return 'rgb('+'{:d}'.format(int(R*255+0.5))+','+'{:d}'.format(int(G*255+0.5))+\
           ','+'{:d}'.format(int(B*255+0.5))+')'

def tri_indices(simplices):
    #simplices is a numpy array defining the simplices of the triangularization
    #returns the lists of indices i, j, k

    return ([triplet[c] for triplet in simplices] for c in range(3))

def plotly_trisurf(x, y, z, simplices, colorscale='Rainbow', plot_edges=False):
    #x, y, z are lists of coordinates of the triangle vertices 
    #simplices are the simplices that define the triangularization;
    #simplices  is a numpy array of shape (no_triangles, 3)
    #insert here the  type check for input data

    points3D=np.vstack((x,y,z)).T
    
    I,J,K=simplices[:,0], simplices[:,1], simplices[:,2]

    triangles=go.Mesh3d(x=x,
                     y=y,
                     z=z,
                     #facecolor=facecolor,
                     i=I,
                     j=J,
                     k=K,
                     name='',
                     intensity=z,
                     colorbar = go.ColorBar(
                        title='z'
                     ),
                     colorscale=colorscale,
                     showscale=True
                    )

    if not plot_edges:# the triangle sides are not plotted 
        return [triangles]
    else:
        #define the lists Xe, Ye, Ze, of x, y, resp z coordinates of edge end points for each triangle
        #None separates data corresponding to two consecutive triangles
        eee = np.concatenate((points3D[simplices], [[[None,None,None]]] * simplices.shape[0]),axis=1)
        Xe, Ye, Ze=eee[:,:,0].flatten(), eee[:,:,1].flatten(), eee[:,:,2].flatten()

        #define the lines to be plotted
        lines=go.Scatter3d(x=Xe,
                        y=Ye,
                        z=Ze,
                        mode='lines',
                        line=dict(color= 'rgb(50,50,50)', width=1.5)
               )
        return [triangles, lines]

def plt_mesh(xyz, ijk, mirror=None, colorscale='Rainbow', slicex=None, slicey=None, slicez=None,
             ar=None, onlysurf=False, plot_edges=False, camera=None, **kw):
    if not plot_edges: 
        data = [go.Mesh3d(
            x=xyz[:,0],
            y=xyz[:,1],
            z=xyz[:,2],
            i=ijk[:,0],
            j=ijk[:,1],
            k=ijk[:,2],
            intensity=xyz[:,2],
            colorbar = go.ColorBar(
                title='z'
            ),
            colorscale=colorscale,
            showscale=True
        )]
        if mirror is not None:
            if 'x' in mirror or 'all' in mirror:
                data.append(
                    go.Mesh3d(
                    x=-xyz[:,0],
                    y=xyz[:,1],
                    z=xyz[:,2],
                    i=ijk[:,0],
                    j=ijk[:,1],
                    k=ijk[:,2],
                    intensity=xyz[:,2],
                    colorscale=colorscale,
                    showscale=True
                ))
            if 'y' in mirror or 'all' in mirror:
                data.append(
                    go.Mesh3d(
                    x=xyz[:,0],
                    y=-xyz[:,1],
                    z=xyz[:,2],
                    i=ijk[:,0],
                    j=ijk[:,1],
                    k=ijk[:,2],
                    intensity=xyz[:,2],
                    colorscale=colorscale,
                    showscale=True
                ))
            if 'all' in mirror:
                data.append(
                    go.Mesh3d(
                    x=-xyz[:,0],
                    y=-xyz[:,1],
                    z=xyz[:,2],
                    i=ijk[:,0],
                    j=ijk[:,1],
                    k=ijk[:,2],
                    intensity=xyz[:,2],
                    colorscale=colorscale,
                    showscale=True
                ))
    else:
        data= plotly_trisurf(xyz[:,0], xyz[:,1], xyz[:,2], ijk, colorscale=colorscale, plot_edges=True)
    if ar is None:
        layoutprop = {}
    else:
        layoutprop = {'scene':{'aspectratio':{'x':ar[0],'y':ar[1],'z':ar[2]}}}
    layout=go.Layout(
        xaxis=go.XAxis(
            title='x'
        ),
        yaxis=go.YAxis(
            title='y'
        ),
        **layoutprop
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
    #if slicex is not None:
    fig = go.Figure(data=data, layout=layout)
    if camera is not None:
        camera = dict(eye=dict(x=camera[0], y=camera[1], z=camera[2]))
        fig.update_layout(scene_camera=camera)
    ply.iplot(fig)


def plt_figure_mpl(x,y, xlim=None, ylim=None, xticks=None, yticks=None,
                xticklabels=None, yticklabels=None, name=None, 
                xlabel='', ylabel='', title='', legend=True,
                legendloc='upper right', yrotate=0,
                linestyle='solid', marker='None', color=None,
                aspect=None, 
                text=None, textsize=15, textprops={}, texttransAxes=False, notransAxesforlog=False,
                figsize=None, inset=None, insetlim=None, insetxticks=None, insetyticks=None,
                insetgrid=False,
                titlefont=17, axisfont=16, tickfont=11, legendfont=15,
                figname=None, grid=True, style='seaborn-white',
                minor_ticks=True, hideminorticklabels=False,
                vline=None, vlineprops={},
                mathfont='cm', dpi=200, logscale=False, ret=False, **kw):
    """Plot figure

    Arguments:
        x {array} -- x coords, can be [coords] for single curve, or [coords_c1, coords_c2, ...] for multi curve
        y {array} -- y coordinates corresponding to x

    Keyword Arguments:
        name {array like} -- names of multi curves, useful for multi curves, 
                    ['name_c1', 'name_c2', ...] (default: {None})
        xlabel {str} -- x axis label (default: {''})
        ylabel {str} -- y axis label (default: {''})
        xticks {array} -- xticks positions (default: {None})
        yticks {array} -- yticks positions (default: {None})
        xticklabels {array of str} -- xtick labels, corresponding to xticks (default: {None})
        yticklabels {array of str} -- ytick labels, corresponding to yticks (default: {None})
        title {str} -- title of figure (default: {''})
        legendloc {str or tuple of float} -- position of lengend, 
                    like (1.05,1) will place at upper right out of frame. (default: {'upper right'})
        linestyle {str} -- line style of plot, 'None' for no line (scatter plot). (default: {'solid'})
        marker {str} -- marker style of plot. (default: {'None'})
        inset {array like} -- [x,y, width, height], no if None (default: {None})
        insetlim {array like} -- [xi, xf, yi, yf] of inset. (default: {None})
        text {array} -- text on plot, can be [x,y, 'text'] or [x, y, 'text', fontsize] multi text, 
                    [[x1,y1,'text1'], [x2, y2, 'text2', fontsize], ...] (default: {None})
        textsize {int} -- text font size (default: {15})
        figsize {array} --  (width, height) (default: None)
        figname {str} -- file name for saving the figure, if None, do not save (default: {None})
        grid {bool} -- use grid on plot (default: {True})
        style {str} -- plot style of matplotlib, see 'print(plt.style.available)' (default: {'ggplot'})
        minor_ticks {bool} -- if show minor ticks (mainly for log scale plot) (default: {True})
        hideminorticklabels {bool} -- if hide minor ticklabels (useful for log scale plot) (default: {False})
        vline {float or array like} -- plot vertical lines at x = vline. (default: {None})
        font {str} -- math font (default: {'stix'})
        dpi {int} -- dpi for showing the figure (default: {200})
        **kw -- other kargs for matplotlib.polt or scatter.
    """
    def set_font(ft, defont):
        #if ft is None:
        #    return defont
        if isinstance(ft, (int, float)):
            defont['size'] = ft
            return defont
        elif isinstance(ft, str):
            defont['family'] = ft
            return defont
        else:
            return ft
    # if title_font is None:
    title_font= {'family':'DejaVu Sans', 'size':'20'}
    axis_font = {'family':'DejaVu Sans', 'size':'27'}
    tick_font = {'family':'DejaVu Sans', 'size':'15'}
    legend_font = {'family':'DejaVu Sans', 'size':'17'}
    titlefont = set_font(titlefont, title_font)
    axisfont = set_font(axisfont, axis_font)
    tickfont = set_font(tickfont, tick_font)
    legendfont = set_font(legendfont, legend_font)
    
    if len(np.shape(x)) == 1:
        if len(np.shape(x[0])) == 0:
            x = [x]
            y = [y]
    if len(np.shape(text)) == 1:
        if len(np.shape(text[0])) == 0:
            text = [text]
    #if len(np.shape(scatter)) == 0:
    #    scatter = [scatter] * len(x)
    if len(np.shape(linestyle)) == 0:
        linestyle = [linestyle] * len(x)
    if len(np.shape(marker)) == 0:
        marker = [marker] * len(x)
    if color is None:
        color = ['C' + str(i%10) for i in range(len(x))]
    elif len(np.shape(color)) ==0:
        color = [color]* len(x)
    if name is None:
        name = ['line' + str(i+1) for i in range(len(x))]

    if style is not None:
        plt.style.use(style)
    else:
        mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['mathtext.fontset'] = mathfont
    
    fig, ax = plt.subplots(figsize=figsize)
    
    if len(name) == len(x):
        for i in range(len(x)):
            if not logscale:
                ax.plot(x[i], y[i], label=name[i], linestyle=linestyle[i], marker=marker[i], color=color[i], **kw)
            else:
                ax.loglog(x[i], y[i], label=name[i], linestyle=linestyle[i], marker=marker[i], color=color[i], **kw)
    else:
        print('error: name should have same length as x and y')
        return 0
    if text is not None:
        if not hasattr(notransAxesforlog, '__len__'):
            notransAxesforlog = [notransAxesforlog]*len(text)
        for i in range(len(text)):
            if (logscale and not notransAxesforlog[i]) or texttransAxes:
                #textprops['transform']=ax.transAxes
                transform={'transform':ax.transAxes}
            else:
                transform = {}
            tt = text[i]
            if len(tt) == 3:
                ax.text(tt[0], tt[1], tt[2], fontsize=textsize, **transform, **textprops)
            elif len(tt) == 4:
                ax.text(tt[0], tt[1], tt[2], fontsize=tt[3], **transform, **textprops)
    if legend and len(x) > 1 and name.count(None) != len(name):
        legendfontprop = {} if legendfont is None else {'prop':font_manager.FontProperties(**legendfont)}
        ax.legend(loc=legendloc, **legendfontprop)

    #ax.plot(x,y)
    if not minor_ticks:
        ax.minorticks_off()
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
    if xticks is not None:
        ax.set_xticks(xticks, minor=False)
    if yticks is not None:
        ax.set_yticks(yticks, minor=False)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels, fontdict=None, minor=False)
    if yticklabels is not None:
        ax.set_yticklabels(yticklabels, fontdict=None, minor=False)
    if tickfont is not None:
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontproperties(font_manager.FontProperties(**tickfont))
    axisfontprop = {} if axisfont is None else {'fontproperties':font_manager.FontProperties(**axisfont)}
    ax.set_xlabel(xlabel, **axisfontprop)
    ax.set_ylabel(ylabel, rotation=yrotate, **axisfontprop)
    titlefontprop = {} if titlefont is None else {'fontproperties':font_manager.FontProperties(**titlefont)}
    ax.set_title(title, **titlefontprop)
    if aspect is not None:
        if len(np.shape(aspect)) == 1:
            ylmin, ylmax = ax.get_ylim()
            xlmin, xlmax = ax.get_xlim()
            if logscale:
                ylmin = np.log10(ylmin)
                ylmax = np.log10(ylmax)
                xlmin = np.log10(xlmin)
                xlmax = np.log10(xlmax)
            asp = (xlmax-xlmin)/(ylmax-ylmin) * aspect[1]/aspect[0]
            ax.set_aspect(asp)
        else:
            ax.set_aspect(aspect)
    ax.grid(grid)

    # plot vertical lines
    if vline is not None:
        if hasattr(vline, '__len__'):
            for xps in vline:
                ax.axvline(xps, **vlineprops)

    if (inset is not None) and len(inset) == 4:
        axi = ax.inset_axes(inset)
        for i in range(len(x)):
            if not logscale:
                axi.plot(x[i], y[i], label=name[i], linestyle=linestyle[i], marker=marker[i], color=color[i], **kw)
            else:
                axi.loglog(x[i], y[i], label=name[i], linestyle=linestyle[i], marker=marker[i], color=color[i], **kw)
        if insetlim is not None:
            axi.set_xlim(insetlim[0], insetlim[1])
            axi.set_ylim(insetlim[2], insetlim[3])
        if insetxticks is not None:
            axi.set_xticks(insetxticks)
        if insetyticks is not None:
            axi.set_yticks(insetyticks)
        axi.grid(insetgrid)
        ax.indicate_inset_zoom(axi)
        
    fig.set_dpi(dpi)
    fig.tight_layout()
    if hideminorticklabels:
        plt.setp(ax.get_xminorticklabels(), visible=False)
        plt.setp(ax.get_yminorticklabels(), visible=False)
    plt.show()
    if figname is not None:
        fig.savefig(figname, format='pdf', bbox_inches='tight')
    if ret:
        return fig, ax

def plt_any_mpl(ax, fig, xlim=None, ylim=None, xticks=None, yticks=None,
                xticklabels=None, yticklabels=None, 
                xlabel='', ylabel='', title='', uselegend=False,
                legendloc='upper right', 
                aspect=None,
                text=None, textsize=15, textprops={}, texttransAxes=False, notransAxesforlog=False,
                figsize=None,
                titlefont=None, axisfont=None, tickfont=None, legendfont=None,
                figname=None, grid=True, style='seaborn-white',
                minor_ticks=True, hideminorticklabels=False,
                vline=None, vlineprops={'c':'r', 'ls':'--', 'lw':1},
                mathfont='cm', dpi=200, logscale=False, **kw):
    """Plot any figure

    Arguments:
        ax {object} -- axis object of matplot
        fig {object} -- figure object of matplot

    Keyword Arguments:
        xlabel {str} -- x axis label (default: {''})
        ylabel {str} -- y axis label (default: {''})
        xticks {array} -- xticks positions (default: {None})
        yticks {array} -- yticks positions (default: {None})
        xticklabels {array of str} -- xtick labels, corresponding to xticks (default: {None})
        yticklabels {array of str} -- ytick labels, corresponding to yticks (default: {None})
        title {str} -- title of figure (default: {''})
        legendloc {str or tuple of float} -- position of lengend, 
                    like (1.05,1) will place at upper right out of frame. (default: {'upper right'})
        text {array} -- text on plot, can be [x,y, 'text'] or [x, y, 'text', fontsize] multi text, 
                    [[x1,y1,'text1'], [x2, y2, 'text2', fontsize], ...] (default: {None})
        textsize {int} -- text font size (default: {15})
        figname {str} -- file name for saving the figure, if None, do not save (default: {None})
        grid {bool} -- use grid on plot (default: {True})
        style {str} -- plot style of matplotlib, see 'print(plt.style.available)' (default: {'ggplot'})
        minor_ticks {bool} -- if show minor ticks (mainly for log scale plot) (default: {True})
        hideminorticklabels {bool} -- if hide minor ticklabels (useful for log scale plot) (default: {False})
        vline {float or array like} -- plot vertical lines at x = vline. (default: {None})
        font {str} -- math font (default: {'stix'})
        dpi {int} -- dpi for showing the figure (default: {200})
        **kw -- other kargs for matplotlib.polt or scatter.
    """
    def set_font(ft, defont):
        #if ft is None:
        #    return defont
        if isinstance(ft, (int, float)):
            defont['size'] = ft
            return defont
        elif isinstance(ft, str):
            defont['family'] = ft
            return defont
        else:
            return ft
    title_font= {'family':'DejaVu Sans', 'size':'16'}
    axis_font = {'family':'DejaVu Sans', 'size':'14'}
    tick_font = {'family':'DejaVu Sans', 'size':'13'}
    legend_font = {'family':'DejaVu Sans', 'size':'14'}
    titlefont = set_font(titlefont, title_font)
    axisfont = set_font(axisfont, axis_font)
    tickfont = set_font(tickfont, tick_font)
    legendfont = set_font(legendfont, legend_font)
    
    if style is not None:
        plt.style.use(style)
    else:
        mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams['mathtext.fontset'] = mathfont

    if text is not None:
        if not hasattr(notransAxesforlog, '__len__'):
            notransAxesforlog = [notransAxesforlog]*len(text)
        for i in range(len(text)):
            if (logscale and not notransAxesforlog[i]) or texttransAxes:
                #textprops['transform']=ax.transAxes
                transform={'transform':ax.transAxes}
            else:
                transform = {}
            tt = text[i]
            if len(tt) == 3:
                ax.text(tt[0], tt[1], tt[2], fontsize=textsize, **transform, **textprops)
            elif len(tt) == 4:
                ax.text(tt[0], tt[1], tt[2], fontsize=tt[3], **transform, **textprops)
    if uselegend:
        legendfontprop = {} if legendfont is None else {'prop':font_manager.FontProperties(**legendfont)}
        ax.legend(loc=legendloc, **legendfontprop)
    #ax.plot(x,y)
    if not minor_ticks:
        ax.minorticks_off()
    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])
    if xticks is not None:
        ax.set_xticks(xticks, minor=False)
    if yticks is not None:
        ax.set_yticks(yticks, minor=False)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels, fontdict=None, minor=False)
    if yticklabels is not None:
        ax.set_yticklabels(yticklabels, fontdict=None, minor=False)
    if tickfont is not None:
        for label in (ax.get_xticklabels() + ax.get_yticklabels()):
            label.set_fontproperties(font_manager.FontProperties(**tickfont))
    axisfontprop = {} if axisfont is None else {'fontproperties':font_manager.FontProperties(**axisfont)}
    ax.set_xlabel(xlabel, **axisfontprop)
    ax.set_ylabel(ylabel, **axisfontprop)
    #ax.set_xlabel(xlabel, fontproperties=axisfontprop)
    #ax.set_ylabel(ylabel, fontproperties=axisfontprop)
    titlefontprop = {} if titlefont is None else {'fontproperties':font_manager.FontProperties(**titlefont)}
    ax.set_title(title, **titlefontprop)
    if aspect is not None:
        if len(np.shape(aspect)) == 1:
            ylmin, ylmax = ax.get_ylim()
            xlmin, xlmax = ax.get_xlim()
            if logscale:
                ylmin = np.log10(ylmin)
                ylmax = np.log10(ylmax)
                xlmin = np.log10(xlmin)
                xlmax = np.log10(xlmax)
            asp = (xlmax-xlmin)/(ylmax-ylmin) * aspect[1]/aspect[0]
            ax.set_aspect(asp)
        else:
            ax.set_aspect(aspect)
    ax.grid(grid)

    # plot vertical lines
    if vline is not None:
        if hasattr(vline, '__len__'):
            for xps in vline:
                ax.axvline(xps, **vlineprops)
    fig.set_dpi(dpi)
    fig.tight_layout()
    if hideminorticklabels:
        plt.setp(ax.get_xminorticklabels(), visible=False)
        plt.setp(ax.get_yminorticklabels(), visible=False)
    #plt.show()
    if figname is not None:
        fig.savefig(figname, format='pdf', bbox_inches='tight')

def plt_curv_fit_mpl(x, y, title='', xlabel='$x$', ylabel='$y$', 
    showlegend=True, name=None,
    showslope=True, triscale=1/7, tripos=[3/4, 1/4],
    textpos_adjust=[-0.3,0.5], slopetextsize=12,
    logfit=True, logscale=True, fitdataindex=None, text=None,
    showfitteddata=True,
    ret=False, retfit=False, retall=False, **kw):
    """Plot the line, and fit the line
    
    Arguments:
        x {array} -- x coordinates
        y {array} -- y coordinates
    
    Keyword Arguments:
        title {str} -- title of plot (default: {''})
        xlabel {str} -- xaxis label (default: {'$x$'})
        ylabel {str} -- yaxis label (default: {'$y$'})
        showslope {bool} -- show slope(draw a triangle, and label slope) (default: {True})
        triscale {float} -- scale of slope triangle, length= triscale * x_range (default: {1/8})
        tripos {list} -- triangle position, in x range and y range (default: {[3/4, 1/4]})
        textpos_adjust {list} -- adjust the position slope text, in the length and height of triangle (default: {[-0.1,0.1]})
        slopetextsize {int} -- slope text size (default: {10})
        logfit {bool} -- if fit log(x) and log(y) instead x, y (default: {False})
        logscale {bool} -- if use log scale for plot (default: {False})
        fitdataindex {arrary or iterable} -- fit the data of indices specified in this, None for all (default: {None})
        ret {bool} -- if return the fitted function (default: {False})
    
    Returns:
        function -- fitted function
    """
    x= np.asarray(x)
    y = np.asarray(y)
    if fitdataindex is None:
        fitdataindex = range(len(x))
    
    xf = x[fitdataindex]
    
    if logfit:
        lx = np.log(xf)
        ly = np.log(y[fitdataindex])
        p1 = np.polyfit(lx, ly, 1)
        f = np.poly1d(p1)
        lyf = f(lx)
        yf = np.exp(lyf)
    else:
        lx = xf
        ly = y[fitdataindex]
        p1 = np.polyfit(lx,ly,1)
        f = np.poly1d(p1)
        yf = f(lx)
    
    tri_x = np.array([0,1,1,0])
    tri_y = np.array([0,0,p1[0],0])
    xr = lx[-1] - lx[0]
    yr = ly[-1] - ly[0]
    tri_x = tri_x * triscale * xr + lx[0] + tripos[0] * xr
    tri_y = tri_y * triscale * xr + ly[0] + tripos[1] * yr
    #tri_y[2] = np.exp(f(np.log(tri_x[2])))
    if logfit:
        textslope = [
            np.exp(
                (tri_x[0]+tri_x[2])/2 + textpos_adjust[0]* (tri_x[1] - tri_x[0])), 
                np.exp(
                (tri_y[0]+ tri_y[2])/2 + textpos_adjust[1]* (tri_y[2] - tri_y[1])),
                r'${0:+g}$'.format(np.around(p1[0],3)), slopetextsize]
    else:
        textslope = [
                (tri_x[0]+tri_x[2])/2 + textpos_adjust[0]* (tri_x[1] - tri_x[0]),
                (tri_y[0]+ tri_y[2])/2 + textpos_adjust[1]* (tri_y[2] - tri_y[1]),
                r'${0:+g}$'.format(np.around(p1[0],3)), slopetextsize]
    if text is not None:
        if len(np.shape(text)) == 1:
            if len(np.shape(text[0])) == 0:
                text = [text]
        text.append(textslope)
    else:
        text = textslope
    tri_x = np.exp(tri_x)
    tri_y = np.exp(tri_y)
    #print(textslope)
    
    #name=('data', 'fit: {0}{2}+{1}'.format(np.around(p1[0],3), np.around(p1[1],3), xlabel)), 
    if showfitteddata:
        if name is None:
            if showlegend:
                name=['data', 'fitting', None]
            else:
                name=[None, None, None]
        datax = [x, xf, tri_x]
        datay = [y,yf, tri_y]
        linestyle=['None', 'solid', 'solid']
        marker=['.', 'None', 'None']
        color = ['C0','C1', 'C1']
    else:
        if name is None:
            name=[None, None]
        datax = [x, tri_x]
        datay = [y, tri_y]
        linestyle=['None', 'solid']
        marker=['.', 'None']
        color=['C0', 'C0']
    
    #notransAxesforlog = [False]*len(text)
    notransAxesforlog=True
    
    pltret =  plt_figure_mpl(datax, datay, linestyle=linestyle, marker=marker,
    name=name,
    title=title, xlabel=xlabel, ylabel=ylabel,
    text= text, notransAxesforlog=notransAxesforlog,
    #notransAxesforlog=[False,False,False,True],
    #textsize= slopetextsize,
    logscale=logscale,
    color=color,
    ret=ret,
    **kw)
    if retall:
        return datax, datay, text, f
    elif ret and not retfit:
        return pltret
    elif (not ret) and retfit:
        print(p1)
        return f
    elif ret and retfit:
        return pltret[0], pltret[1], f


def sample_function(func, points, tol=0.01, min_points=16, max_level=16,
                    sample_transform=None, funckw=None):
    """
    Sample a 1D function to given tolerance by adaptive subdivision.

    The result of sampling is a set of points that, if plotted,
    produces a smooth curve with also sharp features of the function
    resolved.

    Parameters
    ----------
    func : callable
        Function func(x) of a single argument. It is assumed to be vectorized.
    points : array-like, 1D
        Initial points to sample, sorted in ascending order.
        These will determine also the bounds of sampling.
    tol : float, optional
        Tolerance to sample to. The condition is roughly that the total
        length of the curve on the (x, y) plane is computed up to this
        tolerance.
    min_point : int, optional
        Minimum number of points to sample.
    max_level : int, optional
        Maximum subdivision depth.
    sample_transform : callable, optional
        Function w = g(x, y). The x-samples are generated so that w
        is sampled.

    Returns
    -------
    x : ndarray
        X-coordinates
    y : ndarray
        Corresponding values of func(x)

    Notes
    -----
    This routine is useful in computing functions that are expensive
    to compute, and have sharp features --- it makes more sense to
    adaptively dedicate more sampling points for the sharp features
    than the smooth parts.

    Examples
    --------
    >>> def func(x):
    ...     '''Function with a sharp peak on a smooth background'''
    ...     a = 0.001
    ...     return x + a**2/(a**2 + x**2)
    ...
    >>> x, y = sample_function(func, [-1, 1], tol=1e-3)

    >>> import matplotlib.pyplot as plt
    >>> xx = np.linspace(-1, 1, 12000)
    >>> plt.plot(xx, func(xx), '-', x, y[0], '.')
    >>> plt.show()

    """
    x, y = _sample_function(func, points, values=None, mask=None, depth=0,
                            tol=tol, min_points=min_points, max_level=max_level,
                            sample_transform=sample_transform, funckw=funckw)
    return x.flatten(), y.flatten()

def _sample_function(func, points, values=None, mask=None, tol=0.05,
                     depth=0, min_points=16, max_level=16,
                     sample_transform=None, funckw=None):
    points = np.array(points)
    if funckw is not None:
        func = lambda x: func(x, **funckw)

    if values is None:
        #values = np.atleast_2d(func(points))
        values = np.atleast_2d(list(recursive_map(func, points)))

    #if mask is None:
    #    mask = Ellipsis

    if depth > max_level:
        # recursion limit
        return points, values

    x_a = points[...,:-1]
    x_b = points[...,1:]
    if mask is not None:
        x_a = x_a[..., mask]
        x_b = x_b[..., mask]
        
    x_c = .5*(x_a + x_b)
    y_c = np.atleast_2d(list(recursive_map(func, x_c)))

    x_2 = np.r_[points, x_c]
    #print(values)
    #print(y_c)
    y_2 = np.r_['-1', values, y_c]
    j = np.argsort(x_2)

    x_2 = x_2[...,j]
    y_2 = y_2[...,j]

    # -- Determine the intervals at which refinement is necessary

    if len(x_2) < min_points:
        mask = np.ones([len(x_2)-1], dtype=bool)
    else:
        # represent the data as a path in N dimensions (scaled to unit box)
        if sample_transform is not None:
            y_2_val = sample_transform(x_2, y_2)
        else:
            y_2_val = y_2

        p = np.r_['0',
                  x_2[None,:],
                  y_2_val.real.reshape(-1, y_2_val.shape[-1]),
                  y_2_val.imag.reshape(-1, y_2_val.shape[-1])
                  ]

        sz = (p.shape[0]-1)//2

        xscale = x_2.ptp(axis=-1)
        #print(y_2_val, y_2_val.shape)
        yscale = abs(y_2_val.ptp(axis=-1)).ravel()
        #print(yscale, 'd')

        p[0] /= xscale
        p[1:sz+1] /= yscale[:,None]
        p[sz+1:]  /= yscale[:,None]

        # compute the length of each line segment in the path
        dp = np.diff(p, axis=-1)
        s = np.sqrt((dp**2).sum(axis=0))
        s_tot = s.sum()

        # compute the angle between consecutive line segments
        dp /= s
        dcos = np.arccos(np.clip((dp[:,1:] * dp[:,:-1]).sum(axis=0), -1, 1))

        # determine where to subdivide: the condition is roughly that
        # the total length of the path (in the scaled data) is computed
        # to accuracy `tol`
        dp_piece = dcos * .5*(s[1:] + s[:-1])
        mask = (dp_piece > tol * s_tot)

        mask = np.r_[mask, False]
        mask[1:] |= mask[:-1].copy()


    # -- Refine, if necessary

    if mask.any():
        return _sample_function(func, x_2, y_2, mask, tol=tol, depth=depth+1,
                                min_points=min_points, max_level=max_level,
                                sample_transform=sample_transform)
    else:
        return x_2, y_2
    
def recursive_map(func, seq):
    """recursive map sequence to func. return a generator. to list with list()
    Example:
    def f(x):
        if x>2:
            return x**2
        else:
            return x
    list(recursive_map(f, [1,2, [3,4]]))
    results: [1,2, [9,16]]
    """
    for item in seq:
        if isinstance(item, Sequence):
            yield type(item)(recursive_map(func, item))
        elif isinstance(item, np.ndarray):
            yield list(recursive_map(func, item))
        else:
            yield func(item)

def plt_func1d(func, r, funckw=None, tol=0.01, 
                min_points=16, max_level=16,
                sample_transform=None, ret=False, **kw):
    """plot a 1d function like func(x)
    Arguments:
        func {callable} -- a callable function
        r {list like} -- plot domain, [x_min, x_max] or a list of x coords.
    Keyword Arguments:
        funckw {dict} -- kargs for func. (default: None)
        tol {float} -- Tolerance to sample to. The condition is roughly that the total
            length of the curve on the (x, y) plane is computed up to this tolerance.
        min_point {int} -- Minimum number of points to sample. (default: 16)
        max_level {int} -- Maximum subdivision depth. (default: 16)
        sample_transform {callable} -- Function w = g(x, y). The x-samples are generated so 
            that w is sampled. (default: None)
        ret {bool} -- return x, y if True. (default: False)
        **kw -- kargs for plt_figure_mpl.
    """
    if hasattr(func, '__len__'):
        if not hasattr(r[0], '__len__'):
            r = [r] * len(func)
    else:
        func = [func]
        r = [r]
    xx, yy = [], []
    for i in range(len(func)):
        x, y = sample_function(func[i], r[i], tol=tol, min_points=min_points, max_level=max_level,
                            sample_transform=sample_transform, funckw=funckw)
        xx.append(x)
        yy.append(y)
    #x, y = sample_function(func, r, tol=tol, min_points=min_points, max_level=max_level,
    #                        sample_transform=sample_transform, funckw=funckw)
    plt_figure_mpl(xx, yy, **kw)
    if ret:
        return xx, yy

def plt_curve(x, y, name =None,title='', xaxis='', yaxis='', ar=(10,1), ret=False):
    """plot curve with plotly, can plot multi curves
    
    Arguments:
        x {array} -- x coords, can be [coords] for single curve, or [coords_c1, coords_c2, ...] for multi curve
        y {array} -- y coords, cooresponding to y
    
    Keyword Arguments:
        name {array like} -- names of multi curves, useful for multi curves, ['name_c1', 'name_c2', ...] (default: {None})
        title {str} -- title of plot (default: {''})
        xaxis {str} -- xaxis label (default: {''})
        yaxis {str} -- yaxis label (default: {''})
        ar {tuple} -- aspect ratio, [x,y] (default: {(10,1)})
        ret {bool} -- if Ture, return the plotly fig object (default: {False})
    """

    if len(np.shape(x)) == 1:
        if len(np.shape(x[0])) == 0:
            x = [x]
            y = [y]
            name = [name]
    if name is None:
        data=[go.Scatter(
            x=x[i],
            y=y[i]
        ) for i in range(len(x))]
    elif len(name) == len(x):
        data=[go.Scatter(
            x=x[i],
            y=y[i],
            name=name[i]
        ) for i in range(len(x))]
    else:
        print('error: name should have same length as x and y')
        return 0
    layout = dict(
        title=title,
        scene=dict(
            aspectratio=dict(
                    x=ar[0],
                    y=ar[1],
                )
            ),
        xaxis=dict(title=xaxis),
        yaxis=dict(title=yaxis)
    )
    fig = go.Figure(data=data, layout=layout)
    ply.iplot(fig)
    if ret:
        return fig

def plt_3d_surface_mpl(x,y,z, xlabel='', ylabel='', zlabel='', **kw):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x,y,z, **kw)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)


def combine_plots_ply(figs, title=None, name=None, xdomian=None, ret=False):
    data = []
    for f in figs:
        for i in range(len(f.data)):
            d = f.data[i]
            if name is not None:
                d['name'] = name[i]
            data.append(d)
    nlayout = f.layout
    if title is not None:
        nlayout['title'] = title
    newf = go.Figure(data=data, layout=nlayout)
    #newf = go.Figure(data=data, layout=f.layout)
    ply.iplot(newf)
    if ret:
        return newf

def plt_curv_fit(x,y, title='', xaxis='x', yaxis='y', ar=(10,1), ret=False):
    p1 = np.polyfit(x,y,1)
    f = np.poly1d(p1)
    yf = f(x)
    plt_curve([x,x], [y,yf],name=('data', 'fit: {0}{2}+{1}'.format(np.around(p1[0],3),np.around(p1[1],3), xaxis)), title=title, xaxis=xaxis, yaxis=yaxis,ar=ar)
    if ret:
        return f

def contour_plot_ply(x,y,z, title='', xlabel='', ylabel='', ret=False):
    data = [go.Contour(
        x=x,
        y=y,
        z=z
    )]
    layout = dict(
        title=title,
        xaxis=dict(title=xlabel),
        yaxis=dict(title=ylabel)
    )
    fig = go.Figure(data=data, layout=layout)
    ply.iplot(fig)
    if ret:
        return fig

def contour_plot(x, y, z, name=None, xlabel='', ylabel='', title='', 
                legendloc='upper right', cmap='coolwarm', levelnum=None,
                xlim = None, ylim = None, xticks = None, yticks=None,
                mirror=False, mirror_a=1, figsize=(8,6),
                style=None, font='cm',
                figname=None, dpi=200, ret=False, **kw):
    """Plot figure

    Arguments:
        x {array} -- x coords, can be [coords] for single curve, or [coords_c1, coords_c2, ...] for multi curve
        y {array} -- y coordinates corresponding to x

    Keyword Arguments:
        name {array like} -- names of multi curves, useful for multi curves, ['name_c1', 'name_c2', ...] (default: {None})
        xlabel {str} -- x axis label (default: {''})
        ylabel {str} -- y axis label (default: {''})
        title {str} -- title of figure (default: {''})
        legendloc {str or tuple of float} -- position of lengend, 
                    like (1.05,1) will place at upper right out of frame. (default: {'upper right'})
        scatter {bool or array of bool} -- plot scatter instead of curve if true for ith data. (default: {False})
        figname {str} -- file name for saving the figure, if None, do not save (default: {None})
        dpi {int} -- dpi for showing the figure (default: {200})
        **kw -- other kargs for matplotlib.polt or scatter.
    """
    
    if style is not None:
        plt.style.use(style)
    else:
        mpl.rcParams.update(mpl.rcParamsDefault)
    plt.rcParams["mathtext.fontset"] = font
    fig, ax = plt.subplots(figsize=figsize)

    #ax.legend(loc=legendloc)
    #ax.plot(x,y)
    if mirror:
        if mirror_a == 1:
            y = np.concatenate((-y[::-1][:-1],y))
            z = np.concatenate((z[::-1][:-1],z))
    if levelnum is not None:
        ax.contour(x,y,z, levelnum, **kw)
        ct = ax.contourf(x,y,z, levelnum, cmap=cmap)
    else:
        ax.contour(x,y,z, **kw)
        ct = ax.contourf(x,y,z, cmap=cmap)
    fig.colorbar(ct, ax=ax)
    if xlim is not None:
        plt.xlim(xlim[0],xlim[1])
        #ax.set_xlim(xlim[0],xlim[1])
    if ylim is not None:
        plt.ylim(ylim[0], ylim[1])
        #ax.set_ylim(ylim[0],ylim[1])
    if xticks is not None:
        plt.xticks(xticks)
    if yticks is not None:
        plt.yticks(yticks)
    ax.set_xlabel(xlabel, fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    ax.set_title(title, fontsize=19)
    fig.set_dpi(dpi)
    fig.tight_layout()
    plt.show()
    if figname is not None:
        fig.savefig(figname, format='pdf', bbox_inches='tight')
    if ret:
        return fig, ax

def find_tension_range(t, ediff, kind='cubic', threshold=0.01, x0=None, use_root=True):
    f = interp1d(t, ediff, kind=kind)
    if x0 is None:
        x0 = (t[0]+t[-2])/2
    minimum = fmin(f, x0)
    min_val = f(minimum[0])
    if use_root:
        return brentq(f, minimum, t[-2])
    return minimum[0]
    #bounds = (lambda x: f(x)-threshold*min_val)

def flatten(lst):
    for x in lst:
        if isinstance(x, list):
            for x in flatten(x):
                yield x
        else:
            yield x
def ak_fft(k, y, x):
    """Get coffcient a[k] for discreate input of y, with constant space
    
    Arguments:
        lst {[type]} -- [description]
    """
    return trapz(y*np.exp(-2j*np.pi*k*np.arange(y.size)/y.size), x)/(x[-1]-x[0])

def car2pol(x,y):
    '''cartesian to polar coordinates'''
    return np.hypot(x,y), np.arctan2(y, x)

def pol2car(r, t):
    '''polar to cartesian coordinates'''
    return r * np.cos(t), r * np.sin(t)

def plt_loglog(x, y, slope=None, xlabel=None, ylabel=None, figname=None,
               figsize=None, marker='.', slopetextsize=10, slopetext=None,
               tripos=[0.6, 0.3], textpos_adjust=[-0.1,0.15], triscale=0.5, 
               ax=None, xticklabels=None, yticklabels=None, xticks=None, yticks=None,
               xminorticklabels_off=None, 
               yminorticklabels_off=None):
    if ax is None:
        kws = {'figsize': figsize} if figsize is not None else {}
        fig, ax = plt.subplots(**kws)
    ax.loglog(x, y, marker=marker)
    
    if slope is not None:
        if slopetext is None:
            print('You need set slope text')
            return
        ylmin, ylmax = ax.get_ylim()
        xlmin, xlmax = ax.get_xlim()
        ylmin, ylmax = np.log(ylmin), np.log(ylmax)
        xlmin, xlmax = np.log(xlmin), np.log(xlmax)
        xr = xlmax - xlmin
        yr = ylmax - ylmin

        if tripos[0] > tripos[1]:
            tri_x = np.array([0,1,1,0])
            tri_y = np.array([0,0,slope,0])
        else:
            tri_x = np.array([0,0,1,0])
            tri_y = np.array([0,slope, slope, 0])

        tri_x = triscale * tri_x
        tri_y = triscale * tri_y
        # print(tri_x, tri_y)
        tri_x = tri_x * triscale * xr + xlmin + tripos[0] * xr
        tri_y = tri_y * triscale * xr + ylmin + tripos[1] * yr
        # print(tri_x, tri_y)
        textslope = [np.exp((tri_x[0]+tri_x[2])/2 + textpos_adjust[0]* (tri_x[2] - tri_x[0])), 
                    np.exp((tri_y[0]+ tri_y[2])/2 + textpos_adjust[1]* (tri_y[2] - tri_y[0])),
                    slopetext, slopetextsize]
        # print(tri_x, tri_y)
        # print(tri_x, tri_y)
        tri_x = np.exp(tri_x)
        tri_y = np.exp(tri_y)
        # print(tri_x, tri_y)
        # print(tri_x, tri_y)
        ax.text(textslope[0], textslope[1], textslope[2], fontsize=slopetextsize)
        ax.loglog(tri_x, tri_y, 'black')
    
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xticks is not None:
        ax.set_xticks(xticks)
    if yticks is not None:
        ax.set_yticks(yticks)
    if xticklabels is not None:
        ax.set_xticklabels(xticklabels)
    if yticklabels is not None:
        ax.set_yticklabels(yticklabels)
    if xminorticklabels_off is not None:
        ax.tick_params(axis='x', which='minor', labelbottom=(not xminorticklabels_off))
    if yminorticklabels_off is not None:
        ax.tick_params(axis='y', which='minor', labelleft=(not yminorticklabels_off))

    if figname is not None:
        fig.savefig(figname)