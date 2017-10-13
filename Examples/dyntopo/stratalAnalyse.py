##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling companion.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
Here we set usefull functions used to analyse stratigraphic sequences from Badlands outputs.
"""

import os
import math
import h5py
import errno
import numpy as np
import pandas as pd
from cmocean import cm
import colorlover as cl
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import xml.etree.ElementTree as ETO
import scipy.ndimage.filters as filters
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage.filters import gaussian_filter

import plotly
from plotly import tools
from plotly.graph_objs import *
plotly.offline.init_notebook_mode()

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

def readSea(seafile):
    """
    Plot sea level curve.
    Parameters
    ----------
    variable: seafile
        Absolute path of the sea-lelve data.
    """

    df=pd.read_csv(seafile, sep=r'\s+',header=None)
    SLtime,sealevel = df[0],df[1]

    return SLtime,sealevel

def viewData(x0 = None, y0 = None, width = 800, height = 400, linesize = 3, color = '#6666FF',
             xlegend = 'xaxis', ylegend = 'yaxis', title = 'view data'):
    """
    Plot multiple data on a graph.
    Parameters
    ----------
    variable: x0, y0
        Data for plot
    variable: width, height
        Figure width and height.
    variable: linesize
        Requested size for the line.
    variable: color
        
    variable: xlegend
        Legend of the x axis.
    variable: ylegend
        Legend of the y axis.
    variable: title
        Title of the graph.
    """
    trace = Scatter(
        x=x0,
        y=y0,
        mode='lines',
        line=dict(
            shape='line',
            color = color,
            width = linesize
        ),
        fill=None
    )

    layout = dict(
            title=title,
            font=dict(size=10),
            width=width,
            height=height,
            showlegend = False,
            xaxis=dict(title=xlegend,
                       ticks='outside',
                       zeroline=False,
                       showline=True,
                       mirror='ticks'),
            yaxis=dict(title=ylegend,
                       ticks='outside',
                       zeroline=False,
                       showline=True,
                       mirror='ticks')
            )

    fig = Figure(data=[trace], layout=layout)
    plotly.offline.iplot(fig)

    return

def buildShore(cs = None, cs_b = None, sealevel = None, sealevel_b = None):
    """
    Calculate the shoreline trajectory (shoreID), the change of accommodation (accom)
    and sedimentation (sed) at shoreline, the end point of each depostional layer (depoend).
    Parameters
    ----------
    variable: cs, cs_b
        The cross-section at time t and previous timestep (t-dt)
    variable: sealevel, sealevel_b
        The value of sea-level at time t and previous timestep (t-dt)
    """
    shoreID = np.amax(np.where(cs.secDep[cs.nz-1]>=sealevel)[0])
    shoreID_b = np.amax(np.where(cs_b.secDep[cs_b.nz-1]>=sealevel_b)[0])
    accom = sealevel - cs.secDep[cs_b.nz-1][shoreID_b]
    sed = cs.secDep[cs.nz-1][shoreID_b] - cs.secDep[cs_b.nz-1][shoreID_b]
    depoend = np.amax(np.where(cs.secTh[cs.nz-1][shoreID_b:len(cs.secTh[0])]>0.001)[0]) + shoreID_b

    return shoreID, accom, sed, depoend

def build_shoreTrajectory(x, y, grad, sl, nbout, cTC='rgba(56,110,164,0.8)', cDRC='rgba(60,165,67,0.8)',
                          cARC='rgba(112,54,127,0.8)', cSTC='rgba(252,149,7,0.8)'):
    """
    Automatic delimitation of shoreline trajectory classes.
    Parameters
    ----------
    variable: x
        display steps
    variable: y
        shoreline position
    variable: grad
        shoreline position gradient
    variable: sl
        sealevel position
    variable: nbout
        number of output steps

    color schema for different classes
    """

    # Find intersection between line zero and the shoreline gradient trajectory
    ids = np.argwhere(np.diff(np.sign(grad - np.zeros(len(grad)))) != 0).reshape(-1) + 0
    # Number of points to consider
    nbclass = len(ids)

    # Check if there are still some points after the last intersection
    final = False
    if ids[-1]<len(grad):
        nbclass += 1
        final = True
    # Build the color list
    STcolors_ST = []

    ci0 = 0
    i0 = 0
    for k in range(nbclass):
        if k == nbclass-1:
            if not final:
                exit
            else:
                i1 = nbout
                ci1 = nbout
                i2 = -1
                sl1 = sl[i0]
                sl2 = sl[-1]
        else:
            i1 = ids[k]
            ci1 = int(x[ids[k]])
            i2 = ids[k]-1
            sl1 = sl[i0]
            sl2 = sl[ids[k]-1]
        if grad[i2] < 0:
            for p in range(ci0,ci1):
                STcolors_ST.append(cTC)
        elif grad[i2] > 0 and sl1 >= sl2:
            for p in range(ci0,ci1):
                STcolors_ST.append(cDRC)
        elif grad[i2] > 0 and sl1 < sl2:
            for p in range(ci0,ci1):
                STcolors_ST.append(cARC)
        else:
            for p in range(ci0,ci1):
                STcolors_ST.append(cSTC)
        if k < nbclass-1:
            i0 = ids[k]
            ci0 = int(x[ids[k]])

    return STcolors_ST

def build_accomSuccession(x, y, grad, nbout, cR='rgba(51,79,217,0.8)', cPA='rgba(252,149,7,0.8)',
                          cAPD='rgba(15,112,2,0.8)'):
    """
    Automatic delimitation of accommodation succession sequence sets.
    Parameters
    ----------
    variable: x
        display steps
    variable: y
        AS curve
    variable: grad
        shoreline position gradient
    variable: nbout
        number of output steps

    color schema for different classes
    """

    # Find intersection between line zero and the AS curve
    ids1 = np.argwhere(np.diff(np.sign(y - np.zeros(len(y)))) != 0).reshape(-1) + 0
    # Find intersection between line zero and the AS gradient
    ids2 = np.argwhere(np.diff(np.sign(grad - np.zeros(len(y)))) != 0).reshape(-1) + 0
    # Combine ids together
    ids = np.concatenate((ids1,ids2))
    ids.sort(kind='mergesort')

    # Number of points to consider
    nbclass = len(ids)

    # Check if there are still some points after the last intersection
    final = False
    if ids[-1]<len(grad):
        nbclass += 1
        final = True

    # Build the color list
    STcolors_AS = []

    ci0 = 0
    i0 = 0
    for k in range(nbclass):
        if k == nbclass-1:
            if not final:
                exit
            else:
                i1 = nbout
                ci1 = nbout
                i2 = -1
        else:
            i1 = ids[k]
            ci1 = int(x[ids[k]])
            i2 = ids[k]-1
        if y[i2-1] >= 0:
            for p in range(ci0,ci1):
                STcolors_AS.append(cR)
        elif y[i2-1] < 0 and grad[i2-1] >= 0:
            for p in range(ci0,ci1):
                STcolors_AS.append(cAPD)
        elif y[i2-1] < 0 and grad[i2-1] < 0:
            for p in range(ci0,ci1):
                STcolors_AS.append(cPA)
        if k < nbclass-1:
            i0 = ids[k]
            ci0 = int(x[ids[k]])

    return STcolors_AS

def depthID(cs = None, sealevel = None, envIDs = None):
    """
    Calculate the position of different depositional environments for Wheeler diagram.
    Parameters
    ----------
    variable: sealevel
        The value of sea level through time.
    variable: envIDs
        Range of water depth of each depostional environment.
    """
    envID = np.zeros(len(envIDs))
    for i in range(len(envIDs)):
        envID[i] = np.amax(np.where((cs.secDep[cs.nz-1]) > (sealevel - envIDs[i]))[0])

    return envID

def viewSection(width = 800, height = 400, cs = None, dnlay = None,
                rangeX = None, rangeY = None, linesize = 3, title = 'Cross section'):
    """
    Plot multiple cross-sections data on a graph.
    Parameters
    ----------
    variable: width
        Figure width.
    variable: height
        Figure height.
    variable: cs
        Cross-sections dataset.
    variable: dnlay
        Layer step to plot the cross-section.
    variable: rangeX, rangeY
        Extent of the cross section plot.
    variable: linesize
        Requested size for the line.
    variable: title
        Title of the graph.
    """
    nlay = len(cs.secDep)
    colors = cl.scales['9']['div']['BrBG']
    hist = cl.interp( colors, nlay )
    colorrgb = cl.to_rgb( hist )

    trace = {}
    data = []

    trace[0] = Scatter(
        x=cs.dist,
        y=cs.secDep[0],
        mode='lines',
        line=dict(
            shape='line',
            width = linesize+2,
            color = 'rgb(0, 0, 0)'
        )
    )
    data.append(trace[0])

    for i in range(1,nlay-1,dnlay):
        trace[i] = Scatter(
            x=cs.dist,
            y=cs.secDep[i],
            mode='lines',
            line=dict(
                shape='line',
                width = linesize,
                color = 'rgb(0,0,0)'
            ),
            opacity=0.5,
            fill='tonexty',
            fillcolor=colorrgb[i]
        )
        data.append(trace[i])

    trace[nlay-1] = Scatter(
        x=cs.dist,
        y=cs.secDep[nlay-1],
        mode='lines',
        line=dict(
            shape='line',
            width = linesize+2,
            color = 'rgb(0, 0, 0)'
        ),
        fill='tonexty',
        fillcolor=colorrgb[nlay-1]
    )
    data.append(trace[nlay-1])

    trace[nlay] = Scatter(
        x=cs.dist,
        y=cs.secDep[0],
        mode='lines',
        line=dict(
            shape='line',
            width = linesize+2,
            color = 'rgb(0, 0, 0)'
        )
    )
    data.append(trace[nlay])

    if rangeX is not None and rangeY is not None:
        layout = dict(
                title=title,
                font=dict(size=10),
                width=width,
                height=height,
                showlegend = False,
                xaxis=dict(title='distance [m]',
                            range=rangeX,
                            ticks='outside',
                            zeroline=False,
                            showline=True,
                            mirror='ticks'),
                yaxis=dict(title='elevation [m]',
                            range=rangeY,
                            ticks='outside',
                            zeroline=False,
                            showline=True,
                            mirror='ticks')
        )
    else:
        layout = dict(
                title=title,
                font=dict(size=10),
                width=width,
                height=height
        )
    fig = Figure(data=data, layout=layout)
    plotly.offline.iplot(fig)

    return


def viewSection_Depo(width = 8, height = 5, cs = None, IPs = None, dnlay = None, color = None,
                      rangeX = None, rangeY = None, linesize = 3, title = 'Cross section'):
    """
    Plot stratal stacking pattern colored by deposition depth.
    Parameters
    ----------
    variable: cs
        Cross-sections dataset.
    variable: dnlay
        Layer step to plot the cross-section.
    variable: colors
        Colors for different ranges of water depth (i.e. depositional environments).
    variable: rangeX, rangeY
        Extent of the cross section plot.
    variable: linesize
        Requested size for the line.
    variable: title
        Title of the graph.
    """
    fig = plt.figure(figsize = (width,height))
    plt.rc("font", size=10)

    ax = fig.add_subplot(111)
    layID = []
    p = 0
    xi00 = cs.dist
    # color = ['limegreen','sandybrown','khaki','lightsage','c','dodgerblue']
    for i in range(0,cs.nz+1,dnlay):
        if i == cs.nz:
            i = cs.nz-1
        layID.append(i)
        if len(layID) > 1:
            for j in range(0,int(np.amax(xi00))):
                if (j<IPs[0][i/dnlay]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color='limegreen')
                elif (j<IPs[1][i/dnlay]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[0])
                elif (j<IPs[2][i/dnlay]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[1])
                elif (j<IPs[3][i/dnlay]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[2])
                elif (j<IPs[4][i/dnlay]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[3])
                elif (j<IPs[5][i/dnlay]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[4])
                else:
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color='teal')
                    plt.fill_between(xi00, cs.secDep[layID[p]], 0, color='white')
        p=p+1
    for i in range(0,cs.nz,dnlay):
        if i>0:
            plt.plot(xi00,cs.secDep[i],'-',color='k',linewidth=0.2)
    plt.plot(xi00,cs.secDep[cs.nz-1],'-',color='k',linewidth=0.7)
    plt.plot(xi00,cs.secDep[0],'-',color='k',linewidth=0.7)
    plt.xlim( rangeX )
    plt.ylim( rangeY )

    return

def viewSection_Depth(cs = None, IPs = None, dnlay = None, color = None,
                      rangeX = None, rangeY = None, linesize = 3, title = 'Cross section'):
    """
    Plot stratal stacking pattern colored by water depth.
    Parameters
    ----------
    variable: cs
        Cross-sections dataset.
    variable: dnlay
        Layer step to plot the cross-section.
    variable: color
        Colors for different ranges of water depth (i.e. depositional environments).
    variable: rangeX, rangeY
        Extent of the cross section plot.
    variable: linesize
        Requested size for the line.
    variable: title
        Title of the graph.
    """
    fig = plt.figure(figsize = (9,5))
    plt.rc("font", size=10)

    ax = fig.add_subplot(111)
    layID = []
    p = 0
    xi00 = cs.dist
    # color = ['limegreen','sandybrown','khaki','lightsage','c','dodgerblue']
    for i in range(0,cs.nz,dnlay):
        if i == cs.nz:
            i = i-1
        layID.append(i)
        if len(layID) > 1:
            for j in range(0,600):
                if (j<IPs[i][0]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[0])
                elif (j<IPs[i][1]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[0])
                elif (j<IPs[i][2]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[1])
                elif (j<IPs[i][3]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[2])
                elif (j<IPs[i][4]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[3])
                elif (j<IPs[i][5]):
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[4])
                else:
                    plt.fill_between([xi00[j],xi00[j+1]], [cs.secDep[layID[p-1]][j], cs.secDep[layID[p-1]][j+1]], color=color[5])
                    plt.fill_between(xi00, strat.secDep[layID[p]], 0, color='white')
        p=p+1
    for i in range(0,cs.nz,dnlay):
        if i>0:
            plt.plot(xi00,cs.secDep[i],'-',color='k',linewidth=0.2)
    plt.plot(xi00,cs.secDep[cs.nz-1],'-',color='k',linewidth=0.7)
    plt.plot(xi00,cs.secDep[0],'-',color='k',linewidth=0.7)
    plt.xlim( rangeX )
    plt.ylim( rangeY )

    return

def viewSectionST(width = 800, height = 400, cs = None, dnlay = None, colors=None,
                  rangeX = None, rangeY = None, linesize = 3, title = 'Cross section'):
    """
    Plot multiple cross-sections colored by system tracts on a graph.
    Parameters
    ----------
    variable: width
        Figure width.
    variable: height
        Figure height.
    variable: cs
        Cross-sections dataset.
    variable: dnlay
        Layer step to plot the cross-section.
    variable: colors
        System tract color scale.
    variable: rangeX, rangeY
        Extent of the cross section plot.
    variable: linesize
        Requested size for the line.
    variable: title
        Title of the graph.
    """
    nlay = len(cs.secDep)

    trace = {}
    data = []

    trace[0] = Scatter(
        x=cs.dist,
        y=cs.secDep[0],
        mode='lines',
        line=dict(
            shape='line',
            width = linesize+2,
            color = 'rgb(0, 0, 0)'
        )
    )
    data.append(trace[0])

    for i in range(1,nlay-1,dnlay):
        trace[i] = Scatter(
            x=cs.dist,
            y=cs.secDep[i],
            mode='lines',
            line=dict(
                shape='line',
                width = linesize,
                color = 'rgb(0,0,0)'
            ),
            opacity=0.5,
            fill='tonexty',
            fillcolor=colors[i]
        )
        data.append(trace[i])

    trace[nlay-1] = Scatter(
        x=cs.dist,
        y=cs.secDep[nlay-1],
        mode='lines',
        line=dict(
            shape='line',
            width = linesize+2,
            color = 'rgb(0, 0, 0)'
        ),
        fill='tonexty',
        fillcolor=colors[nlay-1]
    )
    data.append(trace[nlay-1])

    trace[nlay] = Scatter(
        x=cs.dist,
        y=cs.secDep[0],
        mode='lines',
        line=dict(
            shape='line',
            width = linesize+2,
            color = 'rgb(0, 0, 0)'
        )
    )
    data.append(trace[nlay])

    if rangeX is not None and rangeY is not None:
        layout = dict(
                title=title,
                font=dict(size=10),
                width=width,
                height=height,
                showlegend = False,
                xaxis=dict(title='distance [m]',
                            range=rangeX,
                            ticks='outside',
                            zeroline=False,
                            showline=True,
                            mirror='ticks'),
                yaxis=dict(title='elevation [m]',
                            range=rangeY,
                            ticks='outside',
                            zeroline=False,
                            showline=True,
                            mirror='ticks')
        )
    else:
        layout = dict(
                title=title,
                font=dict(size=10),
                width=width,
                height=height
        )
    fig = Figure(data=data, layout=layout)
    plotly.offline.iplot(fig)

    return

def viewWheeler(width = 800, height = 400, time = None, rangeE = None, shoreID = None, 
                dnlay = None, npts = None, color = None, rangeX = None, rangeY = None, 
                linesize = 3, title = 'Wheeler Diagram', xlegend = 'xaxis', ylegend = 'yaxis'):
    """
    Plot wheeler diagram colored by deposition environment on a graph.
    Parameters
    ----------
    variable: width, height
        Figure width and height.
    variable: time
        Stacking time of each extracted layer.
    variable: rangeE
        Depositional environments extent.
    variable: shoreID
        Shoreline position through time.
    variable: dnlay, npts
        Layer step to plot the Wheeler diagram, number of extracted outfiles.
    variable: color
        Depositional environments color scale.
    variable: rangeX, rangeY
        Extent of the cross section plot.
    variable: linesize
        Requested size for the line.
    variable: title
        Title of the graph.
    variable: xlegend
        Legend of the x axis.
    variable: ylegend
        Legend of the y axis.
    """

    fig = plt.figure(figsize = (width,height))
    plt.rc("font", size=10)
    
    patch_handles = []
    for i, d in enumerate(rangeE):
        patch_handles.append(plt.barh(time,d,color=color[i],align='edge',left=d, height=dnlay/10., edgecolor = "none"))
    
    for j in range(0,npts): 
        plt.axhline(time[j], color='k', linewidth=0.5)
        
    plt.plot(shoreID, time,'ko',markersize=3)
    #
    plt.xlim(rangeX)
    plt.ylim(rangeY)
    plt.xlabel(xlegend)
    plt.ylabel(ylegend)
    #     
    plt.title(title)
    
    return

def getStack(cs = None, posit = None, envIDs = None, color = None, dn = None):
    """
    Compute vertical stacking parameters.
    variable: cs
        Cross-sections dataset.
    variable: posit
        Position of wells.
    variable: envIDs
        Depositional environment.
    variable: colors
        Depositional environments color scale.
    variable: dn
        Number of output steps.
    """
    color_fill = []
    for i in range(len(posit)):
        for j in range(0,cs.nz,dn):
            if ((cs.secElev[j][posit[i]]) > (- envIDs[0])):
                color_fill.append(color[0])
            elif (cs.secElev[j][posit[i]]) > (- envIDs[1]):
                color_fill.append(color[0])
            elif (cs.secElev[j][posit[i]]) > (- envIDs[2]):
                color_fill.append(color[1])
            elif (cs.secElev[j][posit[i]]) > (- envIDs[3]):
                color_fill.append(color[2])
            elif (cs.secElev[j][posit[i]]) > (- envIDs[4]):
                color_fill.append(color[3])
            elif (cs.secElev[j][posit[i]]) > (- envIDs[5]):
                color_fill.append(color[4])
            else:
                color_fill.append(color[5])          
    
    nbout = cs.nz
    layth = []
    for m in range(len(posit)):
        nz = cs.nz-1
        layth.append(cs.secDep[nz][posit[m]])
        for n in range(1,int(nbout/dn)):
            if nz-n*dn >= 0:
                layth.append(-sum(cs.secTh[(nz-n*dn):(nz-(n-1)*dn)])[posit[m]])
    
    colorFill = np.reshape(color_fill, (len(posit), int(nbout/dn)))
    layTh = np.reshape(layth, (len(posit),  int(nbout/dn)))
    
    return colorFill, layTh

def viewStack(width = 800, height = 400, layTh = None, colorFill = None):
    """
    Plot wheeler diagram colored by deposition environment on a graph.
    Parameters
    ----------
    variable: width, height
        Figure width and height.
    variable: layTh
        Layer thickness for each wells.
    variable: colorFill
        Layer environment color for each wells.
    """

    fig = plt.figure(figsize = (width,height))
    plt.rc("font", size=10)
    
    ax = fig.add_axes([0.2,0.06,0.82,0.91])

    data = layTh
    for k in range(len(data)):
        bottom = np.cumsum(data[k], axis=0)
        colors = np.fliplr([colorFill[k]])[0]
        plt.bar(2*k, data[k][0], color = 'w', edgecolor='lightgrey', hatch = '/')
        for j in range(1, data[k].shape[0]):
            plt.bar(2*k, data[k][j], color=colors[j], edgecolor='black', bottom=bottom[j-1])

    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.axes.get_xaxis().set_visible(False)
    ax.tick_params(axis='both', labelsize=8)
    ax.yaxis.set_ticks_position('left')
    
    #plt.xlim(-0.4,10)
    #plt.ylim(-800,0)

    plt.ylabel('Elevation (m)',fontsize=10)
    #plt.yticks(fontsize=10)
    
    return

class stratalSection:
    """
    Class for creating stratigraphic cross-sections from Badlands outputs.
    """

    def __init__(self, folder=None, ncpus=1):
        """
        Initialization function which takes the folder path to Badlands outputs
        and the number of CPUs used to run the simulation.

        Parameters
        ----------
        variable : folder
            Folder path to Badlands outputs.
        variable: ncpus
            Number of CPUs used to run the simulation.
        """

        self.folder = folder
        if not os.path.isdir(folder):
            raise RuntimeError('The given folder cannot be found or the path is incomplete.')

        self.ncpus = ncpus
        if ncpus > 1:
            raise RuntimeError('Multi-processors function not implemented yet!')

        self.x = None
        self.y = None
        self.xi = None
        self.yi = None
        self.dx = None
        self.dist = None
        self.dx = None
        self.nx = None
        self.ny = None
        self.nz = None
        self.dep = None
        self.th = None
        self.elev = None
        self.xsec = None
        self.ysec = None
        self.secTh = []
        self.secDep = []
        self.secElev = []

        self.shoreID = None
        self.accom_shore = None
        self.sed_shore = None
        self.depoend = None

        self.shoreID_gs = None
        self.accom_shore_gs = None
        self.sed_shore_gs = None
        self.depoend_gs = None

        return

    def _buildShoreline(self, cs = None, cs_b = None, sealevel = None, sealevel_b = None):

        shoreID = np.amax(np.where(cs.secDep[cs.nz-1]>=sealevel)[0])
        shoreID_b = np.amax(np.where(cs_b.secDep[cs_b.nz-1]>=sealevel_b)[0])
        accom = sealevel - cs.secDep[cs_b.nz-1][shoreID_b]
        sed = cs.secDep[cs.nz-1][shoreID_b] - cs.secDep[cs_b.nz-1][shoreID_b]
        depoend = np.amax(np.where(cs.secTh[cs.nz-1][shoreID:len(cs.secTh[0])]>0.001)[0]) + shoreID

        return shoreID, accom, sed, depoend

    def buildParameters(self, npts, strat_all, sealevel):

        # Calculate cross-section parameters
        shoreID = np.zeros(npts)
        accom_shore = np.zeros(npts)
        sed_shore = np.zeros(npts)
        depoend = np.zeros(npts)

        shoreID[0], accom_shore[0], sed_shore[0], depoend[0] = self._buildShoreline(cs = strat_all[0],
                                            cs_b = strat_all[0], sealevel = sealevel[0], sealevel_b = sealevel[0])

        for i in range(1,npts):
            shoreID[i], accom_shore[i], sed_shore[i], depoend[i],  = self._buildShoreline(cs = strat_all[i], cs_b = strat_all[i-1], sealevel = sealevel[i], sealevel_b = sealevel[i-1])

        self.shoreID = shoreID
        self.accom_shore = accom_shore
        self.sed_shore = sed_shore
        self.depoend = depoend

        # Gaussian smooth
        self.shoreID_gs = filters.gaussian_filter1d(shoreID,sigma=1)
        self.accom_gs = filters.gaussian_filter1d(accom_shore,sigma=1)
        self.sed_gs = filters.gaussian_filter1d(sed_shore,sigma=1)
        self.depoend_gs = filters.gaussian_filter1d(depoend,sigma=1)

        return

    def loadStratigraphy(self, timestep=0):
        """
        Read the HDF5 file for a given time step.
        Parameters
        ----------
        variable : timestep
            Time step to load.
        """

        for i in range(0, self.ncpus):
            df = h5py.File('%s/sed.time%s.p%s.hdf5'%(self.folder, timestep, i), 'r')
            #print(list(df.keys()))
            coords = np.array((df['/coords']))
            layDepth = np.array((df['/layDepth']))
            layElev = np.array((df['/layElev']))
            layThick = np.array((df['/layThick']))
            if i == 0:
                x, y = np.hsplit(coords, 2)
                dep = layDepth
                elev = layElev
                th = layThick

        self.dx = x[1]-x[0]
        self.x = x
        self.y = y
        self.nx = int((x.max() - x.min())/self.dx+1)
        self.ny = int((y.max() - y.min())/self.dx+1)
        self.nz = dep.shape[1]
        self.xi = np.linspace(x.min(), x.max(), self.nx)
        self.yi = np.linspace(y.min(), y.max(), self.ny)
        self.dep = dep.reshape((self.ny,self.nx,self.nz))
        self.elev = elev.reshape((self.ny,self.nx,self.nz))
        self.th = th.reshape((self.ny,self.nx,self.nz))

        return

    def _cross_section(self, xo, yo, xm, ym, pts):
        """
        Compute cross section coordinates.
        """

        if xm == xo:
            ysec = np.linspace(yo, ym, pts)
            xsec = np.zeros(pts)
            xsec.fill(xo)
        elif ym == yo:
            xsec = np.linspace(xo, xm, pts)
            ysec = np.zeros(pts)
            ysec.fill(yo)
        else:
            a = (ym-yo)/(xm-xo)
            b = yo - a * xo
            xsec = np.linspace(xo, xm, pts)
            ysec = a * xsec + b

        return xsec, ysec

    def buildSection(self, xo = None, yo = None, xm = None, ym = None,
                    pts = 100, gfilter = 5):
        """
        Extract a slice from the 3D data set and compute the stratigraphic layers.
        Parameters
        ----------
        variable: xo, yo
            Lower X,Y coordinates of the cross-section.
        variable: xm, ym
            Upper X,Y coordinates of the cross-section.
        variable: pts
            Number of points to discretise the cross-section.
        variable: gfilter
            Gaussian smoothing filter.
        """

        if xm > self.x.max():
            xm = self.x.max()

        if ym > self.y.max():
            ym = self.y.max()

        if xo < self.x.min():
            xo = self.x.min()

        if yo < self.y.min():
            yo = self.y.min()

        xsec, ysec = self._cross_section(xo, yo, xm, ym, pts)
        self.dist = np.sqrt(( xsec - xo )**2 + ( ysec - yo )**2)
        self.xsec = xsec
        self.ysec = ysec
        for k in range(self.nz):
            # Thick
            rect_B_spline = RectBivariateSpline(self.yi, self.xi, self.th[:,:,k])
            data = rect_B_spline.ev(ysec, xsec)
            secTh = filters.gaussian_filter1d(data,sigma=gfilter)
            secTh[secTh < 0] = 0
            self.secTh.append(secTh)

            # Elev
            rect_B_spline1 = RectBivariateSpline(self.yi, self.xi, self.elev[:,:,k])
            data1 = rect_B_spline1.ev(ysec, xsec)
            secElev = filters.gaussian_filter1d(data1,sigma=gfilter)
            self.secElev.append(secElev)

            # Depth
            rect_B_spline2 = RectBivariateSpline(self.yi, self.xi, self.dep[:,:,k])
            data2 = rect_B_spline2.ev(ysec, xsec)
            secDep = filters.gaussian_filter1d(data2,sigma=gfilter)
            self.secDep.append(secDep)

        # Ensure the spline interpolation does not create underlying layers above upper ones
        topsec = self.secDep[self.nz-1]
        for k in range(self.nz-2,-1,-1):
            secDep = self.secDep[k]
            self.secDep[k] = np.minimum(secDep, topsec)
            topsec = self.secDep[k]

        return
