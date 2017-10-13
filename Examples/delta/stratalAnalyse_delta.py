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
    grad[np.where(grad==0)[0]] = -0.1
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
            i2 = ids[k]
            sl1 = sl[i0]
            sl2 = sl[i2]
        if grad[i2] > 0:
            for p in range(ci0,ci1):
                STcolors_ST.append(cTC)
        elif grad[i2] < 0 and sl1 >= sl2:
            for p in range(ci0,ci1):
                STcolors_ST.append(cDRC)
        elif grad[i2] < 0 and sl1 < sl2:
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

def depthID(cs = None, sealevel = None, envIDs = None, side = 'left'):
    """
    Calculate the position of different depositional environments for Wheeler diagram.
    Parameters
    ----------
    variable: sealevel
        The value of sea level through time.
    variable: envIDs
        Range of water depth of each depostional environment.
    variable: side
        Which side of the cross-section: 'left' or 'right'.
    """
    
    if side == 'left':
        envID = np.zeros(len(envIDs))
        for i in range(len(envIDs)):
            envID[i] = np.amin(np.where((cs.secDep[cs.nz-1]) > (sealevel - envIDs[i]))[0])
    
    if side == 'right':
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
                height_bar = None, npts = None, color = None, rangeX = None, rangeY = None, 
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
    variable: height_bar, npts
        Height of the bar, number of extracted outfiles.
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
    plt.rc("font", size=9)
    
    patch_handles = []
    left = np.zeros(rangeE.shape[1])
    for i, d in enumerate(rangeE):
        patch_handles.append(plt.barh(time,d,color=color[i],align='center',left=left, height=height_bar, edgecolor = "none"))
        left += d
    
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

        # right side 
        self.shoreID_r = None
        self.accom_r = None
        self.sed_r = None
        self.depoend_r = None
        
        # left side
        self.shoreID_l = None
        self.accom_l = None
        self.sed_l = None
        self.depoend_l = None
        
        return
    
    def _buildShoreline(self, cs = None, cs_b = None, sealevel = None, sealevel_b = None, style = 'delta'):
        """
        Calculate the shoreline trajectory (shoreID), the change of accommodation (accom)
        and sedimentation (sed) at shoreline, the end point of each depostional layer (depoend).
        Parameters
        ----------
        variable: cs, cs_b
            The cross-section at time t and previous timestep (t-dt)
        variable: sealevel, sealevel_b
            The value of sea-level at time t and previous timestep (t-dt)
        variable : style
            Model style, can be 'delta' or 'basin'.
        """
        
        if style == 'basin':
            # right side
            shoreID = np.amax(np.where(cs.secDep[cs.nz-1]<=sealevel)[0])
            shoreID_b = np.amax(np.where(cs_b.secDep[cs_b.nz-1]<=sealevel_b)[0])
            accom = sealevel - cs.secDep[cs_b.nz-1][shoreID_b]
            sed = cs.secDep[cs.nz-1][shoreID_b] - cs.secDep[cs_b.nz-1][shoreID_b]
            depoend = 0
        
            # left side
            shoreID1 = np.amin(np.where(cs.secDep[cs.nz-1]<=sealevel)[0])
            shoreID1_b = np.amin(np.where(cs_b.secDep[cs_b.nz-1]<=sealevel_b)[0])
            accom1 = sealevel - cs.secDep[cs_b.nz-1][shoreID1_b]
            sed1 = cs.secDep[cs.nz-1][shoreID1_b] - cs.secDep[cs_b.nz-1][shoreID1_b]
            depoend1 = 0
            
        if style == 'delta':
            # right side
            shoreID = np.amax(np.where(cs.secDep[cs.nz-1]>=sealevel)[0])
            shoreID_b = np.amax(np.where(cs_b.secDep[cs_b.nz-1]>=sealevel_b)[0])
            accom = sealevel - cs.secDep[cs_b.nz-1][shoreID_b]
            sed = cs.secDep[cs.nz-1][shoreID_b] - cs.secDep[cs_b.nz-1][shoreID_b]
            depoend = np.amax(np.where(cs.secTh[cs.nz-1][shoreID:len(cs.secTh[0])]>0.001)[0]) + shoreID
        
            # left side
            shoreID1 = np.amin(np.where(cs.secDep[cs.nz-1]>=sealevel)[0])
            shoreID1_b = np.amin(np.where(cs_b.secDep[cs_b.nz-1]>=sealevel_b)[0])
            accom1 = sealevel - cs.secDep[cs_b.nz-1][shoreID1_b]
            sed1 = cs.secDep[cs.nz-1][shoreID1_b] - cs.secDep[cs_b.nz-1][shoreID1_b]
            depoend1 = np.amin(np.where(cs.secTh[cs.nz-1]>0.001)[0])

        return shoreID, accom, sed, depoend, shoreID1, accom1, sed1, depoend1

    def buildParameters(self, npts = None, strat_all = None, sealevel = None, gfilter = 0, style = 'delta'):
        """
        Calculate the shoreline trajectory, accommodation change, sedimentation change and the endpoint of deposition.
        Parameters
        ----------
        variable : npts
            Number of outputs.
        variable : strat_all
            Strata dataset.
        variable : sealevel
            Sea level positions.
        variable : gfilter
            Gaussian smoothing filter.
        variable : style
            Model style, can be 'delta' or 'basin'.
        """

        shoreID_r = np.zeros(npts)
        accom_r = np.zeros(npts)
        sed_r = np.zeros(npts)
        depoend_r = np.zeros(npts)
        
        shoreID_l = np.zeros(npts)
        accom_l = np.zeros(npts)
        sed_l = np.zeros(npts)
        depoend_l = np.zeros(npts)
        
        # time 0
        shoreID_r[0], accom_r[0], sed_r[0], depoend_r[0], shoreID_l[0], accom_l[0], sed_l[0], depoend_l[0] = self._buildShoreline(cs = strat_all[0], cs_b = strat_all[0], sealevel = sealevel[0], sealevel_b = sealevel[0], style = style)
        # time 1-npts
        for i in range(1,npts):
            shoreID_r[i], accom_r[i], sed_r[i], depoend_r[i], shoreID_l[i], accom_l[i], sed_l[i], depoend_l[i]  = self._buildShoreline(cs = strat_all[i], cs_b = strat_all[i-1], sealevel = sealevel[i], sealevel_b = sealevel[i-1], style = style)
            
        self.shoreID_r = filters.gaussian_filter1d(shoreID_r, sigma=gfilter)
        self.accom_r = filters.gaussian_filter1d(accom_r, sigma=gfilter)
        self.sed_r = filters.gaussian_filter1d(sed_r, sigma=gfilter)
        self.depoend_r = filters.gaussian_filter1d(depoend_r, sigma=gfilter)
        
        self.shoreID_l = filters.gaussian_filter1d(shoreID_l, sigma=gfilter)
        self.accom_l = filters.gaussian_filter1d(accom_l, sigma=gfilter)
        self.sed_l = filters.gaussian_filter1d(sed_l, sigma=gfilter)
        self.depoend_l = filters.gaussian_filter1d(depoend_l, sigma=gfilter)

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
