##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling companion.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
Here we set usefull functions used to analyse morphometrics from Badlands outputs.
"""

import os
import math
import h5py
import errno
import pandas
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
import xml.etree.ElementTree as ETO
from scipy.interpolate import RectBivariateSpline

import plotly
from plotly.graph_objs import *
plotly.offline.init_notebook_mode()

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

class morphoGrid:
    """
    Class for analysing morphometrics from Badlands outputs.
    """

    def __init__(self, folder=None, ncpus=1, bbox=None, dx=None):
        """
        Initialization function which takes the folder path to Badlands outputs
        and the number of CPUs used to run the simulation. It also takes the
        bounding box and discretization value at which one wants to interpolate
        the data.

        Parameters
        ----------
        variable : folder
            Folder path to Badlands outputs.

        variable: ncpus
            Number of CPUs used to run the simulation.

        variable: bbox
            Bounding box extent SW corner and NE corner.

        variable: dx
            Discretisation value in metres.

        """

        self.folder = folder
        if not os.path.isdir(folder):
            raise RuntimeError('The given folder cannot be found or the path is incomplete.')

        self.ncpus = ncpus
        self.x = None
        self.y = None
        self.z = None
        self.discharge = None
        self.logdischarge = None
        self.cumchange = None
        self.dx = None
        self.grad = None
        self.aspect = None
        self.hcurv = None
        self.vcurv = None
        self.Zbc = None
        self.hillshade = None
        self.nx = None
        self.ny = None

        if dx == None:
            raise RuntimeError('Discretization space value is required.')
        self.dx = dx
        self.bbox = bbox

        return

    def loadHDF5(self, timestep=0):
        """
        Read the HDF5 file for a given time step.

        Parameters
        ----------
        variable : timestep
            Time step to load.

        """

        for i in range(0, self.ncpus):
            df = h5py.File('%s/tin.time%s.p%s.hdf5'%(self.folder, timestep, i), 'r')
            coords = np.array((df['/coords']))
            cumdiff = np.array((df['/cumdiff']))
            discharge = np.array((df['/discharge']))
            if i == 0:
                x, y, z = np.hsplit(coords, 3)
                c = cumdiff
                d = discharge
            else:
                c = np.append(c, cumdiff)
                d = np.append(d, discharge)
                x = np.append(x, coords[:,0])
                y = np.append(y, coords[:,1])
                z = np.append(z, coords[:,2])

        if self.bbox == None:
            self.nx = int((x.max() - x.min())/self.dx+1)
            self.ny = int((y.max() - y.min())/self.dx+1)
            self.x = np.linspace(x.min(), x.max(), self.nx)
            self.y = np.linspace(y.min(), y.max(), self.ny)
            self.bbox = np.zeros(4,dtype=float)
            self.bbox[0] = x.min()
            self.bbox[1] = y.min()
            self.bbox[2] = x.max()
            self.bbox[3] = y.max()
        else:
            if self.bbox[0] < x.min():
                self.bbox[0] = x.min()
            if self.bbox[2] > x.max():
                self.bbox[2] = x.max()
            if self.bbox[1] < y.min():
                self.bbox[1] = y.min()
            if self.bbox[3] > y.max():
                self.bbox[3] = y.max()
            self.nx = int((self.bbox[2] - self.bbox[0])/self.dx+1)
            self.ny = int((self.bbox[3] - self.bbox[1])/self.dx+1)
            self.x = np.linspace(self.bbox[0], self.bbox[2], self.nx)
            self.y = np.linspace(self.bbox[1], self.bbox[3], self.ny)

        self.x, self.y = np.meshgrid(self.x, self.y)
        xyi = np.dstack([self.x.flatten(), self.y.flatten()])[0]
        XY = np.column_stack((x,y))
        tree = cKDTree(XY)
        distances, indices = tree.query(xyi, k=3)
        z_vals = z[indices][:,:,0]
        d_vals = d[indices][:,:,0]
        c_vals = c[indices][:,:,0]
        
        zi = np.zeros(len(xyi))
        di = np.zeros(len(xyi))
        ci = np.zeros(len(xyi))
        onIDs = np.where(distances[:,0] > 0)[0]
        zi[onIDs] = np.average(z_vals[onIDs,:],weights=(1./distances[onIDs,:]), axis=1)
        di[onIDs] = np.average(d_vals[onIDs,:],weights=(1./distances[onIDs,:]), axis=1)
        ci[onIDs] = np.average(c_vals[onIDs,:],weights=(1./distances[onIDs,:]), axis=1)

        onIDs = np.where(distances[:,0] == 0)[0]
        
        if len(onIDs) > 0:
            zi[onIDs] = z[indices[onIDs,0],0]
            di[onIDs] = d[indices[onIDs,0],0]
            ci[onIDs] = c[indices[onIDs,0],0]

        self.z = np.reshape(zi,(self.ny,self.nx))
        self.discharge = np.reshape(di,(self.ny,self.nx))
        self.cumchange = np.reshape(ci,(self.ny,self.nx))

        logdis = self.discharge
        IDs = np.where(logdis<1.)
        logdis[IDs] = 1.
        self.logdischarge = logdis

        return

    def _assignBCs(self):
        """
        Pads the boundaries of a grid. Boundary condition pads the boundaries
        with equivalent values to the data margins, e.g. x[-1,1] = x[1,1].
        It creates a grid 2 rows and 2 columns larger than the input.

        """

        self.Zbc = np.zeros((self.ny + 2, self.nx + 2))
        self.Zbc[1:-1,1:-1] = self.z

        # Assign boundary conditions - sides
        self.Zbc[0, 1:-1] = self.z[0, :]
        self.Zbc[-1, 1:-1] = self.z[-1, :]
        self.Zbc[1:-1, 0] = self.z[:, 0]
        self.Zbc[1:-1, -1] = self.z[:,-1]

        # Assign boundary conditions - corners
        self.Zbc[0, 0] = self.z[0, 0]
        self.Zbc[0, -1] = self.z[0, -1]
        self.Zbc[-1, 0] = self.z[-1, 0]
        self.Zbc[-1, -1] = self.z[-1, 0]

        return

    def _calcFiniteSlopes(self):
        """
        Calculate slope with 2nd order/centered difference method.

        """

        self._assignBCs()
        Sx = (self.Zbc[1:-1, 2:] - self.Zbc[1:-1, :-2]) / (2*self.dx)
        Sy = (self.Zbc[2:,1:-1] - self.Zbc[:-2, 1:-1]) / (2*self.dx)

        return Sx, Sy

    def hillShade(self, az=315, altitude=45):
        """
        Creates a shaded relief from a surface raster by considering the
        illumination source angle and shadows.

        Parameters
        ----------
        variable : az
            Azimuth angle of the light source.The azimuth is expressed in positive
            degrees from 0 to 360, measured clockwise from north.
            The default is 315 degrees.

        variable : altitude
            Altitude angle of the light source above the horizon. The altitude is
            expressed in positive degrees, with 0 degrees at the horizon and 90
            degrees directly overhead. The default is 45 degrees.
        """

        # Convert angular measurements to radians
        azRad, elevRad = (360 - az + 90)*np.pi/180, (90 - altitude)*np.pi/180

        # Calculate slope in X and Y directions
        Sx, Sy = self._calcFiniteSlopes()
        #Smag = np.sqrt(Sx**2 + Sy**2)

        # Angle of aspect
        AspectRad = np.arctan2(Sy, Sx)

        # Magnitude of slope in radians
        SmagRad = np.arctan(np.sqrt(Sx**2 + Sy**2))

        self.hillshade = 255.0 * ((np.cos(elevRad) * np.cos(SmagRad)) + \
              (np.sin(elevRad)* np.sin(SmagRad) * np.cos(azRad - AspectRad)))

        return

    def getParams(self):
        """
        Define aspect, gradient and horizontal/vertical curvature using a
        quadratic polynomial method.
        """

        # Assign boundary conditions
        if self.Zbc == None:
            self.Zbc = self._assignBCs()

        # Neighborhood definition
        # z1     z2     z3
        # z4     z5     z6
        # z7     z8     z9
        z1 = self.Zbc[2:, :-2]
        z2 = self.Zbc[2:,1:-1]
        z3 = self.Zbc[2:,2:]
        z4 = self.Zbc[1:-1, :-2]
        z5 = self.Zbc[1:-1,1:-1]
        z6 = self.Zbc[1:-1, 2:]
        z7 = self.Zbc[:-2, :-2]
        z8 = self.Zbc[:-2, 1:-1]
        z9 = self.Zbc[:-2, 2:]

        # Compute coefficient values
        zz = z2+z5
        r = ((z1+z3+z4+z6+z7+z9)-2.*(z2+z5+z8))/(3. * self.dx**2)
        t = ((z1+z2+z3+z7+z8+z9)-2.*(z4+z5+z6))/(3. * self.dx**2)
        s = (z3+z7-z1-z9)/(4. * self.dx**2)
        p = (z3+z6+z9-z1-z4-z7)/(6.*self.dx)
        q = (z1+z2+z3-z7-z8-z9)/(6.*self.dx)
        u = (5.*z1+2.*(z2+z4+z6+z8)-z1-z3-z7-z9)/9.
        #
        with np.errstate(invalid='ignore',divide='ignore'):
            self.grad = np.arctan(np.sqrt(p**2+q**2))
            self.aspect = np.arctan(q/p)
            self.hcurv = -(r*q**2-2.*p*q*s+t*p**2) / \
                ((p**2+q**2)*np.sqrt(1+p**2+q**2))
            self.vcurv = -(r*p**2+2.*p*q*s+t*q**2) /  \
                ((p**2+q**2)*np.sqrt(1+p**2+q**2))

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

        return xsec,ysec


    def viewSection(self, xo = None, yo = None, xm = None, ym = None, pts = 100, vData = None,
                    width = 800, height = 400, color = 'green', linesize = 3,
                    markersize = 5, title = 'Cross section'):
        """
        Extract a slice from the 3D data set and plot required data on a graph.

        Parameters
        ----------

        variable: xo, yo
            Lower X,Y coordinates of the cross-section

        variable: xm, ym
            Upper X,Y coordinates of the cross-section

        variable: pts
            Number of points to discretise the cross-section

        variable: vData
            Dataset to plot.

        variable: width
            Figure width.

        variable: height
            Figure height.

        variable: color
            Color scale.

        variable: linesize, markersize
            Requested size for the line and markers.

        variable: title
            Title of the graph.
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
        rect_B_spline = RectBivariateSpline(self.y[:,0], self.x[0,:], vData)
        datasec = rect_B_spline.ev(ysec, xsec)
        dist = np.sqrt(( xsec - xo )**2 + ( ysec - yo )**2)

        data = Data([
           Scatter(
                x=dist,
                y=datasec,
                mode='lines+markers',
                name="'spline'",
                line=dict(
                    shape='spline',
                    color = color,
                    width = linesize

                ),
                marker = dict(
                    symbol='circle',
                    size = markersize,
                    color = 'white',
                    line = dict(
                        width = 1,
                        color = 'black'
                    )
                )
            )
            ])
        layout = dict(
            title=title,
            width=width,
            height=height
            )

        fig = Figure(data=data, layout=layout)
        plotly.offline.iplot(fig)

        return

    def extractSection(self, xo = None, yo = None, xm = None, ym = None, pts = 100, vData = None,
                    view = True, width = 800, height = 400, color = 'green', linesize = 3,
                    markersize = 5, title = 'Cross section'):
        """
        Extract a slice from the 3D data set and plot required data on a graph.

        Parameters
        ----------

        variable: xo, yo
            Lower X,Y coordinates of the cross-section

        variable: xm, ym
            Upper X,Y coordinates of the cross-section

        variable: pts
            Number of points to discretise the cross-section

        variable: vData
            Dataset to plot.

        variable: view
            Show the section plot.

        variable: width
            Figure width.

        variable: height
            Figure height.

        variable: color
            Color scale.

        variable: linesize, markersize
            Requested size for the line and markers.

        variable: title
            Title of the graph.

        Return:

        variable: dist, datasec
            X, Y values for the profile

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
        rect_B_spline = RectBivariateSpline(self.y[:,0], self.x[0,:], vData)
        datasec = rect_B_spline.ev(ysec, xsec)
        dist = np.sqrt(( xsec - xo )**2 + ( ysec - yo )**2)

        if view:
            data = Data([
               Scatter(
                    x=dist,
                    y=datasec,
                    mode='lines+markers',
                    name="'spline'",
                    line=dict(
                        shape='spline',
                        color = color,
                        width = linesize

                    ),
                    marker = dict(
                        symbol='circle',
                        size = markersize,
                        color = 'white',
                        line = dict(
                            width = 1,
                            color = 'black'
                        )
                    )
                )
                ])
            layout = dict(
                title=title,
                width=width,
                height=height
                )

            fig = Figure(data=data, layout=layout)
            plotly.offline.iplot(fig)

        return dist, datasec

    def profile_mean(self,a):
        return sum(a) / len(a)

    def profile_min(self,a):
        return min(a)

    def profile_max(self,a):
        return max(a)

    def statProfiles(self, pData = None, pDist = None, width = 800, height = 400, color = 'green', linesize = 2,
                    title = 'Section Min, Mean, Max '):
        """
        Plot profile mean, max and min.

        Parameters
        ----------

        variable: pData
            Dataset to plot along Y axis.

        variable: pDist
            Dataset to plot along X axis.

        variable: width
            Figure width.

        variable: height
            Figure height.

        variable: color
            Color scale.

        variable: linesize, markersize
            Requested size for the line and markers.

        variable: title
            Title of the graph.

        Return:

        variable: minZ, meanZ, maxZ
            Y values for the profile (minZ, meanZ, maxZ)
        """

        meanZ = map(self.profile_mean, zip(*pData))
        minZ = map(self.profile_min, zip(*pData))
        maxZ = map(self.profile_max, zip(*pData))

        trace0 = Scatter(
            x=pDist,
            y=maxZ,
            mode='lines',
            line=dict(
                shape='spline',
                width = 0.5,
                color = 'rgb(0, 0, 0)'
            ),
            name='max'
        )

        trace1 = Scatter(
            x=pDist,
            y=minZ,
            mode='lines',
            line=dict(
                shape='spline',
                width = 0.5,
                color = 'rgb(0, 0, 0)'
            ),
            opacity=0.5,
            fill='tonexty',
            fillcolor=color,
            name='min'
        )

        trace2 = Scatter(
            x=pDist,
            y=meanZ,
            mode='lines',
            line=dict(
                shape='spline',
                width = linesize,
                color = 'rgb(0, 0, 0)'
            ),
            name='mean'
        )
        data = [trace0,trace1,trace2]

        layout = dict(
            title=title,
            width=width,
            height=height
        )

        fig = Figure(data=data, layout=layout)
        plotly.offline.iplot(fig)

        return minZ, meanZ, maxZ

    def timeProfiles(self, pData = None, pDist = None, width = 800, height = 400, linesize = 2,
                    title = 'Profile evolution with time'):
        """
        Plot profile mean, max and min.

        Parameters
        ----------

        variable: pData
            Dataset to plot along Y axis.

        variable: pDist
            Dataset to plot along X axis.

        variable: width
            Figure width.

        variable: height
            Figure height.

        variable: color
            Color scale.

        variable: linesize, markersize
            Requested size for the line and markers.

        variable: title
            Title of the graph.

        Return:

        variable: minZ, meanZ, maxZ
            Y values for the profile (minZ, meanZ, maxZ)
        """

        trace = {}
        data = []

        for i in range(0,len(pData)):
            trace[i] = Scatter(
                x=pDist,
                y=pData[i,:],
                mode='lines',
                line=dict(
                    shape='spline',
                    width = linesize,
                    #color = color
                ),
            )
            data.append(trace[i])

        layout = dict(
            title=title,
            width=width,
            height=height
        )

        fig = Figure(data=data, layout=layout)
        plotly.offline.iplot(fig)

    def viewGrid(self, width = 800, height = 800,
                 Dmin = None, Dmax = None, color = None, reverse=False,
                 Data = None, title='Grid'):
        """
        Use Plotly library to visualise a dataset in 2D.

        Parameters
        ----------

        variable: width
            Figure width.

        variable: height
            Figure height.

        variable: Dmin
            Colorbar minimal value.

        variable: Dmax
            Colorbar maximal value.

        variable: color
            Color scale.

        variable: reverse
            Reverse color scale.

        variable: Data
            Dataset to plot.

        variable: title
            Title of the graph.
        """

        if color == None:
            color = 'Picnic'

        data = [
            Heatmap(
                z = Data, colorscale = color,\
                    zmin = Dmin, zmax = Dmax,
                    reversescale=reverse
                )
            ]
        dy = self.bbox[3]-self.bbox[1]
        dx = self.bbox[2]-self.bbox[0]
        if dx>=dy:
            dr = 0.5 * (dx-dy)
            rangeX = [self.bbox[0],self.bbox[2]]
            rangeY = [self.bbox[1]-dr,self.bbox[3]+dr]
        else:
            dr = 0.5 * (dy-dx)
            rangeX = [self.bbox[0]-dr,self.bbox[2]+dr]
            rangeY = [self.bbox[1],self.bbox[3]]

        layout = Layout(
            title=title,
            autosize=True,
            width=width,
            height=height,
            scene=Scene(
                xaxis=XAxis(autorange=False, range=rangeX, nticks=10, \
                            gridcolor='rgb(255, 255, 255)', \
                            gridwidth=2,zerolinecolor='rgb(255, 255, 255)', \
                            zerolinewidth=2),
                yaxis=YAxis(autorange=False, range=rangeY, nticks=10, \
                            gridcolor='rgb(255, 255, 255)', \
                            gridwidth=2,zerolinecolor='rgb(255, 255, 255)', \
                            zerolinewidth=2),
                bgcolor="rgb(244, 244, 248)"
            )
        )

        fig = Figure(data=data, layout=layout)
        plotly.offline.iplot(fig)

        return

    def viewScatter3D(self, width = 800, height = 800, colors='Viridis',
                 dataX = None, dataY = None, dataZ = None, title='Scatter plot'):
        """
        Use Plotly library to visualise a dataset in 3D.

        Parameters
        ----------

        variable: width
            Figure width.

        variable: height
            Figure height.

        variable: colors
            Color scale.

        variable: dataX
            Data for X-axis.

        variable: dataY
            Data for Y-axis.

        variable: dataZ
            Data for Z-axis.

        variable: title
            Title of the graph.
        """

        #trace = {}
        data = []
        #A = np.asarray(dataX) / np.asarray(dataY)
        #A[np.isnan(A)] = 0
        #A[np.isinf(A)] = max(A[A<1000])+1
        trace = Scatter3d(
           x=dataX,
           y=dataY,
           z=dataZ,
           mode='markers',
           marker=dict(
                    size=8,
                    #color=A,
                    #colorscale=colors,
                    opacity=0.8
                )
           )
        data.append(trace)

        layout = dict(
            title=title,
            width=width,
            height=height,
            margin=dict(
                l=0,
                r=0,
                b=0,
                t=0
            ),
            scene=Scene(
                xaxis=XAxis(title='dip'),
                yaxis=YAxis(title='slip'),
                zaxis=ZAxis(title='sed')
            )
        )

        fig = Figure(data=data, layout=layout)
        plotly.offline.iplot(fig)

        return

    def viewScatter(self, width = 800, height = 800,
                 dataX = None, dataY = None, title='Scatter plot'):
        """
        Use Plotly library to visualise a dataset in 2D.

        Parameters
        ----------

        variable: width
            Figure width.

        variable: height
            Figure height.

        variable: dataX
            Data for X-axis.

        variable: dataY
            Data for Y-axis.

        variable: title
            Title of the graph.
        """

        #trace = {}
        data = []

        trace = Scatter(
           x=dataX,
           y=dataY,
           mode='markers',
           )
        data.append(trace)

        layout = dict(
            title=title,
            width=width,
            height=height
        )

        fig = Figure(data=data, layout=layout)
        plotly.offline.iplot(fig)

        return

    def viewSurf(self, width = 800, height = 800,
                 zmin = None, zmax = None, color = None, reverse=False,
                 vData = None, subsample = 1, title='Surface'):
        """
        Use Plotly library to visualise a dataset over a surface in 3D.

        Parameters
        ----------

        variable: width
            Figure width.

        variable: height
            Figure height.

        variable: zmin
            Minimal Z-axis value.

        variable: zmax
            Maximal Z-axis value.

        variable: color
            Color scale.

        variable: reverse
            Reverse color scale.

        variable: vData
            Dataset to plot.

        variable: subsample
            Subsampling data everythin nth value.

        variable: title
            Title of the graph.
        """

        if color == None:
            color = 'YIGnBu'

        if zmin == None:
            zmin = vData.min()

        if zmax == None:
            zmax = vData.max()

        data = Data([
                Surface(
                    x = self.x[::subsample,::subsample],
                    y = self.y[::subsample,::subsample],
                    z = vData[::subsample,::subsample],
                    colorscale = color,
                    reversescale=reverse
                )
            ])

        dy = self.bbox[3]-self.bbox[1]
        dx = self.bbox[2]-self.bbox[0]
        if dx>=dy:
            dr = 0.5 * (dx-dy)
            rangeX = [self.bbox[0],self.bbox[2]]
            rangeY = [self.bbox[1]-dr,self.bbox[3]+dr]
        else:
            dr = 0.5 * (dy-dx)
            rangeX = [self.bbox[0]-dr,self.bbox[2]+dr]
            rangeY = [self.bbox[1],self.bbox[3]]

        layout = Layout(
            title=title,
            autosize=True,
            width=width,
            height=height,
            scene=Scene(
                zaxis=ZAxis(range=[zmin, zmax], \
                            autorange=False,nticks=10, \
                            gridcolor='rgb(255, 255, 255)', \
                            gridwidth=2,zerolinecolor='rgb(255, 255, 255)', \
                            zerolinewidth=2),

                xaxis=XAxis(autorange=False, range=rangeX, \
                            nticks=10, gridcolor='rgb(255, 255, 255)', \
                            gridwidth=2,zerolinecolor='rgb(255, 255, 255)', \
                            zerolinewidth=2),

                yaxis=YAxis(autorange=False, range=rangeY, nticks=10, \
                            gridcolor='rgb(255, 255, 255)', \
                            gridwidth=2,zerolinecolor='rgb(255, 255, 255)', \
                            zerolinewidth=2),

                bgcolor="rgb(244, 244, 248)"
            )
        )

        fig = Figure(data=data, layout=layout)
        plotly.offline.iplot(fig)

        return
