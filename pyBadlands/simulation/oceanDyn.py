##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module encapsulates functions related to hydrodynamic calculation based on SWAN model.
"""

import time
import math
import numpy
import mpi4py.MPI as mpi
from scipy import interpolate
from scipy.spatial import cKDTree
from collections import OrderedDict
from matplotlib import _cntr as cntr
from scipy.interpolate import interpn
from scipy.ndimage.filters import gaussian_filter


import matplotlib.pyplot as plt
import cmocean as cmo
from matplotlib import cm
from pylab import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

import pyBadlands.libUtils.simswan as swan


class oceanDyn():

    def __init__(self, res, xv, yv):
        """
        Constructor.
        """

        self.rank = mpi.COMM_WORLD.rank
        self.size = mpi.COMM_WORLD.size
        self.comm = mpi.COMM_WORLD
        self.fcomm = mpi.COMM_WORLD.py2f()

        self.gravity = 9.81

        self.xv = xv
        self.yv = yv
        self.xi, self.yi = numpy.meshgrid(xv, yv, sparse=False, indexing='ij')
        xx = numpy.ravel(self.xi,order='F')
        yy = numpy.ravel(self.yi,order='F')
        self.xyi = numpy.vstack((xx, yy)).T

        self.res = res

        self.wnx = int(round((xx.max()-xx.min())/self.res+1))
        self.wny = int(round((yy.max()-yy.min())/self.res+1))
        self.wavX = numpy.linspace(xx.min(),xx.max(),self.wnx)
        self.wavY = numpy.linspace(yy.min(),yy.max(),self.wny)

        self.xyTIN = None
        self.tree= None

        return

    def build_tree(self, xyTIN):
        """
        Update wave mesh.

        Parameters
        ----------
        variable : xyTIN
            Numpy float-type array containing the coordinates for each nodes in the TIN (in m)
        """

        # Update TIN grid kdtree for interpolation
        self.xyTIN = xyTIN
        self.tree = cKDTree(xyTIN)
        tindx = xyTIN[1,0] - xyTIN[0,0]
        self.searchpts = max(int(self.res*self.res/(tindx*tindx)),4)

        return

    def swan_init(self, input, z, wID, sl):
        """
        This function initialise swan model.

        Parameters
        ----------

        variable : input
            Simulation input parameters time step.

        variable: z
            Elevation.

        variable: wID
             Wave number ID.

        variable: sl
            Sealevel elevation.
        """

        # Interpolate TIN elevation on wave mesh
        distances, indices = self.tree.query(self.xyi, k=self.searchpts)
        z_vals = z[indices]

        with numpy.errstate(invalid='ignore'):
            nz = numpy.average(z_vals,weights=(1./distances), axis=1)

        onIDs = numpy.where(distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            nz[onIDs] = z[indices[onIDs,0]]

        rZ = numpy.reshape(nz,(self.wnx,self.wny),order='F')

        wl = input.wavelist[wID]
        cl = input.climlist[wID]

        wD1 = 0
        wD2 = 0
        if self.rank == 0:
            wD1 = input.waveWd[wl][cl] + 2.*input.waveWs[wl][cl]*(numpy.random.rand()-0.5)
            wD2 = input.waveWdd[wl][cl] + 2.*input.waveWs[wl][cl]*(numpy.random.rand()-0.5)
        wD1 = self.comm.bcast(wD1, root=0)
        wD2 = self.comm.bcast(wD2, root=0)

        swan.model.init(self.fcomm, input.swanFile, input.swanInfo,
                        input.swanBot, input.swanOut, rZ, input.waveWh[wl][cl],
                        input.waveWp[wl][cl], wD1, input.waveWu[wl][cl],
                        wD2, input.waveSide[wl][cl], self.res, input.waveBase, sl)

        return

    def swan_run(self, input, force, TINz, wID):
         """
         This function run swan model and computes near-bed velocity
         combining cross-shore and long-shore component.

         Parameters
         ----------

         variable : input
            Simulation input parameters time step.

         variable: force
            Forcing conditions.

         variable: TINz
            TIN elevation.

         variable: WID
             Wave number ID.
         """

         force.wavU = []
         force.wavV = []
         force.wavH = []
         force.wavPerc = []

         # Interpolate TIN elevation on wave mesh
         distances, indices = self.tree.query(self.xyi, k=self.searchpts)
         z_vals = TINz[indices]

         with numpy.errstate(divide='ignore'):
             nz = numpy.average(z_vals,weights=(1./distances), axis=1)

         onIDs = numpy.where(distances[:,0] == 0)[0]
         if len(onIDs) > 0:
             nz[onIDs] = TINz[indices[onIDs,0]]

         rZ = numpy.reshape(nz,(self.wnx,self.wny),order='F')

         # Loop through the different wave climates and store swan output information
         force.wclim = input.climNb[input.wavelist[wID]]

         for clim in range(force.wclim):

            # Define next wave regime
            tw = time.clock()
            wID += 1
            wl = input.wavelist[wID]
            cl = input.climlist[wID]

            wD1 = 0
            wD2 = 0
            if self.rank == 0:
                wD1 = input.waveWd[wl][cl] + 2.*input.waveWs[wl][cl]*(numpy.random.rand()-0.5)
                wD2 = input.waveWdd[wl][cl] + 2.*input.waveWs[wl][cl]*(numpy.random.rand()-0.5)
            wD1 = self.comm.bcast(wD1, root=0)
            wD2 = self.comm.bcast(wD2, root=0)

            # Run SWAN model
            wavU, wavD, H = swan.model.run(self.fcomm, rZ, input.waveWh[wl][cl], input.waveWp[wl][cl], wD1,
                                            input.waveWu[wl][cl], wD2, input.waveSide[wl][cl], force.sealevel)

            # Define velocity current
            fU = gaussian_filter(wavU, sigma=1.)
            fH = gaussian_filter(H, sigma=1.)
            U = fU * numpy.cos(wavD)
            V = fU * numpy.sin(wavD)

            # Plot the wave dynamic using matplotlib
            #self._bottomCurrents(U, V, force.sealevel, rZ)

            # Convert velocity from wave mesh to TIN
            cU = interpn( (self.wavX, self.wavY), U, (self.xyTIN), method='linear')
            cV = interpn( (self.wavX, self.wavY), V, (self.xyTIN), method='linear')
            cH = interpn( (self.wavX, self.wavY), fH, (self.xyTIN), method='linear')

            if self.rank == 0:
                print 'Swan model of waves field %d and climatic conditions %d:' %(wl,clim)
                print 'took %0.02f seconds to run.' %(time.clock()-tw)

            landIDs = numpy.where(TINz>=force.sealevel)[0]
            cU[landIDs] = 0.
            cV[landIDs] = 0.
            cH[landIDs] = 0.

            # Save computed velocity for checking
            #self._dumpData(cU,cV,cH)

            # Store percentage of each climate and induced bottom currents velocity
            force.wavU.append(cU)
            force.wavV.append(cV)
            force.wavH.append(cH)
            force.wavPerc.append(input.wavePerc[wl][cl])

         return int(wID)

    def _dumpData(self, U, V, H):
        """
        Visualise wave induced bottom current streamlines based on SWAN wave model.
        """

        tmp = numpy.column_stack((self.xyTIN,U))
        arr = numpy.column_stack((tmp,V))
        arr1 = numpy.column_stack((arr,H))
        numpy.savetxt("wavedata.csv", arr1, delimiter=",")

        return

    def _bottomCurrents(self, U, V, sl, rZ, gauss=0, dens=6., color=cmo.cm.matter, fsize=(12,12)):
        """
        Visualise wave induced bottom current streamlines based on SWAN wave model.
        Parameters
        ----------
        variable : U, V
            Component of bottom velocity along A and Y axis.
        variable : gauss
            Gaussian filter for smoothing bottom current value.
        variable : dens
            Streamline density.
        variable : Color
            color map from cmocean.
        variable : size
            Plot size.
        variable : fname
            Save PNG filename.
        variable : dpi
            Figure resolution.
        variable : r1,r2,c1,c2
            Zoom to a specific area by defining row/col ids.
        """

        fig = plt.figure(figsize = fsize)
        ax = plt.subplot(1, 1, 1)

        speed = numpy.sqrt(U**2+V**2)

        dataExtent = [numpy.amin(self.wavX), numpy.amax(self.wavX),
                    numpy.amin(self.wavY), numpy.amax(self.wavY)]
        # Current map
        if gauss > 0:
            svel = gaussian_filter(speed, gauss)
            im = ax.imshow(svel.T, interpolation = 'bilinear', cmap=color, vmin=0, vmax=0.3, extent=dataExtent, origin='lower')
        else:
            im = ax.imshow(speed.T, interpolation = 'bilinear', cmap=color, vmin=0, vmax=0.3,extent=dataExtent, origin='lower')

        # Colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        plt.colorbar(im, cax=cax)

        # Shoreline contour
        levels = [0]
        CS = ax.contour(rZ.T-sl, levels, colors='k', extent=dataExtent, linewidths=1.2, origin='lower')

        # Streamlines
        lw = 2.*speed.T / speed.max()
        # print (self.wavX[c1:c2]-self.wavX[c1]).shape
        # print lw.shape,U[c1:c2,r1:r2].shape
        strm = ax.streamplot(self.wavX, self.wavY,
                    U.T, V.T, color='k', density=dens, linewidth=0.5, arrowsize=0.5)

        plt.show()

        return
