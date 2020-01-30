##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module defines the **erodibility** and **thickness** of previously defined stratigraphic layers.
"""
import os
import glob
import time
import h5py
import numpy
import pandas
from scipy import interpolate
from scipy.spatial import cKDTree

class eroMesh():
    """
    This class builds the erodibility and thickness of underlying initial stratigraphic layers.

    Args:
        layNb: total number of erosion stratigraphic layers
        eroMap: erodibility map for each erosion stratigraphic layers
        eroVal: erodibility value for each erosion stratigraphic layers
        eroTop: erodibility value for reworked sediment
        thickMap: thickness map for each erosion stratigraphic layers
        thickVal: thickness value for each erosion stratigraphic layers
        xyTIN: numpy float-type array containing the coordinates for each nodes in the TIN (in m)
        regX: numpy array containing the X-coordinates of the regular input grid.
        regY: numpy array containing the Y-coordinates of the regular input grid.
        bPts: boundary points for the TIN.
        ePts: boundary points for the regular grid.
        folder: name of the output folder.
        rfolder: restart folder.
        rstep: restart step.
    """

    def __init__(self, layNb, eroMap, eroVal, eroTop, thickMap, thickVal, xyTIN,
                 regX, regY, bPts, ePts, folder, rfolder=None, rstep=0):

        self.regX = regX
        self.regY = regY
        self.layNb = layNb + 1
        nbPts = len(xyTIN[:,0])
        self.folder = folder

        # Build erosion layers
        # If we restart a simulation
        if rstep > 0:
            if os.path.exists(rfolder):
                folder = rfolder+'/h5/'
            else:
                raise ValueError('The restart folder is missing or the given path is incorrect.')
            df = h5py.File('%s/h5/erolay.time%s.hdf5'%(rfolder, rstep), 'r')
            self.thickness = numpy.array((df['/elayDepth']))
            self.Ke = numpy.array((df['/elayKe']))

            # Get erodibility from erosive layer thicknesses
            self.erodibility = numpy.zeros(nbPts)
            for k in range(self.layNb):
                existIDs = numpy.where(numpy.logical_and(self.thickness[:,k] > 0., self.erodibility[:,k] == 0.))[0]
                self.erodibility[existIDs] = self.Ke[existIDs,k]
                if(len(numpy.where(self.erodibility == 0)[0]) == 0):
                    break

        # Build the underlying erodibility mesh and associated thicknesses
        else:
            # Initial top layer (id=0) is for reworked sediment (freshly deposited)
            self.thickness = numpy.zeros((nbPts,self.layNb), dtype=float)
            self.Ke = numpy.zeros((nbPts,self.layNb), dtype=float)
            self.thickness[:,0] = 0
            self.Ke[:,0] = eroTop

            # Define inside area kdtree
            inTree = cKDTree(xyTIN[bPts:ePts+bPts,:])
            dist, inID = inTree.query(xyTIN[:bPts,:],k=1)
            inID += bPts

            # Loop through the underlying layers
            for l in range(1,self.layNb):
                # Uniform erodibility value
                if eroMap[l-1] == None:
                    self.Ke[:,l] = eroVal[l-1]
                # Erodibility map
                else:
                    eMap = pandas.read_csv(str(eroMap[l-1]), sep=r'\s+', engine='c',
                                       header=None, na_filter=False, dtype=numpy.float, low_memory=False)
                    reMap = numpy.reshape(eMap.values,(len(self.regX), len(self.regY)), order='F')
                    self.Ke[bPts:,l] = interpolate.interpn( (self.regX, self.regY), reMap, xyTIN[bPts:,:], method='nearest')
                    # Assign boundary nodes
                    tmpK = self.Ke[bPts:,l]
                    self.Ke[:bPts,l] = tmpK[inID]

                # Uniform thickness value
                if thickMap[l-1] == None:
                    self.thickness[:,l] = thickVal[l-1]
                # Thickness map
                else:
                    tMap = pandas.read_csv(str(thickMap[l-1]), sep=r'\s+', engine='c',
                                       header=None, na_filter=False, dtype=numpy.float, low_memory=False)
                    rtMap = numpy.reshape(tMap.values,(len(self.regX), len(self.regY)), order='F')
                    self.thickness[bPts:,l] = interpolate.interpn( (self.regX, self.regY), rtMap, xyTIN[bPts:,:], method='linear')
                    # Assign boundary nodes
                    tmpH = self.thickness[bPts:,l]
                    self.thickness[:bPts,l] = tmpH[inID]

            # Define active layer erodibility
            self.erodibility = numpy.zeros(nbPts)
            for l in range(1,self.layNb):
                # Get erodibility coefficients from active underlying layers
                tmpIDs = numpy.where(numpy.logical_and(self.thickness[:,l] > 0., self.erodibility[:] == 0.))[0]
                self.erodibility[tmpIDs] = self.Ke[tmpIDs,l]
                if(len(numpy.where(self.erodibility == 0)[0]) == 0):
                    break

            # Bottom layer is supposed to be infinitely thick
            self.thickness[:,self.layNb-1] += 1.e6

        return

    def getErodibility(self, cumThick):
        """
        Get the erodibility values for the surface based on underlying erosive stratigraphic layer.

        Args:
            cumThick: numpy float-type array containing the cumulative erosion/deposition of the nodes in the TIN
        """

        # Update deposition
        depIDs = numpy.where(cumThick>=0.)[0]
        self.thickness[depIDs,0] += cumThick[depIDs]

        # Update erosion
        eroIDs = numpy.where(cumThick<0.)[0]
        if len(eroIDs) > 0:
            for k in range(self.layNb):
                # Update thickness for remaining layers
                eIDs = numpy.where(numpy.logical_and(self.thickness[eroIDs,k] > 0., self.thickness[eroIDs,k] >= -cumThick[eroIDs]))[0]
                if len(eIDs) > 0:
                    self.thickness[eroIDs[eIDs],k] += cumThick[eroIDs[eIDs]]
                    cumThick[eroIDs[eIDs]] = 0.

                # Nullify eroded layer thicknesses and update erosion values
                eIDs = numpy.where(numpy.logical_and(self.thickness[eroIDs,k] > 0., cumThick[eroIDs] < 0.))[0]
                if len(eIDs) > 0:
                    cumThick[eroIDs[eIDs]] += self.thickness[eroIDs[eIDs],k]
                    self.thickness[eroIDs[eIDs],k] = 0.

                # Ensure non-negative values
                tmpIDs = numpy.where(self.thickness[:,k] < 0.)[0]
                if len(tmpIDs) > 0:
                    self.thickness[tmpIDs,k] = 0.

                if(len(numpy.where(cumThick < 0)[0]) == 0):
                    break

        # Update surface erodibility map
        self.erodibility = numpy.zeros(len(cumThick))
        for k in range(self.layNb):
            # Get erodibility coefficients from active underlying layers
            tmpIDs = numpy.where(numpy.logical_and(self.thickness[:,k] > 0., self.erodibility[:] == 0.))[0]
            self.erodibility[tmpIDs] = self.Ke[tmpIDs,k]

            if(len(numpy.where(self.erodibility == 0)[0]) == 0):
                break

        return

    def write_hdf5_erolay(self, outstep):
        """
        This function writes the HDF5 file containing erosive layers information.

        Args:
            outstep: output time step.
        """

        eh5file = self.folder+'/h5/erolay.time'+str(outstep)+'.hdf5'
        ptsNb = len(self.erodibility)
        with h5py.File(eh5file, "w") as f:

            # Write erosive layers depth
            f.create_dataset('elayDepth',shape=(ptsNb,self.layNb), dtype='float64', compression='gzip')
            f["elayDepth"][:,:self.layNb] = self.thickness

            # Write erodibility for each layers
            f.create_dataset('elayKe',shape=(ptsNb,self.layNb), dtype='float64', compression='gzip')
            f["elayKe"][:,:self.layNb] = self.Ke

        return
