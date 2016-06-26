##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module defines the erodibility and thickness of previously defined stratigraphic layers.
"""
import os
import glob
import time
import h5py
import numpy
import mpi4py.MPI as mpi
from scipy import interpolate
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator

class eroMesh():
    """
    This class builds the erodibility and thickness of underlying initial stratigraphic layers.
    """

    def __init__(self, layNb, eroMap, eroVal, eroTop, thickMap, thickVal, xyTIN,
                 regX, regY, folder, rfolder=None, rstep=0):
        """
        Constructor.

        Parameters
        ----------
        variable: layNb
            Total number of erosion stratigraphic layers

        variable: eroMap
            Erodibility map for each erosion stratigraphic layers

        variable: eroVal
            Erodibility value for each erosion stratigraphic layers

        variable: eroTop
            Erodibility value for reworked sediment

        variable: thickMap
            Thickness map for each erosion stratigraphic layers

        variable: thickVal
            Thickness value for each erosion stratigraphic layers

        variable : xyTIN
            Numpy float-type array containing the coordinates for each nodes in the TIN (in m)

        variable : folder
            Name of the output folder.

        float : regX
            Numpy array containing the X-coordinates of the regular input grid.

        float : regY
            Numpy array containing the Y-coordinates of the regular input grid.

        variable: rfolder, rstep
            Restart folder and step.
        """

        self.regX = regX
        self.regY = regY
        self.layNb = layNb + 1
        nbPts = len(xyTIN[:,0])
        self.folder = folder

        # Build erosion layers
        if rstep > 0:
            if os.path.exists(rfolder):
                folder = rfolder+'/h5/'
            else:
                raise ValueError('The restart folder is missing or the given path is incorrect.')
            df = h5py.File('%s/h5/erolay.time%s.hdf5'%(rfolder, rstep), 'r')
            self.thickness = numpy.array((df['/elayDepth']))
            self.Ke = numpy.array((df['/elayKe']))

            # Get erodibility from erosive layer thicknesses
            self.erodibility = np.zeros(nbPts)
            for k in range(self.layNb):
                existIDs = numpy.where(numpy.logical_and(self.thickness[:,k] > 0., self.erodibility[:,k] == 0.))[0]
                self.erodibility[existIDs] = self.Ke[existIDs,k]
                if(len(numpy.where(self.erodibility == 0)[0]) == 0):
                    break
        else:
            self.thickness = numpy.zeros((nbPts,self.layNb), dtype=float)
            self.Ke = numpy.zeros((nbPts,self.layNb), dtype=float)
            self.thickness[:,0] = 0
            self.Ke[:,0] = eroTop
            for l in range(1,self.layNb):
                if eroMap[l-1] == None:
                    self.Ke[:,l] = eroVal[l-1]
                else:
                    eMap = pandas.read_csv(str(eroMap[l-1]), sep=r'\s+', engine='c',
                                       header=None, na_filter=False, dtype=numpy.float, low_memory=False)
                    reMap = numpy.reshape(eMap.values,(len(self.regX), len(self.regY)), order='F')
                    self.Ke[:,l] = interpolate.interpn( (self.regX, self.regY), reMap, xyTIN, method='nearest')
                if thickMap[l-1] == None:
                    self.thickness[:,l] = thickVal[l-1]
                else:
                    tMap = pandas.read_csv(str(thickMap[l-1]), sep=r'\s+', engine='c',
                                       header=None, na_filter=False, dtype=numpy.float, low_memory=False)
                    rtMap = numpy.reshape(tMap.values,(len(self.regX), len(self.regY)), order='F')
                    self.thickness[:,l] = interpolate.interpn( (self.regX, self.regY), rtMap, xyTIN, method='linear')
            self.erodibility = self.Ke[:,1]
            self.thickness[:,self.layNb-1] += 1.e6

        return

    def getErodibility(self, cumdiff):
        """
        Get the erodibility values for the surface based on underlying erosive stratigraphic layer.

        variable : cumdiff
            Numpy float-type array containing the cumulative erosion/deposition of the nodes in the TIN
        """

        cumThick = numpy.copy(cumdiff)

        # Update deposition
        depIDs = numpy.where(cumThick>=0.)[0]
        self.thickness[depIDs,0] += cumThick[depIDs]

        # Update erosion
        eroIDs = numpy.where(cumThick<0.)[0]
        if len(eroIDs) > 0:
            for k in range(self.layNb):
                # Update thickness for non completely eroded layers
                eIDs = numpy.where(numpy.logical_and(self.thickness[eroIDs,k] > 0., self.thickness[eroIDs,k]+cumThick[eroIDs]>=0.))[0]
                if len(eIDs) > 0:
                    self.thickness[eroIDs[eIDs],k] += cumThick[eroIDs[eIDs]]
                    cumThick[eroIDs[eIDs]] = 0.
                # Nullify completely eroded layer thicknesses and update erosion values
                eIDs = numpy.where(numpy.logical_and(self.thickness[eroIDs,k] > 0., cumThick[eroIDs]<0.))[0]
                if len(eIDs) > 0:
                    cumThick[eroIDs[eIDs]] += self.thickness[eroIDs[eIDs],k]
                    self.thickness[eroIDs[eIDs],k] = 0.
                if(len(numpy.where(cumThick < 0)[0]) == 0):
                    break
                tmpIDs = numpy.where(self.thickness[:,k]<0.)[0]
                self.thickness[tmpIDs,k] = 0.

        # Update surface erodibility from active underlying layers
        self.erodibility = np.zeros(len(cumThick))
        for k in range(self.layNb):
            existIDs = numpy.where(numpy.logical_and(self.thickness[:,k] > 0., self.erodibility[:,k] == 0.))[0]
            self.erodibility[existIDs] = self.Ke[existIDs,k]
            if(len(numpy.where(self.erodibility == 0)[0]) == 0):
                break

        return

    def write_hdf5_erolay(self, outstep):
        """
        This function writes the HDF5 file containing erosive layers information.

        Parameters
        ----------

        variable : outstep
            Output time step.
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        if rank == 0:
            eh5file = self.folder+'/h5/erolay.time'+str(outstep)+'.hdf5'
            ptsNb = len(self.tXY[inIDs,0])
            with h5py.File(eh5file, "w") as f:

                # Write erosive layers depth
                f.create_dataset('elayDepth',shape=(self.ptsNb,self.layNb), dtype='float32', compression='gzip')
                f["elayDepth"][:,:self.layNb] = self.thickness

                # Write erodibility for each layers
                f.create_dataset('elayKe',shape=(self.ptsNb,self.layNb), dtype='float32', compression='gzip')
                f["elayKe"][:,:self.layNb] = self.Ke

        comm.Barrier()

        return
