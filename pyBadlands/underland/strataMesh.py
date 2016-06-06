##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module defines the stratigraphic layers based on the TIN nodes.
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
from pyBadlands.libUtils import FASTloop

class strataMesh():
    """
    This class builds stratigraphic layer on each depositional point of the regular mesh.

    Parameters
    ----------
    variable: stratIn
        Numpy array flaging the presence of a stratigraphic layer for each node

    variable: stratElev
        Numpy array containing the relative elevation of the layer at the time of deposition

    variable: stratThick
        Numpy array containing the thickness of each stratigraphic layer

    variable: stratDepth
        Numpy array containing the current depth of each stratigraphic layer
    """

    def __init__(self, sdx, bbX, bbY, layNb, RowProc, ColProc, xyTIN, folder,
                 h5file, cumdiff=0, rfolder=None, rstep=0):
        """
        Constructor.

        Parameters
        ----------
        variable: sdx
            Discretisation value [m]

        variable: bbX
            Extent of stratal regular grid along X

        variable: bbY
            Extent of stratal regular grid along Y

        variable: layNb
            Total number of stratigraphic layers

        variable: RowProc
            Number of processor partition along the X axis

        variable: ColProc
            Number of processor partition along the Y axis

        variable : xyTIN
            Numpy float-type array containing the coordinates for each nodes in the TIN (in m)

        variable : folder
            Name of the output folder.

        variable: h5file
            First part of the hdf5 file name.

        variable: cumdiff
            Numpy array containing  cumulative erosion/deposition from previous simulation.

        variable: rfolder, rstep
            Restart folder and step.
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        self.ids = None
        self.ptsNb = None
        self.prevload = 0.
        self.tree = None
        self.folder = folder
        self.h5file = h5file
        self.step = 0

        # User defined parameter
        self.dx = sdx

        # Create stratal regular grid
        self.nx = round((bbX[1]-bbX[0])/sdx - 0.5)+1
        self.ny = round((bbY[1]-bbY[0])/sdx - 0.5)+1
        xgrid = numpy.linspace(bbX[0],bbX[1],num=self.nx)
        ygrid = numpy.linspace(bbY[0],bbY[1],num=self.ny)
        xi, yi = numpy.meshgrid(xgrid, ygrid)
        self.xyi = numpy.dstack([xi.flatten(), yi.flatten()])[0]

        # Partition mesh
        self.buildPartition(bbX, bbY, RowProc, ColProc)
        self.ptsNb = len(self.ids)

        if rstep > 0:
            if os.path.exists(rfolder):
                folder = rfolder+'/h5/'
                fileCPU = 'sed.time%s.p*.hdf5'%rstep
                restartncpus = len(glob.glob1(folder,fileCPU))
                if restartncpus == 0:
                    raise ValueError('The requested time step for the restart simulation cannot be found in the restart folder.')
            else:
                raise ValueError('The restart folder is missing or the given path is incorrect.')

            if restartncpus != size:
                raise ValueError('When using the stratal model you need to run the restart simulation with the same number of processors as the previous one.')

            df = h5py.File('%s/h5/sed.time%s.p%s.hdf5'%(rfolder, rstep, rank), 'r')
            layDepth = numpy.array((df['/layDepth']))
            layElev = numpy.array((df['/layElev']))
            layThick = numpy.array((df['/layThick']))
            rstlays = layDepth.shape[1]
            layNb +=  rstlays
            self.step = rstlays

        # Define global stratigraphic dataset
        self.stratIn = numpy.zeros([self.ptsNb],dtype=int)
        self.stratElev = numpy.zeros([self.ptsNb,layNb])
        self.stratThick = numpy.zeros([self.ptsNb,layNb])
        self.stratDepth = numpy.zeros([self.ptsNb,layNb])

        if rstep > 0:
            self.stratDepth[:,:rstlays] = layDepth
            self.stratElev[:,:rstlays] = layElev
            self.stratThick[:,:rstlays] = layThick

        # Define TIN grid kdtree for interpolation
        self.xyTIN = xyTIN
        self.tree = cKDTree(self.xyTIN)
        tindx = self.xyTIN[1,0] - self.xyTIN[0,0]
        self.searchpts = max(int(sdx*sdx/(tindx*tindx)),4)

        if rstep > 0:
            scumload = numpy.zeros(len(self.xyi))
            distances, indices = self.tree.query(self.xyi, k=self.searchpts)

            if len(cumdiff[indices].shape) == 3:
                cum_vals = cumdiff[indices][:,:,0]
            else:
                cum_vals = cumdiff[indices]
            fcum = numpy.average(cum_vals,weights=(1./distances), axis=1)
            onIDs = numpy.where(distances[:,0] == 0)[0]
            if len(onIDs) > 0:
                fcum[onIDs] = cumdiff[indices[onIDs,0]]
            self.prevload = fcum

        return

    def update_stratal_parameters(self, xyTIN, cumdiff, elev):
        """
        Update TIN variables after 3D displacements.

        variable : xyTIN
            Numpy float-type array containing the coordinates for each nodes in the TIN (in m)
        """

        self.xyTIN = xyTIN
        self.tree = cKDTree(self.xyTIN)

        # Update elevation and cumlative change

        return

    def buildStrata(self, elev, cumdiff, sea, rank, write=0, outstep=0):
        """
        Build the stratigraphic layer on the regular grid.

        variable : elev
            Numpy float-type array containing the elevation of the nodes in the TIN

        variable : cumdiff
            Numpy float-type array containing the cumulative erosion/deposition of the nodes in the TIN

        variable : sea
            Sea level elevation

        variable : rank
            Rank of the given processor

        variable : write
            Flag for output generation

        variable : outstep
            Step for output generation
        """

        selev = numpy.zeros(len(self.xyi))
        scumload = numpy.zeros(len(self.xyi))
        distances, indices = self.tree.query(self.xyi, k=self.searchpts)

        if len(elev[indices].shape) == 3:
            elev_vals = elev[indices][:,:,0]
            cum_vals = cumdiff[indices][:,:,0]
        else:
            elev_vals = elev[indices]
            cum_vals = cumdiff[indices]

        felev = numpy.average(elev_vals,weights=(1./distances), axis=1)
        fcum = numpy.average(cum_vals,weights=(1./distances), axis=1)
        onIDs = numpy.where(distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            felev[onIDs] = elev[indices[onIDs,0]]
            fcum[onIDs] = cumdiff[indices[onIDs,0]]
        scumload = fcum - self.prevload
        self.prevload = fcum
        selev = felev
        # Update stratal elevation
        self.stratElev[self.ids,self.step] =  selev[self.ids]-sea

        # Update stratal deposition
        localCum = scumload[self.ids]
        depIDs = numpy.where(localCum>0.)[0]
        self.depoLayer(self.ids[depIDs], localCum)

        # Update stratal erosion
        eroIDs = numpy.where(localCum<0.)[0]
        self.eroLayer(self.ids[eroIDs], localCum)

        if write>0:
            self.layerMesh(selev[self.ids])
            self.write_hdf5_stratal(outstep,rank)

        self.step += 1

        return

    def buildPartition(self, bbX, bbY, RowProc, ColProc):
        """
        Define a partition for the stratal mesh.

        Parameters
        ----------
        variable: ptsNb
            Number of points in the regular grid

        variable: layNb
            Total number of stratigraphic layers
        """

        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # extent of X partition
        Xstart = numpy.zeros( RowProc )
        Xend = numpy.zeros( RowProc )
        for p in range(RowProc):
            if p == 0:
                Xstart[p] = p*self.nx + bbX[0]
                Xend[p] = Xstart[p] + self.nx + self.dx
            else:
                Xstart[p] = p*self.nx + bbX[0] - self.dx
                Xend[p] = Xstart[p] + self.nx + 2 * self.dx
        Xend[RowProc-1] = bbX[1]

        # extent of Y partition
        Ystart = numpy.zeros( ColProc )
        Yend = numpy.zeros( ColProc )
        for p in range(ColProc):
            if p == 0:
                Ystart[p] = p*self.ny + bbY[0]
                Yend[p] = Ystart[p] + self.ny + self.dx
            else:
                Ystart[p] = p*self.ny + bbY[0] - self.dx
                Yend[p] = Ystart[p] + self.ny + 2 * self.dx
        Yend[ColProc-1] = bbY[1]

        # Define partitions ID globally
        size = ColProc * RowProc
        Xst = numpy.zeros( size )
        Xed = numpy.zeros( size )
        Yst = numpy.zeros( size )
        Yed = numpy.zeros( size )
        for q in range(RowProc):
            for p in range(ColProc):
                Xst[q + p * RowProc] = Xstart[q]
                Yst[q + p * RowProc] = Ystart[p]
                Xed[q + p * RowProc] = Xend[q]
                Yed[q + p * RowProc] = Yend[p]

        # Loop over node coordinates and find if they belong to local partition
        # Note: used a Cython/Fython class to increase search loop performance... in libUtils
        partID = FASTloop.part.overlap(self.xyi[:,0],self.xyi[:,1],Xst[rank],
                                        Yst[rank],Xed[rank],Yed[rank])

        # Extract local domain nodes global ID
        self.ids = numpy.where(partID > -1)[0]

        return

    def depoLayer(self, ids, depo):
        """
        Add deposit to current stratigraphic layer.

        Parameters
        ----------
        variable: ids
            Index of points subject to deposition

        variable: depo
            Value of the deposition for the given point [m]
        """

        # Initialise node deposition flag
        tmpIDs = numpy.where(self.stratIn[ids]==0)
        self.stratIn[ids[tmpIDs]] = 1

        # Add deposit to the considered layer time
        self.stratThick[ids,self.step] += depo[ids]

        return

    def eroLayer(self, nids, erosion):
        """
        Erode top stratigraphic layers.

        Parameters
        ----------
        variable: nids
            Index of points subject to erosion

        variable: erosion
            Value of the erosion for the given points [m]
        """

        # Perform erosion on nodes containing stratigraphic layers
        tmpIDs = numpy.where(self.stratIn[nids] == 1)[0]
        if len(tmpIDs) == 0:
            return

        # Update node indices and associated erosion values
        ids = nids[tmpIDs]
        ero = -erosion[ids]

        # Compute cumulative stratal thicknesses
        cumThick = numpy.cumsum(self.stratThick[ids,self.step::-1],axis=1)[:,::-1]

        # Find nodes with no remaining stratigraphic thicknesses
        tmpIDs = numpy.where(ero>=cumThick[:,0])[0]
        self.stratIn[ids[tmpIDs]] = 0
        self.stratThick[ids[tmpIDs],:self.step+1] = 0.

        # Erode remaining stratal layers
        if len(tmpIDs) < len(ids):
            ero[tmpIDs] = 0.

            # Clear all stratigraphy points which are eroded
            cumThick[cumThick < ero.reshape((len(ids),1))] = 0
            mask = (cumThick > 0).astype(int) == 0
            tmpH = self.stratThick[ids,:self.step+1]
            tmpH[mask] = 0
            self.stratThick[ids,:self.step+1] = tmpH

            # Update thickness of top stratigraphic layer
            eroIDs = numpy.bincount(numpy.nonzero(cumThick)[0]) - 1
            eroVals = cumThick[numpy.arange(len(ids)),eroIDs]-ero
            eroVals[tmpIDs] = 0.
            self.stratThick[ids,eroIDs] = eroVals

        return

    def layerMesh(self, topsurf):
        """
        Define stratigraphic layers mesh.

        Parameters
        ----------

        variable: topsurf
            Elevation of the regular surface
        """

        # Clear points with no stratigraphic layer
        tmpIDs = numpy.where(self.stratIn == 0)[0]
        #surf = numpy.tile(topsurf[tmpIDs].transpose(), (1, self.step+1)).reshape(self.step+1,len(tmpIDs)).transpose()
        surf = numpy.array([topsurf[tmpIDs],]*int(self.step+1)).transpose()
        self.stratDepth[tmpIDs,:self.step+1] = surf
        self.stratThick[tmpIDs,:self.step+1] = 0.

        # Find points with stratigraphic layers
        tmpIDs = numpy.where(self.stratIn == 1)[0]
        if len(tmpIDs) == 0:
            return

        # Compute cumulative stratal thicknesses
        cumThick = numpy.cumsum(self.stratThick[tmpIDs,self.step::-1],axis=1)[:,::-1]

        # Updata stratal depth
        #surf = numpy.tile(topsurf[tmpIDs].transpose(), (1, self.step+1)).reshape(self.step+1,len(tmpIDs)).transpose()
        surf = numpy.array([topsurf[tmpIDs],]*int(self.step+1)).transpose()
        self.stratDepth[tmpIDs,:self.step+1] = surf - cumThick

        return

    def write_hdf5_stratal(self, outstep, rank):
        """
        This function writes for each processor the HDF5 file containing sub-surface information.

        Parameters
        ----------

        variable : outstep
            Output time step.

        variable : rank
            ID of the local partition.
        """

        sh5file = self.folder+'/'+self.h5file+str(outstep)+'.p'+str(rank)+'.hdf5'
        with h5py.File(sh5file, "w") as f:

            # Write node coordinates
            f.create_dataset('coords',shape=(self.ptsNb,2), dtype='float32', compression='gzip')
            f["coords"][:,:2] = self.xyi[self.ids]

            # Write stratal layers depth per cells
            f.create_dataset('layDepth',shape=(self.ptsNb,self.step+1), dtype='float32', compression='gzip')
            f["layDepth"][:,:self.step+1] = self.stratDepth[:,:self.step+1]

            # Write stratal layers elevations per cells
            f.create_dataset('layElev',shape=(self.ptsNb,self.step+1), dtype='float32', compression='gzip')
            f["layElev"][:,:self.step+1] = self.stratElev[:,:self.step+1]

            # Write stratal layers thicknesses per cells
            f.create_dataset('layThick',shape=(self.ptsNb,self.step+1), dtype='float32', compression='gzip')
            f["layThick"][:,:self.step+1] = self.stratThick[:,:self.step+1]

        return

'''
class strataMesh():
    """
    This class builds stratigraphic layer on each depositional point of the unstructured mesh.

    Parameters
    ----------
    string : timeID
        Numpy array containing the stratigraphic time layer index

    string : elev
        Numpy array containing the relative elevation of the layer at the time of deposition

    variable: thick
        Numpy array containing the thickness of each stratigraphic layer
    """

    def __init__(self, timeslice, ):
        self.timeID = numpy.array([],dtype=int)
        self.elev = numpy.array([])
        self.thick = numpy.array([])


    def add2layer(self, step, elev, depo):
        """
        Add deposit to top stratigraphic layer.

        Parameters
        ----------
        variable: step
            Time step of the running simulation

        variable: elev
            Initial elevation of the layer deposit relative to sea level [m]

        variable: depo
            Value of the deposition for the given point [m]
        """

        # First time we store some deposits on this node
        if self.timeID.size == 0:
            self.timeID = numpy.append(self.timeID, step)
            self.elev = numpy.append(self.elev, elev)
            self.thick = numpy.append(self.thick, depo)
        else:
            # New stratigraphic layer deposit
            if self.timeID[-1] < step:
                self.timeID = numpy.append(self.timeID, step)
                self.elev = numpy.append(self.elev, elev)
                self.thick = numpy.append(self.thick, depo)
            # Add to existing stratigraphic layer
            else:
                self.thick[-1] += depo

        return

    def take2layer(self, ero):
        """
        Erode top stratigraphic layer.

        Parameters
        ----------
        variable: ero
            Value of the erosion for the given point [m]

        Return
        ----------
        variable: deleteStrat
            Flag to delete completely eroded stratal layer
        """

        deleteStrat = False

        if ero <= self.thick[-1]:
            self.thick[-1] -= ero
        elif ero >= numpy.sum(self.thick):
            deleteStrat = True
            self.timeID = numpy.array([],dtype=int)
            self.elev = numpy.array([])
            self.thick = numpy.array([])
        else:
            cumthick = numpy.cumsum(self.thick[::-1])[::-1]
            eroLay = numpy.where(cumthick<ero)[0]
            self.thick = numpy.delete(self.thick,eroLay)
            self.timeID = numpy.delete(self.timeID,eroLay)
            self.elev = numpy.delete(self.elev,eroLay)
            stillero = ero - cumthick[eroLay[0]]
            if stillero > 0.:
                self.thick[-1] -= ero

        if self.thick.size > 0:
            if self.thick[-1] < 1.e-6:
                size = self.thick.size
                if size == 1:
                    deleteStrat = True
                    self.timeID = numpy.array([],dtype=int)
                    self.elev = numpy.array([])
                    self.thick = numpy.array([])
                else:
                    self.thick = numpy.delete(self.thick,[size-1])
                    self.timeID = numpy.delete(self.timeID,[size-1])
                    self.elev = numpy.delete(self.elev,[size-1])

        return deleteStrat
'''
