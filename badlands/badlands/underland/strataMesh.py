##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module defines the **stratigraphic layers** based on a **regular mesh**.
"""
import os
import glob
import time
import h5py
import numpy
from scipy import interpolate
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator

if 'READTHEDOCS' not in os.environ:
    from badlands import flowalgo

class strataMesh():
    """
    This class builds stratigraphic layer on each depositional point of the regular mesh.

    Args:
        sdx: discretisation value [m]
        bbX: extent of stratal regular grid along X
        bbY: extent of stratal regular grid along Y
        layNb: total number of stratigraphic layers
        xyTIN: numpy float-type array containing the coordinates for each nodes in the TIN (in m)
        folder: name of the output folder.
        h5file: first part of the hdf5 file name.
        poro0: initial porosity.
        poroC: porosity evolution coefficient.
        cumdiff: numpy array containing  cumulative erosion/deposition from previous simulation.
        rfolder: restart folder.
        rstep: restart step.
    """

    def __init__(self, sdx, bbX, bbY, layNb, xyTIN, folder, h5file, poro0, poroC,
                 cumdiff=0, rfolder=None, rstep=0):

        self.ids = None
        self.ptsNb = None
        self.oldload = None
        self.tree = None
        self.folder = folder
        self.h5file = h5file+'.time'
        self.step = 0
        self.upper = None
        self.lower = None
        self.poro0 = poro0
        self.poroC = poroC
        self.xyTIN = xyTIN

        # User defined parameter
        self.dx = sdx

        # Create stratal regular grid
        self.nx = int(round((bbX[1]-bbX[0])/sdx - 0.5)+1)
        self.ny = int(round((bbY[1]-bbY[0])/sdx - 0.5)+1)
        self.xgrid = numpy.linspace(bbX[0],bbX[1],num=self.nx)
        self.ygrid = numpy.linspace(bbY[0],bbY[1],num=self.ny)
        xi, yi = numpy.meshgrid(self.xgrid, self.ygrid)
        self.xyi = numpy.dstack([xi.flatten(), yi.flatten()])[0]

        # Partition mesh
        self.buildPartition(bbX, bbY)
        self.ptsNb = len(self.ids)

        if rstep > 0:
            if os.path.exists(rfolder):
                folder = rfolder+'/h5/'
                fileCPU = 'sed.time%s.hdf5'%(rstep)
                restartncpus = len(glob.glob1(folder,fileCPU))
                if restartncpus == 0:
                    raise ValueError('The requested time step for the restart simulation cannot be found in the restart folder.')
            else:
                raise ValueError('The restart folder is missing or the given path is incorrect.')

            if restartncpus != 1:
                raise ValueError('When using the stratal model you need to run the restart simulation with the same number of processors as the previous one.')

            df = h5py.File('%s/h5/sed.time%s.hdf5'%(rfolder, rstep), 'r')
            layDepth = numpy.array((df['/layDepth']))
            layElev = numpy.array((df['/layElev']))
            layThick = numpy.array((df['/layThick']))
            layPoro = numpy.array((df['/layPoro']))
            rstlays = layDepth.shape[1]
            layNb +=  rstlays
            self.step = rstlays

        # Define global stratigraphic dataset
        self.stratIn = numpy.zeros([self.ptsNb],dtype=int)
        self.stratElev = numpy.zeros([self.ptsNb,layNb])
        self.stratThick = numpy.zeros([self.ptsNb,layNb])
        self.stratPoro = numpy.zeros([self.ptsNb,layNb])
        self.stratDepth = numpy.zeros([self.ptsNb,layNb])
        self.stratPoro.fill(self.poro0)

        if rstep > 0:
            self.stratDepth[:,:rstlays] = layDepth
            self.stratElev[:,:rstlays] = layElev
            self.stratThick[:,:rstlays] = layThick
            self.stratPoro[:,:rstlays] = layPoro

        # Define TIN grid kdtree for interpolation
        self.tree = cKDTree(xyTIN)
        tindx = xyTIN[1,0] - xyTIN[0,0]
        self.searchpts = max(int(sdx*sdx/(tindx*tindx)),4)

        if rstep > 0:
            distances, indices = self.tree.query(self.xyi, k=self.searchpts)
            weights = 1.0 / distances**2
            if len(cumdiff[indices].shape) == 3:
                cum_vals = cumdiff[indices][:,:,0]
            else:
                cum_vals = cumdiff[indices]
            fcum = numpy.average(cum_vals, weights=weights, axis=1)
            onIDs = numpy.where(distances[:,0] == 0)[0]
            if len(onIDs) > 0:
                fcum[onIDs] = cumdiff[indices[onIDs,0]]
            self.oldload = cumdiff

        return

    def update_TIN(self, xyTIN):
        """
        Update TIN mesh after 3D displacements.

        Args:
            xyTIN: numpy float-type array containing the coordinates for each nodes in the TIN (in m)
        """

        # Update TIN grid kdtree for interpolation
        self.tree = cKDTree(xyTIN)
        self.xyTIN = xyTIN

        return

    def move_mesh(self, dispX, dispY, cumdiff, verbose=False):
        """
        Update stratal mesh after 3D displacements.

        Args:
            dispX: numpy float-type array containing X-displacement for each nodes in the stratal mesh
            dispY: numpy float-type array containing Y-displacement for each nodes in the stratal mesh
            cumdiff: numpy float-type array containing the cumulative erosion/deposition of the nodes in the TIN
            verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).
        """

        # Move coordinates
        walltime = time.clock()
        st_time = walltime
        moveXY = numpy.zeros([self.xyi.shape[0],2])
        moveXY[:,0] = self.xyi[:,0] + dispX
        moveXY[:,1] = self.xyi[:,1] + dispY

        # Define point ids in local partition which needs to be send to neighbourhood
        l0 = self.nx
        l1 = l0 + self.nx
        u1 = self.ptsNb - self.nx
        u0 = u1 - self.nx
        shape = (l1-l0, self.step+1)

        if verbose:
            print(" - move stratal mesh ", time.clock() - walltime)

        walltime = time.clock()
        deformXY = moveXY
        deformThick = self.stratThick[:,:self.step+1]
        deformPoro = self.stratPoro[:,:self.step+1]
        deformElev = self.stratElev[:,:self.step+1]
        if verbose:
            print(" - create deformed stratal mesh arrays ", time.clock() - walltime)

        # Build the kd-tree
        walltime = time.clock()
        deformtree = cKDTree(deformXY)
        if verbose:
            print(" - create deformed stratal mesh kd-tree ", time.clock() - walltime)

        walltime = time.clock()
        distances, indices = deformtree.query(self.xyi, k=4)
        if verbose:
            print(" - query stratal mesh kd-tree ", time.clock() - walltime)

        # Compute inverse weighting distance
        walltime = time.clock()
        with numpy.errstate(divide='ignore'):
            w = 1.0 / distances**2
        w3D = w.reshape((len(self.xyi),4,1))
        weights = numpy.tile(w3D, (1,1,self.step+1))

        # Perform interpolation
        tmpIDs = numpy.where(distances[:,0] == 0)[0]
        if len(tmpIDs) > 0:
            self.stratThick[tmpIDs,:self.step+1] = deformThick[indices[tmpIDs,0],:self.step+1]
            self.stratPoro[tmpIDs,:self.step+1]  = deformPoro[indices[tmpIDs,0],:self.step+1]
            self.stratElev[tmpIDs,:self.step+1]  = deformElev[indices[tmpIDs,0],:self.step+1]
            tmpID = numpy.where(distances[:,0] > 0)[0]
            if len(tmpID) > 0:
                self.stratThick[tmpID,:self.step+1] = numpy.average(deformThick[indices[tmpID,:],:],
                                                                    weights=weights[tmpID,:], axis=1)
                self.stratPoro[tmpID,:self.step+1] = numpy.average(deformPoro[indices[tmpID,:],:],
                                                                    weights=weights[tmpID,:], axis=1)
                self.stratElev[tmpID,:self.step+1] = numpy.average(deformElev[indices[tmpID,:],:],
                                                                    weights=weights[tmpID,:], axis=1)

        else:
            self.stratThick[:,:self.step+1] = numpy.average(deformThick[indices,:],weights=weights, axis=1)
            self.stratPoro[:,:self.step+1] = numpy.average(deformPoro[indices,:],weights=weights, axis=1)
            self.stratElev[:,:self.step+1] = numpy.average(deformElev[indices,:],weights=weights, axis=1)

        # Reset depostion flag
        self.stratIn.fill(0)
        tmpID = numpy.where(numpy.amax(self.stratThick[:,:self.step+1], axis=1)>0)[0]
        self.stratIn[tmpID] = 1
        if verbose:
            print(" - perform stratal mesh interpolation ", time.clock() - walltime)

        if verbose:
            print(" - moving stratal mesh function ", time.clock() - st_time)

        self.oldload = numpy.copy(cumdiff)

        return

    def buildStrata(self, elev, cumdiff, sea, boundsPt, write=0, outstep=0):
        """
        Build the stratigraphic layer on the regular grid.

        Args:
            elev: numpy float-type array containing the elevation of the nodes in the TIN
            cumdiff: numpy float-type array containing the cumulative erosion/deposition of the nodes in the TIN
            sea: sea level elevation
            boundPts: number of nodes on the edges of the TIN surface.
            write: flag for output generation
            outstep: sStep for output generation

        Returns:
            - sub_poro - numpy array containing the subsidence induced by porosity change.
        """

        selev = numpy.zeros(len(self.xyi))
        distances, indices = self.tree.query(self.xyi, k=self.searchpts)
        with numpy.errstate(divide='ignore'):
            weights = 1.0 / distances**2

        if self.oldload is not None:
            load_diff = cumdiff - self.oldload
        else:
            load_diff = cumdiff
        if len(elev[indices].shape) == 3:
            elev_vals = elev[indices][:,:,0]
            cum_vals = load_diff[indices][:,:,0]
        else:
            elev_vals = elev[indices]
            cum_vals = load_diff[indices]

        felev = numpy.average(elev_vals,weights=weights, axis=1)
        fcum = numpy.average(cum_vals,weights=weights, axis=1)
        onIDs = numpy.where(distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            felev[onIDs] = elev[indices[onIDs,0]]
            fcum[onIDs] = load_diff[indices[onIDs,0]]
        self.oldload = numpy.copy(cumdiff)
        selev = felev

        # Update stratal elevation
        self.stratElev[self.ids,self.step] =  selev[self.ids]-sea

        localCum = fcum[self.ids]
        # Update stratal erosion
        eroIDs = numpy.where(localCum<0.)[0]
        self.eroLayer(self.ids[eroIDs], localCum)

        # Update stratal deposition
        depIDs = numpy.where(localCum>0.)[0]
        subs = self.depoLayer(self.ids[depIDs], localCum)
        subsi = numpy.reshape(subs,(len(self.ygrid),len(self.xgrid)))
        subs_values = RegularGridInterpolator((self.ygrid, self.xgrid), subsi)
        sub_poro = numpy.zeros(len(self.xyTIN[:,0]))
        sub_poro[boundsPt:] = subs_values((self.xyTIN[boundsPt:,1],self.xyTIN[boundsPt:,0]))
        sub_poro[sub_poro>0.] = 0.

        self.oldload += sub_poro
        if write>0:
            self.layerMesh(selev[self.ids]+subs[self.ids])
            self.write_hdf5_stratal(outstep-1)

        self.step += 1

        return sub_poro

    def buildPartition(self, bbX, bbY):
        """
        Define a partition for the stratal mesh.

        Args:
            bbX: box boundary for X-axis
            bbY: box boundary for Y-axis
        """

        size = 1
        # extent of X partition
        Yst = numpy.zeros( 1 )
        Yed = numpy.zeros( 1 )
        partYID = numpy.zeros( (1,2) )
        nbY = int((self.ny-1) / 1)
        p = 0
        Yst[p] = bbY[0]
        Yed[p] = Yst[p] + nbY*self.dx
        partYID[p,0] = 0
        partYID[p,1] = (nbY+1)*self.nx
        Yed[size-1] = bbY[1]
        partYID[size-1,1] = self.ny*self.nx

        # Get send/receive data ids for each processors
        self.upper = numpy.zeros( (1,2) )
        self.lower = numpy.zeros( (1,2) )
        self.upper[0,0] = partYID[0,1]
        self.upper[0,1] = partYID[0,1] + self.nx
        self.lower[0,0] = partYID[0,0]
        self.lower[0,1] = partYID[0,0] - self.nx

        # Define partitions ID globally
        Xst = numpy.zeros( 1 )
        Xed = numpy.zeros( 1 )
        Xst += bbX[0]
        Xed += bbX[1]

        # Loop over node coordinates and find if they belong to local partition
        # Note: used a Cython/Fython class to increase search loop performance... in utils
        partID = flowalgo.overlap(self.xyi[:,0],self.xyi[:,1],Xst[0],
                                        Yst[0],Xed[0],Yed[0])

        # Extract local domain nodes global ID
        self.ids = numpy.where(partID > -1)[0]

        return

    def depoLayer(self, ids, depo):
        """
        Add deposit to current stratigraphic layer.

        Args:
            ids: index of points subject to deposition
            depo: value of the deposition for the given point (m)

        Returns:
            - subs - numpy array containing the subsidence induced by porosity change.
        """

        # Initialise node deposition flag
        tmpIDs = numpy.where(self.stratIn[ids]==0)
        self.stratIn[ids[tmpIDs]] = 1

        # Add deposit to the considered layer time
        self.stratThick[ids,self.step] += depo[ids]

        # Define porosity values
        subs = numpy.zeros(self.ptsNb)
        if self.poro0 > 0:
            self.stratPoro[ids,self.step] = self.poro0
            cumThick = numpy.cumsum(self.stratThick[ids,self.step::-1],axis=1)[:,::-1]
            poro = self.poro0*numpy.exp(-self.poroC*cumThick/1000.)
            poro[poro<0.1] = 0.1
            #tmpid = numpy.where(self.stratPoro[ids,:self.step+1]<poro)[0]
            tmp1,tmp2 = numpy.where(numpy.logical_and(self.stratPoro[ids,:self.step+1]<poro,self.stratThick[ids,:self.step+1]>0.))
            if len(tmp1)>0:
                 poro[tmp1,tmp2] = self.stratPoro[ids[tmp1],tmp2]
            nh = self.stratThick[ids,:self.step+1]*(1.-self.stratPoro[ids,:self.step+1])/(1.-poro)
            nh = numpy.minimum(self.stratThick[ids,:self.step+1],nh)
            nh[nh<0.] = 0.
            # Update layer thickness
            self.stratThick[ids,:self.step+1] = nh
            # Update porosity
            self.stratPoro[ids,:self.step+1] = poro
            # Subsidence due to porosity change
            cumThick2 = numpy.cumsum(self.stratThick[ids,self.step::-1],axis=1)[:,::-1]
            subs[ids] = (cumThick2-cumThick)[:,0]
            # subs[ids] = numpy.sum(nh-self.stratThick[ids,:self.step+1],axis=1)
            subs[subs>0.] = 0.

        return subs

    def eroLayer(self, nids, erosion):
        """
        Erode top stratigraphic layers.

        Args:
            nids: index of points subject to erosion
            erosion: value of the erosion for the given points (m)
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

        Args:
            topsurf: elevation of the regular surface
        """

        # Clear points with no stratigraphic layer
        tmpIDs = numpy.where(self.stratIn == 0)[0]
        #surf = numpy.tile(topsurf[tmpIDs].transpose(), (1, self.step+1)).reshape(self.step+1,len(tmpIDs)).transpose()
        surf = numpy.array([topsurf[tmpIDs],]*int(self.step+2)).transpose()
        self.stratDepth[tmpIDs,:self.step+2] = surf
        self.stratThick[tmpIDs,:self.step+2] = 0.


        # Find points with stratigraphic layers
        tmpIDs = numpy.where(self.stratIn == 1)[0]
        if len(tmpIDs) == 0:
            return

        # Compute cumulative stratal thicknesses
        cumThick = numpy.cumsum(self.stratThick[tmpIDs,self.step+1::-1],axis=1)[:,::-1]

        # Updata stratal depth
        surf = numpy.array([topsurf[tmpIDs],]*int(self.step+2)).transpose()
        self.stratDepth[tmpIDs,:self.step+2] = surf - cumThick

        return

    def write_hdf5_stratal(self, outstep):
        """
        This function writes for each processor the HDF5 file containing sub-surface information.

        Args:
            outstep: output time step.
        """

        sh5file = self.folder+'/'+self.h5file+str(outstep)+'.hdf5'
        with h5py.File(sh5file, "w") as f:
            # Write node coordinates
            f.create_dataset('coords',shape=(self.ptsNb,2), dtype='float64', compression='gzip')
            f["coords"][:,:2] = self.xyi[self.ids]

            # Write stratal layers depth per cells
            f.create_dataset('layDepth',shape=(self.ptsNb,self.step+2), dtype='float64', compression='gzip')
            f["layDepth"][:,:self.step+2] = self.stratDepth[self.ids,:self.step+2]

            # Write stratal layers elevations per cells
            f.create_dataset('layElev',shape=(self.ptsNb,self.step+2), dtype='float64', compression='gzip')
            f["layElev"][:,:self.step+2] = self.stratElev[self.ids,:self.step+2]

            # Write stratal layers thicknesses per cells
            f.create_dataset('layThick',shape=(self.ptsNb,self.step+2), dtype='float64', compression='gzip')
            f["layThick"][:,:self.step+2] = self.stratThick[self.ids,:self.step+2]

            # Write stratal layers porosity per cells
            f.create_dataset('layPoro',shape=(self.ptsNb,self.step+2), dtype='float64', compression='gzip')
            f["layPoro"][:,:self.step+2] = self.stratPoro[self.ids,:self.step+2]

        return
