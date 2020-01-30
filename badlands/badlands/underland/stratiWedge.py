##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module defines the **stratigraphic layers** based on the **irregular TIN grid**.
"""
import os
import glob
import time
import h5py
import numpy
import pandas
from scipy import interpolate
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator

if 'READTHEDOCS' not in os.environ:
    from badlands import pdalgo

class stratiWedge():
    """
    This class builds stratigraphic layers over time based on **erosion/deposition** values.

    Args:
        layNb: total number of stratigraphic layers
        elay: regular grid layer thicknesses
        xyTIN: numpy float-type array containing the coordinates for each nodes in the TIN (in m)
        bPts: boundary points for the TIN.
        ePts: boundary points for the regular grid.
        thickMap: filename containing initial layer parameters
        activeh: active layer thickness.
        folder: name of the output folder.
        h5file: first part of the hdf5 file name.
        rockNb: number of type of sedimentary rocks represented in the current simulation.
        regX: numpy array containing the X-coordinates of the regular input grid.
        regY: numpy array containing the Y-coordinates of the regular input grid.
        elev: numpy arrays containing the elevation of the TIN nodes.
        cumdiff: numpy array containing cumulative erosion/deposition from previous simulation.
        rockCk: numpy arrays containing the different rock erodibility values.
        rfolder: restart folder.
        rstep: restart step.
    """

    def __init__(self, layNb, elay, xyTIN, bPts, ePts, thickMap, activeh, folder, h5file, rockNb,
                 regX, regY, elev, rockCk, cumdiff=0, rfolder=None, rstep=0):

        # Number of points on the TIN
        self.ptsNb = len(xyTIN)
        self.oldload = None
        self.folder = folder
        self.h5file = h5file
        self.rockNb = rockNb
        self.initlay = elay
        self.alayR = None
        self.activeh = activeh
        self.rockCk = rockCk

        # In case we restart a simulation
        if rstep > 0:
            if os.path.exists(rfolder):
                folder = rfolder+'/h5/'
                fileCPU = 'stratal.time%s.hdf5'%rstep
                restartncpus = len(glob.glob1(folder,fileCPU))
                if restartncpus == 0:
                    raise ValueError('The requested time step for the restart simulation cannot be found in the restart folder.')
            else:
                raise ValueError('The restart folder is missing or the given path is incorrect.')

            if restartncpus != 1:
                raise ValueError('When using the stratal model you need to run the restart simulation with the same number of processors as the previous one.')

            df = h5py.File('%s/h5/stratal.time%s.hdf5'%(rfolder, rstep), 'r')

            paleoDepth = numpy.array((df['/paleoDepth']))
            eroLay =  paleoDepth.shape[1]
            self.step = paleoDepth.shape[1]

            # Elevation at time of deposition (paleo-depth)
            self.paleoDepth = numpy.zeros((self.ptsNb,layNb+eroLay),order='F')
            self.layerThick = numpy.zeros((self.ptsNb,layNb+eroLay),order='F')
            self.paleoDepth[:,:eroLay] = paleoDepth
            # Deposition thickness for each type of sediment
            self.depoThick = numpy.zeros((self.ptsNb,layNb+eroLay,rockNb),order='F')
            for r in range(rockNb):
                self.depoThick[:,:eroLay,r] = numpy.array((df['/depoThickRock'+str(r)]))
            self.layerThick[:,:eroLay] = numpy.sum(self.depoThick[:,:eroLay,:],axis=-1)
            # Define cumulative deposition/erosion from previous simulation
            self.oldload = cumdiff
        else:
            eroLay = elay+1
            self.step = eroLay

            tmpTH = numpy.zeros(self.ptsNb)
            # Elevation at time of deposition (paleo-depth)
            self.paleoDepth = numpy.zeros((self.ptsNb,layNb+eroLay),order='F')
            # Deposition thickness for each type of sediment
            self.depoThick = numpy.zeros((self.ptsNb,layNb+eroLay,rockNb),order='F')
            self.layerThick = numpy.zeros((self.ptsNb,layNb+eroLay),order='F')
            # Rock type array
            rockType = -numpy.ones(self.ptsNb,dtype=int)

            # If predefined layers exists
            if elay > 0:
                # Build the underlying erodibility mesh and associated thicknesses

                # Define inside area kdtree
                inTree = cKDTree(xyTIN[bPts:ePts+bPts,:])
                dist, inID = inTree.query(xyTIN[:bPts,:],k=1)
                inID += bPts

                # Data is stored from top predefined layer to bottom.
                self.paleoDepth[:,eroLay] = elev
                for l in range(1,eroLay):
                    thMap = pandas.read_csv(str(thickMap[l-1]), sep=r'\s+', engine='c',
                                       header=None, na_filter=False, dtype=numpy.float, low_memory=False)
                    # Extract thickness values
                    tmpH = thMap.values[:,0]
                    tH = numpy.reshape(tmpH,(len(regX), len(regY)), order='F')
                    # Nearest neighbours interpolation to extract rock type values
                    tmpS = thMap.values[:,1].astype(int)
                    tS = numpy.reshape(tmpS,(len(regX), len(regY)), order='F')
                    rockType[bPts:] = interpolate.interpn( (regX, regY), tS, xyTIN[bPts:,:], method='nearest')
                    # Linear interpolation to define underlying layers on the TIN
                    tmpTH.fill(0.)
                    tmpTH[bPts:] = interpolate.interpn( (regX, regY), tH, xyTIN[bPts:,:], method='linear')
                    for r in range(rockNb):
                        ids = numpy.where(numpy.logical_and(rockType==r,tmpTH>0.))
                        self.depoThick[ids,eroLay-l,r] = tmpTH[ids]
                        self.paleoDepth[ids,eroLay-l] = self.paleoDepth[ids,eroLay-l+1]-tmpTH[ids]
                        if eroLay-l==1:
                            self.depoThick[ids,0,r] = 1.e6
                            self.paleoDepth[ids,0] = self.paleoDepth[ids,1]-1.e6
                        # Add an infinite rock layer with the same characteristics as the deepest one
                        self.depoThick[:bPts,eroLay-l,r] = self.depoThick[inID,eroLay-l,r]
                        if r == 0:
                            self.paleoDepth[:bPts,eroLay-l] = self.paleoDepth[:bPts,eroLay-l+1]-self.depoThick[:bPts,eroLay-l,r]
                        else:
                            self.paleoDepth[:bPts,eroLay-l] -= self.depoThick[:bPts,eroLay-l,r]
                        if eroLay-l==1:
                            ids = numpy.where(self.depoThick[:bPts,eroLay-l,r]>0)[0]
                            self.depoThick[ids,0,r] = 1.e6
                            self.paleoDepth[ids,0] = self.paleoDepth[ids,1]-1.e6
                    self.layerThick[:,eroLay-l] = numpy.sum(self.depoThick[:,eroLay-l,:],axis=-1)
                    self.layerThick[:,0] = 1.e6
            else:
                # Add an infinite rock layer with the same characteristics as the deepest one
                self.depoThick[:,0,0] = 1.e6
                self.layerThick[:,0] = self.depoThick[:,0,0]
                self.paleoDepth[:,0] = elev
                self.step = 1

    def get_active_layer(self, actlay, verbose=False):
        """
        This function extracts the active layer based on the underlying stratigraphic architecture.

        Args:
            actlay: active layer elevation based on nodes elevation [m]
            verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).
        """

        time0 = time.clock()
        self.alayR = pdalgo.getactlay(actlay, self.layerThick[:,:self.step+1],
                                    self.depoThick[:,:self.step+1,:])
        if verbose:
            print("   - Get active layer composition ", time.clock() - time0)
            time0 = time.clock()

        return

    def update_layers(self, erosion, deposition, elev, verbose=False):
        """
        This function updates the stratigraphic layers based active layer composition.

        Args:
            erosion: value of the erosion for the given points (m)
            deposition: value of the deposition for the given point (m)
            elev: numpy float-type array containing Z coordinates of the local TIN nodes.
            verbose : (bool) when :code:`True`, output additional debug information (default: :code:`False`).
        """

        time0 = time.clock()
        newH, newS = pdalgo.updatestrati(self.depoThick[:,:self.step+1,:], self.layerThick[:,:self.step+1],
                                    erosion, deposition)

        self.depoThick[:,:self.step+1,:] = newS
        self.layerThick[:,:self.step+1] = newH
        self.paleoDepth[:,self.step] = elev

        if verbose:
            print("   - Update erosion/deposition ", time.clock() - time0)

        return

    def write_hdf5_stratigraphy(self, lGIDs, outstep):
        """
        This function writes for each processor the HDF5 file containing sub-surface information.

        Args:
            lGIDs: global node IDs for considered partition.
            outstep: output time step.
        """

        sh5file = self.folder+'/'+self.h5file+str(outstep)+'.hdf5'
        with h5py.File(sh5file, "w") as f:

            # Write stratal layers paeleoelevations per cells
            f.create_dataset('paleoDepth',shape=(len(lGIDs),self.step+1), dtype='float64', compression='gzip')
            f["paleoDepth"][lGIDs,:self.step+1] = self.paleoDepth[lGIDs,:self.step+1]

            # Write stratal layers thicknesses per cells
            for r in range(self.rockNb):
                f.create_dataset('depoThickRock'+str(r),shape=(len(lGIDs),self.step+1), dtype='float64', compression='gzip')
                f['depoThickRock'+str(r)][lGIDs,:self.step+1] = self.depoThick[lGIDs,:self.step+1,r]

        return
