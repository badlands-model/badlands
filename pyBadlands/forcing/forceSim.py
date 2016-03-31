##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module defines several functions used to force Badlands simulation with external
processes related to climate, tectonic and sea level.
"""
import os
import numpy 
import pandas
import triangle
import mpi4py.MPI as mpi
from scipy import interpolate
from scipy.spatial import cKDTree
from scipy.interpolate import griddata

class forceSim: 
    """ 
    This class defines external forcing parameters.
    
    Parameters
    ----------
    string : seafile
        Path to the sea level fluctuation file (if any).
        
    float : sea0
        Relative sea level position in case no sea level curve is provided (default is 0.).
        
    string : MapRain
        Numpy array containing the rain map file names.
        
    float : TimeRain
        Numpy array containing the start and end times for each rain event in years.
        
    float : ValRain
        Value of precipitation rate for each rain event in m/a.
    
    string : MapDisp
        Numpy array containing the cumulative displacement map file names.
        
    float : TimeDisp
        Numpy array containing the start and end times for each displacement period in years.
        
    float : regX
        Numpy array containing the X-coordinates of the regular input grid.
        
    float : regY
        Numpy array containing the Y-coordinates of the regular input grid.
        
    float : Tdisplay
        Display interval (in years).
        
    """
    
    def __init__(self, seafile = None, sea0 = 0., MapRain = None, TimeRain = None, ValRain = None,
                 MapDisp = None, TimeDisp = None, regX = None, regY = None, Tdisplay = 0.):
        
        self.regX = regX
        self.regY = regY
        self.Map_rain = MapRain
        self.rainVal = ValRain
        self.T_rain = TimeRain
        self.next_rain = None
        
        self.Map_disp = MapDisp
        self.T_disp = TimeDisp
        self.next_disp = None
        
        self.sea0 = sea0
        self.seafile = seafile
        self.sealevel = None
        self.seatime = None
        self.seaval = None
        self.seaFunc = None
        
        self.tXY = None
        self.dispX = None
        self.dispY = None
        self.dispZ = None
        self.merge3d = None
        self.time3d = None
        
        self.next_display = None
        self.time_display = Tdisplay
        
        if self.seafile != None:
            self._build_Sea_function()
        
        return
    
    def _build_Sea_function(self):
        """
        Using Pandas library to read the sea level file and define sea level interpolation
        function based on Scipy 1D linear function.
        """
        
        # Read sea level file
        seadata = pandas.read_csv(self.seafile, sep=' ', engine='c', 
                               header=None, na_filter=False,
                               dtype=numpy.float, low_memory=False)
        
        self.seatime = seadata.values[:,0]
        self.seaval = seadata.values[:,1]
        
        self.seaFunc = interpolate.interp1d(self.seatime, self.seaval, kind='linear')
        
        return
    
    def getSea(self, time):
        """
        Computes for a given time the sea level according to input file parameters.
        
        Parameters
        ----------
        float : time
            Requested time for which to compute sea level elevation.
        """
        
        if self.seafile == None:
            self.sealevel = self.sea0        
        else:
            if time < self.seatime.min():
                time = self.seatime.min()
            if time > self.seatime.max():
                time = self.seatime.max()

            self.sealevel = self.seaFunc(time)
        
        return
    
    def load_Rain_map(self, time, tXY):
        """
        Load rain map for a given period and perform interpolation from regular grid to unstructured TIN one.
        
        Parameters
        ----------
        float : time
            Requested time interval rain map to load.

        float : tXY
            Unstructured grid (TIN) XY coordinates.
            
        Return
        ----------
        variable: tinRain
            Numpy array containing the updated rainfall for the local domain.

        """
        
        event = numpy.where( (time - self.T_rain[:,0]) >= 0)[0][0]
        
        if not (time >= self.T_rain[event,0]) and not (time < self.T_rain[event,1]):
            raise ValueError('Problem finding the rain map to load!')
            
        self.next_rain = self.T_rain[event,1]
        
        if self.Map_rain[event] == None:
            tinRain = numpy.zeros(len(tXY[:,0]), dtype=float)
            tinRain = self.rainVal[event]
        else:
            rainMap = pandas.read_csv(str(self.Map_rain[event]), sep=' ', engine='c', header=None, na_filter=False, \
                               dtype=numpy.float, low_memory=False)  
            
            rectRain = numpy.reshape(rainMap.values,(len(self.regX), len(self.regY)),order='F')
            tinRain = interpolate.interpn( (self.regX, self.regY), rectRain, tXY, method='linear')
                  
        return tinRain
    
    def disp_border(self, disp, neighbours, edge_lenght, boundPts):
        """ 
        This function defines the displacement of the TIN edges.

        Parameters
        ----------
        variable : disp
            Numpy arrays containing the internal nodes displacement value.

        variable : neighbours
            Numpy integer-type array containing for each nodes its neigbhours IDs.

        variable : edge_lenght
            Numpy float-type array containing the lenghts to each neighbour.

        variable : boundPts
            Number of nodes on the edges of the TIN surface.

        Return
        ----------
        variable: disps
            Numpy array containing the updated displacements on the edges.
        """
        disp[:boundPts] = 1.e7
        missedPts = []
        for id in range(boundPts):
            ngbhs = neighbours[id,:]
            ids = numpy.where(ngbhs >= boundPts)[0]
            if len(ids) == 1:
                disp[id] = disp[ngbhs[ids]]
            elif len(ids) > 1:
                lselect = edge_lenght[id,ids]
                picked = numpy.argmin(lselect)
                disp[id] = disp[ngbhs[ids[picked]]]
            else:
                missedPts = numpy.append(missedPts,id)

        if len(missedPts) > 0 :
            for p in range(len(missedPts)):
                id = missedPts[p]
                ngbhs = neighbours[id,:]
                ids = numpy.where((disp[ngbhs] < 9.e6) & (ngbhs >= 0))[0]
                if len(ids) == 0:
                    raise ValueError('Error while getting boundary displacement for point ''%d''.' % id)
                lselect = edge_lenght[id,ids]
                picked = numpy.argmin(lselect)
                disp[id] = disp[ngbhs[ids[picked]]]

        return disp
    
    def load_Tecto_map(self, time, tXY):
        """
        Load vertical displacement map for a given period and perform interpolation from regular grid to unstructured TIN one.
        
        Parameters
        ----------
        float : time
            Requested time interval rain map to load.

        float : tXY
            Unstructured grid (TIN) XY coordinates.
            
        Return
        ----------
        variable: tinDisp
            Numpy array containing the updated displacement rate for the local domain.

        """
        
        event = numpy.where( (time - self.T_disp[:,0]) >= 0)[0][0]
        
        if not (time >= self.T_disp[event,0]) and not (time < self.T_disp[event,1]):
            raise ValueError('Problem finding the displacements map to load!')
        
        self.next_disp = self.T_disp[event,1]
        
        if self.Map_disp[event] != None:
            dispMap = pandas.read_csv(str(self.Map_disp[event]), sep=' ', engine='c', header=None, na_filter=False, \
                               dtype=numpy.float, low_memory=False)  
            
            rectDisp = numpy.reshape(dispMap.values,(len(self.regX), len(self.regY)),order='F')
            tinDisp = interpolate.interpn( (self.regX, self.regY), rectDisp, tXY, method='linear')
            dt = (self.T_disp[event,1] - self.T_disp[event,0])
            if dt <= 0:
                raise ValueError('Problem computing the displacements rate for event %d.'%event)
            tinDisp = tinDisp / dt
        else:
            tinRain = numpy.zeros(len(tXY[:,0]), dtype=float)
            
        return tinDisp
    
    def load_Disp_map(self, time, tXY, inIDs):
        """
        Load 3D displacements map for a given period and perform interpolation from regular grid to unstructured TIN one.
        
        Parameters
        ----------
        float : time
            Requested time interval rain map to load.

        float : tXY
            Unstructured grid (TIN) XY coordinates.
            
        integer : inDs
            List of unstructured vertices contained in each partition.
            
        """
        # Initialise MPI communications
        comm = mpi.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
    
        self.tXY = tXY
    
        totPts = len(tXY[:,0])
        dispX = numpy.zeros(totPts, dtype=float)
        dispY = numpy.zeros(totPts, dtype=float)
        dispZ = numpy.zeros(totPts, dtype=float)
        dispX.fill(-1.e6)
        dispY.fill(-1.e6)
        dispZ.fill(-1.e6)
        
        dpXY = tXY[inIDs]
        locPts = len(inIDs)
        event = numpy.where( (time - self.T_disp[:,0]) >= 0)[0][0]
                
        if not (time >= self.T_disp[event,0]) and not (time < self.T_disp[event,1]):
            raise ValueError('Problem finding the 3D displacements map to load!')
            
        if self.time3d > 0.:
            if time + self.time3d > self.T_disp[event,1]:
                self.next_disp = self.T_disp[event,1]
            else:
                self.next_disp = self.time3d + time
        else:
            self.next_disp = self.T_disp[event,1]
            
        if self.Map_disp[event] != None:
            disps = pandas.read_csv(str(self.Map_disp[event]), sep=' ', engine='c', header=None, na_filter=False, \
                               dtype=numpy.float, low_memory=False)  
            
            disprX = numpy.reshape(disps.values[:,0],(len(self.regX), len(self.regY)),order='F')
            disprY = numpy.reshape(disps.values[:,1],(len(self.regX), len(self.regY)),order='F')
            disprZ = numpy.reshape(disps.values[:,2],(len(self.regX), len(self.regY)),order='F')
            
            dispX[inIDs] = interpolate.interpn( (self.regX, self.regY), disprX, dpXY, method='linear')
            dispY[inIDs] = interpolate.interpn( (self.regX, self.regY), disprY, dpXY, method='linear')
            dispZ[inIDs] = interpolate.interpn( (self.regX, self.regY), disprZ, dpXY, method='linear')
        else:
            raise ValueError('Problem finding the 3D displacements map to load!')
        
        comm.Allreduce(MPI.IN_PLACE, dispX, op=MPI.MAX)
        comm.Allreduce(MPI.IN_PLACE, dispY, op=MPI.MAX)
        comm.Allreduce(MPI.IN_PLACE, dispZ, op=MPI.MAX)
        
        if self.time3d > 0.:
            rate = (self.next_disp - time) / (self.T_disp[event,1] - self.T_disp[event,0])
            dispX = dispX * rate
            dispY = dispY * rate
            dispZ = dispZ * rate
            
        self.dispX = dispX
        self.dispY = dispY
        self.dispZ = dispZ
        
        return
                 
    def apply_XY_dispacements(self, area, fixIDs, elev, cum):
        """
        Apply horizontal displacements and check if any points need to be merged.
        
        Parameters
        ----------
        float : area
            Averaged area of the irregular grid delaunay cells.

        integer : fixIDs
            Number of unstructured vertices which needs to stay fix (edges and borders nodes).
            
        float : elev
            Numpy array with elevation of previous TIN nodes.
        
        float : cum
            Numpy array with erosion/deposition values from previous TIN nodes.
            
        Return
        ----------
        variable: tinMesh
            Delaunay mesh generated after displacements.
            
        variable: newelev
            Numpy array containing the updated elevation for the new TIN.
            
        variable: newcum
            Numpy array containing the updated erosion/deposition values for the new TIN.

        """ 
        elev += self.dispZ 
        self.tXY[fixIDs:,0] += self.dispX[fixIDs:]
        self.tXY[fixIDs:,1] += self.dispY[fixIDs:]
        
        tree = cKDTree(self.tXY)
        pairs = tree.query_pairs(self.merge3d)
        
        if len(pairs) > 0:
            pairIDs = numpy.array(list(pairs))
            nonfixIDs = numpy.where(numpy.logical_and( pairIDs[:,0] >= fixIDs, pairIDs[:,1] >= fixIDs))[0]
            
            mXY = numpy.empty(shape=[len(nonfixIDs),2], dtype=float)
            mXY[:,0] = 0.5*(self.tXY[pairIDs[nonfixIDs,0],0] + self.tXY[pairIDs[nonfixIDs,1],0])
            mXY[:,1] = 0.5*(self.tXY[pairIDs[nonfixIDs,0],1] + self.tXY[pairIDs[nonfixIDs,1],1])
            
            mergedIDs = numpy.unique(pairIDs[nonfixIDs,:].flatten())
            
            distances, indices = tree.query(mXY, k=3)
            z_vals = elev[indices]
            cum_vals = cum[indices]
            z_avg = numpy.average(z_vals, weights=(1./distances**2),axis=1)
            cum_avg = numpy.average(cum_vals, weights=(1./distances**2),axis=1)
            
            newXY = numpy.delete(self.tXY, mergedIDs, 0)
            newelev = numpy.delete(elev, mergedIDs, 0)
            newcum = numpy.delete(cum, mergedIDs, 0)
            
            newXY = numpy.concatenate((newXY, mXY), axis=0)
            newelev = numpy.concatenate((newelev, z_avg), axis=0)
            newcum = numpy.concatenate((newcum, cum_avg), axis=0)
        else:
            newXY = self.tXY
            newelev = elev
            newcum = cum
        
        newTIN = triangle.triangulate( dict(vertices=newXY),'Da'+str(area))
        
        if len(newTIN['vertices'][:,0]) > len(newXY[:,0]):
            addPts = newTIN['vertices'][len(newXY[:,0]):,:2]
            dist, ids = tree.query(addPts, k=3)
            zvals = elev[ids]
            cumvals = cum[ids]
            zavg = numpy.average(zvals, weights=(1./dist**2),axis=1)
            cumavg = numpy.average(cumvals, weights=(1./dist**2),axis=1)
            newelev = numpy.concatenate((newelev, zavg), axis=0)
            newcum = numpy.concatenate((newcum, cumavg), axis=0) 
        
        elif len(newTIN['vertices'][:,0]) < len(newXY):
            raise ValueError('Problem building the TIN after 3D displacements.')
            
        return newTIN, newelev, newcum
             
    def update_TIN_data(self, elevation, cumdiff):
        """
        Apply horizontal displacements and check if any points need to be merged.
        
        Parameters
        ----------
        float : elevation
            Elevation of previous TIN nodes.
        
        float : cumdiff
            Erosion/deposition values from previous TIN nodes.
            
        Return
        ----------
        variable: newelev
            Numpy array containing the updated elevation for the new TIN.
            
        variable: newdiff
            Numpy array containing the updated erosion/deposition values for the new TIN.

        """ 
        
        newelev = griddata(self.tXY,elevation,(self.tinMesh['vertices'][:,0], self.tinMesh['vertices'][:,1]), 
                           method='nearest')
        
        newdiff = griddata(self.tXY,cumdiff,(self.tinMesh['vertices'][:,0], self.tinMesh['vertices'][:,1]), 
                           method='nearest')
        
        return newelev, newdiff