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
from scipy import interpolate

class forceSim: 
    """ 
    This class defines external forcing parameters.
    
    Parameters
    ----------
    string : seafile
        Path to the sea level fluctuation file (if any).
        
    float : sea0
        Relative sea level position in case no sea level curve is provided (default is 0.).
        
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
        
        event = numpy.where( (time - self.T_rain[:,0]) >= 0)[0] 
        
        tinRain = numpy.zeros(len(tXY[:,0]), dtype=float)
        
        if not (time >= self.T_rain[event,0]) and not (time < self.T_rain[event,1]):
            raise ValueError('Problem finding the rain map to load!')
            
        self.next_rain = self.T_rain[event,1]
        
        if self.Map_rain[event] == None:
            tinRain = self.rainVal[event]
        else:               
            rainMap = pandas.read_csv(str(self.Map_rain[event][0]), sep=' ', engine='c', header=None, na_filter=False, \
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
        
        event = numpy.where( (time - self.T_disp[:,0]) >= 0)[0] 
        
        tinDisp = numpy.zeros(len(tXY[:,0]), dtype=float)
        
        if not (time >= self.T_disp[event,0]) and not (time < self.T_disp[event,1]):
            raise ValueError('Problem finding the displacement map to load!')
        
        self.next_disp = self.T_disp[event,1]
        
        if self.Map_disp[event] != None:
            dispMap = pandas.read_csv(str(self.Map_disp[event][0]), sep=' ', engine='c', header=None, na_filter=False, \
                               dtype=numpy.float, low_memory=False)  
            
            rectDisp = numpy.reshape(dispMap.values,(len(self.regX), len(self.regY)),order='F')
            tinDisp = interpolate.interpn( (self.regX, self.regY), rectDisp, tXY, method='linear')
            dt = (self.T_disp[event,1] - self.T_disp[event,0])
            if dt <= 0:
                raise ValueError('Problem computing the displacement rate for event %d.'%event)
            tinDisp = tinDisp / dt
                                     
        return tinDisp
    
    def load_Disp_map(self, time, tXY):
        """
        Load 3D displacements map for a given period and perform interpolation from regular grid to unstructured TIN one.
        
        Parameters
        ----------
        float : time
            Requested time interval rain map to load.

        float : tXY
            Unstructured grid (TIN) XY coordinates.
            
        Return
        ----------
        variable: tinDisp
            Numpy array containing the updated 3D displacement rate for the local domain.

        """
        event = numpy.where( (time - self.T_disp[:,0]) >= 0)[0] 
        
        tinDisp = numpy.zeros(len(tXY[:,0]), dtype=float)
        
        return tinDisp
                 
            
