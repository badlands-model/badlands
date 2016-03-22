##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module encapsulates functions related to Badlands stream network computation.
"""

import math
import numpy
import warnings
import mpi4py.MPI as MPI

from pyBadlands.libUtils import SFDalgo as SFD
from pyBadlands.libUtils import FLOWalgo
from pyBadlands.libUtils import FLWnetwork

class flowNetwork: 
    """
    Class for handling flow network computation based on Braun & Willett's 
    algorithm.
    """
    
    def __init__(self):

        '''Initialization.
        '''
        self.base = None
        self.localbase = None
        self.receivers = None
        self.arrdonor = None
        self.delta = None
        self.donors = None
        self.localstack = None
        self.stack = None
        self.partFlow = None
        self.maxdonors = 0
        self.CFL = None
        self.erodibility = None
        self.m = None
        self.n = None
        self.mindt = None
        
        self.sendprocID = None
        self.rcvprocID = None
        self.rcvprocNb = None
        self.rcvsendID = None
        self.localGIDs = None
        self.localNodes = 0
        
        self.discharge = None
        self.localsedflux = None
        self.maxh = None
        self.maxdep = None
        self.diff_flux = None
        self.chi = None
        self.basinID = None
        #self.drainage_area = None
        
        return
    
    def SFD_receivers(self, fillH, elev, neighbours, edges, distances, globalIDs, sea):
        """
        Single Flow Direction function computes downslope flow directions by inspecting the neighborhood 
        elevations around each node. The SFD method assigns a unique flow direction towards the steepest 
        downslope neighbor.
        
        Parameters
        ----------
        variable : fillH
            Numpy array containing the filled elevations from Planchon & Darboux depression-less algorithm.
            
        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.
            
        variable : neighbours
            Numpy integer-type array with the neighbourhood IDs.
            
        variable : edges
            Numpy real-type array with the voronoi edges length for each neighbours of local IDs.
            
        variable : distances
            Numpy real-type array with the distances between each connection in the local TIN.
            
        variable: globalIDs
            Numpy integer-type array containing for local nodes their global IDs.
                        
        variable : sea
            Current elevation of sea level.
        """
       
        # Initialise MPI communications
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        # Call the SFD function from libUtils
        base, receivers, maxh, maxdep, diff_flux = SFD.sfdcompute.directions( fillH, elev, \
            neighbours, edges, distances, globalIDs, sea)
        
        # Send local base level globally
        comm.Allreduce(MPI.IN_PLACE,base,op=MPI.MAX)
        self.base = numpy.where(base >= 0)[0]

        # Send local receivers globally
        comm.Allreduce(MPI.IN_PLACE,receivers,op=MPI.MAX)
        self.receivers = receivers
        
        # Send local maximum deposition globally
        comm.Allreduce(MPI.IN_PLACE,maxh,op=MPI.MAX)
        self.maxh = maxh
        
        # Send local maximum deposition globally
        comm.Allreduce(MPI.IN_PLACE,maxdep,op=MPI.MAX)
        self.maxdep = maxdep
        
        # Send local diffusion flux globally
        comm.Allreduce(MPI.IN_PLACE,diff_flux,op=MPI.MAX)
        self.diff_flux = diff_flux
        
        return 
    
    def _donors_number_array(self):
        """
        Creates an array containing the number of donors for each node.
        """
    
        self.arrdonor = None
        numPts = len(self.receivers)
        self.arrdonor = numpy.zeros(numPts, dtype=int)
        maxID = numpy.max(self.receivers)
        self.arrdonor[:(maxID+1)] = numpy.bincount(self.receivers)
        self.maxdonors = self.arrdonor.max()
        
        return 

    def _delta_array(self):
        """
        Creates the "delta" array, which is a list containing, for each
        node, the array index where that node's donor list begins.
        """

        self.delta = None
        nbdonors = len(self.arrdonor)
        self.delta = numpy.zeros( nbdonors+1, dtype=int)
        self.delta.fill(nbdonors)
        self.delta[-2::-1] -= numpy.cumsum(self.arrdonor[::-1])

        return 
    
    def ordered_node_array(self):
        """
        Creates an array of node IDs that is arranged in order from downstream 
        to upstream.
        """
        
        # Build donors array for each node
        self._donors_number_array()
        
        # Create the delta array
        self._delta_array()
        
        # Using libUtils stack create the ordered node array
        self.donors,lstcks = FLWnetwork.fstack.build(self.localbase,self.receivers,self.delta)

        # Create local stack
        stids = numpy.where(lstcks > -1 )[0]        
        self.localstack = lstcks[stids]
        
        return  
    
    def compute_flow(self, Acell, rain, parallel=False):
        """
        Calculates the drainage area and water discharge at each node.
        
        Parameters
        ----------
        variable : Acell
            Numpy float-type array containing the voronoi area for each nodes (in m2)
            
        variable : rain
            Numpy float-type array containing the precipitation rate for each nodes (in m/a).
            
        variable : parallel
            Boolean to inform if the discharge algorithm need to be ran in parallel (True) or serial (False).
        """
        
        # Initialise MPI communications
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
            
        numPts = len(Acell)
        
        self.discharge = numpy.zeros(numPts, dtype=float) 
        self.discharge[self.stack] = Acell[self.stack]*rain[self.stack]

        # Compute discharge using libUtils
        if(parallel):
            discharge = FLOWalgo.flowcompute.discharge(self.localstack,self.receivers,self.discharge)
            comm.Allreduce(MPI.IN_PLACE,discharge,op=MPI.MAX)
        else:
            discharge = FLOWalgo.flowcompute.discharge(self.stack,self.receivers,self.discharge)
        
        self.discharge = discharge
        
        return     
    
    def compute_parameters(self, xycoords):
        """
        Calculates the catchment IDs and the Chi parameter (Willett 2014).
        
        Parameters
        ----------
        variable : xycoords
            Numpy float-type array containing X, Y coordinates of the TIN nodes.
        """
        
        # Initialise MPI communications
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        # Compute discharge using libUtils
        splexp = self.m / self.n
        chi, basinID = FLOWalgo.flowcompute.parameters(self.localstack,self.receivers,self.discharge,xycoords,splexp)
        comm.Allreduce(MPI.IN_PLACE,chi,op=MPI.MAX)
        comm.Allreduce(MPI.IN_PLACE,basinID,op=MPI.MAX)
        
        self.chi = chi
        self.basinID = basinID
        
        return   
    
    def compute_sedflux(self, Acell, elev, fillH, xycoords, xymin, xymax, diff_flux, dt, sealevel, parallel=False):
        """
        Calculates the sediment flux at each node.
        
        Parameters
        ----------
        variable : Acell
            Numpy float-type array containing the voronoi area for each nodes (in m2)
            
        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.
            
        variable : fillH
            Numpy array containing the filled elevations from Planchon & Darboux depression-less algorithm.
            
        variable : xycoords
            Numpy float-type array containing X, Y coordinates of the TIN nodes.
            
        variable : xymin
            Numpy array containing the minimal XY coordinates of TIN grid (excuding border cells).
            
        variable : xymax
            Numpy array containing the maximal XY coordinates of TIN grid (excuding border cells).
            
        variable : diff_flux
            Numpy arrays representing the fluxes due to linear diffusion equation. 
            
        variable : dt
            Real value corresponding to the maximal stability time step. 
            
        variable : sealevel
            Real value giving the sea-level height at considered time step.
            
        variable : parallel
            Boolean to inform if the discharge algorithm need to be ran in parallel (True) or serial (False).
        """
        
        # Initialise MPI communications
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        # Compute sediment flux using libUtils
        if(parallel):
            sedflux, newdt = FLOWalgo.flowcompute.sedflux(self.localstack,self.receivers,xycoords,\
                     Acell,xymin,xymax,self.maxh,self.maxdep,self.discharge,fillH,elev,diff_flux, \
                     self.erodibility,self.m,self.n,sealevel,dt)
            timestep = numpy.zeros(1)
            timestep[0] = newdt
            comm.Allreduce(MPI.IN_PLACE,timestep,op=MPI.MIN)
            newdt = timestep[0]
            
            comm.Allreduce(MPI.IN_PLACE,sedflux,op=MPI.MAX)
            tempIDs = numpy.where(sedflux < -9.5e5)
            sedflux[tempIDs] = 0.
            newdt = min(self.mindt,newdt)
            sedrate = sedflux 
        else:
            sedflux, newdt = FLOWalgo.flowcompute.sedflux(self.stack,self.receivers,xycoords,\
                     Acell,xymin,xymax,self.maxh,self.maxdep,self.discharge,fillH,elev,diff_flux, \
                     self.erodibility,self.m,self.n,sealevel,dt)
            tempIDs = numpy.where(sedflux < -9.5e5)
            sedflux[tempIDs] = 0.
            newdt = min(self.mindt,newdt)
            sedrate = sedflux 
            
        return newdt,sedrate    
    
    def dt_pstability(self, xy, elev, locIDs):
        """ 
        This pure python function computes the maximal timestep to ensure computation stability 
        of the flow processes. This CFL-like condition is computed using erodibility 
        coefficients, discharges plus elevations and distances between TIN nodes and 
        their respective reveivers. 
        
        Parameters
        ----------
        variable : xy
            Numpy arrays containing the XY coordinates of the TIN nodes.
            
        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.
            
        variable: locIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """
        
        # Initialise MPI communications
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        # Compute elevation difference and distance between receiver and donor
        dz0 = elev[locIDs] - elev[self.receivers[locIDs]]
        dx = xy[locIDs,0] - xy[self.receivers[locIDs],0]
        dy = xy[locIDs,1] - xy[self.receivers[locIDs],1]
        dist0 = numpy.sqrt(dx**2 + dy**2)
        dshg0 = self.discharge[locIDs]
        
        # Limit to positive values of dz and discharge
        tmpID = numpy.where(dz0 > 0.)
        tmpdz = dz0[tmpID]
        tmpdist = dist0[tmpID]
        tmpdisc = dshg0[tmpID]
        tmpID2 = numpy.where( tmpdisc > 0.)
        discharge = tmpdisc[tmpID2]
        dist = tmpdist[tmpID2]
        dz = tmpdz[tmpID2]
        
        # Compute the local value for time stability
        CFL = numpy.zeros(1)
        CFL[0] = numpy.amin(dist / (self.erodibility * (discharge)**(self.m) * \
                ( dz / dist )**(self.n-1.) ))
        
        # Global mimimum value for diffusion stability
        comm.Allreduce(MPI.IN_PLACE,CFL,op=MPI.MIN)
        self.CFL = CFL[0]
        
    def dt_fstability(self, xy, elev, locIDs):
        """ 
        This pyfortran function computes the maximal timestep to ensure computation stability 
        of the flow processes. This CFL-like condition is computed using erodibility 
        coefficients, discharges plus elevations and distances between TIN nodes and 
        their respective reveivers. 
        
        Parameters
        ----------
        variable : xy
            Numpy arrays containing the XY coordinates of the TIN nodes.
            
        variable : elev
            Numpy arrays containing the elevation of the TIN nodes.
            
        variable: locIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """
        
        # Initialise MPI communications
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        # Compute the local value for time stability
        dt = FLOWalgo.flowcompute.flowcfl(locIDs,self.receivers,xy,elev,self.discharge, \
                                      self.erodibility,self.m,self.n)
        
        # Global mimimum value for diffusion stability
        CFL = numpy.zeros(1)
        CFL[0] = dt
        comm.Allreduce(MPI.IN_PLACE,CFL,op=MPI.MIN)
        self.CFL = CFL[0]