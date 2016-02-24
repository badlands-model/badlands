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

import SFDalgo as SFD
import FLOWalgo
import FLWnetwork

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
        
        self.sendprocID = None
        self.rcvprocID = None
        self.rcvprocNb = None
        self.rcvsendID = None
        self.localGIDs = None
        self.localNodes = 0
        
        self.discharge = None
        #self.drainage_area = None
        
        return
    
    def SFD_receivers(self, fillH, neighbours, globalIDs):
        """
        Single Flow Direction function computes downslope flow directions by inspecting the neighborhood 
        elevations around each node. The SFD method assigns a unique flow direction towards the steepest 
        downslope neighbor.
        
        Parameters
        ----------
        variable : fillH
            Numpy array containing the filled elevations from Planchon & Darboux depression-less algorithm.
            
        variable : neighbours
            Numpy integer-type array with the neighbourhood IDs.
            
        variable: globalIDs
            Numpy integer-type array containing for local nodes their global IDs.
        """
       
        # Initialise MPI communications
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()
        
        # Call the SFD function from libUtils
        base, receivers = SFD.sfdcompute.directions( fillH, neighbours, globalIDs)
        
        # Send local base level globally
        comm.Allreduce(MPI.IN_PLACE,base,op=MPI.MAX)
        self.base = numpy.where(base >= 0)[0]

        # Send local receivers globally
        comm.Allreduce(MPI.IN_PLACE,receivers,op=MPI.MAX)
        self.receivers = receivers
        
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
    
    def compute_flow(self, Acell, rain):
        """
        Calculates the drainage area and water discharge at each node.
        
        Parameters
        ----------
        variable : Acell
            Numpy float-type array containing the voronoi area for each nodes (in m2)
            
        variable : rain
            Numpy float-type array containing the precipitation rate for each nodes (in m/a).
        """
        
        if self.stack is None:
            self.ordered_node_array()
            
        numPts = len(self.stack)
        
        #self.drainage_area = numpy.zeros(numPts, dtype=float) + Acell
        self.discharge = numpy.zeros(numPts, dtype=float) 
        self.discharge[self.stack] = Acell[self.stack]*rain[self.stack]

        # Compute drainage area and discharge using libUtils
        self.discharge = FLOWalgo.flowcompute.discharge(self.stack,self.receivers,self.discharge)
        
        return     
    