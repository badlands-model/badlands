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
from pyBadlands.libUtils import ORmodel
from scipy.ndimage.filters import gaussian_filter
from scipy import interpolate
from scipy.spatial import cKDTree

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

    bool : orographic
        Numpy boolean array defining orographic calculation if any.

    float : rbgd
        Numpy array of background precipitation.

    float : rmin
        Numpy array of minimal precipitation.

    float : rmax
        Numpy array of maximal precipitation.

    float : windx
        Numpy array of wind velocity along X.

    float : windy
        Numpy array of wind velocity along Y.

    float : tauc
        Numpy array of time conversion from cloud water to hydrometeors.

    float : tauf
        Numpy array of time for hydrometeor fallout.

    float : nm
        Numpy array of moist stability frequency.

    float : cw
        Numpy array of uplift sensitivity factor.

    float : hw
        Numpy array of depth of the moist layer.

    float : ortime
        Numpy array of rain computation time step.

    string : MapDisp
        Numpy array containing the cumulative displacement map file names.

    float : TimeDisp
        Numpy array containing the start and end times for each displacement period in years.

    float : regX
        Numpy array containing the X-coordinates of the regular input grid.

    float : regY
        Numpy array containing the Y-coordinates of the regular input grid.

    float : rivPos
        Numpy array containing the XY position of the rivers.

    float : rivTime
        Numpy array containing the active time for the rivers.

    float : rivQws
        Numpy array containing the water and sediment discharge for the rivers.

    float : rivWth
        Numpy array containing the width of the rivers.

    int : rivNb
        Number of rivers.

    float : Tdisplay
        Display interval (in years).
    """

    def __init__(self, seafile = None, sea0 = 0., MapRain = None, TimeRain = None, ValRain = None,
                 orographic = None, rbgd = None, rmin = None, rmax = None, windx = None, windy = None,
                 tauc = None, tauf = None, nm = None, cw = None, hw = None, ortime = None, MapDisp = None,
                 TimeDisp = None, regX = None, regY = None, rivPos = None, rivTime = None, rivQws = None,
                 rivWth = None, rivNb = 0, Tdisplay = 0.):

        self.regX = regX
        self.regY = regY
        self.xi, self.yi = numpy.meshgrid(regX, regY, indexing='xy')
        self.xyi = numpy.dstack([self.xi.flatten(), self.yi.flatten()])[0]
        self.tree = None
        self.dx = None

        self.rivNb = rivNb
        self.rivPos = rivPos
        self.rivTime = rivTime
        self.rivQws = rivQws
        self.rivWth = rivWth
        self.rivQs = None
        self.rivQw = None

        self.Map_rain = MapRain
        self.rainVal = ValRain
        self.T_rain = TimeRain
        self.orographic = orographic
        self.rbgd = rbgd
        self.rmin = rmin
        self.rmax = rmax
        self.windx = windx
        self.windy = windy
        self.tauc = tauc
        self.tauf = tauf
        self.nm = nm
        self.cw = cw
        self.hw = hw
        self.ortime = ortime
        self.next_rain = None

        self.Map_disp = MapDisp
        self.T_disp = TimeDisp
        self.injected_disps = None
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
        self.next_flexure = None
        self.next_layer = None
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
        seadata = pandas.read_csv(self.seafile, sep=r'\s+', engine='c',
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

    def getRivers(self, time):
        """
        Finds for a given time the active rivers and allocates corresponding points with
        water and sediment discharge values.

        Parameters
        ----------
        float : time
            Requested time for which to compute sea level elevation.

        Return
        ----------
        variable: rivQw
            Numpy array containing flow discharge from rivers.

        variable: rivQs
            Numpy array containing sediment discharge from rivers.
        """

        self.rivQw = numpy.zeros(len(self.tXY))
        self.rivQs = numpy.zeros(len(self.tXY))

        if self.rivNb > 0:
            active = numpy.where(numpy.logical_and(self.rivTime[:,0] <= time, self.rivTime[:,1] > time))[0]
            rivNb = len(active)
            if rivNb > 0:
                riv_xy = self.rivPos[active,:]
                radius = self.rivWth[active]
                indices = self.tree.query_ball_point(riv_xy, radius)
                for r in range(rivNb):
                    riv_qw = self.rivQws[active[r],0]
                    riv_qs = self.rivQws[active[r],1]
                    if len(indices[r]) > 0:
                        self.rivQw[indices[r]] += riv_qw / len(indices[r])
                        self.rivQs[indices[r]] += riv_qs / len(indices[r])
                    else:
                        distances, ids = self.tree.query(riv_xy[r], k=1)
                        self.rivQw[ids] += riv_qw
                        self.rivQs[ids] += riv_qs

        return

    def update_force_TIN(self, tXY):
        """
        Update TIN variables after 3D displacements.

        variable : tXY
            Numpy float-type array containing the coordinates for each nodes in the TIN (in m2)
        """
        self.tXY = tXY
        self.tree = cKDTree(self.tXY)
        self.dx = self.tXY[1,0] - self.tXY[0,0]

        return

    def get_Rain(self, time, elev, inIDs):
        """
        Get rain value for a given period and perform interpolation from regular grid to unstructured TIN one.

        Parameters
        ----------
        float : time
            Requested time interval rain map to load.

        float : elev
            Unstructured grid (TIN) Z coordinates.

        integer : inDs
            List of unstructured vertices contained in each partition.

        Return
        ----------
        variable: tinRain
            Numpy array containing the updated rainfall for the local domain.
        """

        events = numpy.where( (self.T_rain[:,1] - time) <= 0)[0]
        event = len(events)

        if not (time >= self.T_rain[event,0]) and not (time < self.T_rain[event,1]):
            raise ValueError('Problem finding the rain map to load!')

        if self.orographic[event]:
            tinRain = self.build_OrographicRain_map(event, elev, inIDs)
            self.next_rain = min(time + self.ortime[event], self.T_rain[event,1])
        elif self.Map_rain[event] == None:
            tinRain = numpy.zeros(len(self.tXY[inIDs,0]), dtype=float)
            tinRain = self.rainVal[event]
            self.next_rain = self.T_rain[event,1]
        else:
            rainMap = pandas.read_csv(str(self.Map_rain[event]), sep=r'\s+', engine='c',
                               header=None, na_filter=False, dtype=numpy.float, low_memory=False)

            rectRain = numpy.reshape(rainMap.values,(len(self.regX), len(self.regY)),order='F')
            tinRain = interpolate.interpn( (self.regX, self.regY), rectRain, self.tXY[inIDs,:], method='linear')
            self.next_rain = self.T_rain[event,1]

        return tinRain

    def build_OrographicRain_map(self, event, elev, inIDs):
        """
        Build rain map using SMith & Barstad (2004) model for a given period and perform interpolation from regular grid to
        unstructured TIN one.

        Parameters
        ----------
        float : event
            rain event number.

        float : elev
            Unstructured grid (TIN) Z coordinates.

        integer : inDs
            List of unstructured vertices contained in each partition.

        Return
        ----------
        variable: tinRain
            Numpy array containing the updated rainfall for the local domain.
        """

        # Interpolate elevation on regular grid
        distances, indices = self.tree.query(self.xyi, k=8)
        if len(elev[indices].shape) == 3:
            elev_vals = elev[indices][:,:,0]
        else:
            elev_vals = elev[indices]
        oelev = numpy.average(elev_vals,weights=(1./distances), axis=1)
        onIDs = numpy.where(distances[:,0] == 0)[0]
        if len(onIDs) > 0:
            oelev[onIDs] = elev[indices[onIDs,0]]
        oelev -= self.sealevel
        oelev = oelev.clip(0)
        regZ = numpy.reshape(oelev,(len(self.regX), len(self.regY)),order='F')
        # Use Smith & Barstad model
        rectRain = ORmodel.orographicrain.compute(regZ, self.dx, self.windx[event], self.windy[event],
            self.rmin[event], self.rmax[event], self.rbgd[event], self.nm[event], self.cw[event],
            self.hw[event], self.tauc[event], self.tauf[event])

        # Apply smoothing here
        smthRain = gaussian_filter(rectRain, sigma=3)

        # Interpolate
        tinRain = interpolate.interpn( (self.regX, self.regY), smthRain, self.tXY[inIDs,:], method='linear')

        return tinRain

    def disp_border(self, disp, neighbours, edge_length, boundPts):
        """
        This function defines the displacement of the TIN edges.

        Parameters
        ----------
        variable : disp
            Numpy arrays containing the internal nodes displacement value.

        variable : neighbours
            Numpy integer-type array containing for each nodes its neigbhours IDs.

        variable : edge_length
            Numpy float-type array containing the lengths to each neighbour.

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
                lselect = edge_length[id,ids]
                picked = numpy.argmin(lselect)
                disp[id] = disp[ngbhs[ids[picked]]]
            else:
                missedPts = numpy.append(missedPts,id)

        if len(missedPts) > 0 :
            for p in range(len(missedPts)):
                id = int(missedPts[p])
                ngbhs = neighbours[id,:]
                ids = numpy.where((disp[ngbhs] < 9.e6) & (ngbhs >= 0))[0]
                if len(ids) == 0:
                    raise ValueError('Error while getting boundary displacement for point ''%d''.' % id)
                lselect = edge_length[id,ids]
                picked = numpy.argmin(lselect)
                disp[id] = disp[ngbhs[ids[picked]]]

        return disp

    def load_Tecto_map(self, time, inIDs):
        """
        Load vertical displacement map for a given period and perform interpolation from regular grid to unstructured TIN one.

        Parameters
        ----------
        float : time
            Requested time interval rain map to load.

        integer : inDs
            List of unstructured vertices contained in each partition.

        Return
        ----------
        variable: tinDisp
            Numpy array containing the updated displacement rate for the local domain.
        """

        events = numpy.where( (self.T_disp[:,1] - time) <= 0)[0]
        event = len(events)

        if not (time >= self.T_disp[event,0]) and not (time < self.T_disp[event,1]):
            raise ValueError('Problem finding the displacements map to load!')

        self.next_disp = self.T_disp[event,1]

        if self.injected_disps is not None or self.Map_disp[event] != None:
            if self.injected_disps is not None:
                dispMap = self.injected_disps
            else:
                dispMap = pandas.read_csv(str(self.Map_disp[event]), sep=r'\s+', engine='c', header=None, na_filter=False, \
                                   dtype=numpy.float, low_memory=False).values


            rectDisp = numpy.reshape(dispMap,(len(self.regX), len(self.regY)),order='F')
            tinDisp = interpolate.interpn( (self.regX, self.regY), rectDisp, self.tXY[inIDs,:], method='linear')
            dt = (self.T_disp[event,1] - self.T_disp[event,0])
            if dt <= 0:
                raise ValueError('Problem computing the displacements rate for event %d.'%event)
            tinDisp = tinDisp / dt
        else:
            tinDisp = numpy.zeros(len(self.tXY[:,0]), dtype=float)

        return tinDisp

    def load_Disp_map(self, time, tXY, inIDs, strata=False, sXY=None, insIDs=None):
        """
        Load 3D displacements map for a given period and perform interpolation from regular grid to unstructured TIN one.

        Parameters
        ----------
        float : time
            Requested time interval rain map to load.

        float : tXY
            Unstructured grid (TIN) XY coordinates.

        integer : inIDs
            List of unstructured vertices contained in each partition.

        boolean : strata
            Stratigraphic module flag.

        float : sXY
            Stratigraphic regular grid XY coordinates.

        integer : insIDs
            List of stratigraphic vertices contained in each partition.
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
        dpXY = tXY[inIDs,:]

        if strata:
            totsPts = len(sXY[:,0])
            sdispX = numpy.zeros(totsPts, dtype=float)
            sdispY = numpy.zeros(totsPts, dtype=float)
            sdispZ = numpy.zeros(totsPts, dtype=float)
            dpsXY = sXY[insIDs,:]

        events = numpy.where( (self.T_disp[:,1] - time) <= 0)[0]
        event = len(events)

        if not (time >= self.T_disp[event,0]) and not (time < self.T_disp[event,1]):
            raise ValueError('Problem finding the 3D displacements map to load!')

        if self.time3d > 0.:
            if time + self.time3d > self.T_disp[event,1]:
                self.next_disp = self.T_disp[event,1]
            else:
                self.next_disp = self.time3d + time
        else:
            self.next_disp = self.T_disp[event,1]

        update = False
        if self.injected_disps is not None or self.Map_disp[event] != None:
            dispX.fill(-1.e6)
            dispY.fill(-1.e6)
            dispZ.fill(-1.e6)

            if self.injected_disps is not None:
                dvals = self.injected_disps
            else:
                disps = pandas.read_csv(str(self.Map_disp[event]), sep=r'\s+', engine='c', header=None, na_filter=False, \
                            dtype=numpy.float, low_memory=False)
                dvals = disps.values


            disprX = numpy.reshape(dvals[:,0],(len(self.regX), len(self.regY)),order='F')
            disprY = numpy.reshape(dvals[:,1],(len(self.regX), len(self.regY)),order='F')
            disprZ = numpy.reshape(dvals[:,2],(len(self.regX), len(self.regY)),order='F')
            dispX[inIDs] = interpolate.interpn( (self.regX, self.regY), disprX, dpXY, method='linear')
            dispY[inIDs] = interpolate.interpn( (self.regX, self.regY), disprY, dpXY, method='linear')
            dispZ[inIDs] = interpolate.interpn( (self.regX, self.regY), disprZ, dpXY, method='linear')
            comm.Allreduce(mpi.IN_PLACE, dispX, op=mpi.MAX)
            comm.Allreduce(mpi.IN_PLACE, dispY, op=mpi.MAX)
            comm.Allreduce(mpi.IN_PLACE, dispZ, op=mpi.MAX)
            update = True

            if strata:
                sdispX.fill(-1.e6)
                sdispY.fill(-1.e6)
                sdispX[insIDs] = interpolate.interpn( (self.regX, self.regY), disprX, dpsXY, method='linear')
                sdispY[insIDs] = interpolate.interpn( (self.regX, self.regY), disprY, dpsXY, method='linear')
                comm.Allreduce(mpi.IN_PLACE, sdispX, op=mpi.MAX)
                comm.Allreduce(mpi.IN_PLACE, sdispY, op=mpi.MAX)

        if self.time3d > 0. and (self.injected_disps is not None or self.Map_disp[event] != None):
            rate = (self.next_disp - time) / (self.T_disp[event,1] - self.T_disp[event,0])
            assert rate > 0
            dispX = dispX * rate
            dispY = dispY * rate
            dispZ = dispZ * rate
            if strata:
                sdispX = sdispX * rate
                sdispY = sdispY * rate

        self.dispX = dispX
        self.dispY = dispY
        self.dispZ = dispZ

        if strata:
            return update, sdispX, sdispY
        else:
            return update

    def apply_XY_dispacements(self, area, fixIDs, telev, tcum, tflex=None, scum=None,
                                      Te=None, Ke=None, flexure=0, strat=0, ero=0):
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

        float : tcum
            Numpy array with erosion/deposition values from previous TIN nodes.

        float : tflex
            Numpy array with cumulative flexural values from previous TIN nodes.

        float : scum
            Numpy array with erosion/deposition used for stratal mesh.

        float : Te
            Numpy array with thickness used for erosional mesh.

        float : Ke
            Numpy array with erodibility used for erosional mesh.

        integer : flexure
            Integer flagging flexural isostasy.

        integer : strat
            Integer flagging stratigraphic mesh model.

        integer : ero
            Integer flagging erosional mesh model.

        Return
        ----------
        variable: tinMesh
            Delaunay mesh generated after displacements.

        variable: newelev
            Numpy array containing the updated elevation for the new TIN.

        variable: newcum
            Numpy array containing the updated erosion/deposition values for the new TIN.

        variable: newcumf
            Numpy array containing the updated cumulative flexural values for the new TIN.

        variable: newscum
            Numpy array containing the updated erosion/deposition values used in the stratal mesh.

        variable: newKe
            Numpy array containing the updated erodibility values for the new TIN.

        variable: newTe
            Numpy array containing the updated thickness values used in the erosional mesh.

        """

        # Apply displacements to TIN points (excluding boundary points)
        telev += self.dispZ
        tXY = numpy.copy(self.tXY)
        tXY[fixIDs:,0] += self.dispX[fixIDs:]
        tXY[fixIDs:,1] += self.dispY[fixIDs:]

        # Specify inside simulation area parameters
        dx = tXY[1,0] - tXY[0,0]
        minX = min(tXY[:fixIDs,0]) + dx
        minY = min(tXY[:fixIDs,1]) + dx
        maxX = max(tXY[:fixIDs,0]) - dx
        maxY = max(tXY[:fixIDs,1]) - dx

        # Find points outside the area
        xID = numpy.where(numpy.logical_or( tXY[fixIDs:,0] <= minX, tXY[fixIDs:,0] >= maxX))[0]
        yID = numpy.where(numpy.logical_or( tXY[fixIDs:,1] <= minY, tXY[fixIDs:,1] >= maxY))[0]
        tIDs = numpy.concatenate((xID, yID), axis=0)
        tID = numpy.unique(tIDs)
        tID += fixIDs

        # Delete outside domain points if any
        if len(tID) > 0:
            self.tXY = numpy.delete(self.tXY, tID, 0)
            elev = numpy.delete(telev, tID, 0)
            cum = numpy.delete(tcum, tID, 0)
            if flexure == 1:
                cumf = numpy.delete(tflex, tID, 0)
            if strat == 1:
                stcum = numpy.delete(scum, tID, 0)
            if ero == 1:
                lay = Ke.shape()[1]
                mKe = numpy.zeros((len(cum),lay))
                mTe = numpy.zeros((len(cum),lay))
                for k in range(lay):
                    mKe[:,k] = numpy.delete(Ke[:,k], tID, 0)
                    mTe[:,k] = numpy.delete(Te[:,k], tID, 0)
        else:
            self.tXY = tXY
            elev = telev
            cum = tcum
            if flexure == 1:
                cumf = tflex
            if strat == 1:
                stcum = scum
            if ero == 1:
                lay = Ke.shape()[1]
                mKe = Ke
                mTe = Te

        # Create KDTree with deformed points and find points which needs to be merged
        tree = cKDTree(self.tXY)
        pairs = tree.query_pairs(self.merge3d)
        # For points which require merging define a new point and
        # interpolate parameters based on merged points
        if len(pairs) > 0:
            pairIDs = numpy.array(list(pairs))
            nonfixIDs = numpy.where(numpy.logical_and( pairIDs[:,0] >= fixIDs, pairIDs[:,1] >= fixIDs))[0]
            mXY = numpy.empty(shape=[len(nonfixIDs),2], dtype=float)
            mXY[:,0] = 0.5*(self.tXY[pairIDs[nonfixIDs,0],0] + self.tXY[pairIDs[nonfixIDs,1],0])
            mXY[:,1] = 0.5*(self.tXY[pairIDs[nonfixIDs,0],1] + self.tXY[pairIDs[nonfixIDs,1],1])
            mergedIDs = numpy.unique(pairIDs[nonfixIDs,:].flatten())
            distances, indices = tree.query(mXY, k=3)
            weights = 1.0 / distances**2
            onIDs = numpy.where(distances[:,0] == 0)[0]
            if len(onIDs) > 0:
                raise ValueError('Problem: IDs after merging is on previous vertex position.')
            if len(elev[indices].shape) == 3:
                z_vals = elev[indices][:,:,0]
                cum_vals = cum[indices][:,:,0]
                if flexure == 1:
                    cumf_vals = cumf[indices][:,:,0]
                if strat == 1:
                    scum_vals = stcum[indices][:,:,0]
            else:
                z_vals = elev[indices]
                cum_vals = cum[indices]
                if flexure == 1:
                    cumf_vals = cumf[indices]
                if strat == 1:
                    scum_vals = stcum[indices]
            z_avg = numpy.average(z_vals, weights=weights,axis=1)
            cum_avg = numpy.average(cum_vals, weights=weights,axis=1)
            if flexure == 1:
                cumf_avg = numpy.average(cumf_vals, weights=weights,axis=1)
            if strat == 1:
                scum_avg = numpy.average(scum_vals, weights=weights,axis=1)
            if ero == 1:
                Te_avg = numpy.zeros((len(z_avg),lay))
                Ke_avg = numpy.zeros((len(z_avg),lay))
                for k in range(lay):
                    if len(mTe[indices,k].shape) == 3:
                        Te_vals = mTe[indices,k][:,:,0]
                        Ke_vals = mKe[indices,k][:,:,0]
                    else:
                        Te_vals = mTe[indices,k]
                        Ke_vals = mKe[indices,k]
                    Te_avg[:,k] = numpy.average(Te_vals, weights=weights,axis=1)
                    Ke_avg[:,k] = Ke_vals[indices[0]]
            # Delete points that have been merged
            newXY = numpy.delete(self.tXY, mergedIDs, 0)
            newelev = numpy.delete(elev, mergedIDs, 0)
            newcum = numpy.delete(cum, mergedIDs, 0)
            if flexure == 1:
                newcumf = numpy.delete(cumf, mergedIDs, 0)
            if strat == 1:
                newscum = numpy.delete(stcum, mergedIDs, 0)
            if ero == 1:
                newKe = numpy.zeros((len(newcum),lay))
                newTe = numpy.zeros((len(newcum),lay))
                for k in range(lay):
                    newKe[:,k] = numpy.delete(mKe[:,k], mergedIDs, 0)
                    newTe[:,k] = numpy.delete(mTe[:,k], mergedIDs, 0)
            # Add new points to the deformed TIN
            newXY = numpy.concatenate((newXY, mXY), axis=0)
            newelev = numpy.concatenate((newelev, z_avg), axis=0)
            newcum = numpy.concatenate((newcum, cum_avg), axis=0)
            if flexure == 1:
                newcumf = numpy.concatenate((newcumf, cumf_avg), axis=0)
            if strat == 1:
                newscum = numpy.concatenate((newscum, scum_avg), axis=0)
            if ero == 1:
                nKe = numpy.zeros((len(newcum),lay))
                nTe = numpy.zeros((len(newcum),lay))
                for k in range(lay):
                    nKe[:,k] = numpy.concatenate((newKe[:,k], Ke_avg[:,k]), axis=0)
                    nTe[:,k] = numpy.concatenate((newTe[:,k], Te_avg[:,k]), axis=0)
        else:
            newXY = self.tXY
            newelev = elev
            newcum = cum
            if flexure == 1:
                newcumf = cumf
            if strat == 1:
                newscum = stcum
            if ero == 1:
                nKe = mKe
                nTe = mTe

        # Based on new points build the triangulation
        newTIN = triangle.triangulate( dict(vertices=newXY),'Da'+str(area))
        # If some points have been added during the triangulation update the TIN
        # and interpolate neighbouring paramters to these new points
        if len(newTIN['vertices'][:,0]) > len(newXY[:,0]):
            addPts = newTIN['vertices'][len(newXY[:,0]):,:2]
            dist, ids = tree.query(addPts, k=3)
            weights = 1.0/dist**2
            if len(elev[ids].shape) == 3:
                zvals = elev[ids][:,:,0]
                cumvals = cum[ids][:,:,0]
                if flexure == 1:
                    cumfvals = cumf[ids][:,:,0]
                if strat == 1:
                    scumvals = stcum[ids][:,:,0]
            else:
                zvals = elev[ids]
                cumvals = cum[ids]
                if flexure == 1:
                    cumfvals = cumf[ids]
                if strat == 1:
                    scumvals = stcum[ids]
            if ero == 1:
                Tevals = numpy.zeros((len(zvals),lay))
                Kevals = numpy.zeros((len(zvals),lay))
                for k in range(lay):
                    if len(mTe[ids,k].shape) == 3:
                        Tevals[:,k] = mTe[ids,k][:,:,0]
                        Kevals[:,k] = mKe[ids,k][:,:,0]
                    else:
                        Tevals[:,k] = mTe[ids,k]
                        Kevals[:,k] = mKe[ids,k]
            zavg = numpy.average(zvals, weights=weights,axis=1)
            cumavg = numpy.average(cumvals, weights=weights,axis=1)
            if flexure == 1:
                cumfavg = numpy.average(cumfvals, weights=weights,axis=1)
            if strat == 1:
                scumavg = numpy.average(scumvals, weights=weights,axis=1)
            if ero == 1:
                Teavg = numpy.zeros((len(zavg),lay))
                Keavg = numpy.zeros((len(zavg),lay))
                for k in range(lay):
                    Teavg[:,k] = numpy.average(Tevals[:,k], weights=weights,axis=1)
                    Keavg[:,k] = Kevals[indices[0],k]
            newelev = numpy.concatenate((newelev, zavg), axis=0)
            newcum = numpy.concatenate((newcum, cumavg), axis=0)
            if flexure == 1:
                newcumf = numpy.concatenate((newcumf, cumfavg), axis=0)
            if strat == 1:
                newscum = numpy.concatenate((newscum, scumavg), axis=0)
            if ero == 1:
                newTe = numpy.zeros((len(newelev),lay))
                newKe = numpy.zeros((len(newelev),lay))
                for k in range(lay):
                    newTe[:,k] = numpy.concatenate((nTe[:,k], Teavg[:,k]), axis=0)
                    newKe[:,k] = numpy.concatenate((nKe[:,k], Keavg[:,k]), axis=0)
        elif len(newTIN['vertices'][:,0]) < len(newXY):
            raise ValueError('Problem building the TIN after 3D displacements.')

        if flexure == 0:
            newcumf = None
        if strat == 0:
            newscum = None
        if ero == 0:
            newKe = None
            newTe = None

        return newTIN, newelev, newcum, newcumf, newscum, newKe, newTe
