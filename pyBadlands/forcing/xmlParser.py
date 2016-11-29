##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module encapsulates parsing functions of Badlands XmL input file.
"""
import os
import glob
import numpy
import shutil
import xml.etree.ElementTree as ET

class xmlParser:
    """
    This class defines XmL input file variables.

    Parameters
    ----------
    string : inputfile
        The XmL input file name.
    """

    def __init__(self, inputfile = None, makeUniqueOutputDir=True):
        """
        If makeUniqueOutputDir is set, we create a uniquely-named directory for
        the output. If it's clear, we blindly accept what's in the XML file.
        """

        if inputfile==None:
            raise RuntimeError('XmL input file name must be defined to run a Badlands simulation.')
        if not os.path.isfile(inputfile):
            raise RuntimeError('The XmL input file name cannot be found in your path.')
        self.inputfile = inputfile

        self.demfile = None
        self.btype = 'slope'
        self.perc_dep = 0.
        self.slp_cr = 0.
        self.fillmax = 1.
        self.Afactor = 1
        self.nopit = 0

        self.restart = False
        self.rForlder = None
        self.rStep = 0
        self.tStart = None
        self.tEnd = None
        self.tDisplay = None
        self.minDT = 1.
        self.maxiDT = 1.e6

        self.stratdx = 0.
        self.laytime = 0.
        self.region = 0
        self.llcXY = None
        self.urcXY = None

        self.seapos = 0.
        self.sealimit = 100.
        self.seafile = None

        self.disp3d = False
        self.tectNb = None
        self.tectTime = None
        self.tectFile = None
        self.merge3d = None
        self.time3d = None

        self.riverNb = None
        self.riverTime = None
        self.riverPos = None
        self.riverQws = None
        self.riverWidth = None

        self.rainNb = None
        self.rainVal = None
        self.rainTime = None
        self.oroRain = False
        self.orographic = None
        self.ortime = None
        self.rbgd = None
        self.rmin = None
        self.rmax = None
        self.windx = None
        self.windy = None
        self.tauc = 1000.
        self.tauf = 1000.
        self.nm = 0.005
        self.cw = 0.005
        self.hw = 3000.

        self.depo = 1
        self.SPLm = 0.5
        self.SPLn = 1.
        self.SPLero = 0.
        self.maxDT = None
        self.alluvial = 0.
        self.bedrock = 0.
        self.esmooth = None
        self.dsmooth = None

        self.spl = False
        self.capacity = False
        self.filter = False
        self.Hillslope = False
        self.nHillslope = False

        self.CDa = 0.
        self.CDm = 0.
        self.Sc = 0.
        self.makeUniqueOutputDir = makeUniqueOutputDir

        self.outDir = None
        self.sh5file = 'h5/sed'

        self.th5file = 'h5/tin.time'
        self.txmffile = 'xmf/tin.time'
        self.txdmffile = 'tin.series.xdmf'

        self.fh5file = 'h5/flow.time'
        self.fxmffile = 'xmf/flow.time'
        self.fxdmffile = 'flow.series.xdmf'

        self.flexure = False
        self.ftime = None
        self.fnx = None
        self.fny = None
        self.dmantle = None
        self.dsediment = None
        self.youndMod = None
        self.elasticH = None
        self.elasticGrid = None
        self.flexbounds = []

        self.erolays = None
        self.eroMap = None
        self.eroVal = None
        self.thickMap = None
        self.thickVal = None

        self._get_XmL_Data()

        return

    def _get_XmL_Data(self):
        """
        Main function used to parse the XmL input file.
        """

        # Load XmL input file
        tree = ET.parse(self.inputfile)
        root = tree.getroot()

        # Extract grid structure information
        grid = None
        grid = root.find('grid')
        if grid is not None:
            element = grid.find('demfile')
            if element is not None:
                self.demfile = element.text
                if not os.path.isfile(self.demfile):
                    raise ValueError('DEM file is missing or the given path is incorrect.')
            element = None
            element = grid.find('boundary')
            if element is not None:
                self.btype = element.text
                if self.btype != 'slope' and self.btype != 'flat' and self.btype != 'wall':
                    raise ValueError('Error in the definition of the grid structure: Boundary type is either: flat, slope or wall')
            else:
                self.btype = 'slope'
            element = None
            element = grid.find('resfactor')
            if element is not None:
                self.Afactor = int(element.text)
                if self.Afactor < 1:
                    self.Afactor = 1
            else:
                self.Afactor = 1
            element = None
            element = grid.find('nopit')
            if element is not None:
                self.nopit = int(element.text)
                if self.nopit < 1:
                    self.nopit = 0
                else:
                    self.nopit = 1
            else:
                self.nopit = 0
        else:
            raise ValueError('Error in the XmL file: grid structure definition is required!')

        # Extract time structure information
        time = None
        time = root.find('time')
        if time is not None:
            element = None
            element = time.find('restart')
            if element is not None:
                self.restart = True
            if self.restart:
                for rst in time.iter('restart'):
                    element = None
                    element = rst.find('rfolder')
                    if element is not None:
                        self.rfolder = element.text
                    else:
                        raise ValueError('Error the restart folder name needs to be defined for the restart function to work')
                    element = None
                    element = rst.find('rstep')
                    if element is not None:
                        self.rstep = int(element.text)
                    else:
                        raise ValueError('Error the restart step number needs to be defined for the restart function to work')
            element = None
            element = time.find('start')
            if element is not None:
                self.tStart = float(element.text)
            else:
                raise ValueError('Error in the definition of the simulation time: start time declaration is required')
            element = None
            element = time.find('end')
            if element is not None:
                self.tEnd = float(element.text)
            else:
                raise ValueError('Error in the definition of the simulation time: end time declaration is required')
            if self.tStart > self.tEnd:
                raise ValueError('Error in the definition of the simulation time: start time is greater than end time!')
            element = None
            element = time.find('display')
            if element is not None:
                self.tDisplay = float(element.text)
            else:
                raise ValueError('Error in the definition of the simulation time: display time declaration is required')
            if (self.tEnd - self.tStart) % self.tDisplay != 0:
                raise ValueError('Error in the definition of the simulation time: display time needs to be a multiple of simulation time.')
            element = None
            element = time.find('mindt')
            if element is not None:
                self.minDT = float(element.text)
            else:
                self.minDT = 1.
            element = None
            element = time.find('maxdt')
            if element is not None:
                self.maxiDT = float(element.text)
            else:
                self.maxiDT = self.tDisplay
        else:
            raise ValueError('Error in the XmL file: time structure definition is required!')

        # Extract stratigraphic structure information
        strat = None
        strat = root.find('strata')
        if strat is not None:
            element = None
            element = strat.find('stratdx')
            if element is not None:
                self.stratdx = float(element.text)
            else:
                self.stratdx = 0.
            element = None
            element = strat.find('laytime')
            if element is not None:
                self.laytime = float(element.text)
            else:
                self.laytime = self.tDisplay
            if self.laytime >  self.tDisplay:
                 self.laytime = self.tDisplay
            if self.tDisplay % self.laytime != 0:
                raise ValueError('Error in the XmL file: stratal layer interval needs to be an exact multiple of the display interval!')

            element = None
            element = strat.find('region')
            if element is not None:
                self.region = int(element.text)
                self.llcXY = numpy.empty((self.region,2))
                self.urcXY = numpy.empty((self.region,2))
            element = None
            id = 0
            for dom in strat.iter('domain'):
                element = None
                element = dom.find('llcx')
                self.llcXY[id,0] = float(element.text)
                element = None
                element = dom.find('llcy')
                self.llcXY[id,1] = float(element.text)
                element = None
                element = dom.find('urcx')
                self.urcXY[id,0] = float(element.text)
                element = None
                element = dom.find('urcy')
                self.urcXY[id,1] = float(element.text)
                id += 1
            if id != self.region:
                raise ValueError('Number of region %d does not match with the number of declared region parameters %d.' %(self.region,id))

        # Extract sea-level structure information
        sea = None
        sea = root.find('sea')
        if sea is not None:
            element = None
            element = sea.find('position')
            if element is not None:
                self.seapos = float(element.text)
            else:
                self.seapos = 0.
            element = None
            element = sea.find('limit')
            if element is not None:
                self.sealimit = float(element.text)
            else:
                self.sealimit = 100.
            element = None
            element = sea.find('curve')
            if element is not None:
                self.seafile = element.text
                if not os.path.isfile(self.seafile):
                    raise ValueError('Sea level file is missing or the given path is incorrect.')
            else:
                self.seafile = None
        else:
            self.seapos = 0.
            self.sealimit = 100.
            self.seafile = None

        # Extract Tectonic structure information
        tecto = None
        tecto = root.find('tectonic')
        if tecto is not None:
            element = None
            element = tecto.find('disp3d')
            if element is not None:
                tmp3d = int(element.text)
                if tmp3d == 0:
                    self.disp3d = False
                else:
                    self.disp3d = True
            else:
                self.disp3d = False
            element = None
            element = tecto.find('merge3d')
            if element is not None:
                self.merge3d = float(element.text)
            else:
                self.merge3d = 0.
            element = None
            element = tecto.find('time3d')
            if element is not None:
                self.time3d = float(element.text)
            else:
                self.time3d = 0.
            element = None
            element = tecto.find('events')
            if element is not None:
                tmpNb = int(element.text)
                if tmpNb == 0:
                    raise ValueError('The number of tectonic events needs to be at least 1 or comment the tectonic structure.')

            else:
                raise ValueError('The number of tectonic events needs to be defined.')
            tmpFile = numpy.empty(tmpNb,dtype=object)
            tmpTime = numpy.empty((tmpNb,2))
            id = 0
            for disp in tecto.iter('disp'):
                element = None
                element = disp.find('dstart')
                if element is not None:
                    tmpTime[id,0] = float(element.text)
                else:
                    raise ValueError('Displacement event %d is missing start time argument.'%id)
                element = None
                element = disp.find('dend')
                if element is not None:
                    tmpTime[id,1] = float(element.text)
                else:
                    raise ValueError('Displacement event %d is missing end time argument.'%id)
                if tmpTime[id,0] >= tmpTime[id,1]:
                    raise ValueError('Displacement event %d start and end time values are not properly defined.'%id)
                if id > 0:
                    if tmpTime[id,0] < tmpTime[id-1,1]:
                        raise ValueError('Displacement event %d start time needs to be >= than displacement event %d end time.'%(id,id-1))
                element = None
                element = disp.find('dfile')
                if element is not None:
                    tmpFile[id] = element.text
                    if not os.path.isfile(tmpFile[id]):
                        raise ValueError('Displacement file %s is missing or the given path is incorrect.'%(tmpFile[id]))
                else:
                    raise ValueError('Displacement event %d is missing file argument.'%id)
                id += 1
            if id != tmpNb:
                raise ValueError('Number of events %d does not match with the number of declared displacement parameters %d.' %(tmpNb,id))

            # Create continuous displacement series
            self.tectNb = tmpNb
            if tmpTime[0,0] > self.tStart:
                self.tectNb += 1
            for id in range(1,tmpNb):
                if tmpTime[id,0] > tmpTime[id-1,1]:
                    self.tectNb += 1
            if tmpTime[tmpNb-1,1] < self.tEnd:
                self.tectNb += 1
            self.tectFile = numpy.empty(self.tectNb,dtype=object)
            self.tectTime = numpy.empty((self.tectNb,2))
            id = 0
            if tmpTime[id,0] > self.tStart:
                self.tectFile[id] = None
                self.tectTime[id,0] = self.tStart
                self.tectTime[id,1] = tmpTime[0,0]
                id += 1
            self.tectFile[id] = tmpFile[0]
            self.tectTime[id,:] = tmpTime[0,:]
            id += 1
            for p in range(1,tmpNb):
                if tmpTime[p,0] > tmpTime[p-1,1]:
                    self.tectFile[id] = None
                    self.tectTime[id,0] = tmpTime[p-1,1]
                    self.tectTime[id,1] = tmpTime[p,0]
                    id += 1
                self.tectFile[id] = tmpFile[p]
                self.tectTime[id,:] = tmpTime[p,:]
                id += 1
            if tmpTime[tmpNb-1,1] < self.tEnd:
                self.tectFile[id] = None
                self.tectTime[id,0] = tmpTime[tmpNb-1,1]
                self.tectTime[id,1] = self.tEnd
        else:
            self.tectNb = 1
            self.tectTime = numpy.empty((self.tectNb,2))
            self.tectTime[0,0] = self.tEnd + 1.e5
            self.tectTime[0,1] = self.tEnd + 2.e5
            self.tectFile = numpy.empty((self.tectNb),dtype=object)
            self.tectFile = None

        # Extract Rivers structure information
        rivers = None
        rivers = root.find('rivers')
        if rivers is not None:
            element = None
            element = rivers.find('riverNb')
            if element is not None:
                self.riverNb = int(element.text)
            else:
                raise ValueError('The number of rivers needs to be defined.')
            self.riverTime = numpy.empty((self.riverNb,2))
            self.riverPos = numpy.empty((self.riverNb,2))
            self.riverQws = numpy.empty((self.riverNb,2))
            self.riverWidth = numpy.empty((self.riverNb))
            id = 0
            for riv in rivers.iter('river'):
                if id >= self.riverNb:
                    raise ValueError('The number of rivers does not match the value of riverNb.')
                element = None
                element = riv.find('rstart')
                if element is not None:
                    self.riverTime[id,0] = float(element.text)
                else:
                    raise ValueError('River %d is missing start time argument.'%id)
                element = None
                element = riv.find('rend')
                if element is not None:
                    self.riverTime[id,1] = float(element.text)
                else:
                    raise ValueError('River %d is missing end time argument.'%id)
                element = None
                element = riv.find('rposX')
                if element is not None:
                    self.riverPos[id,0] = float(element.text)
                else:
                    raise ValueError('River %d is missing X position.'%id)
                element = None
                element = riv.find('rposY')
                if element is not None:
                    self.riverPos[id,1] = float(element.text)
                else:
                    raise ValueError('River %d is missing Y position.'%id)
                element = None
                element = riv.find('rwidth')
                if element is not None:
                    self.riverWidth[id] = float(element.text)
                else:
                    raise ValueError('River %d is missing width.'%id)
                element = None
                element = riv.find('rQw')
                if element is not None:
                    self.riverQws[id,0] = float(element.text)
                else:
                    raise ValueError('River %d is missing water discharge.'%id)
                element = None
                element = riv.find('rQs')
                if element is not None:
                    self.riverQws[id,1] = float(element.text)
                else:
                    raise ValueError('River %d is missing sediment discharge.'%id)
                id += 1
        else:
            self.riverNb = 0

        # Extract Precipitation structure information
        precip = None
        precip = root.find('precipitation')
        if precip is not None:
            element = None
            element = precip.find('climates')
            if element is not None:
                tmpNb = int(element.text)
            else:
                raise ValueError('The number of climatic events needs to be defined.')
            tmpVal = numpy.empty(tmpNb)
            tmpMap = numpy.empty(tmpNb,dtype=object)
            tmpOro = numpy.empty(tmpNb,dtype=bool)
            tmpTime = numpy.empty((tmpNb,2))
            tmpoTime = numpy.empty(tmpNb)
            tmprbgd = numpy.empty(tmpNb)
            tmprmin = numpy.empty(tmpNb)
            tmprmax = numpy.empty(tmpNb)
            tmpwdx = numpy.empty(tmpNb)
            tmpwdy = numpy.empty(tmpNb)
            tmptauc = numpy.empty(tmpNb)
            tmptauf = numpy.empty(tmpNb)
            tmpnm = numpy.empty(tmpNb)
            tmpcw = numpy.empty(tmpNb)
            tmphw = numpy.empty(tmpNb)
            id = 0
            for clim in precip.iter('rain'):
                if id >= tmpNb:
                    raise ValueError('The number of climatic events does not match the number of defined climates.')
                element = None
                element = clim.find('rstart')
                if element is not None:
                    tmpTime[id,0] = float(element.text)
                else:
                    raise ValueError('Rain climate %d is missing start time argument.'%id)
                element = None
                element = clim.find('rend')
                if element is not None:
                    tmpTime[id,1] = float(element.text)
                else:
                    raise ValueError('Rain climate %d is missing end time argument.'%id)
                if tmpTime[id,0] >= tmpTime[id,1]:
                    raise ValueError('Rain climate %d start and end time values are not properly defined.'%id)
                if id > 0:
                    if tmpTime[id,0] < tmpTime[id-1,1]:
                        raise ValueError('Rain climate %d start time needs to be >= than rain climate %d end time.'%(id,id-1))
                element = None
                element = clim.find('map')
                if element is not None:
                    tmpMap[id] = element.text
                    if not os.path.isfile(tmpMap[id]):
                        raise ValueError('Rain map file %s is missing or the given path is incorrect.'%(tmpMap[id]))
                else:
                    tmpMap[id] = None
                element = None
                element = clim.find('rval')
                if element is not None:
                    tmpVal[id] = float(element.text)
                else:
                    tmpVal[id] = 0.
                element = None
                element = clim.find('ortime')
                if element is not None:
                    tmpoTime[id] = float(element.text)
                    tmpOro[id] = True
                    self.oroRain = True
                else:
                    tmpoTime[id] = 0.
                    tmpOro[id] = False
                element = None
                element = clim.find('rbgd')
                if element is not None:
                    tmprbgd[id] = float(element.text)
                else:
                    tmprbgd[id] = 0.
                element = None
                element = clim.find('rmin')
                if element is not None:
                    tmprmin[id] = float(element.text)
                else:
                    tmprmin[id] = 0.
                element = None
                element = clim.find('rmax')
                if element is not None:
                    tmprmax[id] = float(element.text)
                else:
                    tmprmax[id] = 0.
                element = None
                element = clim.find('windx')
                if element is not None:
                    tmpwdx[id] = float(element.text)
                else:
                    tmpwdx[id] = 0.
                element = None
                element = clim.find('windy')
                if element is not None:
                    tmpwdy[id] = float(element.text)
                else:
                    tmpwdy[id] = 0.
                element = None
                element = clim.find('tauc')
                if element is not None:
                    tmptauc[id] = float(element.text)
                else:
                    tmptauc[id] = 1000.
                element = None
                element = clim.find('tauf')
                if element is not None:
                    tmptauf[id] = float(element.text)
                else:
                    tmptauf[id] = 1000.
                element = None
                element = clim.find('nm')
                if element is not None:
                    tmpnm[id] = float(element.text)
                else:
                    tmpnm[id] = 0.005
                element = None
                element = clim.find('cw')
                if element is not None:
                    tmpcw[id] = float(element.text)
                else:
                    tmpcw[id] = 0.005
                element = None
                element = clim.find('hw')
                if element is not None:
                    tmphw[id] = float(element.text)
                else:
                    tmphw[id] = 3000.
                id += 1
            if id != tmpNb:
                raise ValueError('Number of climates %d does not match with the number of declared rain parameters %d.' %(tmpNb,id))

            # Create continuous precipitation series
            self.rainNb = tmpNb
            if tmpTime[0,0] > self.tStart:
                self.rainNb += 1
            for id in range(1,tmpNb):
                if tmpTime[id,0] > tmpTime[id-1,1]:
                    self.rainNb += 1
            if tmpTime[tmpNb-1,1] < self.tEnd:
                self.rainNb += 1
            self.rainVal = numpy.empty(self.rainNb)
            self.rainMap = numpy.empty(self.rainNb,dtype=object)
            self.rainTime = numpy.empty((self.rainNb,2))
            self.orographic = numpy.empty(self.rainNb,dtype=bool)
            self.ortime = numpy.empty(self.rainNb)
            self.rbgd = numpy.empty(self.rainNb)
            self.rmin = numpy.empty(self.rainNb)
            self.rmax = numpy.empty(self.rainNb)
            self.windx = numpy.empty(self.rainNb)
            self.windy = numpy.empty(self.rainNb)
            self.tauc = numpy.empty(self.rainNb)
            self.tauf = numpy.empty(self.rainNb)
            self.nm = numpy.empty(self.rainNb)
            self.cw = numpy.empty(self.rainNb)
            self.hw = numpy.empty(self.rainNb)

            id = 0
            if tmpTime[id,0] > self.tStart:
                self.rainMap[id] = None
                self.orographic[id] = False
                self.rainVal[id] = 0.
                self.rbgd[id] = 0.
                self.rmin[id] = 0.
                self.rmax[id] = 0.
                self.windx[id] = 0.
                self.windy[id] = 0.
                self.tauc[id] = 1000.
                self.tauf[id] = 1000.
                self.nm[id] = 0.005
                self.cw[id] = 0.005
                self.hw[id] = 3000.
                self.rainTime[id,0] = self.tStart
                self.rainTime[id,1] = tmpTime[0,0]
                self.ortime[id] = tmpTime[0,0] - self.tStart
                id += 1
            self.rainMap[id] = tmpMap[0]
            self.rainTime[id,:] = tmpTime[0,:]
            self.rainVal[id] = tmpVal[0]
            self.orographic[id] = tmpOro[0]
            self.ortime[id] = tmpoTime[0]
            self.rbgd[id] = tmprbgd[0]
            self.rmin[id] = tmprmin[0]
            self.rmax[id] = tmprmax[0]
            self.windx[id] = tmpwdx[0]
            self.windy[id] = tmpwdy[0]
            self.tauc[id] = tmptauc[0]
            self.tauf[id] = tmptauf[0]
            self.nm[id] = tmpnm[0]
            self.cw[id] = tmpcw[0]
            self.hw[id] = tmphw[0]

            id += 1
            for p in range(1,tmpNb):
                if tmpTime[p,0] > tmpTime[p-1,1]:
                    self.rainMap[id] = None
                    self.rainVal[id] = 0.
                    self.rainTime[id,0] = tmpTime[p-1,1]
                    self.rainTime[id,1] = tmpTime[p,0]
                    self.rbgd[id] = 0.
                    self.rmin[id] = 0.
                    self.rmax[id] = 0.
                    self.windx[id] = 0.
                    self.windy[id] = 0.
                    self.tauc[id] = 1000.
                    self.tauf[id] = 1000.
                    self.nm[id] = 0.005
                    self.cw[id] = 0.005
                    self.hw[id] = 3000.
                    self.ortime[id] = tmpTime[p,0] - tmpTime[p-1,1]
                    id += 1
                self.rainMap[id] = tmpMap[p]
                self.rainTime[id,:] = tmpTime[p,:]
                self.rainVal[id] = tmpVal[p]
                self.orographic[id] = tmpOro[p]
                self.ortime[id] = tmpoTime[p]
                self.rbgd[id] = tmprbgd[p]
                self.rmin[id] = tmprmin[p]
                self.rmax[id] = tmprmax[p]
                self.windx[id] = tmpwdx[p]
                self.windy[id] = tmpwdy[p]
                self.tauc[id] = tmptauc[p]
                self.tauf[id] = tmptauf[p]
                self.nm[id] = tmpnm[p]
                self.cw[id] = tmpcw[p]
                self.hw[id] = tmphw[p]
                id += 1
            if tmpTime[tmpNb-1,1] < self.tEnd:
                self.rainMap[id] = None
                self.rainVal[id] = 0.
                self.rainTime[id,0] = tmpTime[tmpNb-1,1]
                self.rainTime[id,1] = self.tEnd
                self.orographic[id] = False
                self.rainVal[id] = 0.
                self.rbgd[id] = 0.
                self.rmin[id] = 0.
                self.rmax[id] = 0.
                self.windx[id] = 0.
                self.windy[id] = 0.
                self.tauc[id] = 1000.
                self.tauf[id] = 1000.
                self.nm[id] = 0.005
                self.cw[id] = 0.005
                self.hw[id] = 3000.
                self.ortime[id] = self.tEnd - tmpTime[tmpNb-1,1]
        else:
            self.rainNb = 1
            self.rainVal = numpy.empty(self.rainNb)
            self.rainTime = numpy.empty((self.rainNb,2))
            self.rainMap = numpy.empty((self.rainNb),dtype=object)
            self.orographic = numpy.empty(self.rainNb,dtype=bool)
            self.ortime = numpy.empty(self.rainNb)
            self.rbgd = numpy.empty(self.rainNb)
            self.rmin = numpy.empty(self.rainNb)
            self.rmax = numpy.empty(self.rainNb)
            self.windx = numpy.empty(self.rainNb)
            self.windy = numpy.empty(self.rainNb)
            self.tauc = numpy.empty(self.rainNb)
            self.tauf = numpy.empty(self.rainNb)
            self.nm = numpy.empty(self.rainNb)
            self.cw = numpy.empty(self.rainNb)
            self.hw = numpy.empty(self.rainNb)
            self.rainVal[0] = 0.
            self.rainTime[0,0] = self.tStart
            self.rainTime[0,1] = self.tEnd
            self.rainMap[0] = None
            self.orographic[0] = False
            self.rainVal[0] = 0.
            self.rbgd[0] = 0.
            self.rmin[0] = 0.
            self.rmax[0] = 0.
            self.windx[0] = 0.
            self.windy[0] = 0.
            self.tauc[0] = 1000.
            self.tauf[0] = 1000.
            self.nm[0] = 0.005
            self.cw[0] = 0.005
            self.hw[0] = 3000.
            self.ortime[0] = self.tEnd - self.tStart

        # Extract Stream Power Law structure parameters
        spl = None
        spl = root.find('sp_law')
        if spl is not None:
            self.spl = True
            element = None
            element = spl.find('dep')
            if element is not None:
                self.depo = int(element.text)
            else:
                self.depo = 1
            element = None
            element = spl.find('slp_cr')
            if element is not None:
                self.slp_cr = float(element.text)
            else:
                self.slp_cr = 0.
            element = None
            element = spl.find('perc_dep')
            if element is not None:
                self.perc_dep = float(element.text)
            else:
                self.perc_dep = 0.
            element = None
            element = spl.find('fillmax')
            if element is not None:
                self.fillmax = float(element.text)
            else:
                self.fillmax = 1.
            element = None
            element = spl.find('m')
            if element is not None:
                self.SPLm = float(element.text)
            else:
                self.SPLm = 0.5
            element = None
            element = spl.find('n')
            if element is not None:
                self.SPLn = float(element.text)
            else:
                self.SPLn = 1.
            element = None
            element = spl.find('erodibility')
            if element is not None:
                self.SPLero = float(element.text)
            else:
                self.SPLero = 0.
            self.maxDT = None
            self.alluvial = 0.
            self.bedrock = 0.
            self.esmooth = 0.
            self.dsmooth = 0.
        else:
            self.depo = 0
            self.SPLm = 1.
            self.SPLn = 1.
            self.SPLero = 0.

        # Extract Transport Capacity model parameters
        tc = None
        tc = root.find('tc_law')
        if tc is not None:
            self.capacity = True
            if self.spl :
                raise ValueError('Only one of the sediment transport law can be defined.')
            self.depo = 1
            element = None
            element = tc.find('m')
            if element is not None:
                self.SPLm = float(element.text)
            else:
                self.SPLm = 1.
            element = None
            element = tc.find('n')
            if element is not None:
                self.SPLn = float(element.text)
            else:
                self.SPLn = 1.
            element = None
            element = tc.find('transport')
            if element is not None:
                self.SPLero = float(element.text)
            else:
                self.SPLero = 0.
            element = None
            element = tc.find('maxdt')
            if element is not None:
                self.maxDT = float(element.text)
            else:
                self.maxDT = None
            element = None
            element = tc.find('ascale')
            if element is not None:
                self.alluvial = float(element.text)
            else:
                self.alluvial = 0.
            element = None
            element = tc.find('bscale')
            if element is not None:
                self.bedrock = float(element.text)
            else:
                self.bedrock = 0.
        else:
            if not self.spl:
                self.depo = 0
                self.SPLm = 1.
                self.SPLn = 1.
                self.SPLero = 0.

        # Extract linear and nonlinear slope diffusion structure parameters
        creep = None
        creep = root.find('creep')
        if creep is not None:
            element = None
            element = creep.find('caerial')
            if element is not None:
                self.CDa = float(element.text)
            else:
                self.CDa = 0.
            element = None
            element = creep.find('cmarine')
            if element is not None:
                self.CDm = float(element.text)
            else:
                self.CDm = 0.
            element = None
            element = creep.find('cslp')
            if element is not None:
                self.Sc = float(element.text)
                if self.Sc >= 1.:
                    self.nHillslope = True
                else:
                    self.Sc = 0.
                    self.Hillslope = True
            else:
                self.Sc = 0.
                self.Hillslope = True
        else:
            self.CDa = 0.
            self.CDm = 0.
            self.Sc = 0.


        # Loading variable erodibility layers
        erostruct = None
        erostruct = root.find('erocoeff')
        if erostruct is not None:
            element = None
            element = erostruct.find('erolayers')
            if element is not None:
                self.erolays = int(element.text)
                if self.erolays == 0:
                    tmpNb = 1
                    self.eroMap = numpy.empty(1,dtype=object)
                    self.eroVal = numpy.empty(1,dtype=object)
                    self.thickMap = numpy.empty(1,dtype=object)
                    self.thickVal = numpy.empty(1,dtype=object)
                else:
                    tmpNb = self.erolays
                    self.eroMap = numpy.empty(self.erolays,dtype=object)
                    self.eroVal = numpy.empty(self.erolays,dtype=object)
                    self.thickMap = numpy.empty(self.erolays,dtype=object)
                    self.thickVal = numpy.empty(self.erolays,dtype=object)
            else:
                tmpNb = 1
                self.erolays = 0
                self.eroMap = numpy.empty(1,dtype=object)
                self.eroVal = numpy.empty(1,dtype=object)
                self.thickMap = numpy.empty(1,dtype=object)
                self.thickVal = numpy.empty(1,dtype=object)
            id = 0
            for elay in erostruct.iter('erolay'):
                if id >= tmpNb:
                    raise ValueError('The number of erodibility layers events does not match the number of defined layers.')
                element = None
                element = elay.find('erocst')
                if element is not None:
                    self.eroVal[id] = float(element.text)
                else:
                    self.eroVal[id] = None
                element = None
                element = elay.find('eromap')
                if element is not None:
                    self.eroMap[id] = element.text
                    if not os.path.isfile(self.eroMap[id]):
                        raise ValueError('Erodibility map file %s is missing or the given path is incorrect.'%(self.eroMap[id]))
                else:
                    self.eroMap[id] = None
                element = None
                element = elay.find('thcst')
                if element is not None:
                    self.thickVal[id] = float(element.text)
                else:
                    self.thickVal[id] = None
                element = None
                element = elay.find('thmap')
                if element is not None:
                    self.thickMap[id] = element.text
                    if not os.path.isfile(self.thickMap[id]):
                        raise ValueError('Thickness map file %s is missing or the given path is incorrect.'%(self.thickMap[id]))
                else:
                    self.thickMap[id] = None
                id += 1
        else:
            self.erolays = None
            self.eroMap = None
            self.eroVal = None
            self.thickMap = None
            self.thickVal = None

        # Flexural isostasy parameters
        flex = None
        flex = root.find('flexure')
        if flex is not None:
            self.flexure = True
            element = None
            element = flex.find('ftime')
            if element is not None:
                self.ftime = float(element.text)
            else:
                self.fnx = 0
            self.ftime = min(self.ftime,self.tDisplay)
            element = None
            element = flex.find('fnx')
            if element is not None:
                self.fnx = int(element.text)
            else:
                self.fnx = 0
            element = None
            element = flex.find('fny')
            if element is not None:
                self.fny = int(element.text)
            else:
                self.fny = 0
            element = None
            element = flex.find('dmantle')
            if element is not None:
                self.dmantle = float(element.text)
            else:
                self.dmantle = 3300.
            element = None
            element = flex.find('dsediment')
            if element is not None:
                self.dsediment = float(element.text)
            else:
                self.dsediment = 2500.
            element = None
            element = flex.find('youngMod')
            if element is not None:
                self.youngMod = float(element.text)
            else:
                self.youngMod = 65E9
            element = None
            element = flex.find('elasticH')
            if element is not None:
                self.elasticH = float(element.text)
            else:
                self.elasticH = None
            element = None
            element = flex.find('elasticGrid')
            if element is not None:
                self.elasticGrid = element.text
                if not os.path.isfile(self.elasticGrid):
                    raise ValueError('Elastic grid file is missing or the given path is incorrect.')
            else:
                self.elasticGrid = None
            element = None
            element = flex.find('boundary_W')
            if element is not None:
                self.flexbounds.append(element.text)
            else:
                raise ValueError('West boundary condition for flexure is not defined')
            element = None
            element = flex.find('boundary_E')
            if element is not None:
                self.flexbounds.append(element.text)
            else:
                raise ValueError('East boundary condition for flexure is not defined')
            element = None
            element = flex.find('boundary_S')
            if element is not None:
                self.flexbounds.append(element.text)
            else:
                raise ValueError('South boundary condition for flexure is not defined')
            element = None
            element = flex.find('boundary_N')
            if element is not None:
                self.flexbounds.append(element.text)
            else:
                raise ValueError('North boundary condition for flexure is not defined')

        # Extract Gaussian Filter structure parameters
        filter = None
        filter = root.find('filter')
        if filter is not None:
            self.filter = True
            if not self.capacity:
                element = None
                element = filter.find('gtime')
                if element is not None:
                    self.maxDT = float(element.text)
                else:
                    raise ValueError('Gaussian filter time step needs to be defined.')
            elif self.maxDT is None:
                 raise ValueError('tc_law maxdt element needs to be defined.')
            element = None
            element = filter.find('esmooth')
            if element is not None:
                self.esmooth = float(element.text)
            else:
                self.esmooth = 0.
            element = None
            element = filter.find('dsmooth')
            if element is not None:
                self.dsmooth = float(element.text)
            else:
                self.dsmooth = 0.
        else:
            self.esmooth = None
            self.dsmooth = None

        # Get output directory
        out = None
        out = root.find('outfolder')
        if out is not None:
            self.outDir = out.text
        else:
            self.outDir = os.getcwd()+'/out'

        if self.makeUniqueOutputDir:
            if os.path.exists(self.outDir):
                self.outDir += '_'+str(len(glob.glob(self.outDir+str('*')))-1)

            os.makedirs(self.outDir)
            os.makedirs(self.outDir+'/h5')
            os.makedirs(self.outDir+'/xmf')
            shutil.copy(self.inputfile,self.outDir)

        return
