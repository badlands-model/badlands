##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##

import os
import glob
import numpy
import shutil
from decimal import *
import xml.etree.ElementTree as ET

class xmlParser:
    """
    This class parses the XmL input file required to run **badlands**.

    Args:
        inputfile : (str) this is a string containing the XML input file.
        makeUniqueOutputDir : (boolean) uniquely-named directory for the output (default: True)
    """
    def __init__(self, inputfile = None, makeUniqueOutputDir=True):

        if inputfile==None:
            raise RuntimeError('XmL input file name must be defined to run a Badlands simulation.')
        if not os.path.isfile(inputfile):
            raise RuntimeError('The XmL input file name cannot be found in your path.')
        self.inputfile = inputfile

        self.demfile = None
        self.btype = 'slope'
        self.perc_dep = 0.5
        self.slp_cr = 0.
        self.fillmax = 200.
        self.Afactor = 1
        self.nopit = 0
        self.udw = 0
        self.searef = None
        self.poro0 = 0.
        self.poroC = 0.47

        self.restart = False
        self.rfolder = None
        self.rstep = 0
        self.tStart = None
        self.tEnd = None
        self.tDisplay = None
        self.tmesh = 1
        self.mesh = None
        self.minDT = 1.
        self.maxDT = 1.e6

        self.stratdx = 0.
        self.laytime = 0.

        self.seapos = 0.
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
        self.riverRck = None

        self.rainNb = None
        self.rainVal = None
        self.rainTime = None
        self.oroRain = False
        self.orographic = None
        self.orographiclin = None
        self.rzmax = None
        self.ortime = None
        self.rbgd = None
        self.rmin = None
        self.rmax = None
        self.rzmax = None
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
        self.diffnb = 5
        self.diffprop = 0.9
        self.propa = None
        self.propb = None
        self.slpDiffprop = False
        self.spl = False
        self.deepbasin = -10000.
        self.denscrit = 20000.

        self.incisiontype = 0
        self.mp = 0.
        self.mt = 0.
        self.nt = 0.
        self.kt = 0.
        self.kw = 0.
        self.b = 0.
        self.bedslptype = 0

        self.Hillslope = False
        self.CDa = 0.
        self.CDm = 0.
        self.Sc = 0.
        self.Sfail = 0.
        self.Cfail = 0.
        self.CDr = 0.
        self.makeUniqueOutputDir = makeUniqueOutputDir

        self.outDir = None
        self.sh5file = 'h5/sed'
        self.strath5file = 'h5/stratal.time'

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
        self.elasticA1 = None
        self.elasticA2 = None
        self.flexbounds = []

        self.erolays = None
        self.eroMap = None
        self.eroVal = None
        self.thickMap = None
        self.thickVal = None

        self.rockNb = 0
        self.rockCk = None
        self.actlay = 50.
        self.initlayers = 0
        self.layersData = None
        self.laytime = 0.

        self.waveOn = False
        self.waveSed = False
        self.tWave = None
        self.resW = None
        self.waveBase = 20.
        self.d50 = 0.0001
        self.tsteps = 1000
        self.dsteps = 1000
        self.waveNb = 0
        self.climNb = 0
        self.waveTime = None
        self.wavePerc = None
        self.waveWu = None
        self.waveWh = None
        self.waveWp = None
        self.waveWd = None
        self.waveWdd = None
        self.waveWs = None
        self.wCd = None
        self.wCe = None
        self.wEro = None
        self.waveSide = None
        self.wavelist = None
        self.climlist = None

        self.swanFile = None
        self.swanInfo = None
        self.swanBot = None
        self.swanOut = None

        self.carb = False
        self.carbonate = False
        self.carbDepth = None
        self.carbSed = None
        self.carbWave = None
        self.islandPerim = 0.
        self.coastdist = 0.
        self.baseMap = None
        self.tCarb = None
        self.carbNb = None
        self.carbValSp1 = None
        self.carbValSp2 = None
        self.carbTime = None

        self.carbonate2 = False
        self.carbDepth2 = None
        self.carbSed2 = None
        self.carbWave2 = None
        self.islandPerim2 = 0.
        self.coastdist2 = 0.

        self.pelagic = False
        self.pelGrowth = 0.
        self.pelDepth = None

        self._get_XmL_Data()

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
                if self.btype != 'slope' and self.btype != 'flat' and self.btype != 'wall' \
                   and self.btype != 'fixed' and self.btype != 'outlet' and self.btype != 'wall1':
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
            element = None
            element = grid.find('udw')
            if element is not None:
                self.udw = int(element.text)
                if self.udw < 1:
                    self.udw = 0
                else:
                    self.udw = 1
            else:
                self.udw = 0
            element = None
            element = grid.find('searef')
            if element is not None:
                self.searef = eval(element.text)
                if not isinstance(self.searef, tuple):
                    raise ValueError("""searef must be a python tuple (x,y)""")
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
            element = time.find('meshout')
            if element is not None:
                self.tmesh = int(element.text)
            else:
                self.tmesh = 1
            element = None
            element = time.find('mindt')
            if element is not None:
                self.minDT = float(element.text)
            else:
                self.minDT = 1.
            element = None
            element = time.find('maxdt')
            if element is not None:
                self.maxDT = float(element.text)
            else:
                self.maxDT = self.tDisplay
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
            element = strat.find('poroC')
            if element is not None:
                self.poroC = float(element.text)
            element = None
            element = strat.find('poro0')
            if element is not None:
                self.poro0 = float(element.text)
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
            element = sea.find('curve')
            if element is not None:
                self.seafile = element.text
                if not os.path.isfile(self.seafile):
                    raise ValueError('Sea level file is missing or the given path is incorrect.')
            else:
                self.seafile = None
        else:
            self.seapos = 0.
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
                if self.time3d < self.minDT:
                    raise ValueError('The value of time3d cannot be lower than mindt.')
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
            self.riverRck = numpy.empty((self.riverNb),int)
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
                element = riv.find('rQw')
                if element is not None:
                    self.riverQws[id,0] = float(element.text)
                else:
                    raise ValueError('River %d is missing mean annual discharge.'%id)
                element = None
                element = riv.find('rQs')
                if element is not None:
                    # Convert from Mt/a to kg/a
                    qs = float(element.text)*1.e9
                else:
                    raise ValueError('River %d is missing river annual load.'%id)
                element = None
                element = riv.find('rhoS')
                if element is not None:
                    rhoS = float(element.text)
                else:
                    rhoS = 2650.
                # Convert from kg/a to m3/a
                self.riverQws[id,1] = qs/rhoS
                element = None
                element = riv.find('rck')
                if element is not None:
                    self.riverRck[id] = int(element.text)
                else:
                    self.riverRck[id] = 0

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
            tmpOroLin = numpy.empty(tmpNb,dtype=bool)
            tmpTime = numpy.empty((tmpNb,2))
            tmpoTime = numpy.empty(tmpNb)
            tmprbgd = numpy.empty(tmpNb)
            tmprmin = numpy.empty(tmpNb)
            tmprmax = numpy.empty(tmpNb)
            tmprzmax = numpy.empty(tmpNb)
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
                element = clim.find('rzmax')
                if element is not None:
                    tmprzmax[id] = float(element.text)
                    tmpOroLin[id] = True
                    self.oroRain = True
                else:
                    tmprzmax[id] = 0.
                    tmpOroLin[id] = False
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
            self.orographiclin = numpy.empty(self.rainNb,dtype=bool)
            self.ortime = numpy.empty(self.rainNb)
            self.rbgd = numpy.empty(self.rainNb)
            self.rmin = numpy.empty(self.rainNb)
            self.rmax = numpy.empty(self.rainNb)
            self.rzmax = numpy.empty(self.rainNb)
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
                self.orographiclin[id] = False
                self.rainVal[id] = 0.
                self.rbgd[id] = 0.
                self.rmin[id] = 0.
                self.rmax[id] = 0.
                self.rzmax[id] = 0.
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
            self.orographiclin[id] = tmpOroLin[0]
            self.ortime[id] = tmpoTime[0]
            self.rbgd[id] = tmprbgd[0]
            self.rmin[id] = tmprmin[0]
            self.rmax[id] = tmprmax[0]
            self.rzmax[id] = tmprzmax[0]
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
                    self.rzmax[id] = 0.
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
                self.orographiclin[id] = tmpOroLin[p]
                self.ortime[id] = tmpoTime[p]
                self.rbgd[id] = tmprbgd[p]
                self.rmin[id] = tmprmin[p]
                self.rmax[id] = tmprmax[p]
                self.rzmax[id] = tmprzmax[p]
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
                self.orographiclin[id] = False
                self.rainVal[id] = 0.
                self.rbgd[id] = 0.
                self.rmin[id] = 0.
                self.rmax[id] = 0.
                self.rzmax[id] = 0.
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
            self.orographiclin = numpy.empty(self.rainNb,dtype=bool)
            self.ortime = numpy.empty(self.rainNb)
            self.rbgd = numpy.empty(self.rainNb)
            self.rmin = numpy.empty(self.rainNb)
            self.rmax = numpy.empty(self.rainNb)
            self.rzmax = numpy.empty(self.rainNb)
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
            self.orographiclin[0] = False
            self.rbgd[0] = 0.
            self.rmin[0] = 0.
            self.rmax[0] = 0.
            self.rzmax[0] = 0.
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
                if self.slp_cr == 0.:
                    self.perc_dep = 0.
                else:
                    self.perc_dep = 0.5
            element = None
            element = spl.find('fillmax')
            if element is not None:
                self.fillmax = float(element.text)
            else:
                self.fillmax = 200.
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
            element = None
            element = spl.find('diffnb')
            if element is not None:
                self.diffnb = int(element.text)
            else:
                self.diffnb = 5
            element = None
            element = spl.find('dens_cr')
            if element is not None:
                self.denscrit = float(element.text)
            element = None
            element = spl.find('deepbasin')
            if element is not None:
                self.deepbasin = float(element.text)
            element = None
            element = spl.find('diffprop')
            if element is not None:
                self.diffprop = float(element.text)
                if self.diffprop <= 0 or self.diffprop >= 1:
                    raise ValueError('Proportion of marine sediment deposited on downstream nodes needs to be range between ]0,1[')
            else:
                self.diffprop = 0.9
            element = None
            element = spl.find('propa')
            if element is not None:
                self.propa = float(element.text)
                self.slpDiffprop = True
            else:
                self.propa = None
            element = None
            element = spl.find('propb')
            if element is not None:
                self.propb = float(element.text)
            else:
                self.propb = None
            dpropList = [self.propa,self.propb]
            if dpropList.count(None) == 1:
                raise ValueError('Both <propa> and <propb> must be assigned to implement slope-dependent diffprop.')
            if dpropList.count(None) == 2: # Assign them to real vals for input into pdalgo
                self.propa = 0.
                self.propb = 0.
        else:
            self.depo = 0
            self.SPLm = 1.
            self.SPLn = 1.
            self.SPLero = 0.
            self.diffnb = 5

        # Extract flux-dependent function structure parameters
        erof = None
        erof = root.find('sedfluxfunction')
        if erof is not None:
            element = None
            element = erof.find('modeltype')
            if element is not None:
                self.incisiontype = int(element.text)
            else:
                self.incisiontype = 0
            element = None
            element = erof.find('mt')
            if element is not None:
                self.mt = element.text
            else:
                self.mt = 0
            element = None
            element = erof.find('nt')
            if element is not None:
                self.nt = element.text
            else:
                self.nt = 0
            element = None
            element = erof.find('kt')
            if element is not None:
                self.kt = element.text
            else:
                self.kt = 0
            element = None
            element = erof.find('kw')
            if element is not None:
                self.kw = element.text
            else:
                self.kw = 0
            element = None
            element = erof.find('b')
            if element is not None:
                self.b = element.text
            else:
                self.b = 0
            element = None
            element = erof.find('mp')
            if element is not None:
                self.mp = float(element.text)
            else:
                self.mp = 0
            element = None
            element = erof.find('bedslp')
            if element is not None:
                self.bedslptype = int(element.text)
            else:
                self.bedslptype = 0

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
                if self.Sc < 0.:
                    self.Sc = 0.
            else:
                self.Sc = 0.
            element = None
            element = creep.find('sfail')
            if element is not None:
                self.Sfail = float(element.text)
                if self.Sfail < 0.:
                    self.Sfail = 0.
            else:
                self.Sfail = 0.
            element = None
            element = creep.find('cfail')
            if element is not None:
                self.Cfail = float(element.text)
                if self.Cfail < 0.:
                    self.Cfail = 0.
            else:
                self.Cfail = 0.
            element = None
            element = creep.find('criver')
            if element is not None:
                self.CDr = float(element.text)
            else:
                self.CDr = 0.
            self.Hillslope = True
        else:
            self.CDa = 0.
            self.CDm = 0.
            self.CDr = 0.
            self.Sc = 0.

        # Loading erodibility layers
        erost = None
        erost = root.find('erocoeffs')
        if erost is not None:
            element = None
            element = erost.find('actlay')
            if element is not None:
                self.actlay = float(element.text)
            else:
                self.actlay = 50.
            id = 0
            element = None
            element = erost.find('rocktype')
            if element is not None:
                self.rockNb = int(element.text)
            else:
                self.rockNb = 0
            element = None
            element = erost.find('laytime')
            if element is not None:
                self.laytime = float(element.text)
            else:
                self.laytime = self.tDisplay
            if self.laytime >  self.tDisplay:
                 self.laytime = self.tDisplay
            if self.tDisplay % self.laytime != 0:
                raise ValueError('Error in the XmL file: stratal layer interval needs to be an exact multiple of the display interval!')
            id = 0
            self.rockCk = numpy.zeros(self.rockNb)
            for rcktyp in erost.iter('rockero'):
                if id >= self.rockNb:
                    raise ValueError('The number of rock erodibilities values does not match the number of defined rock types.')
                element = None
                element = rcktyp.find('erorock')
                if element is not None:
                    self.rockCk[id] = float(element.text)
                else:
                    self.rockCk[id] = None
                id += 1
            element = None
            element = erost.find('erolayers')
            if element is not None:
                self.initlayers = int(element.text)
            else:
                self.initlayers = 0
            id = 0
            self.layersData = numpy.empty(self.initlayers,dtype=object)
            for rcklay in erost.iter('erolay'):
                if id >= self.initlayers:
                    raise ValueError('The number of erodibility layer maps does not match the number of defined layers.')
                element = None
                element = rcklay.find('laymap')
                if element is not None:
                    self.layersData[id] = element.text
                else:
                    self.layersData[id] = None
                id += 1
        else:
            self.rockNb = 0
            self.actlay = 50.
            self.rockCk = None
            self.initlayers = 0
            self.layersData = None

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
                    self.eroVal = numpy.zeros(1,dtype=object)
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
            element = flex.find('elasticA1')
            if element is not None:
                self.elasticA1 = float(element.text)
            else:
                self.elasticA1 = None
            element = None
            element = flex.find('elasticA2')
            if element is not None:
                self.elasticA2 = float(element.text)
            else:
                self.elasticA2 = None
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

        # Extract global wave field parameters
        wavefield = None
        wavefield = root.find('waveglobal')
        if wavefield is not None:
            element = None
            element = wavefield.find('wmodel')
            if element is not None:
                model = int(element.text)
                if model == 1:
                    self.waveOn = True
                    self.waveSed = False
                else:
                    self.waveOn = False
                    self.waveSed = True
            else:
                self.waveOn = False
                self.waveSed = True
            element = None
            element = wavefield.find('twave')
            if element is not None:
                self.tWave = float(element.text)
                if self.waveOn and Decimal(self.tEnd - self.tStart) % Decimal(self.tWave) != 0.:
                    raise ValueError('Error in the definition of the simulation time: wave interval needs to be a multiple of simulation time.')
            else:
                raise ValueError('Error in the definition of the simulation time: wave interval is required')
            if self.tWave > self.maxDT:
                self.tWave = self.maxDT
            element = None
            element = wavefield.find('wres')
            if element is not None:
                self.resW = float(element.text)
            else:
                raise ValueError('Error the wave grid resolution needs to be defined')
            element = None
            element = wavefield.find('wCd')
            if element is not None:
                self.wCd = float(element.text)
            else:
                self.wCd = 50.
            element = None
            element = wavefield.find('wCe')
            if element is not None:
                self.wCe = float(element.text)
                if self.wCe > 1:
                    self.wCe = 1.
                if self.wCe < 0:
                    self.wCe = 0.1
            else:
                self.wCe = 0.5
            element = None
            element = wavefield.find('wEro')
            if element is not None:
                self.wEro = float(element.text)
                if self.wEro < 0:
                    self.wEro = -self.wEro
            else:
                self.wEro = 0.5
            element = None
            element = wavefield.find('wbase')
            if element is not None:
                self.waveBase = float(element.text)
            else:
                self.waveBase = 20.
            if self.waveOn:
                self.waveBase = 100000.
            element = None
            element = wavefield.find('events')
            if element is not None:
                self.waveNb = int(element.text)
            else:
                raise ValueError('The number of wave temporal events needs to be defined.')
            element = None
            element = wavefield.find('d50')
            if element is not None:
                self.d50 = float(element.text)
            else:
                self.d50 = 0.0001
            element = None
            element = wavefield.find('tsteps')
            if element is not None:
                self.tsteps = int(element.text)
            else:
                self.tsteps = 1000
            element = None
            element = wavefield.find('dsteps')
            if element is not None:
                self.dsteps = int(element.text)
            else:
                self.dsteps = 1000
        else:
            self.waveNb = 0
            self.tWave = self.tEnd - self.tStart + 10000.

        # Extract wave field structure information
        if self.waveNb > 0:
            tmpNb = self.waveNb
            if self.waveOn:
                self.waveWu = []
                self.waveWp = []
                self.waveWdd = []
                self.waveWs = []
                self.waveSide = []
            self.waveWh = []
            self.wavePerc = []
            self.waveWd = []

            self.waveTime = numpy.empty((tmpNb,2))
            self.climNb = numpy.empty(tmpNb, dtype=int)
            w = 0

            for wavedata in root.iter('wave'):
                if w >= tmpNb:
                    raise ValueError('Wave event number above number defined in global wave structure.')
                if wavedata is not None:
                    element = None
                    element = wavedata.find('start')
                    if element is not None:
                        self.waveTime[w,0] = float(element.text)
                        # if w > 0 and self.waveTime[w,0] != self.waveTime[w-1,1]:
                        #     raise ValueError('The start time of the wave field %d needs to match the end time of previous wave data.'%w)
                        # if w == 0 and self.waveTime[w,0] != self.tStart:
                        #     raise ValueError('The start time of the first wave field needs to match the simulation start time.')
                    else:
                        raise ValueError('Wave event %d is missing start time argument.'%w)
                    element = None
                    element = wavedata.find('end')
                    if element is not None:
                        self.waveTime[w,1] = float(element.text)
                    else:
                        raise ValueError('Wave event %d is missing end time argument.'%w)
                    if self.waveTime[w,0] >= self.waveTime[w,1]:
                        raise ValueError('Wave event %d start and end time values are not properly defined.'%w)
                    element = None
                    element = wavedata.find('climNb')
                    if element is not None:
                        self.climNb[w] = int(element.text)
                    else:
                        raise ValueError('Wave event %d is missing climatic wave number argument.'%w)

                    if self.waveOn and Decimal(self.waveTime[w,1]-self.waveTime[w,0]) % Decimal(self.tWave) != 0.:
                        raise ValueError('Wave event %d duration need to be a multiple of the wave interval.'%w)

                    if self.waveSed and Decimal(self.waveTime[w,1]-self.waveTime[w,0]) % Decimal(self.maxDT) != 0.:
                        raise ValueError('Wave event %d duration need to be a multiple of the wave interval.'%w)

                    if self.waveOn:
                        listWu = []
                        listWp = []
                        listWs = []
                        listWdd = []
                        listWside = []
                        listBreak = []
                    listPerc = []
                    listWd = []
                    listWh = []

                    id = 0
                    sumPerc = 0.
                    for clim in wavedata.iter('climate'):
                        if id >= self.climNb[w]:
                            raise ValueError('The number of climatic events does not match the number of defined climates.')
                        element = None
                        element = clim.find('perc')
                        if element is not None:
                            sumPerc += float(element.text)
                            if sumPerc > 1:
                                raise ValueError('Sum of wave event %d percentage is higher than 1.'%w)
                            listPerc.append(float(element.text))
                            if listPerc[id] < 0:
                                raise ValueError('Wave event %d percentage cannot be negative.'%w)
                        else:
                            raise ValueError('Wave event %d is missing percentage argument.'%w)
                        element = None
                        element = clim.find('hs')
                        if element is not None:
                            listWh.append(float(element.text))
                            if listWh[id] < 0:
                                raise ValueError('Wave event %d significant wave height cannot be negative.'%w)
                        else:
                            listWh.append(0.001)
                        element = None
                        element = clim.find('dir')
                        if element is not None:
                            listWd.append(float(element.text))
                            if listWd[id] < 0:
                                raise ValueError('Wave event %d wave direction needs to be set between 0 and 360.'%w)
                            if listWd[id] > 360:
                                raise ValueError('Wave event %d wave direction needs to be set between 0 and 360.'%w)
                        else:
                            listWd.append(0.)

                        if self.waveOn:
                            element = None
                            element = clim.find('per')
                            if element is not None:
                                listWp.append(float(element.text))
                                if listWp[id] <= 0:
                                    raise ValueError('Wave event %d wave period cannot be negative.'%w)
                            else:
                                listWp.append(30.)
                            element = None
                            element = clim.find('windv')
                            if element is not None:
                                listWu.append(float(element.text))
                                if listWu[id] < 0:
                                    raise ValueError('Wave event %d wind velocity cannot be negative.'%w)
                            else:
                                listWu.append(0.)
                            element = None
                            element = clim.find('wdir')
                            if element is not None:
                                listWdd.append(float(element.text))
                                if listWdd[id] < 0:
                                    raise ValueError('Wave event %d wind direction needs to be set between 0 and 360.'%w)
                                if listWdd[id] > 360:
                                    raise ValueError('Wave event %d wind direction needs to be set between 0 and 360.'%w)
                            else:
                                listWdd.append(0.)
                            element = None
                            element = clim.find('spread')
                            if element is not None:
                                listWs.append(float(element.text))
                                if listWs[id] < 0:
                                    raise ValueError('Wave event %d spreading angle needs to be set between 0 and 360.'%w)
                                if listWs[id] > 360:
                                    raise ValueError('Wave event %d spreading angle needs to be set between 0 and 360.'%w)
                            else:
                                listWs.append(0.)
                            element = None
                            element = clim.find('side')
                            if element is not None:
                                listWside.append(int(element.text))
                                if listWside[id] > 8:
                                    raise ValueError('Wave boundary side is between 1 and 8.')
                            else:
                                listWside[id]=1
                        id += 1

                    w += 1
                    self.wavePerc.append(listPerc)
                    self.waveWd.append(listWd)
                    self.waveWh.append(listWh)
                    if self.waveOn:
                        self.waveWu.append(listWu)
                        self.waveWdd.append(listWdd)
                        self.waveWs.append(listWs)
                        self.waveWp.append(listWp)
                        self.waveSide.append(listWside)
                else:
                    raise ValueError('Wave event %d is missing.'%w)

        # Construct a list of climatic events for swan/wavesed model
        if self.waveOn:
            self.wavelist = []
            self.climlist = []
            twsteps = numpy.arange(self.tStart,self.tEnd,self.tWave)
            for t in range(len(twsteps)):
                c = -1
                # Find the wave field active during the time interval
                for k in range(self.waveNb):
                    if self.waveTime[k,0] <= twsteps[t] and self.waveTime[k,1] >= twsteps[t]:
                        c = k
                # Extract the wave climate for the considered time interval
                for p in range(self.climNb[c]):
                    self.wavelist.append(c)
                    self.climlist.append(p)

            # Add a fake final wave field and climate
            self.wavelist.append(self.wavelist[-1])
            self.climlist.append(self.climlist[-1])

            # Create swan model repository and files
            os.makedirs(self.outDir+'/swan')
            self.swanFile = numpy.array(self.outDir+'/swan/swan.swn')
            self.swanInfo = numpy.array(self.outDir+'/swan/swanInfo.swn')
            self.swanBot = numpy.array(self.outDir+'/swan/swan.bot')
            self.swanOut = numpy.array(self.outDir+'/swan/swan.csv')

        # Carbonate parameters
        carbp = None
        carbp = root.find('carb')
        if carbp is not None:
            self.carb = True
            element = None
            element = carbp.find('baseMap')
            if element is not None:
                self.baseMap = element.text
                if not os.path.isfile(self.baseMap):
                    raise ValueError('Basement map file for species1 growth is missing or the given path is incorrect.')
            element = None
            element = carbp.find('tcarb')
            if element is not None:
                self.tCarb = float(element.text)
                if Decimal(self.tEnd - self.tStart) % Decimal(self.tCarb) != 0.:
                    raise ValueError('Error in the definition of the simulation time: carbonate interval needs to be a multiple of simulation time.')
            else:
                raise ValueError('Error in the definition of the simulation time: carbonate interval is required')
            element = None
            element = carbp.find('growth_events')
            if element is not None:
                tmpNb = int(element.text)
            else:
                raise ValueError('The number of carbonate growth events needs to be defined.')
            if self.tCarb > self.maxDT:
                self.tCarb = self.maxDT
            if self.tCarb > self.tWave:
                self.tCarb = self.tWave
            else:
                self.tWave = self.tCarb
                self.maxDT = self.tCarb

            tmpValSp1 = numpy.empty(tmpNb)
            tmpValSp2 = numpy.empty(tmpNb)
            tmpTime = numpy.empty((tmpNb,2))

            id = 0

            for clim in carbp.iter('event'):
                if id >= tmpNb:
                    raise ValueError('The number of reef growth events does not match the number of defined events.')
                element = None
                element = clim.find('growth_sp1')
                if element is not None:
                    tmpValSp1[id] = float(element.text)
                else:
                    tmpValSp1[id] = 0.
                element = None
                element = clim.find('growth_sp2')
                if element is not None:
                    tmpValSp2[id] = float(element.text)
                else:
                    tmpValSp2[id] = 0.
                element = None
                element = clim.find('gstart')
                if element is not None:
                    tmpTime[id,0] = float(element.text)
                else:
                    raise ValueError('Reef growth event %d is missing start time argument.'% int(id+1))
                element = None
                element = clim.find('gend')
                if element is not None:
                    tmpTime[id,1] = float(element.text)
                else:
                    raise ValueError('Reef growth event %d is missing end time argument.'% int(id+1))
                if tmpTime[id,0] >= tmpTime[id,1]:
                    raise ValueError('Reef growth event %d start and end time values are not properly defined.' % int(id+1))
                if id > 0:
                    if tmpTime[id,0] < tmpTime[id-1,1]:
                        raise ValueError('Reef growth event %d start time needs to be >= than rain climate %d end time.'%(int(id+1),id))

                id += 1

            if id != tmpNb:
                raise ValueError('Reef growth event parameter %d does not match with the number of declared reef growth events %d.' %(tmpNb,id))

            # Create continuous reef growth series
            self.carbNb = tmpNb
            if tmpTime[0,0] > self.tStart:
                self.carbNb += 1
            for id in range(1,tmpNb):
                if tmpTime[id,0] > tmpTime[id-1,1]:
                    self.carbNb += 1
            if tmpTime[tmpNb-1,1] < self.tEnd:
                self.carbNb += 1
            self.carbValSp1 = numpy.empty(self.carbNb)
            self.carbValSp2 = numpy.empty(self.carbNb)
            self.carbTime = numpy.empty((self.carbNb,2))

            id = 0

            if tmpTime[id,0] > self.tStart:
                self.carbValSp1[id] = 0.
                self.carbValSp2[id] = 0.
                self.carbTime[id,0] = self.tStart
                self.carbTime[id,1] = tmpTime[0,0]
                id += 1
            self.carbTime[id,:] = tmpTime[0,:]
            self.carbValSp1[id] = tmpValSp1[0]
            self.carbValSp2[id] = tmpValSp2[0]

            id += 1
            for p in range(1,tmpNb):
                if tmpTime[p,0] > tmpTime[p-1,1]:
                    self.carbValSp1[id] = 0.
                    self.carbValSp2[id] = 0.
                    self.carbTime[id,0] = tmpTime[p-1,1]
                    self.carbTime[id,1] = tmpTime[p,0]
                    id += 1
                self.carbTime[id,:] = tmpTime[p,:]
                self.carbValSp1[id] = tmpValSp1[p]
                self.carbValSp2[id] = tmpValSp2[p]
                id += 1
            if tmpTime[tmpNb-1,1] < self.tEnd:
                self.carbValSp1[id] = 0.
                self.carbValSp2[id] = 0.
                self.carbTime[id,0] = tmpTime[tmpNb-1,1]
                self.carbTime[id,1] = self.tEnd
        else:
            self.carbNb = 1
            self.carbValSp1 = numpy.empty(self.carbNb)
            self.carbValSp2 = numpy.empty(self.carbNb)
            self.carbTime = numpy.empty((self.carbNb,2))
            self.carbValSp1[0] = 0.
            self.carbValSp2[0] = 0.
            self.carbTime[0,0] = self.tStart
            self.carbTime[0,1] = self.tEnd


        # Species 1 class
        carb = None
        carb = root.find('species1')
        if carb is not None:
            self.carbonate = True
            element = None
            element = carb.find('depthControl')
            if element is not None:
                self.carbDepth = element.text
                if not os.path.isfile(self.carbDepth):
                    raise ValueError('Species1 depth control file is missing or the given path is incorrect.')
            else:
                self.carbDepth = None
            element = None
            element = carb.find('waveControl')
            if element is not None:
                self.carbWave = element.text
                if not os.path.isfile(self.carbWave):
                    raise ValueError('Species1 wave control file is missing or the given path is incorrect.')
            else:
                self.carbWave = None
            element = None
            element = carb.find('sedControl')
            if element is not None:
                self.carbSed = element.text
                if not os.path.isfile(self.carbSed):
                    raise ValueError('Species1 sedimentation control file is missing or the given path is incorrect.')
            else:
                self.carbSed = None
            element = None
            element = carb.find('isld')
            if element is not None:
                self.islandPerim = float(element.text)
            else:
                self.islandPerim = 0.
            element = None
            element = carb.find('dist')
            if element is not None:
                self.coastdist = float(element.text)
            else:
                self.coastdist = 0.

        # Species 2 class
        carb2 = None
        carb2 = root.find('species2')
        if carb2 is not None:
            self.carbonate2 = True
            element = None
            element = carb2.find('depthControl')
            if element is not None:
                self.carbDepth2 = element.text
                if not os.path.isfile(self.carbDepth2):
                    raise ValueError('Species2 depth control file is missing or the given path is incorrect.')
            else:
                self.carbDepth2 = None
            element = None
            element = carb2.find('waveControl')
            if element is not None:
                self.carbWave2 = element.text
                if not os.path.isfile(self.carbWave2):
                    raise ValueError('Species2 wave control file is missing or the given path is incorrect.')
            else:
                self.carbWave2 = None
            element = None
            element = carb2.find('sedControl')
            if element is not None:
                self.carbSed2 = element.text
                if not os.path.isfile(self.carbSed2):
                    raise ValueError('Species2 sedimentation control file is missing or the given path is incorrect.')
            else:
                self.carbSed2 = None
            element = None
            element = carb2.find('isld')
            if element is not None:
                self.islandPerim2 = float(element.text)
            else:
                self.islandPerim2 = 0.
            element = carb2.find('dist')
            element = None
            if element is not None:
                self.coastdist2 = float(element.text)
            else:
                self.coastdist2 = 0.

        # Pelagic class
        pelagic = None
        pelagic = root.find('pelagic')
        if pelagic is not None:
            self.pelagic = True
            element = None
            element = pelagic.find('growth')
            if element is not None:
                self.pelGrowth = float(element.text)
            else:
                self.pelGrowth = 0.
            element = None
            element = pelagic.find('depthControl')
            if element is not None:
                self.pelDepth = element.text
                if not os.path.isfile(self.pelDepth):
                    raise ValueError('Pelagic depth control file is missing or the given path is incorrect.')
            else:
                raise ValueError('Depth control file on pelagic deposition is required.')

        return
