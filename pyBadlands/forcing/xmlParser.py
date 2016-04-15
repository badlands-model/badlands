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
        self.fillmax = 1.
        self.Afactor = 1

        self.tStart = None
        self.tEnd = None
        self.tDisplay = None
        self.minDT = 1.

        self.seapos = 0.
        self.sealimit = 100.
        self.seafile = None

        self.disp3d = False
        self.tectNb = None
        self.tectTime = None
        self.tectFile = None
        self.merge3d = None
        self.time3d = None

        self.rainNb = None
        self.rainVal = None
        self.rainTime = None
        self.rainMap = None

        self.depo = 1
        self.SPLm = 0.5
        self.SPLn = 1.
        self.SPLero = 0.
        self.maxDT = None
        self.alluvial = 0.
        self.bedrock = 0.

        self.CDa = 0.
        self.CDm = 0.
        self.makeUniqueOutputDir = makeUniqueOutputDir

        self.outDir = None
        self.th5file = 'h5/tin.time'
        self.txmffile = 'xmf/tin.time'
        self.txdmffile = 'tin.series.xdmf'

        self.fh5file = 'h5/flow.time'
        self.fxmffile = 'xmf/flow.time'
        self.fxdmffile = 'flow.series.xdmf'

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
            element = None
            element = grid.find('demfile')
            if element is not None:
                self.demfile = element.text
                if not os.path.isfile(self.demfile):
                    raise ValueError('DEM file is missing or the given path is incorrect.')
            else:
                raise ValueError('Error in the definition of the grid structure: DEM file definition is required!')
            element = None
            element = grid.find('boundary')
            if element is not None:
                self.btype = element.text
                if self.btype != 'slope' and self.btype != 'flat' and self.btype != 'wall':
                    raise ValueError('Error in the definition of the grid structure: Boundary type is either: flat, slope or wall')
            else:
                self.btype = 'slope'
            element = None
            element = grid.find('fillmax')
            if element is not None:
                self.fillmax = float(element.text)
            else:
                self.fillmax = 1.
            element = None
            element = grid.find('resfactor')
            if element is not None:
                self.Afactor = int(element.text)
                if self.Afactor < 1:
                    self.Afactor = 1
            else:
                self.Afactor = 1
        else:
            raise ValueError('Error in the XmL file: grid structure definition is required!')

        # Extract time structure information
        time = None
        time = root.find('time')
        if time is not None:
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
                raise ValueErself.ror('Error in the definition of the simulation time: end time declaration is required')
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
        else:
            raise ValueError('Error in the XmL file: time structure definition is required!')

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
            tmpTime = numpy.empty((tmpNb,2))
            id = 0
            for clim in precip.iter('rain'):
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
                    tmpVal[id] = element.text
                else:
                    tmpVal[id] = 0.
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
            id = 0
            if tmpTime[id,0] > self.tStart:
                self.rainMap[id] = None
                self.rainVal[id] = 0.
                self.rainTime[id,0] = self.tStart
                self.rainTime[id,1] = tmpTime[0,0]
                id += 1
            self.rainMap[id] = tmpMap[0]
            self.rainTime[id,:] = tmpTime[0,:]
            self.rainVal[id] = tmpVal[0]
            id += 1
            for p in range(1,tmpNb):
                if tmpTime[p,0] > tmpTime[p-1,1]:
                    self.rainMap[id] = None
                    self.rainVal[id] = 0.
                    self.rainTime[id,0] = tmpTime[p-1,1]
                    self.rainTime[id,1] = tmpTime[p,0]
                    id += 1
                self.rainMap[id] = tmpMap[p]
                self.rainTime[id,:] = tmpTime[p,:]
                self.rainVal[id] = tmpVal[p]
                id += 1
            if tmpTime[tmpNb-1,1] < self.tEnd:
                self.rainMap[id] = None
                self.rainVal[id] = 0.
                self.rainTime[id,0] = tmpTime[tmpNb-1,1]
                self.rainTime[id,1] = self.tEnd
        else:
            self.rainNb = 1
            self.rainVal = numpy.empty(self.rainNb)
            self.rainTime = numpy.empty((self.rainNb,2))
            self.rainMap = numpy.empty((self.rainNb),dtype=object)
            self.rainVal[0] = 0.
            self.rainTime[0,0] = self.tStart
            self.rainTime[0,1] = self.tEnd
            self.rainMap = None


        # Extract Stream Power Law structure parameters
        spl = None
        spl = root.find('spl')
        if spl is not None:
            element = None
            element = spl.find('dep')
            if element is not None:
                self.depo = int(element.text)
            else:
                self.depo = 1
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
            element = spl.find('maxdt')
            if element is not None:
                self.maxDT = float(element.text)
            else:
                self.maxDT = None
            element = None
            element = spl.find('ascale')
            if element is not None:
                self.alluvial = float(element.text)
            else:
                self.alluvial = 0.
            element = None
            element = spl.find('bscale')
            if element is not None:
                self.bedrock = float(element.text)
            else:
                self.bedrock = 0.
        else:
            self.depo = 0
            self.SPLm = 1.
            self.SPLn = 1.
            self.SPLero = 0.

        # Extract Linear Slope Diffusion structure parameters
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
        else:
            self.CDa = 0.
            self.CDm = 0.

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

        return
