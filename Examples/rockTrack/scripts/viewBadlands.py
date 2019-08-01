##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling companion.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
Here we set usefull functions used to analyse stratigraphic sequences from Badlands outputs.
"""

import warnings
warnings.filterwarnings('ignore')
warnings.simplefilter(action = "ignore", category = FutureWarning)

import re
import lavavu
import h5py
import numpy

def loadStep(folder, step,timev=True):
    """ Loading a Badlands output.

    Parameters
    ----------
    variable: folder
        Folder name of Badlands output
    variable: step
        Time step to load
    """

    # Load output files
    flow = h5py.File(folder+'/h5/flow.time'+str(step)+'.p0.hdf5', 'r')
    tin = h5py.File(folder+'/h5/tin.time'+str(step)+'.p0.hdf5', 'r')

    # Get the sea level history and depositional time for each stratigraphic layer
    with open(folder+'/xmf/tin.time'+str(step)+'.xmf') as fp:
        line = fp.readline()
        while line:
            if 'Time' in line:
                val = [float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
                time = val[0]
            if 'Function' in line:
                val = [float(s) for s in re.findall(r'-?\d+\.?\d*', line)]
                sea = val[2]
            line = fp.readline()
    if timev:
        print 'Rendering for time step '+str(step)+': '+str(time)+' years'

    return tin,flow,sea

def view1Step(tin, flow, sea, scaleZ=20,  maxZ=100, maxED=20,flowlines=False):
    """Visualise a given time step of Badlands output.

    Parameters
    ----------
    variable: tin
        TIN surface
    variable: flow
        Flow network
    variable: sea
        Sea level position
    variable: scaleZ
        Maximum vertical exageration
    variable: maxZ
        Maximum elevation (integer)
    variable: maxED
        Maximum erosion/deposition (integer)
    variable: flowlines
        Boolean to plot flow lines
    """
    lv = lavavu.Viewer(border=False, axis=False, background="grey80")

    #Plot the triangles
    tris = lv.triangles("elevation")
    tris2 = lv.triangles("ero/dep")
    tris3 = lv.triangles("discharge")

    #Load the vertices
    verts = numpy.array(tin["coords"])
    tris.vertices(verts)
    tris.indices(tin["cells"], offset=1)

    tris2.vertices(verts)
    tris2.indices(tin["cells"], offset=1)

    tris3.vertices(verts)
    tris3.indices(tin["cells"], offset=1)

    # Surface map
    tris.values(verts[:,2]-sea, 'height')
    cm = tris.colourmap("world", range=[-maxZ,maxZ])
    cb = tris.colourbar()

    # Erosion depostion map
    tris2.values(tin['cumdiff'], 'erosion/deposition')
    cm2 = tris2.colourmap("polar", range=[-maxED,maxED])
    cb2 = tris2.colourbar()

    # Discharge map
    tris3.values(tin['discharge'], 'discharge')
    cm3 = tris3.colourmap("abyss", reverse=True, log=True)
    cb3 = tris3.colourbar()

    # Flow lines
    if flowlines:
        flines = lv.lines("flow")
        flines.vertices(flow['coords'])
        flines.indices(flow['connect'], offset=1)
        #flines.values(flow['discharge'], 'discharge')
        flines.values(flow['chi'], 'chi')
        #flines.values(flow['sedload'], 'sedload')
        flines.values(flow['basin'], 'basin')
        flines.colourmap()
        flines.colourbar()

    lv.translation(-851.394, 4214.401, -51083.648)
    lv.rotation(-54.936, -14.238, -2.516)
    lv.scale('z', 5)

    lv.control.Panel()
    lv.control.ObjectList()
    lv.control.Range(command='scale z', range=(1,scaleZ), step=1., value=3)
    lv.control.Range(command='background', range=(0,1), step=0.1, value=1)

    if flowlines:
        flines.control.List(options=['chi', 'basin'], property='colourby', label='Flow colour field')
        flines.control.List(options=['spectral', 'cool', 'diverge', 'coolwarm'], property='colourmap', label='Flow colourmap', command="reload")
        flines.control.Range('linewidth', range=(0.5,7.5))

    lv.control.show()

    return

def viewTime(folder,steps=100,it=2,scaleZ=20, maxZ=100,maxED=20,flowlines=False):
    """Visualise sequence of time steps of Badlands output.

    Parameters
    ----------
    variable: folder
        Folder name of Badlands output
    variable: steps
        Time steps to load
    variable: it
        Output interval
    variable: scaleZ
        Maximum vertical exageration
    variable: maxZ
        Maximum elevation (integer)
    variable: maxED
        Maximum erosion/deposition (integer)
    variable: flowlines
        Boolean to plot flow lines
    """

    lv = lavavu.Viewer(border=False, axis=False, background="grey80")

    for s in range(0,steps,it):
        values = []

        tin,flow,sea = loadStep(folder, s, timev=False)

        #Add a new time step
        lv.addstep()

        #Plot the triangles
        if s == 0:
            tris = lv.triangles("elevation")
            tris2 = lv.triangles("ero/dep")
            tris3 = lv.triangles("discharge")

        # Load the vertices
        verts = numpy.array(tin["coords"])
        tris.vertices(verts)
        tris.indices(tin["cells"], offset=1)

        tris2.vertices(verts)
        tris2.indices(tin["cells"], offset=1)

        tris3.vertices(verts)
        tris3.indices(tin["cells"], offset=1)

        # Surface map
        tris.values(verts[:,2]-sea, 'height')
        if s == 0:
            cm = tris.colourmap("world", range=[-maxZ,maxZ])
            cb = tris.colourbar()

        # Erosion depostion map
        tris2.values(tin['cumdiff'], 'erosion/deposition')
        if s == 0:
            cm2 = tris2.colourmap("polar", range=[-maxED,maxED])
            cb2 = tris2.colourbar()

        # Discharge map
        tris3.values(tin['discharge'], 'discharge')
        if s == 0:
            cm3 = tris3.colourmap("abyss", reverse=True, log=True)
            cb3 = tris3.colourbar()

        # Flow lines
        if flowlines:
            if s == 0:
                flines = lv.lines("flow")

            flines.vertices(flow['coords'])
            flines.indices(flow['connect'], offset=1)

            #flines.values(flow['discharge'], 'discharge')
            flines.values(flow['chi'], 'chi')
            #flines.values(flow['sedload'], 'sedload')
            flines.values(flow['basin'], 'basinID')

            if s == 0:
                flines.colourmap()
                flines.colourbar()

    lv.control.Panel()
    lv.control.TimeStepper()
    lv.control.Range(command='scale z', range=(1,scaleZ), step=1., value=3)
    lv.control.Range(command='background', range=(0,1), step=0.1, value=1)

    if flowlines:
        #flines.control.List(options=['discharge', 'chi', 'sedload', 'basin'], property='colourby', label='Flow colour field')
        flines.control.List(options=['chi', 'basinID'], property='colourby', label='Flow colour field')
        flines.control.List(options=['spectral', 'cool', 'diverge', 'coolwarm'], property='colourmap', label='Flow colourmap', command="reload")
        flines.control.Range('linewidth', range=(0.5,7.5))

    lv.control.ObjectList()
    lv.control.show()

    return
