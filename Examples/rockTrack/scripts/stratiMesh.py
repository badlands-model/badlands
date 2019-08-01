##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling companion.      ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
Here we define the main functions to create the TIN based stratigraphic mesh from Badlands outputs.

Example of calling method:

mesh = stratiMesh(folder='output2', xdmfName = 'stratal_series', ncpus=1, layperstep=2, dispTime=5000.)
mesh.outputSteps(startTime=240000.,endTime=245000.)

"""

import os
import math
import h5py
import errno
import numpy as np

import warnings
warnings.simplefilter(action = "ignore", category = FutureWarning)

class stratiMesh:
    """
    Class for creating irregular stratigraphic mesh from Badlands outputs.
    """

    def __init__(self, folder=None, xdmfName = 'stratal_series', ncpus=1, layperstep=1, dispTime=None):
        """
        Initialization function which takes the folder path to Badlands outputs
        and the number of CPUs used to run the simulation.
        Parameters
        ----------
        variable : folder
            Folder path to Badlands outputs.
        variable : xdmfName
            Name of Badlands stratigraphic grid outputs.
        variable: ncpus
            Number of CPUs used to run the simulation.
        variable: layperstep
            Number of layers created between each display
            interval (obtained from the XmL input file).
        variable: dispTime
            Time interval in years used to display badlands outputs.
        """

        self.folder = folder
        if not os.path.isdir(folder):
            raise RuntimeError('The given folder cannot be found or the path is incomplete.')

        self.ncpus = ncpus

        self.x = None
        self.y = None
        self.elev = None
        self.dep = None
        self.th = None
        self.timestep = 0
        self.layperstep = 0
        self.laynb = 0
        self.startStep = None
        self.endStep = None
        self.rockNb = None

        # Assign file names
        self.h5TIN = 'h5/tin.time'
        self.h5Strat = 'h5/stratal.time'
        self.h5Strati = 'h5/strati.time'
        self.xmffile = 'xmf/stratal.time'
        self.xdmfName = xdmfName+'.xdmf'
        self.dispTime = dispTime
        self.tnow = None

        return

    def _loadTIN(self, step, rank):
        """
        Load TIN grid to extract cells connectivity and vertices position.
        Parameters
        ----------
        variable : step
            Specific step at which the TIN variables will be read.
        variable: rank
            TIN file for the considered CPU.
        """

        h5file = self.folder+'/'+self.h5TIN+str(step)+'.p'+str(rank)+'.hdf5'
        df = h5py.File(h5file, 'r')
        coords = np.array((df['/coords']))
        cells = np.array((df['/cells']),dtype=int)

        return coords, cells

    def _loadStrati(self, step, rank):
        """
        Load stratigraphic dataset.
        Parameters
        ----------
        variable : step
            Specific step at which the TIN variables will be read.
        variable: rank
            Stratigraphic grid for the considered CPU.
        """

        h5file = self.folder+'/'+self.h5Strat+str(step)+'.p'+str(rank)+'.hdf5'
        df = h5py.File(h5file, 'r')
        rockNb = len(df.keys())-1
        paleoH = np.array((df['/paleoDepth']))
        rockProp = np.zeros((paleoH.shape[0],paleoH.shape[1],rockNb))

        for r in range(rockNb):
            rockProp[:,:,r] = np.array((df['/depoThickRock'+str(r)]))

        return rockProp, paleoH

    def _write_hdf5(self, xt, yt, zt, cellt, layID, layH, layD, propR,
                    clayID, clayD, clayH, cpropR, step, rank):
        """
        Write the HDF5 file containing the stratigraphic mesh variables.
        Parameters
        ----------
        variable : xt
            X-axis coordinates of the vertices.
        variable : yt
            Y-axis coordinates of the vertices.
        variable : zt
            Z-axis coordinates of the vertices.
        variable : cellt
            Wedge cells connectivity.
        variable : layID
            ID of each layer.
        variable : layH
            Thickness of each layer.
        variable : propR
            Proportion of each rock type in each layer.
        variable : layD
            Paleo depth informs about the elevation at time of deposition.
        variable : clayID
            Cell ID of each layer.
        variable : clayH
            Cell thickness of each layer.
        variable : cpropR
            Cell proportion of each rock type in each layer.
        variable : clayD
            Cell paleo depth informs about the elevation at time of deposition.
        variable: rank
            TIN file for the considered CPU.
        """

        h5file = self.folder+'/'+self.h5Strati+str(step)+'.p'+str(rank)+'.hdf5'
        with h5py.File(h5file, "w") as f:

            # Write node coordinates and elevation
            f.create_dataset('coords',shape=(len(xt),3), dtype='float32', compression='gzip')
            f["coords"][:,0] = xt
            f["coords"][:,1] = yt
            f["coords"][:,2] = zt

            f.create_dataset('cells',shape=(len(cellt[:,0]),6), dtype='int32', compression='gzip')
            f["cells"][:,:] = cellt

            f.create_dataset('layID',shape=(len(xt),1), dtype='int32', compression='gzip')
            f["layID"][:,0] = layID

            f.create_dataset('clayID',shape=(len(cellt[:,0]),1), dtype='int32', compression='gzip')
            f["clayID"][:,0] = clayID

            f.create_dataset('layH',shape=(len(xt),1), dtype='float32', compression='gzip')
            f["layH"][:,0] = layH

            f.create_dataset('clayH',shape=(len(cellt[:,0]),1), dtype='float32', compression='gzip')
            f["clayH"][:,0] = clayH

            f.create_dataset('layD',shape=(len(xt),1), dtype='float32', compression='gzip')
            f["layD"][:,0] = layD

            f.create_dataset('clayD',shape=(len(cellt[:,0]),1), dtype='float32', compression='gzip')
            f["clayD"][:,0] = clayD

            for r in range(self.rockNb):
                attribute = 'propR'+str(r)
                f.create_dataset(attribute,shape=(len(xt),1), dtype='float32', compression='gzip')
                f[attribute][:,0] = propR[:,r]

                attribute = 'cpropR'+str(r)
                f.create_dataset(attribute,shape=(len(cellt[:,0]),1), dtype='float32', compression='gzip')
                f[attribute][:,0] = cpropR[:,r]

        return

    def _write_xmf(self, step, elems, nodes):
        """
        Write the XMF file which load and read the hdf5 parameters files at any given step.
        Parameters
        ----------
        variable : step
            Specific step at which the TIN variables will be read.
        variable: elems
            Number of wedges elements per processor mesh.
        variable: nodes
            Number of irregular points per processor mesh
        """

        xmf_file = self.folder+'/'+self.xmffile+str(step)+'.xmf'
        f= open(str(xmf_file),'w')

        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        f.write(' <Domain>\n')
        f.write('    <Grid GridType="Collection" CollectionType="Spatial">\n')
        f.write('      <Time Type="Single" Value="%s"/>\n'%self.tnow)

        for p in range(self.ncpus):
            pfile = self.h5Strati+str(step)+'.p'+str(p)+'.hdf5'
            f.write('      <Grid Name="Block.%s">\n' %(str(p)))
            f.write('         <Topology Type="Wedge" NumberOfElements="%d" BaseOffset="1">\n'%elems[p])
            f.write('          <DataItem Format="HDF" DataType="Int" ')
            f.write('Dimensions="%d 6">%s:/cells</DataItem>\n'%(elems[p],pfile))
            f.write('         </Topology>\n')

            f.write('         <Geometry Type="XYZ">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 3">%s:/coords</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Geometry>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="layer ID">\n')
            f.write('          <DataItem Format="HDF" NumberType="Int" ')
            f.write('Dimensions="%d 1">%s:/layID</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Cell" Name="layer ID">\n')
            f.write('          <DataItem Format="HDF" NumberType="Int" ')
            f.write('Dimensions="%d 1">%s:/clayID</DataItem>\n'%(elems[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="paleo-depth">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/layD</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Cell" Name="paleo-depth">\n')
            f.write('          <DataItem Format="HDF" NumberType="float" ')
            f.write('Dimensions="%d 1">%s:/clayD</DataItem>\n'%(elems[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="layer th">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/layH</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Cell" Name="layer th">\n')
            f.write('          <DataItem Format="HDF" NumberType="float" ')
            f.write('Dimensions="%d 1">%s:/clayH</DataItem>\n'%(elems[p],pfile))
            f.write('         </Attribute>\n')

            for r in range(self.rockNb):
                f.write('         <Attribute Type="Scalar" Center="Node" Name="prop rock%d">\n'%r)
                f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
                f.write('Dimensions="%d 1">%s:/propR%d</DataItem>\n'%(nodes[p],pfile,r))
                f.write('         </Attribute>\n')

                f.write('         <Attribute Type="Scalar" Center="Cell" Name="prop rock%d">\n'%r)
                f.write('          <DataItem Format="HDF" NumberType="Float"  ')
                f.write('Dimensions="%d 1">%s:/cpropR%d</DataItem>\n'%(elems[p],pfile,r))
                f.write('         </Attribute>\n')

            f.write('      </Grid>\n')

        f.write('    </Grid>\n')
        f.write(' </Domain>\n')
        f.write('</Xdmf>\n')
        f.close()

        return

    def _write_xdmf(self):
        """
        Write the XDMF file which load and read the XMF parameters files for the requested steps.
        """

        f= open(self.folder+'/'+self.xdmfName,'w')

        f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
        f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
        f.write(' <Domain>\n')
        f.write('    <Grid GridType="Collection" CollectionType="Temporal">\n')

        for p in range(self.startStep,self.endStep+1):
            xfile = self.xmffile+str(p)+'.xmf'
            f.write('      <xi:include href="%s" xpointer="xpointer(//Xdmf/Domain/Grid)"/>\n' %xfile)

        f.write('    </Grid>\n')
        f.write(' </Domain>\n')
        f.write('</Xdmf>\n')
        f.close()

        return

    def outputSteps(self, startTime=0, endTime=5000):
        """
        Define the steps that need to be visualise.
        Parameters
        ----------
        variable : startTime
            First Badlands output time to visualise.
        variable: endStep
            Last Badlands output time to visualise.
        """

        self.startTime = startTime
        self.endTime = endTime
        self.tnow = startTime

        self.startStep = int(startTime/self.dispTime)
        self.endStep = int(endTime/self.dispTime)
        if not os.path.isdir(self.folder):
            raise RuntimeError('The given folder cannot be found or the path is incomplete.')

        assert self.startStep<=self.endStep, 'ERROR: End step lower than Start step.'

        for s in range(self.startStep,self.endStep+1):
            print 'Process layers at time [in years]: ',self.tnow
            ptsnb = []
            cellnb = []
            for i in range(0, self.ncpus):

                # Load TIN grid for specific time step
                coords, cells = self._loadTIN(s,i)
                x, y, z = np.hsplit(coords, 3)

                # Load STRATI dataset
                rockTH, paleoH = self._loadStrati(s,i)
                tmpRock = rockTH[:,1:]
                tmpPaleo = paleoH[:,1:]

                # Define dimensions
                ptsNb = len(x)
                layNb = paleoH.shape[1]
                rockNb = rockTH.shape[2]
                if self.rockNb is None:
                    self.rockNb = rockNb

                # Build attributes:
                # Layer number attribute
                #layID = np.array([np.arange(layNb+1),]*ptsNb,dtype=int)
                layID = np.array([np.arange(layNb),]*ptsNb,dtype=int)
                ltmp = layID.flatten(order='F')

                # Thickness of each layer
                layTH = np.sum(tmpRock,axis=2)
                htmp = layTH.flatten(order='F')
                htmp = np.concatenate((np.zeros(ptsNb),htmp), axis=0)

                # Paleo-depth of each layer
                dtmp = tmpPaleo.flatten(order='F')
                dtmp = np.concatenate((dtmp,z[:,0]), axis=0)

                # Elevation of each layer
                cumH = np.cumsum(layTH[:,::-1],axis=1)[:,::-1]
                layZ = z - cumH
                # Add the top surface to the layer elevation record
                ztmp = layZ.flatten(order='F')
                ztmp = np.concatenate((ztmp, z[:,0]), axis=0)

                # Proportion of each rock type
                r,c = np.where(layTH>0)

                rockProp = np.zeros((tmpRock.shape))
                for k in range(rockNb):
                    rockProp[r,c,k] = np.divide(tmpRock[r,c,k],layTH[r,c])
                r1,c1 = np.where(layTH<=0)

                for cc in range(1,layNb-1):
                    tt = np.where(layTH[:,cc]<=0)
                    for k in range(rockNb):
                        rockProp[tt,cc,k] = rockProp[tt,cc-1,k]

                ptmp = rockProp.reshape((ptsNb*(layNb-1),rockNb),order='F')
                ptmp = np.concatenate((np.zeros((ptsNb,rockNb)),ptmp), axis=0)

                # Creation of each layer coordinates
                xtmp = x[:,0]
                ytmp = y[:,0]
                ctmp = cells
                oldcells = ctmp

                # Cell based rock property attributes
                cellP = np.zeros((len(cells),rockNb))

                # Cell layer index
                cellI = np.array([np.arange(1,layNb),]*len(cells),dtype=int)
                cellItmp = cellI.flatten(order='F')

                for l in range(1,layNb):
                    xtmp = np.concatenate((xtmp, x[:,0]), axis=0)
                    ytmp = np.concatenate((ytmp, y[:,0]), axis=0)

                    # Creation of each layer elements
                    newcells = oldcells+len(x)
                    celltmp = np.concatenate((oldcells, newcells), axis=1)
                    cellD = np.sum(dtmp[newcells-1],axis=1)/3.
                    cellH = np.sum(htmp[newcells-1],axis=1)/3.
                    for k in range(rockNb):
                        cellP[:,k] = np.sum(ptmp[newcells-1,k],axis=1)/3.
                    oldcells = newcells
                    if l > 1:
                        ctmp = np.concatenate((ctmp, celltmp), axis=0)
                        cellDtmp = np.concatenate((cellDtmp, cellD), axis=0)
                        cellHtmp = np.concatenate((cellHtmp, cellH), axis=0)
                        cellPtmp = np.concatenate((cellPtmp, cellP), axis=0)
                    else:
                        ctmp = np.copy(celltmp)
                        cellDtmp = np.copy(cellD)
                        cellHtmp = np.copy(cellH)
                        cellPtmp = np.copy(cellP)
                    cellP.fill(0.)

                # Create the HDF5 file
                self._write_hdf5(xtmp, ytmp, ztmp,
                                 ctmp, ltmp, htmp, dtmp, ptmp,
                                 cellItmp, cellDtmp, cellHtmp,
                                 cellPtmp, s, i)

                cellnb.append(len(ctmp))
                ptsnb.append(len(xtmp))

            # Create the XMF file
            self._write_xmf(s, cellnb, ptsnb)
            self.tnow += self.dispTime

        # Create the XDMF file
        self._write_xdmf()

        return
