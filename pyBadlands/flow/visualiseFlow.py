##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module exports the flow network with associated parameters based on hdf5.
"""

import time
import h5py
import numpy
import xml.etree.ElementTree as ETO

def output_Polylines(outPts, rcvIDs, visXlim, visYlim, coordXY):
    """
    This function defines the connectivity array for visualising flow network.

    Parameters
    ----------
    outPts
        Numpy integer-type array containing the output node IDs.

    inIDs
        Numpy array containing the receiver node IDs.

    visXlim
        Numpy array containing the extent of visualisation grid.

    visYlim
        Numpy array containing the extent of visualisation grid.

    coordXY
        Numpy float-type array containing X, Y coordinates of the local nodes.

    Returns
    -------
    flowIDs
        Numpy integer-type array containing the output node IDs for the flow network.

    polyline
        Numpy 2D integer-type array containing the connectivity IDs for each polyline.
    """

    flowIDs = numpy.unique(numpy.concatenate((rcvIDs,outPts)))

    # For every element in rcvIDs array, find the index in flowIDs
    assert( len(numpy.setdiff1d(rcvIDs, flowIDs)) == 0)
    sort = numpy.argsort(flowIDs)
    pos = numpy.searchsorted(flowIDs[sort], rcvIDs)
    lrcvIDs = sort[pos]

    # For every element in outPts array, find the index in flowIDs
    assert( len(numpy.setdiff1d(outPts, flowIDs)) == 0)
    sort1 = numpy.argsort(flowIDs)
    pos1 = numpy.searchsorted(flowIDs[sort1], outPts)
    visIDs = sort1[pos1]

    # Define connectivity array
    connect = numpy.zeros((len(flowIDs),2), dtype=int)
    connect[:,0] = numpy.arange(len(flowIDs))
    connect[:,1] = numpy.arange(len(flowIDs))

    # Find outside domain nodes
    coords = coordXY[flowIDs,:2]
    border = numpy.where((coords[:,0] <= visXlim[0]) | (coords[:,0] >= visXlim[1]) |
                (coords[:,1] <= visYlim[0]) | (coords[:,1] >= visYlim[1]) )[0]

    # For every element in lrcvIDs array, find the index in border
    sortB = numpy.argsort(lrcvIDs)
    posB = numpy.searchsorted(lrcvIDs[sortB], border)
    bIDs = sortB[posB]
    lrcvIDs[bIDs] = -1

    # Trim the connectivity array
    connect[border,0] = -1
    connect[visIDs,1] = lrcvIDs
    connect += 1

    # Define polyline connectivity array
    lineID = numpy.where(connect[:,0] != connect[:,1])[0]
    line = connect[lineID,:2]
    lineIDs = numpy.where((line[:,1] > 0) & (line[:,0] > 0))[0]
    #polyline = line[lineIDs,:2]

    return flowIDs, line[lineIDs,:2]

def write_hdf5(folder, h5file, step, coords, elevation, discharge, chi,
               sedload, basin, connect, rank):
    """
    This function writes for each processor the HDF5 file containing flow network information.

    Parameters
    ----------
    folder
        Name of the output folder.

    h5file
        First part of the hdf5 file name.

    step
        Output visualisation step.

    coords
        Numpy float-type array containing X, Y coordinates of the local TIN nodes.

    elevation
        Numpy float-type array containing Z coordinates of the local TIN nodes.

    discharge
        Numpy float-type array containing the discharge values of the local TIN.

    chi
        Numpy float-type array containing the chi values of the local TIN.

    sedload
        Numpy float-type array containing the sedload values in kg/s.

    basin
        Numpy integer-type array containing the basin IDs values of the local TIN.

    connect
        Numpy 2D integer-type array containing the local nodes IDs for each connected network.

    rank
        ID of the local partition.
    """

    h5file = folder+'/'+h5file+str(step)+'.p'+str(rank)+'.hdf5'
    with h5py.File(h5file, "w") as f:

        # Write node coordinates and elevation
        f.create_dataset('coords',shape=(len(elevation),3), dtype='float32', compression='gzip')
        f["coords"][:,:2] = coords
        f["coords"][:,2] = elevation

        f.create_dataset('connect',shape=(len(connect[:,0]),2), dtype='int32', compression='gzip')
        f["connect"][:,:2] = connect

        f.create_dataset('basin',shape=(len(basin), 1), dtype='int32', compression='gzip')
        f["basin"][:,0] = basin

        f.create_dataset('chi',shape=(len(chi), 1), dtype='float32', compression='gzip')
        f["chi"][:,0] = chi

        f.create_dataset('sedload',shape=(len(sedload), 1), dtype='float32', compression='gzip')
        f["sedload"][:,0] = sedload

        f.create_dataset('discharge',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["discharge"][:,0] = discharge/3.154e7

def _write_xdmf(folder, xdmffile, xmffile, step):
    """
    This function writes the XDmF file which is calling the XmF file.

    Parameters
    ----------
    folder
        Name of the output folder.

    xdmffile
        XDmF file name.

    xmffile
        First part of the XmF file name.

    step
        Output visualisation step.
    """

    f= open(folder+'/'+xdmffile,'w')

    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
    f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
    f.write(' <Domain>\n')
    f.write('    <Grid GridType="Collection" CollectionType="Temporal">\n')

    for p in range(step+1):
        xfile = xmffile+str(p)+'.xmf'
        f.write('      <xi:include href="%s" xpointer="xpointer(//Xdmf/Domain/Grid)"/>\n' %xfile)

    f.write('    </Grid>\n')
    f.write(' </Domain>\n')
    f.write('</Xdmf>\n')
    f.close()

def write_xmf(folder, xmffile, xdmffile, step, time, elems, nodes, h5file, size):
    """
    This function writes the XmF file which is calling each HFD5 file.

    Parameters
    ----------

    folder
        Name of the output folder.

    xmffile
        First part of the XmF file name.

    step
        Output visualisation step.

    time
        Simulation time.

    elems
        Numpy integer-type array containing the number of elements of each local partition.

    nodes
        Numpy integer-type array containing the number of nodes of each local partition.

    h5file
        First part of the hdf5 file name.

    size
        Number of partitions.
    """

    xmf_file = folder+'/'+xmffile+str(step)+'.xmf'
    f= open(str(xmf_file),'w')

    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
    f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
    f.write(' <Domain>\n')
    f.write('    <Grid GridType="Collection" CollectionType="Spatial">\n')
    f.write('      <Time Type="Single" Value="%s"/>\n'%time)

    for p in range(size):
        pfile = h5file+str(step)+'.p'+str(p)+'.hdf5'
        f.write('      <Grid Name="Block.%s">\n' %(str(p)))
        f.write('         <Topology Type="Polyline" NodesPerElement="2" ')
        f.write('NumberOfElements="%d" BaseOffset="1">\n'%elems[p])
        f.write('          <DataItem Format="HDF" DataType="Int" ')
        f.write('Dimensions="%d 2">%s:/connect</DataItem>\n'%(elems[p],pfile))
        f.write('         </Topology>\n')

        f.write('         <Geometry Type="XYZ">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 3">%s:/coords</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Geometry>\n')

        f.write('         <Attribute Type="Scalar" Center="Node" Name="BasinID">\n')
        f.write('          <DataItem Format="HDF" NumberType="Integer" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/basin</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

        f.write('         <Attribute Type="Scalar" Center="Node" Name="Discharge [m3/s]">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/discharge</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

        f.write('         <Attribute Type="Scalar" Center="Node" Name="Chi">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/chi</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

        f.write('         <Attribute Type="Scalar" Center="Node" Name="sedload [m3/s]">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/sedload</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

        f.write('      </Grid>\n')


    f.write('    </Grid>\n')
    f.write(' </Domain>\n')
    f.write('</Xdmf>\n')
    f.close()

    _write_xdmf(folder, xdmffile, xmffile, step)
