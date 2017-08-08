##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
This module exports the TIN surface with associated parameters based on hdf5.
"""

import h5py
import numpy
import xml.etree.ElementTree as ETO

def output_cellsIDs(lGIDs, inIDs, visXlim, visYlim, coords, cells):
    """
    This function defines the cells used for visualising the TIN surface.

    Parameters
    ----------
    lGIDs
        Numpy integer-type array filled with the global vertex IDs for each local grid located
        within the partition (including those on the edges).

    inIDs
        Numpy integer-type array filled with the global vertex IDs for each local grid located
        within the partition (not those on the edges).

    visXlim, visYlim
        Numpy array containing the extent of visualisation grid.

    coords
        Numpy float-type array containing X, Y coordinates of the local TIN nodes.

    cells
        Numpy integer-type array filled with the global cell IDs.

    Returns
    -------
    outPts
        Numpy integer-type array containing the output node IDs.

    cells
        Numpy integer-type array containing the output cell IDs.
    """

    # Find non-overlapping vertices in each local TIN
    findInside = numpy.in1d(lGIDs, inIDs)
    inside = numpy.where(findInside == True)

    # Find cells containing only the inside indices
    inCell = numpy.in1d(cells, inside).reshape(cells.shape)
    intCell = 1*inCell
    sumCell = intCell.sum(axis=1)
    localCell = numpy.where(sumCell>0)[0]
    outcell = cells[localCell]

    # Get the non-border points IDs
    notBorder = numpy.where((coords[lGIDs,0] >= visXlim[0]) & (coords[lGIDs,0] <= visXlim[1]) &
                         (coords[lGIDs,1] >= visYlim[0]) & (coords[lGIDs,1] <= visYlim[1]) )[0]
    notBorderIDs = numpy.zeros(len(lGIDs),dtype=int)
    notBorderIDs.fill(-1)
    notBorderIDs[notBorder] = lGIDs[notBorder]
    findInside2 = numpy.in1d(lGIDs, notBorderIDs)
    inside2 = numpy.where(findInside2 == True)

    inCell2 = numpy.in1d(outcell, inside2).reshape(outcell.shape)
    intCell2 = 1*inCell2
    sumCell2 = intCell2.sum(axis=1)
    localCell2 = numpy.where(sumCell2>2)[0]

    # Get inside nodes
    rav = numpy.ravel(outcell[localCell2])
    allInside = numpy.unique(rav)

    return lGIDs[allInside], outcell[localCell2] + 1

def write_hdf5(folder, h5file, step, coords, elevation, rain, discharge, cumdiff, cumhill,
               cells, rank, rainOn, eroOn, erodibility, area, waveOn, waveH, waveU, waveV,
               rockOn, prop):
    """
    This function writes for each processor the HDF5 file containing surface information.

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

    sealevel
        Level of sea water

    elevation
        Numpy float-type array containing Z coordinates of the local TIN nodes.

    rain
        Numpy float-type array containing rain value of the local TIN nodes.

    discharge
        Numpy float-type array containing the discharge values of the local TIN.

    cumdiff
        Numpy float-type array containing the cumulative elevation changes values of the local TIN.

    cumhill
        Numpy float-type array containing the cumulative elevation changes values for hillslope of the local TIN.

    erodibility
        Numpy float-type array containing the top surface erodibility values of the local TIN.

    cells
        Numpy integer-type array filled with the global cell IDs.

    rank
        ID of the local partition.

    rainOn
        Boolean for orographic precipitation.

    eroOn
        Boolean for erodibility values.

    waveOn
        Boolean for wave activity.

    waveH
        Average wave height.

    waveU,waveV
        Average wave velocity.
    """

    h5file = folder+'/'+h5file+str(step)+'.p'+str(rank)+'.hdf5'
    with h5py.File(h5file, "w") as f:

        # Write node coordinates and elevation
        f.create_dataset('coords',shape=(len(elevation),3), dtype='float32', compression='gzip')
        f["coords"][:,:2] = coords
        f["coords"][:,2] = elevation

        f.create_dataset('cells',shape=(len(cells[:,0]),3), dtype='int32', compression='gzip')
        f["cells"][:,:] = cells

        f.create_dataset('discharge',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["discharge"][:,0] = discharge

        if rainOn:
            f.create_dataset('precipitation',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["precipitation"][:,0] = rain

        if eroOn:
                f.create_dataset('erodibility',shape=(len(discharge), 1), dtype='float32', compression='gzip')
                f["erodibility"][:,0] = erodibility

        f.create_dataset('cumdiff',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["cumdiff"][:,0] = cumdiff

        f.create_dataset('cumhill',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["cumhill"][:,0] = cumhill

        f.create_dataset('area',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["area"][:,0] = area

        if waveOn:
            f.create_dataset('waveH',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["waveH"][:,0] = waveH
            f.create_dataset('waveU',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["waveU"][:,0] = waveU
            f.create_dataset('waveV',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["waveV"][:,0] = waveV

        if rockOn:
            f.create_dataset('depClastic',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["depClastic"][:,0] = prop[:,0]
            f.create_dataset('depSpecies1',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["depSpecies1"][:,0] = prop[:,1]
            f.create_dataset('depSpecies2',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["depSpecies2"][:,0] = prop[:,2]
            f.create_dataset('depPelagic',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["depPelagic"][:,0] = prop[:,3]

def write_hdf5_flexure(folder, h5file, step, coords, elevation, rain, discharge, cumdiff, cumhill,
                       cumflex, cells, rank, rainOn, eroOn, erodibility, area, waveOn, waveH,  waveU, waveV,
                       rockOn, prop):
    """
    This function writes for each processor the HDF5 file containing surface information.

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

    rain
        Numpy float-type array containing rain value of the local TIN nodes.

    discharge
        Numpy float-type array containing the discharge values of the local TIN.

    cumdiff
        Numpy float-type array containing the cumulative elevation changes values of the local TIN.

    cumhill
        Numpy float-type array containing the cumulative elevation changes values for hillslope of the local TIN.

    erodibility
        Numpy float-type array containing the top surface erodibility values of the local TIN.

    cumflex
        Numpy float-type array containing the cumulative flexural changes values of the local TIN.

    cells
        Numpy integer-type array filled with the global cell IDs.

    rank
        ID of the local partition.

    rainOn
        Boolean for orographic precipitation.

    eroOn
        Boolean for erodibility values.

    waveOn
        Boolean for wave activity.

    waveH
        Average wave height.

    waveU,waveV
        Average wave velocity.
    """

    h5file = folder+'/'+h5file+str(step)+'.p'+str(rank)+'.hdf5'
    with h5py.File(h5file, "w") as f:

        # Write node coordinates and elevation
        f.create_dataset('coords',shape=(len(elevation),3), dtype='float32', compression='gzip')
        f["coords"][:,:2] = coords
        f["coords"][:,2] = elevation

        f.create_dataset('cells',shape=(len(cells[:,0]),3), dtype='int32', compression='gzip')
        f["cells"][:,:] = cells

        f.create_dataset('discharge',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["discharge"][:,0] = discharge

        if rainOn:
            f.create_dataset('precipitation',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["precipitation"][:,0] = rain

        if eroOn:
            f.create_dataset('erodibility',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["erodibility"][:,0] = erodibility

        f.create_dataset('cumdiff',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["cumdiff"][:,0] = cumdiff

        f.create_dataset('cumhill',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["cumhill"][:,0] = cumhill

        f.create_dataset('cumflex',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["cumflex"][:,0] = cumflex

        f.create_dataset('area',shape=(len(discharge), 1), dtype='float32', compression='gzip')
        f["area"][:,0] = area

        if waveOn:
            f.create_dataset('waveH',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["waveH"][:,0] = waveH
            f.create_dataset('waveU',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["waveU"][:,0] = waveU
            f.create_dataset('waveV',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["waveV"][:,0] = waveV

        if rockOn:
            f.create_dataset('depClastic',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["depClastic"][:,0] = prop[:,0]
            f.create_dataset('depSpecies1',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["depSpecies1"][:,0] = prop[:,1]
            f.create_dataset('depSpecies2',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["depSpecies2"][:,0] = prop[:,2]
            f.create_dataset('depPelagic',shape=(len(discharge), 1), dtype='float32', compression='gzip')
            f["depPelagic"][:,0] = prop[:,3]

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

    return

def write_xmf(folder, xmffile, xdmffile, step, t, elems, nodes, h5file, sealevel, size,
              flexOn, rainOn, eroOn, waveOn, rockOn):
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

    t
        Simulation time.

    elems
        Numpy integer-type array containing the number of cells of each local partition.

    nodes
        Numpy integer-type array containing the number of nodes of each local partition.

    h5file
        First part of the hdf5 file name.

    sealevel
        Sealevel elevation.

    size
        Number of partitions.

    flexOn
        Boolean for flexural isostasy.

    rainOn
        Boolean for orographic precipitation.

    eroOn
        Boolean for erodibility values.

    waveOn
        Boolean for wave activity.
    """

    xmf_file = folder+'/'+xmffile+str(step)+'.xmf'
    f= open(str(xmf_file),'w')

    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
    f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
    f.write(' <Domain>\n')
    f.write('    <Grid GridType="Collection" CollectionType="Spatial">\n')
    f.write('      <Time Type="Single" Value="%s"/>\n'%t)

    for p in range(size):
        pfile = h5file+str(step)+'.p'+str(p)+'.hdf5'
        f.write('      <Grid Name="Block.%s">\n' %(str(p)))
        f.write('         <Topology Type="Triangle" NumberOfElements="%d" BaseOffset="1">\n'%elems[p])
        f.write('          <DataItem Format="HDF" DataType="Int" ')
        f.write('Dimensions="%d 3">%s:/cells</DataItem>\n'%(elems[p],pfile))
        f.write('         </Topology>\n')

        f.write('         <Geometry Type="XYZ">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 3">%s:/coords</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Geometry>\n')

        f.write('         <Attribute Type="Scalar" Center="Node" Name="Discharge">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/discharge</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

        f.write('         <Attribute Type="Scalar" Center="Node" Name="EroDep">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/cumdiff</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

        f.write('         <Attribute Type="Scalar" Center="Node" Name="EroDep hillslope">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/cumhill</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

        if flexOn:
            f.write('         <Attribute Type="Scalar" Center="Node" Name="Cumflex">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/cumflex</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

        if rainOn:
            f.write('         <Attribute Type="Scalar" Center="Node" Name="Rain">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/precipitation</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

        if eroOn:
            f.write('         <Attribute Type="Scalar" Center="Node" Name="Ke">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/erodibility</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

        if waveOn:
            f.write('         <Attribute Type="Scalar" Center="Node" Name="waveH">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/waveH</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="waveU">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/waveU</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="waveV">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/waveV</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

        if rockOn:
            f.write('         <Attribute Type="Scalar" Center="Node" Name="depClastic">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/depClastic</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="depSpecies1">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/depSpecies1</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="depSpecies2">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/depSpecies2</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

            f.write('         <Attribute Type="Scalar" Center="Node" Name="depPelagic">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/depPelagic</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

        f.write('         <Attribute Type="Scalar" Center="Node" Name="Sealevel">\n')
        f.write('          <DataItem ItemType="Function" Function="$0 * 0.00000000001 + %f" Dimensions="%d 1">\n'%(sealevel,nodes[p]))
        f.write('           <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/cumdiff</DataItem>\n'%(nodes[p],pfile))
        f.write('          </DataItem>\n')
        f.write('         </Attribute>\n')

        f.write('      </Grid>\n')


    f.write('    </Grid>\n')
    f.write(' </Domain>\n')
    f.write('</Xdmf>\n')
    f.close()

    _write_xdmf(folder, xdmffile, xmffile, step)
