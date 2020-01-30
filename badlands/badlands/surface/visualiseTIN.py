##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
##                                                                                   ##
##  This file forms part of the Badlands surface processes modelling application.    ##
##                                                                                   ##
##  For full license and copyright information, please refer to the LICENSE.md file  ##
##  located at the project root, or contact the authors.                             ##
##                                                                                   ##
##~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~#~##
"""
The set of functions below are used to export the TIN mesh and the associated mesh variables
computed with **badlands**.


**Badlands** uses `hdf5`_ library for outputs generation.
To read the **hdf5** files, the code creates a *XML schema* (**.xmf** files) describing
how the data is stored in the **hdf5** files.

Then `xdmf`_ files provides support to visualise temporal evolution of the output with
applications such as **Paraview** or Visit.

.. _`hdf5`: https://www.hdfgroup.org/HDF5/
.. _`xdmf`: http://www.xdmf.org/

.. image:: img/output.png
   :scale: 50 %
   :alt: Badlands output folder
   :align: center

Above is a typical structure that you will find in **badlands** output folder:

- **h5** folder contains the **hdf5** data, all the information computed by the model are stored in these files. You will have at least the *tin* (surface) and *flow* (stream network) dataset and also the *sed* (stratigraphy) data if the stratal structure is computed in your simulation.
- **xmf** folder contains the XML files used to read the **hdf5** files contained in the **h5** folder.
- the **.xml** input file used to build this specific model.
- two **.xdmf** files for the surface (**tin_series.xdmf**) and the flow network (**flow_series.xdmf**) that read the **xmf** files through time.


Important:
    The TIN outputs are **hdf5** files.
"""

import h5py
import numpy
import xml.etree.ElementTree as ETO

def output_cellsIDs(lGIDs, inIDs, visXlim, visYlim, coords, cells):
    """
    This function defines the cells used for visualising the TIN surface.

    Args:
        lGIDs: numpy integer-type local nodes array filled with global node IDs including partition edges.
        inIDs: numpy integer-type local nodes array filled with global node IDs excluding partition edges.
        visXlim: numpy array containing the X-axis extent of visualisation grid.
        visYlim: numpy array containing the Y-axis extent of visualisation grid.
        coords: 2D numpy float-type array containing X, Y coordinates of the local TIN nodes.
        cells: numpy integer-type array filled with global cell IDs.

    Returns
    -------
    outPts
        numpy integer-type array containing the output node IDs.
    cells
        numpy integer-type array containing the output cell IDs.
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

def write_hdf5(folder, h5file, step, coords, elevation, rain, discharge, cumdiff, cumhill, cumfail,
               cells, rainOn, eroOn, erodibility, area, waveOn, waveH, waveS, wavediff,
               rockOn, prop, sealevel):
    """
    This function writes for each processor the **hdf5** file containing surface information.

    Args:
        folder: name of the output folder.
        h5file: first part of the hdf5 file name.
        step: output visualisation step.
        coords: numpy float-type array containing X, Y coordinates of the local TIN nodes.
        sealevel: level of sea water
        elevation: numpy float-type array containing Z coordinates of the local TIN nodes.
        rain: numpy float-type array containing rain value of the local TIN nodes.
        discharge: numpy float-type array containing the discharge values of the local TIN.
        cumdiff: numpy float-type array containing the cumulative elevation changes values of the local TIN.
        cumhill: numpy float-type array containing the cumulative elevation changes values for hillslope of the local TIN.
        wavediff: numpy float-type array containing the cumulative elevation changes values for wave of the local TIN.
        erodibility: numpy float-type array containing the top surface erodibility values of the local TIN.
        cells: numpy integer-type array filled with the global cell IDs.
        rainOn: boolean for orographic precipitation.
        eroOn: boolean for erodibility values.
        waveOn: boolean for wave activity.
        waveH: average wave height.
        waveS: average shear stress on bed.
    """

    h5file = folder+'/'+h5file+str(step)+'.hdf5'
    with h5py.File(h5file, "w") as f:

        # Write node coordinates and elevation
        f.create_dataset('coords',shape=(len(elevation),3), dtype='float64', compression='gzip')
        f["coords"][:,:2] = coords
        f["coords"][:,2] = elevation

        f.create_dataset('cells',shape=(len(cells[:,0]),3), dtype='int32', compression='gzip')
        f["cells"][:,:] = cells

        f.create_dataset('discharge',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        discharge[elevation<sealevel] = 0.
        discharge[discharge<1.] = 1.
        f["discharge"][:,0] = discharge

        if rainOn:
            f.create_dataset('precipitation',shape=(len(discharge), 1), dtype='float64', compression='gzip')
            f["precipitation"][:,0] = rain

        if eroOn:
                f.create_dataset('erodibility',shape=(len(discharge), 1), dtype='float64', compression='gzip')
                f["erodibility"][:,0] = erodibility

        f.create_dataset('cumdiff',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        f["cumdiff"][:,0] = cumdiff

        f.create_dataset('cumhill',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        f["cumhill"][:,0] = cumhill

        f.create_dataset('cumfail',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        f["cumfail"][:,0] = cumfail

        f.create_dataset('area',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        f["area"][:,0] = area

        if waveOn:
            f.create_dataset('waveH',shape=(len(discharge), 1), dtype='float64', compression='gzip')
            f["waveH"][:,0] = waveH
            f.create_dataset('waveS',shape=(len(discharge), 1), dtype='float64', compression='gzip')
            f["waveS"][:,0] = waveS
            f.create_dataset('cumwave',shape=(len(discharge), 1), dtype='float64', compression='gzip')
            f["cumwave"][:,0] = wavediff

        if rockOn:
            nbSed = prop.shape[1]
            for k in range(nbSed):
                if k == 0:
                    f.create_dataset('depClastic',shape=(len(discharge), 1), dtype='float64', compression='gzip')
                    f["depClastic"][:,0] = prop[:,k]
                elif k == 1:
                    f.create_dataset('depSpecies1',shape=(len(discharge), 1), dtype='float64', compression='gzip')
                    f["depSpecies1"][:,0] = prop[:,k]
                elif k == 2:
                    f.create_dataset('depSpecies2',shape=(len(discharge), 1), dtype='float64', compression='gzip')
                    f["depSpecies2"][:,0] = prop[:,k]

def write_hdf5_flexure(folder, h5file, step, coords, elevation, rain, discharge, cumdiff, cumhill, cumfail,
                       cumflex, cells, rainOn, eroOn, erodibility, area, waveOn, waveH,  waveS, wavediff,
                       rockOn, prop, sealevel):
    """
    This function writes for each processor the **hdf5** file containing surface information
    with flexural isostasy turned on.

    Args:
        folder: name of the output folder.
        h5file: first part of the hdf5 file name.
        step: output visualisation step.
        coords: numpy float-type array containing X, Y coordinates of the local TIN nodes.
        sealevel: level of sea water
        elevation: numpy float-type array containing Z coordinates of the local TIN nodes.
        rain: numpy float-type array containing rain value of the local TIN nodes.
        discharge: numpy float-type array containing the discharge values of the local TIN.
        cumdiff: numpy float-type array containing the cumulative elevation changes values of the local TIN.
        cumhill: numpy float-type array containing the cumulative elevation changes values for hillslope of the local TIN.
        wavediff: numpy float-type array containing the cumulative elevation changes values for wave of the local TIN.
        erodibility: numpy float-type array containing the top surface erodibility values of the local TIN.
        cumflex: numpy float-type array containing the cumulative flexural changes values of the local TIN.
        cells: numpy integer-type array filled with the global cell IDs.
        rainOn: boolean for orographic precipitation.
        eroOn: boolean for erodibility values.
        waveOn: boolean for wave activity.
        waveH: average wave height.
        waveS: average shear stress on bed.


    Note:
        Flexural isostasy is obtained from the gFlex_ package available on Github!

        Wickert, A. D. (2016), Open-source modular solutions for flexural isostasy: gFlex v1.0,
        Geosci. Model Dev., 9(3), 997–1017, `doi:10.5194/gmd-9-997-2016`_.

    .. _`doi:10.5194/gmd-9-997-2016`:  https://doi.org/10.5194/gmd-9-997-2016
    .. _gflex: https://github.com/awickert/gFlex

    """

    h5file = folder+'/'+h5file+str(step)+'.hdf5'
    with h5py.File(h5file, "w") as f:

        # Write node coordinates and elevation
        f.create_dataset('coords',shape=(len(elevation),3), dtype='float64', compression='gzip')
        f["coords"][:,:2] = coords
        f["coords"][:,2] = elevation

        f.create_dataset('cells',shape=(len(cells[:,0]),3), dtype='int32', compression='gzip')
        f["cells"][:,:] = cells

        f.create_dataset('discharge',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        discharge[elevation<sealevel] = 0.
        discharge[discharge<1.] = 1.
        f["discharge"][:,0] = discharge

        if rainOn:
            f.create_dataset('precipitation',shape=(len(discharge), 1), dtype='float64', compression='gzip')
            f["precipitation"][:,0] = rain

        if eroOn:
            f.create_dataset('erodibility',shape=(len(discharge), 1), dtype='float64', compression='gzip')
            f["erodibility"][:,0] = erodibility

        f.create_dataset('cumdiff',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        f["cumdiff"][:,0] = cumdiff

        f.create_dataset('cumhill',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        f["cumhill"][:,0] = cumhill

        f.create_dataset('cumfail',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        f["cumfail"][:,0] = cumfail

        f.create_dataset('cumflex',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        f["cumflex"][:,0] = cumflex

        f.create_dataset('area',shape=(len(discharge), 1), dtype='float64', compression='gzip')
        f["area"][:,0] = area

        if waveOn:
            f.create_dataset('waveH',shape=(len(discharge), 1), dtype='float64', compression='gzip')
            f["waveH"][:,0] = waveH
            f.create_dataset('waveS',shape=(len(discharge), 1), dtype='float64', compression='gzip')
            f["waveS"][:,0] = waveS
            f.create_dataset('cumwave',shape=(len(discharge), 1), dtype='float64', compression='gzip')
            f["cumwave"][:,0] = wavediff

        if rockOn:
            nbSed = prop.shape[1]
            for k in range(nbSed):
                if k == 0:
                    f.create_dataset('depClastic',shape=(len(discharge), 1), dtype='float64', compression='gzip')
                    f["depClastic"][:,0] = prop[:,k]
                elif k == 1:
                    f.create_dataset('depSpecies1',shape=(len(discharge), 1), dtype='float64', compression='gzip')
                    f["depSpecies1"][:,0] = prop[:,k]
                elif k == 2:
                    f.create_dataset('depSpecies2',shape=(len(discharge), 1), dtype='float64', compression='gzip')
                    f["depSpecies2"][:,0] = prop[:,k]

def _write_xdmf(folder, xdmffile, xmffile, step):
    """
    This function writes the **XDmF** file which is calling the **XmF** file.

    Args:
        folder: name of the output folder.
        xdmffile: XDmF file name.
        xmffile: first part of the XmF file name.
        step: output visualisation step.
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

def write_xmf(folder, xmffile, xdmffile, step, t, elems, nodes, h5file, sealevel,
              flexOn, rainOn, eroOn, waveOn, rockOn, nbSed):
    """
    This function writes the **XmF** file which is calling each **hdf5** file.

    Args:
        folder: name of the output folder.
        xmffile: first part of the XmF file name.
        xdmffile: XDmF file name.
        step: output visualisation step.
        t: simulation time.
        elems: numpy integer-type array containing the number of cells of each local partition.
        nodes: numpy integer-type array containing the number of nodes of each local partition.
        h5file: first part of the hdf5 file name.
        sealevel: level of sea water
        flexOn: boolean for flexural isostasy.
        rainOn: boolean for orographic precipitation.
        eroOn: boolean for erodibility values.
        waveOn: boolean for wave activity.
        rockOn: boolean for multiple rock types.
        nbSed: number of sediment types considered.
    """

    xmf_file = folder+'/'+xmffile+str(step)+'.xmf'
    f= open(str(xmf_file),'w')

    f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">\n')
    f.write('<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">\n')
    f.write(' <Domain>\n')
    f.write('    <Grid GridType="Collection" CollectionType="Spatial">\n')
    f.write('      <Time Type="Single" Value="%s"/>\n'%t)

    p = 0
    pfile = h5file+str(step)+'.hdf5'
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

    f.write('         <Attribute Type="Scalar" Center="Node" Name="EroDep failure">\n')
    f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
    f.write('Dimensions="%d 1">%s:/cumfail</DataItem>\n'%(nodes[p],pfile))
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

        f.write('         <Attribute Type="Scalar" Center="Node" Name="waveS">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/waveS</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

        f.write('         <Attribute Type="Scalar" Center="Node" Name="EroDep wave">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/cumwave</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

    if rockOn:
        f.write('         <Attribute Type="Scalar" Center="Node" Name="depClastic">\n')
        f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
        f.write('Dimensions="%d 1">%s:/depClastic</DataItem>\n'%(nodes[p],pfile))
        f.write('         </Attribute>\n')

        if nbSed > 1:
            f.write('         <Attribute Type="Scalar" Center="Node" Name="depSpecies1">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/depSpecies1</DataItem>\n'%(nodes[p],pfile))
            f.write('         </Attribute>\n')

        if nbSed > 2:
            f.write('         <Attribute Type="Scalar" Center="Node" Name="depSpecies2">\n')
            f.write('          <DataItem Format="HDF" NumberType="Float" Precision="4" ')
            f.write('Dimensions="%d 1">%s:/depSpecies2</DataItem>\n'%(nodes[p],pfile))
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

    return
