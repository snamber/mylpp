#!/usr/bin/python
"""
A simple routine to load in a LIGGGHTS hybrid dump file containing
contact and contact force data and convert into a .vtk unstructured
grid which can be used to visualise the force network.

Contributing author: Mark Bentley, Space Research Institute, 
Austrian Academy of Sciences, mark.bentley@oeaw.ac.at

evtk is used to write binary VTK files:
https://bitbucket.org/pauloh/pyevtk

The pizza.py bdump command is used to handle LIGGGHTS dump files and
therefore PYTHONPATH must include the pizza/src location.

"""

try: from evtk.vtk  import VtkFile, VtkGroup, VtkUnstructuredGrid
except: print "Need package evtk for dump2force, https://bitbucket.org/pauloh/pyevtk"

try: from bdump import bdump
except:  print "Need bdump moldule for dump2force"

try:
    import numpy as np
    oldnumeric = False
except:
    import Numeric as np
    oldnumeric = True

import sys, os


# TODO: use a try/except here to check for missing modules, and report
# TODO: ask for timestep or timestep range as input
# TODO  make it a class....
# TODO: check for periodic boundaries and avoid connecting points that map across them


# Check for command line arguments
if len(sys.argv) != 2:
        sys.exit('Usage: dump2forcenetwork.py <filename>, where filename is typically dump.<runname>')
        
elif len(sys.argv) == 2: # we have one input param, that should be parsed as a filename
    filename = str(sys.argv[1])
    if not os.path.isfile(filename):
        sys.exit('File ' + filename + ' does not exist!')

splitname = filename.split('.')

if len(splitname) == 2 and splitname[0].lower() == 'dump':
    fileprefix = splitname[1]
else:
    fileprefix = splitname[0]

inputpath = os.path.abspath(filename)
inputdir = os.path.split(inputpath)[0]

# create a sub-directory for the output .vtu files
outputdir = os.path.join(inputdir,fileprefix)
try:
    os.mkdir(outputdir)
except:
    pass

# Read in the dump file - since we can have many contacts (i.e. >> nparticles)
# and many timesteps I will deal with one timestep at a time in memory,
# write to the appropriate .vtu file for a single timestep, then move on.

forcedata = bdump(filename,0)

groupfile = fileprefix
groupfile = os.path.join(inputdir,groupfile)
groupfile = VtkGroup(groupfile)

fileindex = 0
timestep = forcedata.next()

# check that we have the right number of colums (11)
if forcedata.snaps[fileindex].natoms !=0 and len(forcedata.snaps[0].atoms[0]) != 11:
    print "Error - dump file requires all parameters from a compute pair/gran/local id pos force (11 in total)"
    sys.exit

# loop through available timesteps
# 
while timestep >= 0:

    # data are stored as pos1 (3) pos2 (3) id1 id2 force (3) for 11 columns
    # assign names to atom columns (1-N)
    forcedata.map(1,"x1",2,"y1",3,"z1",4,"x2",5,"y2",6,"z2",7,"id1",8,"id2",9,"fx",10,"fy",11,"fz")

    # check for contact data (some timesteps may have no particles in contact)
    #
    # NB. if one loads two datasets into ParaView with defined timesteps, but in which
    # one datasets has some missing, data for the previous timestep are still displayed - 
    # this means that it is better here to generate "empty" files for these timesteps.

    if forcedata.snaps[fileindex].natoms == 0:   
        vtufile = fileprefix+'_'+str(timestep)+'.vtu'
        vtufile = os.path.join(outputdir,vtufile)
        vtuwrite = file(vtufile,'w')
        vtuwrite.write("""<?xml version="1.0"?>
<VTKFile byte_order="LittleEndian" version="0.1" type="UnstructuredGrid">
<UnstructuredGrid>
        <Piece NumberOfCells="0" NumberOfPoints="0">
                <Cells>
                        <DataArray NumberOfComponents="1" offset="0" type="Int64" Name="connectivity"/>
                        <DataArray NumberOfComponents="1" offset="0" type="Int64" Name="offsets"/>
                        <DataArray NumberOfComponents="1" offset="0" type="Int64" Name="types"/>

                </Cells>
        </Piece>
</UnstructuredGrid>
</VTKFile>""")
        
    else:

        # number of cells = number of interactions (i.e. entries in the dump file)
        ncells = len(forcedata.snaps[fileindex].atoms)

        # extract the IDs as an array of integers
        id1 = np.array(forcedata.snaps[fileindex].atoms[:,forcedata.names["id1"]],dtype=long)
        id2 = np.array(forcedata.snaps[fileindex].atoms[:,forcedata.names["id2"]],dtype=long)

        # and convert to lists
        id1 = id1.tolist()
        id2 = id2.tolist()

        # concatenate into a single list
        ids = []
        ids = id1[:]
        ids.extend(id2)

        # convert to a set and back to remove duplicates, then sort
        ids = list(set(ids))
        ids.sort()

        # number of points = number of unique IDs (particles)
        npoints = len(ids)

        # create empty arrays to hold x,y,z data
        x = np.zeros( npoints, dtype=np.float64)
        y = np.zeros( npoints, dtype=np.float64)
        z = np.zeros( npoints, dtype=np.float64)

        print 'Timestep:',str(timestep),'npoints=',str(npoints),'ncells=',str(ncells)

        # Point data = location of each unique particle
        #
        # The order of this data is important since we use the position of each particle
        # in this list to reference particle connectivity! We will use the order of the 
        # sorted ids array to determine this.

        counter = 0   
        for id in ids:
            if id in id1:
                index = id1.index(id)
                xtemp,ytemp,ztemp = forcedata.snaps[fileindex].atoms[index,forcedata.names["x1"]], \
                        forcedata.snaps[fileindex].atoms[index,forcedata.names["y1"]], \
                        forcedata.snaps[fileindex].atoms[index,forcedata.names["z1"]]
            else:
                index = id2.index(id)
                xtemp,ytemp,ztemp = forcedata.snaps[fileindex].atoms[index,forcedata.names["x2"]], \
                        forcedata.snaps[fileindex].atoms[index,forcedata.names["y2"]], \
                        forcedata.snaps[fileindex].atoms[index,forcedata.names["z2"]]
            
            x[counter]=xtemp
            y[counter]=ytemp
            z[counter]=ztemp           
            counter += 1

        # Now create the connectivity list - this corresponds to pairs of IDs, but referencing
        # the order of the ids array, so now we loop through 0..ncells and have to connect 
        # id1 and id2, so I need to see where in ids these correspond to

        # create an empty array to hold particle pairs
        connections = np.zeros( 2*ncells, dtype=int )

        for pair in range(ncells):
            connections[2*pair],connections[2*pair+1] = ids.index(id1[pair]),ids.index(id2[pair])
            
        # The offset array is simply generated from 2*(1..ncells)
        offset=(np.arange(ncells,dtype=int)+1)*2

        # The type array is simply ncells x 3 (i.e. a VTKLine type)
        celltype = np.ones(ncells,dtype=int)*3

        # Finally we need force data for each cell
        force = np.sqrt( np.array(forcedata.snaps[fileindex].atoms[:,forcedata.names["fx"]],dtype=np.float64)**2 + \
                         np.array(forcedata.snaps[fileindex].atoms[:,forcedata.names["fy"]],dtype=np.float64)**2 + \
                         np.array(forcedata.snaps[fileindex].atoms[:,forcedata.names["fz"]],dtype=np.float64)**2 )

        # Now we have enough data to create the file:
        # Points - (x,y,z) (npoints)
        # Cells
        #   Connectivity - connections (ncells,2)
        #   Offset - offset (ncells)
        #   type - celltype (ncells)
        # Celldata - force (ncells)

        # create a VTK unstructured grid (.vtu) file
        vtufile = fileprefix+'_'+str(timestep)
        vtufile = os.path.join(outputdir,vtufile)
        w = VtkFile(vtufile, VtkUnstructuredGrid)
        vtufile += '.vtu'

        w.openGrid()
        w.openPiece(npoints=npoints, ncells=ncells)    

        # Set up Points data
        w.openElement("Points")
        w.addData("points", (x,y,z) )
        w.closeElement("Points")

        # Set up Cell data
        w.openElement("Cells")
        w.addData("connectivity", connections )

        w.addData("offsets", offset)
        w.addData("types", celltype)
        w.closeElement("Cells")

        # Set up force data
        w.openData("Cell", scalars = "force")
        w.addData("force", force)
        w.closeData("Cell")

        # Wrap up
        w.closePiece()
        w.closeGrid()

        # Append binary data
        w.appendData( (x,y,z) )
        w.appendData(connections).appendData(offset).appendData(celltype)
        w.appendData(data = force)
        w.save()

    # Add this file to the group of all timesteps
    groupfile.addFile(filepath = os.path.relpath(vtufile,inputdir), sim_time = timestep)

    fileindex += 1
    timestep = forcedata.next()

# end of main loop - close group file
groupfile.save()

