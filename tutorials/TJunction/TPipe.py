#http://www.salome-platform.org/forum/forum_10/377327504

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'/home/sg5516/OpenFOAM_new/sokratiag-4.1/run/uBend_buoyant')
print(sys.path)
###
### GEOM component
###
#
import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS
#
#
geompy = geomBuilder.New(theStudy)


import smesh

#import geompy, smesh, math

f=.1
raza_cil_mare=2.0*f
raza_cil_mic= 1.7*f
inalt_cil_mare= 15.0*f
inalt_cil_mic= 7.0*f
dim_cutie_mare= 1*f
#dim_cutie_mica= 1*f
angle = math.radians(60)
#dim_plan_mare= 30.0*f
#dim_plan_mic= inalt_cil_mic*2.0
#locatie_cil_mic= 3.0*f
#locatie_cil_mic_z= -raza_cil_mare
#locatie_cutie = dim_cutie_mare/2.0
#transl_cutie= -locatie_cutie

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

def Profile (dim_box):
    # make a hexagonal face
    Points = []
    for i in range (6):
        y = dim_box*math.cos(angle*i)
        z = dim_box*math.sin(angle*i)
        point = geompy.MakeVertex (0.0,y,z)
        Points.append(point)
    Wire_1 = geompy.MakePolyline( Points, True )
    Face_1 = geompy.MakeFaceWires([Wire_1], 1)
    return Face_1

# Make a half of a hexagonal horizontal prism + cut a half of its bases by 45 degree
# the prism
hex_profile   = Profile( dim_cutie_mare )
hex_tool_hor  = geompy.MakePrismVecH( hex_profile, OX, inalt_cil_mare )
# the cutting tool
cut_pln_norm  = geompy.MakeVectorDXDYDZ( -1, 0, 1 )
cut_plane     = geompy.MakePlane( O, cut_pln_norm, 10*f )
cut_box       = geompy.MakePrismVecH( cut_plane, cut_pln_norm, 10*f )
# cut
hex_tool_hor  = geompy.MakeCut( hex_tool_hor, cut_box )
geompy.addToStudy( hex_tool_hor, "hex_tool_hor")

# Make a half of a vertival hexagonal prism
hex_vert_base = geompy.GetInPlace( hex_tool_hor, cut_plane )
hex_tool_vert = geompy.MakePrismVecH( hex_vert_base, OZ, 10*f )
geompy.addToStudy( hex_tool_vert, "hex_tool_vert")

# get a Z upper coord of the hex_profile
hex_points    = geompy.SubShapeAllSorted( hex_profile, geompy.ShapeType["VERTEX"])
Z = 0
for p in hex_points:
    Z = max( Z, geompy.PointCoordinates( p )[2] )

# Make plate tools from the top and bottom faces of the horizontal hexagonal prism
p_up          = geompy.MakeVertex( 2*dim_cutie_mare, 0, Z )
p_down        = geompy.MakeVertex( 2*dim_cutie_mare, 0, -Z )
plate_top     = geompy.GetFaceNearPoint( hex_tool_hor, p_up )
plate_bottom  = geompy.GetFaceNearPoint( hex_tool_hor, p_down )
plate_tool_top= geompy.MakePrismVecH( plate_top,    OZ, 10*f )
plate_tool_bot= geompy.MakePrismVecH( plate_bottom, OZ, -10*f )
geompy.addToStudy( plate_tool_top, "plate_tool_top")
geompy.addToStudy( plate_tool_bot, "plate_tool_bot")

# 45 degree cutting plane passing through (0,0,0)
p1            = geompy.MakeVertex( 0, -f*2, 0 )
p2            = geompy.MakeVertex( 0, +f*2, 0 )
line          = geompy.MakeLineTwoPnt( p1, p2 )
cut_plane2    = geompy.MakePrismVecH( line, cut_pln_norm, 10*f )

# Make remaining tools for Partition
tools_right   = geompy.MakeCompound([ hex_tool_hor, hex_tool_vert,
                                      plate_tool_top, plate_tool_bot, cut_plane2 ])
tools_left    = geompy.MakeRotation( tools_right, OZ, math.pi )
xoz_plane     = geompy.MakePlane( O, OY, 100*f )
yoz_plane     = geompy.MakePlane( O, OX, 100*f )
xoy_plane     = geompy.MakePlane( O, OZ, 100*f )
geompy.addToStudy( tools_left, "tools_left")

# Two cilinders
cyl_hor_start = geompy.MakeVertex( -inalt_cil_mare* 0.5, 0, 0 )
cyl_hor       = geompy.MakeCylinder( cyl_hor_start, OX, raza_cil_mare, 1.5* inalt_cil_mare)
cyl_vert      = geompy.MakeCylinder( O, OZ, raza_cil_mic,  inalt_cil_mic)
geompy.addToStudy( cyl_hor, "cyl_hor_init")
geompy.addToStudy( cyl_vert, "cyl_vert_init")

# Make a surface separating the two cylinders
cyl_section   = geompy.MakeSection( cyl_hor, cyl_vert )
section_edges = geompy.SubShapeAll( cyl_section, geompy.ShapeType["EDGE"])
cyl_boundary  = geompy.MakePrismVecH2Ways( section_edges[0], OY, 10*f )
geompy.addToStudy( cyl_boundary, "cyl_boundary")

# Get needed parts of the cylinders
cyl_vert_part = geompy.MakePartition( [cyl_vert], [ cyl_boundary] )
cyl_vert      = geompy.GetShapesNearPoint( cyl_vert_part,
                                           geompy.MakeVertex( 0, 0, 3*f),
                                           geompy.ShapeType["SOLID"])
cyl_hor       = geompy.MakeCut( cyl_hor, cyl_vert )
geompy.addToStudy( cyl_hor, "cyl_hor")
geompy.addToStudy( cyl_vert, "cyl_vert")

# Final Partition
t_pipe = geompy.MakePartition( [ cyl_vert, cyl_hor ],
                               [ tools_left, tools_right, xoz_plane, yoz_plane, xoy_plane ])

geompy.addToStudy( t_pipe, "t_pipe")


# Meshing

mesh = smesh.Mesh( t_pipe )
algo1D = mesh.Segment()
algo1D.NumberOfSegments( 10 )
mesh.Quadrangle()
mesh.Hexahedron()

# set local nb of segments
segLen = f / 5.
minNbSeg = 3
oppEdgesList = geompy.Propagate( t_pipe )
for oppEdges in oppEdgesList:
    edges = geompy.SubShapeAll( oppEdges, geompy.ShapeType["EDGE"])
    avgLen = 0
    for e in edges:
        avgLen += geompy.BasicProperties( e )[0]
    avgLen /= len(edges)
    nbSeg = int (max( minNbSeg, avgLen / segLen ))
    #print nbSeg
    geompy.addToStudyInFather( t_pipe, oppEdges, "edges_nbSeg=%s"%nbSeg )
    algo1D = mesh.Segment( oppEdges )
    algo1D.NumberOfSegments( nbSeg, UseExisting=True )

mesh.Compute()
