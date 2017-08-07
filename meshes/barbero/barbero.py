# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v7.6.0 with dump python functionality
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'/home/guido/codes/sputnik/meshes/barbero')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

geomObj_1 = geompy.MakeVertex(0, 0, 0)
geomObj_2 = geompy.MakeVertex(2, 0, 0)
geomObj_3 = geompy.MakeVertex(1.414213562373095, 1.414213562373095, 0)
geomObj_4 = geompy.MakeVertex(3, 0, 0)
geomObj_5 = geompy.MakeVertex(3, 3, 0)
geomObj_6 = geompy.MakeVertex(0, 0, 1)
Arc01 = geompy.MakeArcCenter(geomObj_1, geomObj_2, geomObj_3,False)
Line01 = geompy.MakeLineTwoPnt(geomObj_2, geomObj_4)
Line02 = geompy.MakeLineTwoPnt(geomObj_4, geomObj_5)
Line03 = geompy.MakeLineTwoPnt(geomObj_5, geomObj_3)
Line04 = geompy.MakeLineTwoPnt(geomObj_1, geomObj_6)
Wire01 = geompy.MakeWire([Arc01, Line01, Line02, Line03], 1e-07)
Face01 = geompy.MakeFaceWires([Wire01], 1)
geomObj_7 = geompy.MakeFaceHW(3, 3, 1)
geomObj_8 = geompy.MakeTranslation(geomObj_7, 4.5, 1.5, 0)
Mirror_1 = geompy.MakeMirrorByAxis(Face01, Line03)
Shell_3 = geompy.MakeShell([Face01, Mirror_1, geomObj_8])
geomObj_9 = geompy.MakeMirrorByAxis(Shell_3, Line04)
geomObj_10 = geompy.MakeTranslation(geomObj_9, 6, 6, 0)
Shell_4 = geompy.MakeShell([Shell_3, geomObj_10])
Extrusion_1 = geompy.MakePrismVecH(Shell_4, Line04, 1)
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
X0 = geompy.CreateGroup(Extrusion_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(X0, [141, 107])
X1 = geompy.CreateGroup(Extrusion_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(X1, [52, 124])
Y0 = geompy.CreateGroup(Extrusion_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Y0, [83, 117])
Y1 = geompy.CreateGroup(Extrusion_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Y1, [28, 134])
Z0 = geompy.CreateGroup(Extrusion_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Z0, [56, 111, 87, 32, 145, 128])
Z1 = geompy.CreateGroup(Extrusion_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(Z1, [113, 147, 34, 130, 89, 58])
IRON = geompy.CreateGroup(Extrusion_1, geompy.ShapeType["SOLID"])
geompy.UnionIDs(IRON, [60, 91, 132, 2, 36, 115])
geompy.addToStudy( Arc01, 'Arc01' )
geompy.addToStudy( Line01, 'Line01' )
geompy.addToStudy( Line02, 'Line02' )
geompy.addToStudy( Line03, 'Line03' )
geompy.addToStudy( Line04, 'Line04' )
geompy.addToStudy( Wire01, 'Wire01' )
geompy.addToStudy( Face01, 'Face01' )
geompy.addToStudy( Mirror_1, 'Mirror_1' )
geompy.addToStudy( Shell_3, 'Shell_3' )
geompy.addToStudy( Shell_4, 'Shell_4' )
geompy.addToStudy( Extrusion_1, 'Extrusion_1' )
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

geompy.addToStudyInFather( Extrusion_1, X0, 'X0' )
geompy.addToStudyInFather( Extrusion_1, X1, 'X1' )
geompy.addToStudyInFather( Extrusion_1, Y0, 'Y0' )
geompy.addToStudyInFather( Extrusion_1, Y1, 'Y1' )
geompy.addToStudyInFather( Extrusion_1, Z0, 'Z0' )
geompy.addToStudyInFather( Extrusion_1, Z1, 'Z1' )
geompy.addToStudyInFather( Extrusion_1, IRON, 'IRON' )

###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Mesh01 = smesh.Mesh(Extrusion_1)
Regular_1D = Mesh01.Segment()
Nb_Segments_1 = Regular_1D.NumberOfSegments(8)
Nb_Segments_1.SetDistrType( 0 )
Quadrangle_2D = Mesh01.Quadrangle(algo=smeshBuilder.QUADRANGLE)
Hexa_3D = Mesh01.Hexahedron(algo=smeshBuilder.Hexa)
isDone = Mesh01.Compute()
Mesh01.ExportCGNS( r'/home/guido/codes/sputnik/meshes/barbero/papu.cgns', 1, Mesh01)
Mesh01.ExportDAT( r'/home/guido/codes/sputnik/meshes/barbero/papu01.dat' )
X1_1 = Mesh01.GroupOnGeom(X1,'X1',SMESH.FACE)
X1_1.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
X0_1 = Mesh01.GroupOnGeom(X0,'X0',SMESH.FACE)
Y0_1 = Mesh01.GroupOnGeom(Y0,'Y0',SMESH.FACE)
Y1_1 = Mesh01.GroupOnGeom(Y1,'Y1',SMESH.FACE)
Z0_1 = Mesh01.GroupOnGeom(Z0,'Z0',SMESH.FACE)
Z1_1 = Mesh01.GroupOnGeom(Z1,'Z1',SMESH.FACE)
IRON_1 = Mesh01.GroupOnGeom(IRON,'IRON',SMESH.VOLUME)
IRON_1.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
Mesh01.ExportDAT( r'/home/guido/codes/sputnik/meshes/barbero/papu01.dat' )


## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Hexa_3D.GetAlgorithm(), 'Hexa_3D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(Nb_Segments_1, 'Nb. Segments_1')
smesh.SetName(X0_1, 'X0')
smesh.SetName(X1_1, 'X1')
smesh.SetName(Y0_1, 'Y0')
smesh.SetName(Y1_1, 'Y1')
smesh.SetName(Z0_1, 'Z0')
smesh.SetName(Z1_1, 'Z1')
smesh.SetName(Mesh01.GetMesh(), 'Mesh01')
smesh.SetName(IRON_1, 'IRON')


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
