from paraview.simple import *
import numpy as np
import scipy as sc
import scipy
import os
import vtk
import sys
import vtk2numpy_m as vn


#----------------------------------------------------------------# 

fin   = "direct/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
KEYs  = PVTU.GetPointDataInformation().keys()
print KEYs
KEYs  = PVTU.GetCellDataInformation().keys()
print KEYs
Times = np.array(PVTU.TimestepValues)

PlotOverLine1 = PlotOverLine( Input=PVTU, guiName="PlotOverLine1", Source="High Resolution Line Source" )
PlotOverLine1.Source.Point2 = [30.0, 30.0, 0.0]
PlotOverLine1.Source.Point1 = [30.0, 00.0, 0.0]
PlotOverLine1.Source.Resolution = 1000
PlotOverLine1.UpdatePipeline()

key   = "stress"
stress = vn.getPtsData( PlotOverLine1, key )

#key   = "residual"
#Res   = getPtsData( PlotOverLine1, key )
#
#key   = "displ"
#Displ = getPtsData( PlotOverLine1, key )
#
#key   = "energy_interp"
#Ene_i = getPtsData( PlotOverLine1, key )
#
#key   = "energy"
#Ene   = getPtsData( PlotOverLine1, key )
#
#key = "arc_length"
#leng = getPtsData( PlotOverLine1, key )
#
#aux = np.zeros( (leng.shape[0],5) )
#aux[:,4] = Ene_i[:]
#aux[:,3] = Ene[:]
#aux[:,2] = stress[:,0]
#aux[:,1] = Displ[:,0]
#aux[:,0] = leng 
#
#np.savetxt("direct.dat", aux )
#Ene_int_direct = np.trapz(Ene_i,leng)
#force_int_direct = np.trapz(stress[:,0],leng)
#
#fin   = "unif_strains/macro_t_1.pvtu"
#PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
#Times = np.array(PVTU.TimestepValues)
#
#PlotOverLine2 = PlotOverLine( Input=PVTU, guiName="PlotOverLine2", Source="High Resolution Line Source" )
#PlotOverLine2.Source.Point2 = [30.0, 30.0, 0.0]
#PlotOverLine2.Source.Point1 = [30.0 ,0.0, 0.0]
#PlotOverLine2.Source.Resolution = 1000
#PlotOverLine2.UpdatePipeline()
#
#key   = "stress"
#stress = getPtsData( PlotOverLine2, key )
#
#key   = "residual"
#Res   = getPtsData( PlotOverLine2, key )
#
#key   = "displ"
#Displ = getPtsData( PlotOverLine2, key )
#
#key   = "energy"
#Ene   = getPtsData( PlotOverLine2, key )
#
#key = "arc_length"
#leng = getPtsData( PlotOverLine2, key )
#
#aux = np.zeros( (leng.shape[0],4) )
#aux[:,3] = Ene[:]
#aux[:,2] = stress[:,0]
#aux[:,1] = Displ[:,0]
#aux[:,0] = leng 
#
#np.savetxt("unif_strains.dat", aux )
#Ene_int_unif_strains = np.trapz(Ene,leng)
#force_int_unif_strains = np.trapz(stress[:,0],leng)
#
#fin   = "taylor/macro_t_1.pvtu"
#PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
#Times = np.array(PVTU.TimestepValues)
#
#PlotOverLine3 = PlotOverLine( Input=PVTU, guiName="PlotOverLine3", Source="High Resolution Line Source" )
#PlotOverLine3.Source.Point2 = [30.0, 30.0, 0.0]
#PlotOverLine3.Source.Point1 = [30.0 ,0.0, 0.0]
#PlotOverLine3.Source.Resolution = 1000
#PlotOverLine3.UpdatePipeline()
#
#key   = "stress"
#stress = getPtsData( PlotOverLine3, key )
#
#key   = "residual"
#Res   = getPtsData( PlotOverLine3, key )
#
#key   = "displ"
#Displ = getPtsData( PlotOverLine3, key )
#
#key   = "energy"
#Ene   = getPtsData( PlotOverLine3, key )
#
#key = "arc_length"
#leng = getPtsData( PlotOverLine3, key )
#
#aux = np.zeros( (leng.shape[0],4) )
#aux[:,3] = Ene[:]
#aux[:,2] = stress[:,0]
#aux[:,1] = Displ[:,0]
#aux[:,0] = leng 
#
#np.savetxt("taylor.dat", aux )
#
#force_int_taylor = np.trapz(stress[:,0],leng)
#print "force_int_taylor       = ", force_int_taylor       ,"K_taylor       = ", force_int_taylor       / Displ[0,0]
#print "force_int_unif_strains = ", force_int_unif_strains ,"K_unif_strains = ", force_int_unif_strains / Displ[0,0]
#print "force_int_taylor       = ", force_int_direct       ,"K_taylor       = ", force_int_direct       / Displ[0,0]
#
