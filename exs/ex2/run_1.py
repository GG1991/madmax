from paraview.simple import *
import numpy as np
import scipy as sc
import scipy
import os
import vtk
import sys
import vtk2numpy as vn

#----------------------------------------------------------------# 

fin   = "run_1/direct/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
KEYs  = PVTU.GetPointDataInformation().keys()
KEYs  = PVTU.GetCellDataInformation().keys()

PlotOverLine1 = PlotOverLine( Input=PVTU, guiName="PlotOverLine", Source="High Resolution Line Source" )
PlotOverLine1.Source.Point1 = [30.0, 00.0, 0.0]
PlotOverLine1.Source.Point2 = [30.0, 30.0, 0.0]
PlotOverLine1.Source.Resolution = 1000
PlotOverLine1.UpdatePipeline()

key    = "stress"
stress = vn.getPtsData( PlotOverLine1, key )

key    = "displ"
Displ  = vn.getPtsData( PlotOverLine1, key )

key    = "arc_length"
leng   = vn.getPtsData( PlotOverLine1, key )

disp0 = Displ[0,0]

fx_direct = np.trapz(stress[:,0],leng)

#----------------------------------------------------------------# 

fin   = "run_1/taylor_1/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
Times = np.array(PVTU.TimestepValues)

PlotOverLine1 = PlotOverLine( Input=PVTU, guiName="PlotOverLine", Source="High Resolution Line Source" )
PlotOverLine1.Source.Point1 = [30.0, 00.0, 0.0]
PlotOverLine1.Source.Point2 = [30.0, 30.0, 0.0]
PlotOverLine1.Source.Resolution = 1000
PlotOverLine1.UpdatePipeline()

key    = "stress"
stress = vn.getPtsData( PlotOverLine1, key )

key    = "displ"
Displ  = vn.getPtsData( PlotOverLine1, key )

key    = "arc_length"
leng   = vn.getPtsData( PlotOverLine1, key )

fx_taylor_1 = np.trapz(stress[:,0],leng)

#----------------------------------------------------------------# 

fin   = "run_1/taylor_2/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
Times = np.array(PVTU.TimestepValues)

PlotOverLine1 = PlotOverLine( Input=PVTU, guiName="PlotOverLine", Source="High Resolution Line Source" )
PlotOverLine1.Source.Point1 = [30.0, 00.0, 0.0]
PlotOverLine1.Source.Point2 = [30.0, 30.0, 0.0]
PlotOverLine1.Source.Resolution = 1000
PlotOverLine1.UpdatePipeline()

key    = "stress"
stress = vn.getPtsData( PlotOverLine1, key )

key    = "displ"
Displ  = vn.getPtsData( PlotOverLine1, key )

key    = "arc_length"
leng   = vn.getPtsData( PlotOverLine1, key )

fx_taylor_2 = np.trapz(stress[:,0],leng)

#----------------------------------------------------------------# 

fin   = "run_1/unifst_1/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)

PlotOverLine1 = PlotOverLine( Input=PVTU, guiName="PlotOverLine", Source="High Resolution Line Source" )
PlotOverLine1.Source.Point1 = [30.0, 00.0, 0.0]
PlotOverLine1.Source.Point2 = [30.0, 30.0, 0.0]
PlotOverLine1.Source.Resolution = 1000
PlotOverLine1.UpdatePipeline()

key    = "stress"
stress = vn.getPtsData( PlotOverLine1, key )

key    = "displ"
Displ  = vn.getPtsData( PlotOverLine1, key )

key    = "arc_length"
leng   = vn.getPtsData( PlotOverLine1, key )

fx_unifst_1 = np.trapz(stress[:,0],leng)

#----------------------------------------------------------------# 

fin   = "run_1/unifst_2/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)

PlotOverLine1 = PlotOverLine( Input=PVTU, guiName="PlotOverLine", Source="High Resolution Line Source" )
PlotOverLine1.Source.Point1 = [30.0, 00.0, 0.0]
PlotOverLine1.Source.Point2 = [30.0, 30.0, 0.0]
PlotOverLine1.Source.Resolution = 1000
PlotOverLine1.UpdatePipeline()

key    = "stress"
stress = vn.getPtsData( PlotOverLine1, key )

key    = "displ"
Displ  = vn.getPtsData( PlotOverLine1, key )

key    = "arc_length"
leng   = vn.getPtsData( PlotOverLine1, key )

fx_unifst_2 = np.trapz(stress[:,0],leng)

#----------------------------------------------------------------# 

aux = np.zeros( (2,6) )

aux[0,0] = 0.0  ;aux[0,1] = 0.0      ;aux[0,2] = 0.0        ;aux[0,3] = 0.0        ;
aux[1,0] = disp0;aux[1,1] = fx_direct;aux[1,2] = fx_unifst_1;aux[1,3] = fx_unifst_2;

aux[0,4] = 0.0        ;aux[0,5] = 0.0        ;
aux[1,4] = fx_taylor_1;aux[1,5] = fx_taylor_2;

np.savetxt("run_1/run_1.dat", aux )

kd  = fx_direct/disp0  
ku1 = fx_unifst_1/disp0
ku2 = fx_unifst_2/disp0
kt1 = fx_taylor_1/disp0  
kt2 = fx_taylor_2/disp0  

print "k_direct     = ", kd  ,"error =",(kd - kd)/kd *100,"%"
print "k_unifst_1   = ", ku1 ,"error =",(kd - ku1)/kd *100,"%"
print "k_unifst_2   = ", ku2 ,"error =",(kd - ku2)/kd *100,"%"   
print "k_taylor_1   = ", kt1 ,"error =",(kd - kt1)/kd *100,"%"     
print "k_taylor_2   = ", kt2 ,"error =",(kd - kt2)/kd *100,"%"     

