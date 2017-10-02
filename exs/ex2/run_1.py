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

fin   = "run_1/taylor_s/macro_t_1.pvtu"
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

fx_taylor_s = np.trapz(stress[:,0],leng)

#----------------------------------------------------------------# 

fin   = "run_1/taylor_p/macro_t_1.pvtu"
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

fx_taylor_p = np.trapz(stress[:,0],leng)

#----------------------------------------------------------------# 

fin   = "run_1/unifst/macro_t_1.pvtu"
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

fx_unifst = np.trapz(stress[:,0],leng)

aux = np.zeros( (2,6) )

aux[0,0] = 0.0  ;aux[0,1] = 0.0      ;aux[0,2] = 0.0        ;aux[0,3] = 0.0 ;
aux[1,0] = disp0;aux[1,1] = fx_direct;aux[1,2] = fx_unifst  ;aux[1,3] = 0.0 ; 

aux[0,4] = 0.0        ;aux[0,5] = 0.0        ;
aux[1,4] = fx_taylor_s;aux[1,5] = fx_taylor_p;

np.savetxt("run_1/run_1.dat", aux )

kd  = fx_direct  /disp0  
ku1 = fx_unifst/disp0
kt1 = fx_taylor_s/disp0  
kt2 = fx_taylor_p/disp0  

print "force_direct     = ", fx_direct  
print "force_unifst     = ", fx_unifst
print "force_taylor_s   = ", fx_taylor_s     
print "force_taylor_p   = ", fx_taylor_p     

print "k_direct     = ", kd  ,"error =",(kd - kd)/kd *100,"%"
print "k_unifst     = ", ku1 ,"error =",(kd - ku1)/kd *100,"%"
print "k_taylor_s   = ", kt1 ,"error =",(kd - kt1)/kd *100,"%"     
print "k_taylor_p   = ", kt2 ,"error =",(kd - kt2)/kd *100,"%"     

