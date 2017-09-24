from paraview.simple import *
import numpy as np
import scipy as sc
import scipy
import os
import vtk
import sys
import vtk2numpy as vn

#----------------------------------------------------------------# 

aux = np.zeros( (2,4) )

fin   = "run_1/direct/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
KEYs  = PVTU.GetPointDataInformation().keys()
KEYs  = PVTU.GetCellDataInformation().keys()

PlotOverLine1 = PlotOverLine( Input=PVTU, guiName="PlotOverLine1", Source="High Resolution Line Source" )
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

fin   = "run_1/unifst/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)

PlotOverLine2 = PlotOverLine( Input=PVTU, guiName="PlotOverLine2", Source="High Resolution Line Source" )
PlotOverLine2.Source.Point1 = [30.0, 00.0, 0.0]
PlotOverLine2.Source.Point2 = [30.0, 30.0, 0.0]
PlotOverLine2.Source.Resolution = 1000
PlotOverLine2.UpdatePipeline()

key   = "stress"
stress = vn.getPtsData( PlotOverLine2, key )

key   = "displ"
Displ = vn.getPtsData( PlotOverLine2, key )

key = "arc_length"
leng = vn.getPtsData( PlotOverLine2, key )

fx_unifst = np.trapz(stress[:,0],leng)

#----------------------------------------------------------------# 

fin   = "run_1/taylor/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
Times = np.array(PVTU.TimestepValues)

PlotOverLine3 = PlotOverLine( Input=PVTU, guiName="PlotOverLine3", Source="High Resolution Line Source" )
PlotOverLine3.Source.Point1 = [30.0, 00.0, 0.0]
PlotOverLine3.Source.Point2 = [30.0, 30.0, 0.0]
PlotOverLine3.Source.Resolution = 1000
PlotOverLine3.UpdatePipeline()

key   = "stress"
stress = vn.getPtsData( PlotOverLine3, key )

key   = "displ"
Displ = vn.getPtsData( PlotOverLine3, key )

key = "arc_length"
leng = vn.getPtsData( PlotOverLine3, key )

fx_taylor = np.trapz(stress[:,0],leng)

aux[0,0] = 0.0      ; aux[0,1] = 0.0      ; aux[0,2] = 0.0      ; aux[0,3] = 0.0
aux[1,0] = disp0    ; aux[1,1] = fx_direct; aux[1,2] = fx_unifst; aux[1,3] = fx_taylor

np.savetxt("run_1/run_1.dat", aux )

print "k_direct = ", fx_direct/disp0       
print "k_unifst = ", fx_unifst/disp0 
print "k_taylor = ", fx_taylor/disp0       

