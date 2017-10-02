from paraview.simple import *
import numpy as np
import scipy as sc
import scipy
import os
import vtk
import sys
import vtk2numpy as vn


#----------------------------------------------------------------# 

fin   = 'run_1/direct/macro_t_1.pvtu'
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
KEYs  = PVTU.GetPointDataInformation().keys()
KEYs  = PVTU.GetCellDataInformation().keys()

PlotOverLine1 = PlotOverLine(Input=PVTU, guiName="PlotOverLine1", Source="High Resolution Line Source")
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

force = np.trapz(stress[:,0],leng)
print "force direct", force, disp0

kd    = np.trapz(stress[:,0],leng) /disp0

#----------------------------------------------------------------# 

fx = np.zeros( (19,3) )

for i in range(2,21):
    
    print "case ", i
    fin   = 'run_2/homog_'+`i`+'/macro_t_1.pvtu'
    PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
    KEYs  = PVTU.GetPointDataInformation().keys()
    KEYs  = PVTU.GetCellDataInformation().keys()
    
    PlotOverLine1 = PlotOverLine(Input=PVTU, guiName="PlotOverLine1", Source="High Resolution Line Source")
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

    force = np.trapz(stress[:,0],leng)
    print "force ", force, disp0
    
    fx[i-2,0] = i
    fx[i-2,1] = (np.trapz(stress[:,0],leng)/disp0 - kd )/ kd
    fx[i-2,2] = np.trapz(stress[:,0],leng)/disp0
	
np.savetxt("run_2/run_2.dat", fx)
