from paraview.simple import *
import numpy as np
import scipy as sc
import scipy
import vtk
import vtk2numpy as vn
import os
import sys

num_files = 6
displ_x = np.zeros( (num_files,1) )
displ_y = np.zeros( (num_files,1) )
force_x = np.zeros( (num_files,1) )
force_y = np.zeros( (num_files,1) )
file_matrix = np.zeros((num_files,4))

for i in range(0,num_files):

  file_in = "force_x/macro_t_"+str(i)+".pvtu"
  print "check for file "+file_in,
  try:
     f = open(file_in,'rb')
     print " ok"
  except IOError:
    print " error, not found"
    sys.exit()

  PVTU  = XMLPartitionedUnstructuredGridReader(FileName=file_in)
  KEYs  = PVTU.GetPointDataInformation().keys()
  KEYs  = PVTU.GetCellDataInformation().keys()

  PlotOverLine1 = PlotOverLine(Input=PVTU, guiName="PlotOverLine1", Source="High Resolution Line Source")
  PlotOverLine1.Source.Point1 = [30.0, 00.0, 0.0]
  PlotOverLine1.Source.Point2 = [30.0, 30.0, 0.0]
  PlotOverLine1.Source.Resolution = 1000
  PlotOverLine1.UpdatePipeline()

  stress_pvtu = vn.getPtsData(PlotOverLine1, "stress")
  displ_pvtu = vn.getPtsData(PlotOverLine1, "displ")
  leng_pvtu = vn.getPtsData(PlotOverLine1, "arc_length")

  displ_x[i] = displ_pvtu[0,0]
  force_x[i] = np.trapz(stress_pvtu[:,0],leng_pvtu)

for i in range(0,num_files):

  file_in = "force_y/macro_t_"+str(i)+".pvtu"
  print "check for file "+file_in,
  try:
     f = open(file_in,'rb')
     print " ok"
  except IOError:
    print " error, not found"
    sys.exit()

  PVTU  = XMLPartitionedUnstructuredGridReader(FileName=file_in)
  KEYs  = PVTU.GetPointDataInformation().keys()
  KEYs  = PVTU.GetCellDataInformation().keys()

  PlotOverLine1 = PlotOverLine(Input=PVTU, guiName="PlotOverLine1", Source="High Resolution Line Source")
  PlotOverLine1.Source.Point1 = [30.0, 00.0, 0.0]
  PlotOverLine1.Source.Point2 = [30.0, 30.0, 0.0]
  PlotOverLine1.Source.Resolution = 1000
  PlotOverLine1.UpdatePipeline()

  stress_pvtu = vn.getPtsData(PlotOverLine1, "stress")
  displ_pvtu = vn.getPtsData(PlotOverLine1, "displ")
  leng_pvtu = vn.getPtsData(PlotOverLine1, "arc_length")

  displ_y[i] = displ_pvtu[0,1]
  force_y[i] = np.trapz(stress_pvtu[:,2],leng_pvtu)

for i in range(0,num_files):
  file_matrix[i] = [ displ_x[i], displ_y[i], force_x[i], abs(force_y[i]) ]

np.savetxt("direct.dat", file_matrix)
