from paraview.simple import *
import numpy as np
import scipy as sc
import scipy
import os
import vtk
import sys

#----------------------------------------------------------------# 
def getPtsData( _Obj, _Prop=None, _fname=None, _T=-1.0):
    from paraview.numpy_support import vtk_to_numpy

    if(_Prop==None):
      print _Obj.GetPointDataInformation().keys()
      sys.exit()

    if _T==None: 
      _Obj.SMProxy.UpdatePipeline(    )
    else:
      _Obj.SMProxy.UpdatePipeline( _T )
    _Obj.UpdatePipelineInformation()

    B = servermanager.Fetch( _Obj )
    if( B.GetClassName()== 'vtkMultiBlockDataSet' ): B = extract_block(B)[0]
    if( B.GetClassName()== 'vtkMultiBlockDataSet' ): B = extract_block(B)[0]

    Data = vtk_to_numpy( B.GetPointData().GetArray(_Prop) )

    if not _fname==None:
      time = B.GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
      print "  |_Time: %f" % ( time )

      nout = "%s_%s" % (_fname, _Prop)
      Fout = open(nout+".dat", "w")

      Fout.close()

      print "  |_'%s' " % ( Fout.name )

    return Data

#----------------------------------------------------------------# 
def getPtsCoord( _Obj, _T ):
    from paraview.numpy_support import vtk_to_numpy

    _Obj.SMProxy.UpdatePipeline( _T )
    _Obj.UpdatePipelineInformation()

    GetOutput = servermanager.Fetch( _Obj )
    if( GetOutput.GetClassName()== 'vtkMultiBlockDataSet' ): GetOutput = extract_block(GetOutput)[0]
    if( GetOutput.GetClassName()== 'vtkMultiBlockDataSet' ): GetOutput = extract_block(GetOutput)[0]

    vtk_n_pts = GetOutput.GetNumberOfPoints()
    Pts       = []
    for idx in range(vtk_n_pts):
      pt    = GetOutput.GetPoint(idx)
      n_pt  = len(pt)
      coord = []
      for j in range(n_pt): coord.append( pt[j] )
      Pts.append( coord )

    return np.array(Pts)
#----------------------------------------------------------------# 
def conv_macro():

    mesh_size=[5, 10, 15, 20, 25, 30, 35, 40]
    integrals = np.zeros( (len(mesh_size),2) )
    
    for method in methods:
      n = 0
      for i in mesh_size:
         print 'case ', i,'x',i
         os.system('m4 -DNx_m4_m4='+`i`+' -Dhomo_m4='+`method`+' run_multiscale_analysis.sh.m4 > run_multiscale_analysis.sh')
         os.system("chmod 755 run_multiscale_analysis.sh ")
         os.system("./run_multiscale_analysis.sh > macro.out")
         fin   = "macro_t_1.pvtu"
         PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
         Times = np.array(PVTU.TimestepValues)
      
         PlotOverLine3 = PlotOverLine( Input=PVTU, guiName="PlotOverLine3", Source="High Resolution Line Source" )
         PlotOverLine3.Source.Point2 = [30.0, 15.0, 0.0]
         PlotOverLine3.Source.Point1 = [0.0 , 15.0, 0.0]
         PlotOverLine3.Source.Resolution = 1000
         PlotOverLine3.UpdatePipeline()
      
         key   = "energy"
         Ene   = getPtsData( PlotOverLine3, key )
      
         key = "arc_length"
         leng = getPtsData( PlotOverLine3, key )
      
         Ene_int = np.trapz(Ene,leng)
         integrals[n,0] = i
         integrals[n,1] = Ene_int
         n += 1
    
      np.savetxt('integrals_'+method+'.dat', integrals )

#----------------------------------------------------------------# 
def conv_rve():

    rve_size=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    integrals = np.zeros( (len(rve_size),2) )
     
    n = 0
    for i in rve_size:
      print 'case ', i
      os.system('m4 -DNx_m4_m4=11 -Dhomo_m4=homo_unif_strains run_multiscale_analysis.sh.m4 > run_multiscale_analysis.sh')
      os.system("chmod 755 run_multiscale_analysis.sh ")
      os.system('m4 -DNint_m4='+`i`+' ../../meshes/cube_fiber/cube_fiber_2d_analysis.geo.m4 > cube_fiber_2d_analysis.geo')
      os.system('gmsh -2 cube_fiber_2d_analysis.geo')
      os.system('mv cube_fiber_2d_analysis.msh cube_fiber_2d.msh')

      os.system("./run_multiscale_analysis.sh > macro.out")
      fin   = "macro_t_1.pvtu"
      PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
      Times = np.array(PVTU.TimestepValues)
  
      PlotOverLine3 = PlotOverLine( Input=PVTU, guiName="PlotOverLine3", Source="High Resolution Line Source" )
      PlotOverLine3.Source.Point2 = [30.0, 15.0, 0.0]
      PlotOverLine3.Source.Point1 = [0.0 , 15.0, 0.0]
      PlotOverLine3.Source.Resolution = 1000
      PlotOverLine3.UpdatePipeline()
  
      key   = "energy"
      Ene   = getPtsData( PlotOverLine3, key )
  
      key = "arc_length"
      leng = getPtsData( PlotOverLine3, key )
  
      Ene_int = np.trapz(Ene,leng)
      integrals[n,0] = i
      integrals[n,1] = Ene_int
      n += 1
  
    np.savetxt('integrals_rve_convergence.dat', integrals )

#----------------------------------------------------------------# 
def rve_shape_analysis():

    rve_size=[1, 2, 3, 4, 5, 6, 7, 8]
    integrals = np.zeros( (len(rve_size),2) )
     
    n = 0
    for i in rve_size:
      print 'case ', i
      os.system('m4 -DNx_m4_m4=21 -Dhomo_m4=homo_unif_strains run_multiscale_analysis.sh.m4 > run_multiscale_analysis.sh')
      os.system("chmod 755 run_multiscale_analysis.sh ")
      os.system('cp ../../meshes/cube_fiber/struct_fiber_2d_'+`i`+'_'+`i`+'.msh cube_fiber_2d.msh')

      os.system('./run_multiscale_analysis.sh > macro_'+`i`+'.out')
      fin   = "macro_t_1.pvtu"
      PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
      Times = np.array(PVTU.TimestepValues)
  
      PlotOverLine3 = PlotOverLine( Input=PVTU, guiName="PlotOverLine3", Source="High Resolution Line Source" )
      PlotOverLine3.Source.Point2 = [30.0, 15.0, 0.0]
      PlotOverLine3.Source.Point1 = [0.0 , 15.0, 0.0]
      PlotOverLine3.Source.Resolution = 1000
      PlotOverLine3.UpdatePipeline()
  
      key   = "energy"
      Ene   = getPtsData( PlotOverLine3, key )
  
      key = "arc_length"
      leng = getPtsData( PlotOverLine3, key )
  
      aux = np.zeros( (leng.shape[0],2) )
      aux[:,1] = Ene[:]
      aux[:,0] = leng 

      np.savetxt('rve_shape_analysis_'+`i`+'.dat', aux )
      n += 1

#----------------------------------------------------------------# 
def rve_shape_analysis_1():

    rve_size=[1,2,3]
    integrals = np.zeros( (len(rve_size),2) )
     
    n = 0
    for i in rve_size:
      print 'case ', i
      os.system('m4 -Dlc_m4_m4='+`i*1.5`+' -DN_m4_m4='+`i*71`+' -Dnx_m4='+`i`+' run_multiscale_analysis_struc_micro.sh.m4 > run_multiscale_analysis.sh')
      os.system("chmod 755 run_multiscale_analysis.sh ")

      os.system('./run_multiscale_analysis.sh > macro_'+`i`+'.out')
      fin   = "macro_t_1.pvtu"
      PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
      Times = np.array(PVTU.TimestepValues)
  
      PlotOverLine3 = PlotOverLine( Input=PVTU, guiName="PlotOverLine3", Source="High Resolution Line Source" )
      PlotOverLine3.Source.Point2 = [30.0, 15.0, 0.0]
      PlotOverLine3.Source.Point1 = [0.0 , 15.0, 0.0]
      PlotOverLine3.Source.Resolution = 1000
      PlotOverLine3.UpdatePipeline()
  
      key   = "energy"
      Ene   = getPtsData( PlotOverLine3, key )
  
      key = "arc_length"
      leng = getPtsData( PlotOverLine3, key )
  
      aux = np.zeros( (leng.shape[0],2) )
      aux[:,1] = Ene[:]
      aux[:,0] = leng 

      np.savetxt('rve_shape_analysis_1_'+`i`+'.dat', aux )
      n += 1
  
#----------------------------------------------------------------# 
fin   = "direct/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
KEYs  = PVTU.GetPointDataInformation().keys()
print KEYs
KEYs  = PVTU.GetCellDataInformation().keys()
print KEYs
Times = np.array(PVTU.TimestepValues)

PlotOverLine1 = PlotOverLine( Input=PVTU, guiName="PlotOverLine1", Source="High Resolution Line Source" )
PlotOverLine1.Source.Point2 = [30.0, 15.0, 0.0]
PlotOverLine1.Source.Point1 = [0.0 , 15.0, 0.0]
PlotOverLine1.Source.Resolution = 1000
PlotOverLine1.UpdatePipeline()

key   = "residual"
Res   = getPtsData( PlotOverLine1, key )

key   = "displ"
Displ = getPtsData( PlotOverLine1, key )

key     = "energy_interp"
Ene_i   = getPtsData( PlotOverLine1, key )

key   = "energy"
Ene   = getPtsData( PlotOverLine1, key )

key = "arc_length"
leng = getPtsData( PlotOverLine1, key )

aux = np.zeros( (leng.shape[0],5) )
aux[:,4] = Ene_i[:]
aux[:,3] = Ene[:]
aux[:,2] = Res[:,0]
aux[:,1] = Displ[:,0]
aux[:,0] = leng 

np.savetxt("direct.dat", aux )
Ene_int_direct = np.trapz(Ene_i,leng)

fin   = "unif_strains/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
Times = np.array(PVTU.TimestepValues)

PlotOverLine2 = PlotOverLine( Input=PVTU, guiName="PlotOverLine2", Source="High Resolution Line Source" )
PlotOverLine2.Source.Point2 = [30.0, 15.0, 0.0]
PlotOverLine2.Source.Point1 = [0.0 , 15.0, 0.0]
PlotOverLine2.Source.Resolution = 1000
PlotOverLine2.UpdatePipeline()

key   = "residual"
Res   = getPtsData( PlotOverLine2, key )

key   = "displ"
Displ = getPtsData( PlotOverLine2, key )

key   = "energy"
Ene   = getPtsData( PlotOverLine2, key )

key = "arc_length"
leng = getPtsData( PlotOverLine2, key )

aux = np.zeros( (leng.shape[0],4) )
aux[:,3] = Ene[:]
aux[:,2] = Res[:,0]
aux[:,1] = Displ[:,0]
aux[:,0] = leng 

np.savetxt("unif_strains.dat", aux )
Ene_int_unif_strains = np.trapz(Ene,leng)

fin   = "taylor/macro_t_1.pvtu"
PVTU  = XMLPartitionedUnstructuredGridReader(FileName=fin)
Times = np.array(PVTU.TimestepValues)

PlotOverLine3 = PlotOverLine( Input=PVTU, guiName="PlotOverLine3", Source="High Resolution Line Source" )
PlotOverLine3.Source.Point2 = [30.0, 15.0, 0.0]
PlotOverLine3.Source.Point1 = [0.0 , 15.0, 0.0]
PlotOverLine3.Source.Resolution = 1000
PlotOverLine3.UpdatePipeline()

key   = "residual"
Res   = getPtsData( PlotOverLine3, key )

key   = "displ"
Displ = getPtsData( PlotOverLine3, key )

key   = "energy"
Ene   = getPtsData( PlotOverLine3, key )

key = "arc_length"
leng = getPtsData( PlotOverLine3, key )

aux = np.zeros( (leng.shape[0],4) )
aux[:,3] = Ene[:]
aux[:,2] = Res[:,0]
aux[:,1] = Displ[:,0]
aux[:,0] = leng 

np.savetxt("taylor.dat", aux )

Ene_int_taylor = np.trapz(Ene,leng)
print "Ene_int_taylor       = ", Ene_int_taylor
print "Ene_int_unif_strains = ", Ene_int_unif_strains
print "Ene_int_taylor       = ", Ene_int_direct


#conv_macro()
#conv_rve()
#rve_shape_analysis()
rve_shape_analysis_1()
