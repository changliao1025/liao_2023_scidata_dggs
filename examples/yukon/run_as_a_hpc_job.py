
import os,  stat

from pathlib import Path
from os.path import realpath

from shutil import copy2
import numpy as np
import logging
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(format='%(asctime)s %(message)s')
logging.warning('is the time pyhexwatershed simulation started.')


from pyhexwatershed.pyhexwatershed_read_model_configuration_file import pyhexwatershed_read_model_configuration_file
from pyflowline.mesh.dggrid.create_dggrid_mesh import dggrid_find_resolution_by_index
from pyflowline.change_json_key_value import change_json_key_value

sMesh_type = 'dggrid'
iCase_index = 1
dResolution_meter=5000
iFlag_create_job = 0
iFlag_visualization = 1
aExtent_full = None

dLongitude_outlet_degree=-117
dLatitude_outlet_degree=42
pProjection_map = None
sDate='20240102'
sPath = str( Path().resolve() )
iFlag_option = 1
sWorkspace_data = realpath( sPath +  '/data/yukon' )
sWorkspace_input =  str(Path(sWorkspace_data)  /  'input')
sWorkspace_output=  '/compyfs/liao313/04model/pyhexwatershed/yukon'
if not os.path.exists(sWorkspace_output):
    os.makedirs(sWorkspace_output)


dLongitude_outlet_degree= -164.47594
dLatitude_outlet_degree= 63.04269

iMesh_type = 5

#set up dggrid resolution level 
aResolution_index = [10, 11, 12, 13, 14]
nCase = len(aResolution_index)
aCase_index = [10, 11, 12, 13, 14] #aResolution_index
sDggrid_type = 'ISEA3H'
#generate a bash job script
if iFlag_create_job ==1:
    sFilename_job = sWorkspace_output + '/' + sDate  + 'submit.bash'
    ofs = open(sFilename_job, 'w')
    sLine  = '#!/bin/bash' + '\n'
    ofs.write(sLine)

sFilename_configuration_in = realpath( sWorkspace_input +  '/pyhexwatershed_yukon_dggrid.json' )
sFilename_wbd_boundary = '/qfs/people/liao313/data/hexwatershed/yukon/vector/hydrology/boundary.geojson'
sFilename_dem_arctic = '/qfs/people/liao313/data/hexwatershed/yukon/raster/dem/hyd_ar_dem_30s/hyd_ar_dem_30s.tif'
    
if os.path.isfile(sFilename_configuration_in):
    print(sFilename_configuration_in)
else:
    print('This configuration file does not exist: ', sFilename_configuration_in )
    exit()
    
#mpas mesh only has one resolution
iFlag_stream_burning_topology = 1 
iFlag_use_mesh_dem = 0
iFlag_elevation_profile = 0

aExtent = [-60.6,-59.2  ,-3.6, -2.5]
aExtent = None

sFilename_basins_in= '/qfs/people/liao313/workspace/python/liao_2023_scidata_dggs/data/yukon/input/pyflowline_yukon_basins.json'
for iCase in range(4, 5, 1):    
    iResolution_index = aResolution_index[iCase]
    sResolution = "{:0d}".format(iResolution_index)
    iCase_index = aCase_index[iCase]
    dResolution = dggrid_find_resolution_by_index(sDggrid_type, iResolution_index)
    print(dResolution)   

    oPyhexwatershed = pyhexwatershed_read_model_configuration_file(sFilename_configuration_in,
                    iCase_index_in=iCase_index,iFlag_stream_burning_topology_in=iFlag_stream_burning_topology,
                    iFlag_use_mesh_dem_in=iFlag_use_mesh_dem,
                    iFlag_elevation_profile_in=iFlag_elevation_profile,
                    iResolution_index_in = iResolution_index, 
                    sDggrid_type_in=sDggrid_type,
                    sDate_in= sDate, sMesh_type_in= sMesh_type) 
    
    #copy to the output directory
    sWorkspace_output = oPyhexwatershed.sWorkspace_output
    sFilename_configuration_copy= os.path.join( sWorkspace_output, 'pyhexwatershed_configuration_copy.json' )
    #copy the configuration file to the output directory
    copy2(sFilename_configuration_in, sFilename_configuration_copy)
    #also copy the basin configuration file to the output directory
    sFilename_configuration_basins_copy = os.path.join( sWorkspace_output, 'pyflowline_configuration_basins_copy.json' )    
    copy2(sFilename_basins_in, sFilename_configuration_basins_copy)
    sFilename_basins = sFilename_configuration_basins_copy

    #now use the new configuration file to create class object
    sFilename_configuration = sFilename_configuration_copy
    change_json_key_value(sFilename_configuration, 'sFilename_mesh_boundary', sFilename_wbd_boundary) 
    change_json_key_value(sFilename_configuration, 'sFilename_basins', sFilename_basins) #individual basin configuration file
    sFilename_elevation=  sWorkspace_input +  "/" + r'dem.tif'
    change_json_key_value(sFilename_configuration, 'sFilename_dem', sFilename_dem_arctic)
    sPath_dummy = os.path.join( sWorkspace_input, 'dggrid' +  sResolution)
    sFilename_flowline  = os.path.join( sPath_dummy, 'river_networks.geojson' )
    change_json_key_value(sFilename_basins, 'sFilename_flowline_filter', sFilename_flowline, iFlag_basin_in=1) #user provided flowline

    #a minimal length of 5 grid cells for small river
    oPyhexwatershed.pPyFlowline.aBasin[0].dThreshold_small_river = dResolution * 5 
    oPyhexwatershed.pPyFlowline.pyflowline_change_model_parameter('dLongitude_outlet_degree', dLongitude_outlet_degree, iFlag_basin_in= 1)
    oPyhexwatershed.pPyFlowline.pyflowline_change_model_parameter('dLatitude_outlet_degree', dLatitude_outlet_degree, iFlag_basin_in= 1)
    oPyhexwatershed.pPyFlowline.pyflowline_change_model_parameter('sFilename_flowline_filter', sFilename_flowline, iFlag_basin_in= 1)

    if iFlag_create_job == 1:
        oPyhexwatershed._pyhexwatershed_create_hpc_job()
        print(iCase_index)
        sLine  = 'cd ' + oPyhexwatershed.sWorkspace_output + '\n'
        ofs.write(sLine)
        sLine  = 'sbatch submit.job' + '\n'
        ofs.write(sLine)
    else:
        #oPyhexwatershed.pyhexwatershed_export() #for testing  
        pass

    if iFlag_visualization == 1:
        pBasin_hexwatershed = oPyhexwatershed.aBasin[0]
        sWorkspace_output_basin = pBasin_hexwatershed.sWorkspace_output_basin

        sFilename = os.path.join( oPyhexwatershed.sWorkspace_output_hexwatershed, 'hexwatershed.mp4' )
        #oPyhexwatershed._animate(sFilename,  iFlag_type_in = 3,  iFont_size_in= 14)
        #exit()

        #polyline
        sFilename = os.path.join( sWorkspace_output_basin, 'flow_direction.png' )
        #oPyhexwatershed.plot( sVariable_in = 'flow_direction', sFilename_output_in = sFilename, iFont_size_in= 14,iFlag_title_in=1)

        #polygon    
             
        sFilename = os.path.join( sWorkspace_output_basin, 'surface_elevation.png' )    
        #oPyhexwatershed.plot( sVariable_in = 'elevation', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=2500,iFlag_title_in=1, iFlag_colorbar_in = 1)     

        sFilename = os.path.join( sWorkspace_output_basin, 'surface_slope.png' )        
        #oPyhexwatershed.plot( sVariable_in = 'slope', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=0.1, iFlag_title_in=1,iFlag_colorbar_in = 1 )

        sFilename = os.path.join( sWorkspace_output_basin, 'drainage_area.png' )
        #oPyhexwatershed.plot( sVariable_in = 'drainage_area',  sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=8.5E11, iFlag_title_in=1, iFlag_colorbar_in = 1,iFlag_scientific_notation_colorbar_in=1 )

        sFilename = os.path.join( sWorkspace_output_basin, 'travel_distance.png' )
        #oPyhexwatershed.plot( sVariable_in = 'travel_distance', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=3.4E6 ,iFlag_title_in=1,iFlag_colorbar_in=1,iFlag_scientific_notation_colorbar_in=1)
        #mixed
        sFilename = os.path.join( sWorkspace_output_basin, 'flow_direction_w_mesh.png' )
        #oPyhexwatershed.plot( sVariable_in = 'flow_direction_with_mesh', sFilename_output_in = sFilename)  

        sFilename = os.path.join( sWorkspace_output_basin, 'flow_direction_w_observation.png' )
        oPyhexwatershed.plot( sVariable_in = 'flow_direction_with_observation',  sFilename_output_in = sFilename,
                             iFont_size_in= 14,iFlag_title_in=1, 
                             iFlag_openstreetmap_in=0,
                             aExtent_in=aExtent)
  

    iCase_index = iCase_index + 1

if iFlag_create_job ==1:
    ofs.close()
    os.chmod(sFilename_job, stat.S_IREAD | stat.S_IWRITE | stat.S_IXUSR)   
