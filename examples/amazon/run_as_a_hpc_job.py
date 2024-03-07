
import os,  stat

from pathlib import Path
from os.path import realpath
 #import copy2
from shutil import copy2


import logging
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(format='%(asctime)s %(message)s')
logging.warning('is the time pyhexwatershed simulation started.')


from pyhexwatershed.pyhexwatershed_read_model_configuration_file import pyhexwatershed_read_model_configuration_file
from pyflowline.mesh.dggrid.create_dggrid_mesh import dggrid_find_resolution_by_index
from pyhexwatershed.change_json_key_value import change_json_key_value

sMesh_type = 'dggrid'
iCase_index = 1
dResolution_meter=5000
iFlag_create_job = 0
iFlag_visualization = 0
aExtent_full = None

#-49.47644,1.08260 from location capture
dLongitude_outlet_degree=-49.47644
dLatitude_outlet_degree=1.08260

#for case 12
#dLongitude_outlet_degree= -50.58226
#dLatitude_outlet_degree=-0.10195

#case 13
#dLongitude_outlet_degree= -50.54811
#dLatitude_outlet_degree=-0.07100


pProjection_map = None
sDate='20240101'
sPath = str( Path().resolve() )
iFlag_option = 1
sWorkspace_data = realpath( sPath +  '/data/amazon' )
sWorkspace_input =  str(Path(sWorkspace_data)  /  'input')
sWorkspace_output=  '/compyfs/liao313/04model/pyhexwatershed/amazon'

iMesh_type = 5

#set up dggrid resolution level 
aResolution_index = [10, 11, 12, 13, 14]
nCase = len(aResolution_index)
aCase_index = [10, 11, 12, 13, 14] #aResolution_index
sDggrid_type = 'ISEA3H'
#generate a bash job script
if iFlag_create_job == 1:
    sFilename_job = sWorkspace_output + '/' + sDate  + 'submit.bash'
    ofs = open(sFilename_job, 'w')
    sLine  = '#!/bin/bash' + '\n'
    ofs.write(sLine)

sFilename_configuration_in = realpath( sWorkspace_input +  '/pyhexwatershed_amazon_dggrid.json' )

sFilename_wbd_boundary = '/qfs/people/liao313/workspace/python/liao_2023_scidata_dggs/data/amazon/input/mesh_boundary_new.geojson'
#sFilename_wbd_boundary = '/qfs/people/liao313/workspace/python/liao_2023_scidata_dggs/data/amazon/input/mesh_boundary_buffer.geojson'
sFilename_basins = '/qfs/people/liao313/workspace/python/liao_2023_scidata_dggs/data/amazon/input/pyflowline_amazon_basins.json'
#sFilename_dem = ''
    
if os.path.isfile(sFilename_configuration_in):
    print(sFilename_configuration_in)
else:
    print('This configuration file does not exist: ', sFilename_configuration_in )
    exit()
    
#mpas mesh only has one resolution
iFlag_stream_burning_topology = 1 
iFlag_use_mesh_dem = 0
iFlag_elevation_profile = 0

aExtent = [-60.6, -59.2, -3.6, -2.5]
aExtent = None

for iCase in range(1, 2, 1):    
    iResolution_index = aResolution_index[iCase]
    sResolution = "{:0d}".format(iResolution_index)
    iCase_index = aCase_index[iCase]
    dResolution = dggrid_find_resolution_by_index(sDggrid_type, iResolution_index)
    print(dResolution)   
    sPath_dummy = os.path.join( sWorkspace_input, 'dggrid' +  sResolution)
    sFilename_flowline = os.path.join(  sPath_dummy, 'river_networks_extend.geojson' )    
    oPyhexwatershed = pyhexwatershed_read_model_configuration_file(sFilename_configuration_in,
                    iCase_index_in=iCase_index,iFlag_stream_burning_topology_in=iFlag_stream_burning_topology,
                    iFlag_use_mesh_dem_in=iFlag_use_mesh_dem,
                    iFlag_elevation_profile_in=iFlag_elevation_profile,
                    iResolution_index_in = iResolution_index, 
                    sDggrid_type_in=sDggrid_type,
                    sDate_in = sDate, sMesh_type_in= sMesh_type)  

    #copy to the output directory
    sWorkspace_output = oPyhexwatershed.sWorkspace_output
    sFilename_configuration_copy= os.path.join( sWorkspace_output, 'pyhexwatershed_configuration_copy.json' )
    #copy the configuration file to the output directory
    copy2(sFilename_configuration_in, sFilename_configuration_copy)
    #also copy the basin configuration file to the output directory
    #sFilename_configuration_basins_in = oPyhexwatershed.sFilename_basins
    sFilename_configuration_basins_copy = os.path.join( sWorkspace_output, 'pyflowline_configuration_basins_copy.json' )    
    copy2(sFilename_basins, sFilename_configuration_basins_copy)

    #now use the new configuration file to create class object
    sFilename_configuration = sFilename_configuration_copy
    change_json_key_value(sFilename_configuration, 'sFilename_mesh_boundary', sFilename_wbd_boundary) 
    change_json_key_value(sFilename_configuration, 'sFilename_basins', sFilename_configuration_basins_copy) #individual basin configuration file    
    change_json_key_value(sFilename_configuration_basins_copy, 'sFilename_flowline_filter', sFilename_flowline, iFlag_basin_in=1) #user provided flowline
    oPyhexwatershed = None
    oPyhexwatershed = pyhexwatershed_read_model_configuration_file(sFilename_configuration,
                    iCase_index_in=iCase_index,iFlag_stream_burning_topology_in=iFlag_stream_burning_topology,
                    iFlag_use_mesh_dem_in=iFlag_use_mesh_dem,
                    iFlag_elevation_profile_in=iFlag_elevation_profile,
                    iResolution_index_in = iResolution_index, 
                    sDggrid_type_in=sDggrid_type,
                    sDate_in= sDate, sMesh_type_in = sMesh_type)  

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
        oPyhexwatershed.pyhexwatershed_export() #for testing  
        pass

    if iFlag_visualization == 1:
        pBasin_hexwatershed = oPyhexwatershed.aBasin[0]
        sWorkspace_output_basin = pBasin_hexwatershed.sWorkspace_output_basin

        #animation
        sFilename = os.path.join( oPyhexwatershed.sWorkspace_output_hexwatershed, 'hexwatershed.mp4' )
        #oPyhexwatershed._animate(sFilename,  iFlag_type_in = 3,  iFont_size_in= 14)
        #exit()

        #polyline
        sFilename = os.path.join( sWorkspace_output_basin, 'flow_direction.png' )
        oPyhexwatershed.plot( sVariable_in = 'flow_direction', sFilename_output_in = sFilename, iFont_size_in= 14,iFlag_title_in=1)

        #polygon    
       
        sFilename = os.path.join(  sWorkspace_output_basin, 'surface_elevation.png' )    
        oPyhexwatershed.plot( sVariable_in = 'elevation', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=5000,iFlag_title_in=1, iFlag_colorbar_in = 1)     

        sFilename = os.path.join(  sWorkspace_output_basin, 'surface_slope.png' )        
        oPyhexwatershed.plot( sVariable_in = 'slope', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=0.1, iFlag_title_in=1,iFlag_colorbar_in = 1 )

        sFilename = os.path.join( sWorkspace_output_basin, 'drainage_area.png' )
        oPyhexwatershed.plot( sVariable_in = 'drainage_area',  sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=6.0E12, iFlag_title_in=1, iFlag_colorbar_in = 0,iFlag_scientific_notation_colorbar_in=1 )

        sFilename = os.path.join(  sWorkspace_output_basin, 'travel_distance.png' )
        oPyhexwatershed.plot( sVariable_in = 'travel_distance', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=5.8E6 ,iFlag_title_in=1,iFlag_colorbar_in=0,iFlag_scientific_notation_colorbar_in=1)
        #mixed
        sFilename = os.path.join( sWorkspace_output_basin, 'flow_direction_w_mesh.png' )
        #oPyhexwatershed.plot( sVariable_in = 'flow_direction_with_mesh', sFilename_output_in = sFilename)  

        sFilename = os.path.join(  sWorkspace_output_basin, 'flow_direction_w_observation.png' )
        #sFilename = None
        oPyhexwatershed.plot( sVariable_in = 'flow_direction_with_observation', sFilename_output_in = sFilename,
                             iFont_size_in = 14,iFlag_title_in=1, 
                             iFlag_openstreetmap_in=0,
                             aExtent_in=aExtent)
  

    iCase_index = iCase_index + 1

if iFlag_create_job ==1:
    ofs.close()
    os.chmod(sFilename_job, stat.S_IREAD | stat.S_IWRITE | stat.S_IXUSR)   
