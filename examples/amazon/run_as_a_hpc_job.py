
import os,  stat

from pathlib import Path
from os.path import realpath

import numpy as np
import logging
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(format='%(asctime)s %(message)s')
logging.warning('is the time pyhexwatershed simulation started.')


from pyhexwatershed.pyhexwatershed_read_model_configuration_file import pyhexwatershed_read_model_configuration_file
from pyflowline.mesh.dggrid.create_dggrid_mesh import dggrid_find_resolution_by_index

sMesh_type = 'dggrid'
iCase_index = 1
dResolution_meter=5000
iFlag_create_job=0
iFlag_visualization =1
aExtent_full = None
dLongitude_outlet_degree=-117
dLatitude_outlet_degree=42
pProjection_map = None
sDate='20230801'
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
if iFlag_create_job ==1:
    sFilename_job = sWorkspace_output + '/' + sDate  + 'submit.bash'
    ofs = open(sFilename_job, 'w')
    sLine  = '#!/bin/bash' + '\n'
    ofs.write(sLine)

sFilename_configuration_in = realpath( sPath +  '/examples/amazon/pyhexwatershed_amazon_dggrid.json' )

    
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
#aExtent = None
for iCase in range(0, 4, 1):    
    iResolution_index = aResolution_index[iCase]
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

    #a minimal length of 5 grid cells for small river
    oPyhexwatershed.pPyFlowline.aBasin[0].dThreshold_small_river = dResolution * 10 

    if iFlag_create_job == 1:
        oPyhexwatershed._create_hpc_job()
        print(iCase_index)
        sLine  = 'cd ' + oPyhexwatershed.sWorkspace_output + '\n'
        ofs.write(sLine)
        sLine  = 'sbatch submit.job' + '\n'
        ofs.write(sLine)
    else:
        #oPyhexwatershed.export() #for testing  
        pass

    if iFlag_visualization == 1:
        pBasin_hexwatershed = oPyhexwatershed.aBasin[0]
        sWorkspace_output_basin = pBasin_hexwatershed.sWorkspace_output_basin

        #polyline
        sFilename = os.path.join( sWorkspace_output_basin, 'flow_direction.png' )
        #oPyhexwatershed.plot( sVariable_in = 'flow_direction', sFilename_output_in = sFilename, iFont_size_in= 14,iFlag_title_in=1)

        #polygon    

        sFilename = os.path.join(  sWorkspace_output_basin, 'area.png' )    
        #oPyhexwatershed.plot( sVariable_in = 'area', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=None,iFlag_title_in=1, iFlag_colorbar_in = 1,iFlag_scientific_notation_colorbar_in=1)     

        sFilename = os.path.join(  sWorkspace_output_basin, 'surface_elevation.png' )    
        #oPyhexwatershed.plot( sVariable_in = 'elevation', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=5000,iFlag_title_in=0, iFlag_colorbar_in = 0)     

        sFilename = os.path.join(  sWorkspace_output_basin, 'surface_slope.png' )        
        #oPyhexwatershed.plot( sVariable_in = 'slope', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=0.1, iFlag_title_in=0,iFlag_colorbar_in = 0 )

        sFilename = os.path.join( sWorkspace_output_basin, 'drainage_area.png' )
        #oPyhexwatershed.plot( sVariable_in = 'drainage_area',  sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=6.0E12, iFlag_title_in=0, iFlag_colorbar_in = 0,iFlag_scientific_notation_colorbar_in=0 )

        sFilename = os.path.join(  sWorkspace_output_basin, 'travel_distance.png' )
        #oPyhexwatershed.plot( sVariable_in = 'travel_distance', sFilename_output_in = sFilename, iFont_size_in= 14, dData_min_in=0, dData_max_in=5.8E6 ,iFlag_title_in=0,iFlag_colorbar_in=0,iFlag_scientific_notation_colorbar_in=0)
        #mixed
        sFilename = os.path.join( sWorkspace_output_basin, 'flow_direction_w_mesh.png' )
        #oPyhexwatershed.plot( sVariable_in = 'flow_direction_with_mesh', sFilename_output_in = sFilename)  

        sFilename = os.path.join(  sWorkspace_output_basin, 'flow_direction_w_observation_manaus.png' )
        #sFilename = None
        oPyhexwatershed.plot( sVariable_in = 'flow_direction_with_observation',  sFilename_output_in = sFilename,
                             iFont_size_in= 14,iFlag_title_in=0, 
                             iFlag_openstreetmap_in=1,
                             aExtent_in=aExtent)
  

    iCase_index = iCase_index + 1

if iFlag_create_job ==1:
    ofs.close()
    os.chmod(sFilename_job, stat.S_IREAD | stat.S_IWRITE | stat.S_IXUSR)   
