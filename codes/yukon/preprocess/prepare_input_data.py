
import os,  stat

from pathlib import Path
from os.path import realpath
#this is the pre-processing of the data
#basically the input needs at least (1) boundary, (2) dem


#import two functions 

from codes.shared.prepare_regional_river_networks import prepare_regional_river_networks
from codes.shared.prepare_regional_dem import prepare_regional_dem

sPath = str( Path().resolve() )
iFlag_option = 1
sWorkspace_data = realpath( sPath +  '/data/yukon' )
sWorkspace_input =  str(Path(sWorkspace_data)  /  'input')
sWorkspace_output=  '/compyfs/liao313/04model/pyhexwatershed/yukon'


#will copy this later
#sFilename_wbd_boundary = sWorkspace_data + '/boundary.geojson' 
sFilename_wbd_boundary = '/qfs/people/liao313/data/hexwatershed/yukon/vector/hydrology/boundary.geojson'

sFilename_dem_arctic = '/qfs/people/liao313/data/hexwatershed/yukon/raster/dem/hyd_ar_dem_30s/hyd_ar_dem_30s.tif'

sFilename_dem = sWorkspace_input + '/dem.tif'

#clip the large dem using the boundary
#the river dataset is simulated by the REACH library
sFolder_river_networks = '/compyfs/liao313/04model/reach'
aResolution = [30, 18, 10, 6]
#the large river networks 
for iResolution in range(10, 14,1):
    #convert resoluton to string
    sResolution = "{:0d}".format(iResolution) 
    dResolution = aResolution[iResolution-10]
    sResolution1 = "{:0d}".format(dResolution)   

    sFilename_rivernetworks = sFolder_river_networks + '/dggrid'  + sResolution + '/filtered_' + sResolution1 + '/filtered_' + sResolution1 + 'km.shp'
    sFilename_watershed_boundary = sFilename_wbd_boundary

    #make a folder the resolution
    sFolder_resolution = sWorkspace_input + '/dggrid' + sResolution
    if not os.path.exists(sFolder_resolution):
        os.makedirs(sFolder_resolution)

    sFilename_river_networks_out = sFolder_resolution + '/river_networks.geojson'

    prepare_regional_river_networks(sFilename_rivernetworks, sFilename_watershed_boundary, sFilename_river_networks_out)

#prepare_regional_dem(sFilename_dem_arctic,sFilename_wbd_boundary, sFilename_dem) 

print('Done')
