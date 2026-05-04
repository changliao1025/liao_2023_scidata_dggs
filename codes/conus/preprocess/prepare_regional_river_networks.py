
import os,  stat
from osgeo import ogr, osr
from pathlib import Path
from os.path import realpath
#this is the pre-processing of the data
#basically the input needs at least (1) boundary, (2) dem

from pyearth.system.define_global_variables import *
#import two functions

from codes.shared.prepare_regional_river_networks import prepare_regional_river_networks
#from codes.shared.prepare_regional_dem import prepare_regional_dem

sPath = str( Path().resolve() )
iFlag_option = 1
sWorkspace_data = realpath( sPath +  '/data/' )
sWorkspace_input =  str(Path(sWorkspace_data)  /  'conus')
Path(sWorkspace_input).mkdir(parents=True, exist_ok=True)
sWorkspace_output=  '/compyfs/liao313/04model/pyhexwatershed/conus'
Path(sWorkspace_output).mkdir(parents=True, exist_ok=True)

pDriver_geojson = ogr.GetDriverByName('GeoJSON')
#create a geojson file for the boundary of the conus

#set the boundary of the conus
aBound_conus = [-125, 24, -66, 50]
sFilename_wbd_boundary = sWorkspace_input + '/boundary.geojson'
#create a polygon geojson file for the boundary
if os.path.exists(sFilename_wbd_boundary):
    os.remove(sFilename_wbd_boundary)
pDataSource = pDriver_geojson.CreateDataSource(sFilename_wbd_boundary)
pSrs = osr.SpatialReference()
pSrs.ImportFromEPSG(4326)
pLayer = pDataSource.CreateLayer('boundary', geom_type=ogr.wkbPolygon, srs=pSrs)
feature = ogr.Feature(pLayer.GetLayerDefn())
ring = ogr.Geometry(ogr.wkbLinearRing)
ring.AddPoint(aBound_conus[0], aBound_conus[1])
ring.AddPoint(aBound_conus[2], aBound_conus[1])
ring.AddPoint(aBound_conus[2], aBound_conus[3])
ring.AddPoint(aBound_conus[0], aBound_conus[3])
ring.AddPoint(aBound_conus[0], aBound_conus[1])
polygon = ogr.Geometry(ogr.wkbPolygon)
polygon.AddGeometry(ring)
feature.SetGeometry(polygon)
pLayer.CreateFeature(feature)
pDataSource = None
pSrs = None

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



print('Done')
