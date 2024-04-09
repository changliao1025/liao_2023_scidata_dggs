import os
import numpy as np
from pathlib import Path
from os.path import realpath
import matplotlib as mpl
#import gdal
from osgeo import gdal, ogr, osr
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyflowline.external.tinyr.tinyr.tinyr import RTree
from pyhexwatershed.pyhexwatershed_read_model_configuration_file import pyhexwatershed_read_model_configuration_file
from pyearth.gis.gdal.read.raster.gdal_read_geotiff_file import gdal_read_geotiff_file
from pyearth.visual.scatter.scatter_plot_multiple_data import scatter_plot_multiple_data
#read the level 14 dggrid file

sFilename = '/compyfs/liao313/04model/pyhexwatershed/yukon/pyhexwatershed20240102014/pyflowline/dggrid.geojson'

#define the gdal driver to read the geojson file
pDriver_geojson = ogr.GetDriverByName('GeoJSON')

#read the geojson file using gdal
pDataSource = pDriver_geojson.Open(sFilename, 0)
pLayer_mesh = pDataSource.GetLayer()
#get the feature count
iFeatureCount = pLayer_mesh.GetFeatureCount()
print(iFeatureCount)


aLongitude = np.zeros(iFeatureCount)
aLatitude = np.zeros(iFeatureCount)

for i in range(iFeatureCount):
    pFeature_mesh = pLayer_mesh.GetFeature(i)    
    pGeometry_mesh = pFeature_mesh.GetGeometryRef()   
    dLon = pFeature_mesh.GetField("longitude")
    dLat = pFeature_mesh.GetField("latitude")        
    dArea = pFeature_mesh.GetField("area")
    if i==0:
        print(np.sqrt(dArea)/1000)

    aLongitude[i] = dLon
    aLatitude[i] = dLat


#build a rtrees index



#generate n random points

iNumber = 500
aRandom = np.random.randint(0, iFeatureCount, iNumber)

#read the dem data in the dem format
sFilename_dem_in = '/qfs/people/liao313/workspace/python/liao_2023_scidata_dggs/data/yukon/input/dem.tif'
#pDataset_elevation = gdal.Open(sFilename_dem_in, gdal.GA_ReadOnly)
dummy = gdal_read_geotiff_file(sFilename_dem_in)
aDem_in= dummy['dataOut']
dPixelWidth = dummy['pixelWidth']                        
pPixelHeight = dummy['pixelHeight']
dOriginX = dummy['originX']
dOriginY = dummy['originY']
nrow = dummy['nrow']
ncolumn = dummy['ncolumn']
dMissing_value= dummy['missingValue']
#pGeotransform = dummy['geotransform']
pProjection = dummy['projection']
pSpatialRef_target = dummy['spatialReference']
aVariable_obs = np.full(iNumber, np.nan)
pSrs = osr.SpatialReference()  
pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
for j in range(iNumber):    
    aX_target = aLongitude[aRandom[j]]
    aY_target = aLatitude[aRandom[j]]
    #project the coordinates to the dem projection
    #aX_target, aY_target = reproject_coordinates(x, y, pSrs, pSpatialRef_target)
    #convert the coordinates to row and column index
    lRow  = int((aY_target - dOriginY)/pPixelHeight)
    lColumn = int((aX_target - dOriginX)/abs(dPixelWidth))
    dElevation = aDem_in[lRow, lColumn]
    #use x and y to get column and row index from the dem data

    aVariable_obs[j] = dElevation


sPath_parent = str(Path(__file__).parents[3]) # data is located two dir's up
sPath_data = realpath( sPath_parent +  '/data/yukon' )
sWorkspace_input =  str(Path(sPath_data)  /  'input')
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/yukon'

sFilename_configuration_in = realpath( sPath_parent +  '/data/yukon/input/pyhexwatershed_yukon_dggrid.json' )

aResolution_index = [10, 11, 12, 13]
aCase_index = aResolution_index
nCase = len(aResolution_index)
aData_x = list()
aData_y = list()
aLabel_legend=list()

#mpas mesh only has one resolution
iFlag_stream_burning_topology = 1 
iFlag_use_mesh_dem = 0
iFlag_elevation_profile = 0
sMesh_type='dggrid'
sDggrid_type = 'ISEA3H'
sDate= '20240102'
for iCase in range(0, 4, 1):   
    iResolution_index = aResolution_index[iCase]
    iCase_index = aCase_index[iCase]  
    aLabel_legend.append(   'ISEA3H Level ' +  "{:0d}".format(iCase_index)   )
    oPyhexwatershed = pyhexwatershed_read_model_configuration_file(sFilename_configuration_in,
                    iCase_index_in=iCase_index,iFlag_stream_burning_topology_in=iFlag_stream_burning_topology,                   
                    iResolution_index_in = iResolution_index, 
                    sDggrid_type_in=sDggrid_type,
                    sDate_in= sDate, sMesh_type_in= sMesh_type)  

    #read the output file
    pBasin_hexwatershed = oPyhexwatershed.aBasin[0]
    sWorkspace_output_basin = pBasin_hexwatershed.sWorkspace_output_basin
    print(sWorkspace_output_basin)
    sFilename = os.path.join(  sWorkspace_output_basin, 'variable_polygon.geojson' ) 
    pDataSource_variable = pDriver_geojson.Open(sFilename, 0)
    pLayer_variable = pDataSource_variable.GetLayer()
    #get the feature count
    iFeatureCount = pLayer_variable.GetFeatureCount()   
    print(iFeatureCount)
    interleaved = True
    index_mesh = RTree(interleaved=interleaved, max_cap=5, min_cap=2)
    aVariable = list()
    for i in range(iFeatureCount):
        pFeature_variable = pLayer_variable.GetFeature(i)    
        #get the center of the cell
        pGeometry_variable = pFeature_variable.GetGeometryRef()
        aCoords_gcs = get_geometry_coordinates(pGeometry_variable)      
        dElevation = pFeature_variable.GetField("elevation")
        if i ==0:
            dArea = pFeature_variable.GetField("area")
            print( np.sqrt(dArea)/1000  )

        left, right, bottom, top= pGeometry_variable.GetEnvelope()   
        pBound= (left, bottom, right, top)
        index_mesh.insert(i, pBound)  #  
        aVariable.append(dElevation)
  

    aVariable_sample = np.full(iNumber, np.nan)
    for j in range(iNumber):        
        x = aLongitude[aRandom[j]]
        y = aLatitude[aRandom[j]]
        left= x - 1E-5
        right= x + 1E-5
        bottom= y-1E-5
        top=    y+1E-5
        pBound= (left, bottom, right, top)      
        aIntersect = list(index_mesh.search(pBound))
        for k in aIntersect:
            dVariable = aVariable[k]    
            aVariable_sample[j] = dVariable
            pass

    
    aDatax0 = np.asarray(aVariable_obs)
    aDatay0 = np.asarray(aVariable_sample)

    aData_x.append(aDatax0)
    aData_y.append(aDatay0)


sFilename_out = sPath_parent + '/' + 'figures' + '/yukon/' + 'scatterplot_elevation_multiple.png'
aColor = np.full(nCase, None, dtype=object)
aMarker= np.full(nCase, None, dtype=object)
aSize = np.full(nCase, mpl.rcParams['lines.markersize'] ** 2.5 , dtype=object)

nmesh=nCase
aColor= create_diverge_rgb_color_hex(nmesh)
aMarker = [ '.','o','+','x','^']

#exit()
scatter_plot_multiple_data(aData_x, 
                      aData_y,
                      sFilename_out,  
                      iFlag_miniplot_in = None,
                      iFlag_scientific_notation_x_in=0,
                      iFlag_scientific_notation_y_in=0,
                      iSize_x_in = None, 
                      iSize_y_in = None,  
                      iDPI_in = None ,
                      iFlag_log_x_in = 0,
                      iFlag_log_y_in = 0,             
                      dMin_x_in=0,
                      dMax_x_in=2500,
                      dMin_y_in=0,
                      dMax_y_in=2500,     
                      dSpace_x_in = None, 
                      dSpace_y_in = None, 
                      sFormat_x_in =None,
                      sFormat_y_in =None,
                      sLabel_x_in ='HydroSHEDS, m',
                      sLabel_y_in = 'HexWatershed, m' ,                    
                      aColor_in=aColor,
                      aMarker_in=aMarker,
                      aSize_in = aSize,
                      aLabel_legend_in = aLabel_legend,
                      sTitle_in = 'Surface elevation')






