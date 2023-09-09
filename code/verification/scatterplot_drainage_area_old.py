import os, sys, stat
from pathlib import Path
from os.path import realpath
import json
import numpy as np
from osgeo import  osr, gdal, ogr
import matplotlib as mpl
import rtree
from pyearth.visual.scatter.scatter_plot_data import scatter_plot_data
from pyearth.toolbox.reader.text_reader_string import text_reader_string
from pyearth.visual.scatter.scatter_plot_multiple_data import scatter_plot_multiple_data
from pyflowline.external.pyearth.gis.gdal.gdal_functions  import meter_to_degree, degree_to_meter
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex
from pyhexwatershed.pyhexwatershed_read_model_configuration_file import pyhexwatershed_read_model_configuration_file
from pyflowline.formats.read_flowline import read_flowline_geojson
from pyflowline.mesh.dggrid.create_dggrid_mesh import dggrid_find_resolution_by_index
from pyflowline.external.pyearth.gis.gdal.gdal_functions import gdal_read_geotiff_file, reproject_coordinates
sPath_parent = str(Path(__file__).parents[2]) # data is located two dir's up
print(sPath_parent)
sPath_data = realpath( sPath_parent +  '/data/amazon' )
sWorkspace_input =  str(Path(sPath_data)  /  'input')
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/amazon'
aResolution_index = [10, 11, 12, 13]
nCase = len(aResolution_index)
aCase_index = aResolution_index
sDate='20230801'
sMesh_type = 'dggrid'
sFilename_configuration_in = realpath( sPath_parent +  '/examples/amazon/pyhexwatershed_amazon_dggrid.json' )

   
# we define a dictionnary with months that we'll use later
case_dict = dict()

for i in range(1, nCase +1):
    case_dict[i] = 'Case ' +  "{:0d}".format(i) 

print(case_dict)   

aDistance_niws =list()
aData_x=list()
aData_y  =list()
aLabel_legend =list()

#read the point location using gdal 
pSrs = osr.SpatialReference()  
pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
sFilename_confluence = '/qfs/people/liao313/data/hexwatershed/amazon/vector/amztrbmth.geojson'
pDriver = ogr.GetDriverByName('GeoJSON')
pDataSource = pDriver.Open(sFilename_confluence, 0)
pLayer = pDataSource.GetLayer()
nPoint = pLayer.GetFeatureCount()
aLongitude = np.full(nPoint, None, dtype=float)
aLatitude = np.full(nPoint, None, dtype=float)
for i in range(nPoint):
    pFeature = pLayer.GetFeature(i)
    pGeometry = pFeature.GetGeometryRef()
    #get x and y of the point
    dLongitude = pGeometry.GetX()
    dLatitude = pGeometry.GetY()
    aLongitude[i] = dLongitude
    aLatitude[i] = dLatitude

pDataSource = None

#the flow accumenation is in the same order as the point

sFilename_flow_accumulation = '/compyfs/liao313/00raw/hydrology/amazon/CD06_CAMREX_1086/data/amzfloacc.tif'
pRaster = gdal.Open(sFilename_flow_accumulation)
aFlow_accumulation_in, dPixelWidth, dOriginX, dOriginY, \
            nrow, ncolumn,dMissing_value, pGeotransform, pProjection,  pSpatialRef_target = gdal_read_geotiff_file(sFilename_flow_accumulation)

dResolution_degree_in= 0.005
dLatitude_mean = np.mean(aLatitude)
dResolution_flow_accumulation = degree_to_meter( dLatitude_mean ,dResolution_degree_in)

dX_left=dOriginX
dX_right = dOriginX + ncolumn * dPixelWidth
dY_top = dOriginY
dY_bot = dOriginY - nrow * dPixelWidth
aFlow_accumulation = np.full(nPoint, None, dtype=float)
for i in range(nPoint):
    x1 = aLongitude[i]
    y1 = aLatitude[i]
    dX_out,dY_out = reproject_coordinates(x1,y1, pSrs,pSpatialRef_target)   
    dDummy1 = (dX_out - dX_left) / dPixelWidth
    lColumn_index = int(dDummy1)
    dDummy2 = (dY_top - dY_out) / dPixelWidth
    lRow_index = int(dDummy2)
    if lColumn_index >= ncolumn or lColumn_index < 0 \
        or lRow_index >= nrow or lRow_index < 0:        
        #this pixel is out of bound            
        continue
    else:         
        dFlow_accumulation = aFlow_accumulation_in[lRow_index, lColumn_index]     
        aFlow_accumulation[i] = dFlow_accumulation * dResolution_flow_accumulation * dResolution_flow_accumulation / 1.0E6
       
#read the model resolution 
sDggrid_type = 'ISEA3H'
aDrainage_area = np.full(( nCase,nPoint), None, dtype=float)
for iCase in range(0, nCase ):
    #read 
    iCase_index = aCase_index[iCase]  
    iResolution_index = aResolution_index[iCase]
    dResolution = dggrid_find_resolution_by_index(sDggrid_type, iResolution_index)
    print(dResolution)     

    oPyhexwatershed = pyhexwatershed_read_model_configuration_file(sFilename_configuration_in,
                    iCase_index_in=iCase_index,
                    sDate_in= sDate, sMesh_type_in= sMesh_type)  

    #read the variable polygon
    sFilename_polygon = oPyhexwatershed.aBasin[0].sFilename_variable_polygon
    pDataSource = pDriver.Open(sFilename_polygon, 0)
    pLayer_base = pDataSource.GetLayer()
    nFeature_base = pLayer_base.GetFeatureCount()
    index_base = rtree.index.Index()

    for i in range(nFeature_base):
        lID = i 
        pFeature_base = pLayer_base.GetFeature(i)
        pGeometry_base = pFeature_base.GetGeometryRef()    
        left, right, bottom, top= pGeometry_base.GetEnvelope()   
        pBound= (left, bottom, right, top)
        index_base.insert(lID, pBound)  #
        
    #find the confluence point
    for i in range(nPoint):       
        #convert meter to degree
        dLatitude_mean=     aLatitude[i]
        dResolution_degree = meter_to_degree(dResolution, dLatitude_mean) * 2
        pBound= (aLongitude[i] - dResolution_degree, 
                 aLatitude[i] - dResolution_degree, 
                 aLongitude[i] + dResolution_degree, 
                 aLatitude[i] + dResolution_degree)
        aIntersect = list(index_base.intersection(pBound))
        aDummy = list()
        for k in aIntersect:
            pFeature_base = pLayer_base.GetFeature(k)            
            dDrainage_base = pFeature_base.GetField("drainage_area")
            aDummy.append(dDrainage_base/ 1.0E6)
        
        #find the closese to the observation
        dFlow_accumulation = aFlow_accumulation[i]
        #calculate the difference
        aDummy = np.asarray(aDummy)
        aDummy0 = np.abs(aDummy - dFlow_accumulation)
        #find th minimum and its index
        iIndex = np.argmin(aDummy0)
        aDrainage_area[iCase, i] = aDummy[iIndex] 
    
    #close the file
    pDataSource = None


sFilename_out = sPath_parent + '/' + 'figures' + '/' + 'scatterplot_drainage_area_multiple.png'

#5 cases + 1 obs

aColor = np.full(nCase, None, dtype=object)
aMarker= np.full(nCase, None, dtype=object)
aSize = np.full(nCase, mpl.rcParams['lines.markersize'] ** 2 , dtype=object)

nmesh=nCase
aColor= create_diverge_rgb_color_hex(nmesh)
aMarker = [ '.','o','+','x','^']

aFlow_accumulation = np.array(aFlow_accumulation) 
aData_x = list()
aData_y = list()

for i in range(nCase):
    iCase_index = aCase_index[i]  
    aLabel_legend.append(   'ISEA3H Level ' +  "{:0d}".format(iCase_index)   )
    aDatax0 = list()
    aDatay0 = list()
    for j in range(nPoint):
        aDatax0.append(aFlow_accumulation[j])
        aDatay0.append(aDrainage_area[i,j])
    
    aDatax0 = np.asarray(aDatax0)
    aDatay0 = np.asarray(aDatay0)
    aData_x.append(aDatax0)
    aData_y.append(aDatay0)


scatter_plot_multiple_data(aData_x, 
                      aData_y,
                      sFilename_out,  
                      iFlag_scientific_notation_x_in=1,
                      iFlag_scientific_notation_y_in=1,
                      iSize_x_in = None, 
                      iSize_y_in = None,  
                      iDPI_in = None ,
                      iFlag_log_x_in = None,
                      iFlag_log_y_in = None,
                      dMin_x_in = 0, 
                      dMax_x_in = np.max(aFlow_accumulation), 
                      dMin_y_in = 0, 
                      dMax_y_in = np.max(aFlow_accumulation), 
                      dSpace_x_in = None, 
                      dSpace_y_in = None, 
                      sFormat_x_in =None,
                      sFormat_y_in =None,
                      sLabel_x_in ='LBA-ECO data, $km^{2}$',
                      sLabel_y_in = 'HexWatershed, $km^{2}$)' , 
                      aColor_in=aColor,
                      aMarker_in=aMarker,
                      aSize_in = aSize,
                      aLabel_legend_in = aLabel_legend,
                      sTitle_in = 'Drainage area')

print('finished')