
from pathlib import Path
from os.path import realpath
import numpy as np
from osgeo import   ogr
import matplotlib as mpl

from pyearth.visual.scatter.scatter_plot_multiple_data import scatter_plot_multiple_data
from pyflowline.external.pyearth.gis.gdal.gdal_functions  import  degree_to_meter
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex


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


aDistance_niws =list()
aData_x=list()
aData_y  =list()
aLabel_legend =list()
nPoint = 6
dResolution_degree_in= 0.005
dLatitude_mean = np.mean(aLatitude)
dResolution_flow_accumulation = degree_to_meter( dLatitude_mean ,dResolution_degree_in)

aLabel_point= [ 'Amapa Para','Ilha Urucuricaia' , 'Santarem','Linha Rio Madeira', 'Manaus' ,'Leticia']
aFlow_accumulation = np.array([  108.720894,  143.61356, 517.4572, 1135.5105, 1274.3236, 2837.8542]) 
aDrainage_area = np.full((nCase, nPoint), None, dtype=object)

aDrainage_area[0,:] = [    294233.5625,   323656.90625, 676737.375,   1265326,      1422093.125, 2899839]
aDrainage_area[1,:] = [    306896.09375,  323911.03125, 730466.6875,  1402879.25,      1573705.375, 3091742.5]
aDrainage_area[2,:] = [    294067.875,    288789.28125, 696353.375,   1481350.875,       1609921.875, 3217139.75]
aDrainage_area[3,:] = [    313485.8125,   284133.90625, 736129.125,   1466595.75,       1639031, 3329568]

aDrainage_area = aDrainage_area / 1.0E3

sFilename_out = sPath_parent + '/' + 'figures' + '/' + 'scatterplot_travel_distance_multiple.png'

#5 cases + 1 obs

aColor = np.full(nCase, None, dtype=object)
aMarker= np.full(nCase, None, dtype=object)
aSize = np.full(nCase, mpl.rcParams['lines.markersize'] ** 2.5 , dtype=object)

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
                      iFlag_scientific_notation_x_in=0,
                      iFlag_scientific_notation_y_in=0,
                      iSize_x_in = None, 
                      iSize_y_in = None,  
                      iDPI_in = None ,
                      iFlag_log_x_in = 0,
                      iFlag_log_y_in = 0,
                      dMin_x_in = 0, 
                      dMax_x_in = 3000, 
                      dMin_y_in = 0, 
                      dMax_y_in = 3000, 
                      dSpace_x_in = None, 
                      dSpace_y_in = None, 
                      sFormat_x_in =None,
                      sFormat_y_in =None,
                      sLabel_x_in ='LBA-ECO data, km',
                      sLabel_y_in = 'HexWatershed, km' , 
                      aLabel_point_in = aLabel_point,
                      aColor_in=aColor,
                      aMarker_in=aMarker,
                      aSize_in = aSize,
                      aLabel_legend_in = aLabel_legend,
                      sTitle_in = 'Travel distance')

print('finished')