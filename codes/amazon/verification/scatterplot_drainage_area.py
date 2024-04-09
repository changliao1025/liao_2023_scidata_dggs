
from pathlib import Path
from os.path import realpath

import numpy as np
from osgeo import  ogr
import matplotlib as mpl

from pyearth.visual.scatter.scatter_plot_multiple_data import scatter_plot_multiple_data
from pyearth.gis.spatialref.conversion_between_degree_and_meter  import  degree_to_meter
from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex

sPath_parent = str(Path(__file__).parents[3]) # data is located two dir's up
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
nPoint = 7
dResolution_degree_in= 0.005 #0.00833333 
dLatitude_mean = np.mean(aLatitude)
dResolution_flow_accumulation = degree_to_meter( dLatitude_mean ,dResolution_degree_in)

aLabel_point= ['River mouth' , 'Amapa Para','Ilha Urucuricaia' , 'Santarem','Linha Rio Madeira', 'Manaus' ,'Leticia']
aFlow_accumulation = np.array([19653284, 167545,  140895, 1633301, 4580531, 2322894, 318989]) *  dResolution_flow_accumulation * dResolution_flow_accumulation / 1.0E6
aDrainage_area = np.full((nCase, nPoint), None, dtype=object)

aDrainage_area[0,:] = [ 5864467464192,   58870071296,   48481243136, 472691146752,   1344485982208,      670078992384, 83110592512]
aDrainage_area[1,:] = [ 5840830988288,   56849965056,   42998198272, 478462967808,   1458763333632,      679890255872, 106485432320]
aDrainage_area[2,:] = [ 5839470460928,   63006318592,   34917969920, 479425757184,  1350647152640,       681910730752, 107158781952]
aDrainage_area[3,:] = [ 5844697088000,   63871963136,   33154451456, 487281885184,  1357410861056,       693967388672, 109371195392]

aDrainage_area = aDrainage_area / 1.0E6

sFilename_out = sPath_parent + '/' + 'figures' + '/' + 'amazon'+ '/'+ 'scatterplot_drainage_area_multiple.png'

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

    #conver to log10
    aDatax0 = np.log10(aDatax0)
    aDatay0 = np.log10(aDatay0)
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
                      iFlag_log_x_in = 1,
                      iFlag_log_y_in = 1,
                      dMin_x_in = 4, 
                      dMax_x_in = 7, 
                      dMin_y_in = 4, 
                      dMax_y_in = 7, 
                      dSpace_x_in = None, 
                      dSpace_y_in = None, 
                      sFormat_x_in =None,
                      sFormat_y_in =None,
                      sLabel_x_in ='LBA-ECO data, $\mathrm{km}^{2}$',
                      sLabel_y_in = 'HexWatershed, $\mathrm{km}^{2}$' , 
                      aLabel_point_in = aLabel_point,
                      aColor_in=aColor,
                      aMarker_in=aMarker,
                      aSize_in = aSize,
                      aLabel_legend_in = aLabel_legend,
                      sTitle_in = 'Drainage area')

print('finished')