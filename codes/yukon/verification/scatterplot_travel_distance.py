
from pathlib import Path
from os.path import realpath
import numpy as np
from osgeo import   ogr
import matplotlib as mpl

from pyearth.visual.scatter.scatter_plot_multiple_data import scatter_plot_multiple_data

from pyearth.visual.color.create_diverge_rgb_color_hex import create_diverge_rgb_color_hex


sPath_parent = str(Path(__file__).parents[3]) # data is located two dir's up
print(sPath_parent)
sPath_data = realpath( sPath_parent +  '/data/yukon' )
sWorkspace_input =  str(Path(sPath_data)  /  'input')
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/yukon'
aResolution_index = [10, 11, 12, 13]
nCase = len(aResolution_index)
aCase_index = aResolution_index
sDate='20230801'
sMesh_type = 'dggrid'
sFilename_configuration_in = realpath( sPath_parent +  '/examples/yukon/pyhexwatershed_yukon_dggrid.json' )

aDistance_niws =list()
aData_x=list()
aData_y  =list()
aLabel_legend =list()
nPoint = 5
dResolution_degree_in= 0.005


aLabel_point= [ 'Porcupine River','Tanana River' , 'Koyukuk River','Stewart River', 'Pelly River' ]
aFlow_distance = np.array([  1850200,  1157880, 809110, 2258190, 2434300 ]) 
aTravel_distance = np.full((nCase, nPoint), None, dtype=object)

aTravel_distance[0,:] = [    1885185.125,   1293173.625, 923087.25,   2268588.5,       2418954]
aTravel_distance[1,:] = [    1972409.375,  1240093.875, 877802.3125,  2276337,       2416650.5]
aTravel_distance[2,:] = [    2045560.625,    1314244.25, 931740.25,   2462972,       2658606]
aTravel_distance[3,:] = [    2051919.125,   1285394.5, 896079.1875,   2500803.5,      2684039.75]

aTravel_distance = aTravel_distance / 1.0E3

sFilename_out = sPath_parent + '/' + 'figures' + '/yukon/' + 'scatterplot_travel_distance_multiple.png'

#5 cases + 1 obs

aColor = np.full(nCase, None, dtype=object)
aMarker= np.full(nCase, None, dtype=object)
aSize = np.full(nCase, mpl.rcParams['lines.markersize'] ** 2.5 , dtype=object)

nmesh=nCase
aColor= create_diverge_rgb_color_hex(nmesh)
aMarker = [ '.','o','+','x','^']

aFlow_distance = np.array(aFlow_distance) / 1000
aData_x = list()
aData_y = list()

for i in range(nCase):
    iCase_index = aCase_index[i]  
    aLabel_legend.append(   'ISEA3H Level ' +  "{:0d}".format(iCase_index)   )
    aDatax0 = list()
    aDatay0 = list()
    for j in range(nPoint):
        aDatax0.append(aFlow_distance[j])
        aDatay0.append(aTravel_distance[i,j])
    
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
                      sLabel_x_in ='HydroSHEDS data, km',
                      sLabel_y_in = 'HexWatershed, km' , 
                      aLabel_point_in = aLabel_point,
                      aColor_in=aColor,
                      aMarker_in=aMarker,
                      aSize_in = aSize,
                      aLabel_legend_in = aLabel_legend,
                      sTitle_in = 'Travel distance')

print('finished')