
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
nPoint = 6
dResolution_degree_in= 0.004166666666666670078
dLatitude_mean = 65.0
dResolution_flow_accumulation = degree_to_meter( dLatitude_mean ,dResolution_degree_in)


#In Canada, four principal tributaries feed the Yukon River: the Teslin River, the Pelly River, the White River and the Stewart River. In Alaska, the major tributaries are the Porcupine, Tanana and Koyukuk rivers.

aLabel_point= ['River mouth' , 'Porcupine River','Tanana River' , 'Koyukuk River','Stewart River', 'Pelly River' ]
aFlow_accumulation = np.array([9015232, 745848, 1209024, 940227, 536791, 510609]) *  dResolution_flow_accumulation * dResolution_flow_accumulation / 1.0E6

aFlow_accumulation = np.array([83282872, 6071528, 11381646, 8060644, 5110990, 5063349]) *  1.0E-2  # hectares  to km


aDrainage_area = np.full((nCase, nPoint), None, dtype=object)

aDrainage_area[0,:] = [ 744532568313.1809,   56272845384.47094,   103888327151.29428,  22509093797.867867,   62332992315.61696,   865736003.0348423]
aDrainage_area[1,:] = [ 794744915588.0951,   56561379327.308846,   118317170188.92104, 69258779666.25288,   49635494680.54951,      48481180603.46261]
aDrainage_area[2,:] = [ 806095575756.8727,   60986233110.11447,   112256993488.4198, 73972240687.48392,  50212639376.62372,       51847916813.36182]
aDrainage_area[3,:] = [ 817927258510.18,   62493248098.30496,   113411297256.824, 80128585553.19221,  50372956611.87265,       49827864117.934814]

aDrainage_area = aDrainage_area / 1.0E6

sFilename_out = sPath_parent + '/' + 'figures' + '/' + 'yukon'+ '/'+ 'scatterplot_drainage_area_multiple.png'

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
                      dMax_x_in = 6, 
                      dMin_y_in = 4, 
                      dMax_y_in = 6, 
                      dSpace_x_in = None, 
                      dSpace_y_in = None, 
                      sFormat_x_in =None,
                      sFormat_y_in =None,
                      sLabel_x_in ='HydroSHEDS data, $\mathrm{km}^{2}$',
                      sLabel_y_in = 'HexWatershed, $\mathrm{km}^{2}$' , 
                      aLabel_point_in = aLabel_point,
                      aColor_in=aColor,
                      aMarker_in=aMarker,
                      aSize_in = aSize,
                      aLabel_legend_in = aLabel_legend,
                      sTitle_in = 'Drainage area')

print('finished')