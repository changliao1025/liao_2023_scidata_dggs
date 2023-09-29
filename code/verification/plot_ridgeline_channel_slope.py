import os, sys, stat
from pathlib import Path
from os.path import realpath
import json
import numpy as np

from pyearth.visual.ridgeplot.ridgeplot_data_density import ridgeplot_data_density
from pyhexwatershed.pyhexwatershed_read_model_configuration_file import pyhexwatershed_read_model_configuration_file



# getting the data
sPath_parent = str(Path(__file__).parents[3]) # data is located two dir's up
print(sPath_parent)
sPath_data = realpath( sPath_parent +  '/data/susquehanna' )
sWorkspace_input =  str(Path(sPath_data)  /  'input')
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/susquehanna'
nCase  = 14
sDate='20230701'

# we define a dictionnary with months that we'll use later
case_dict = dict()

for i in range(1, nCase +1):
    case_dict[i] = 'Case ' +  "{:0d}".format(i) 

print(case_dict)

aData = list()
for iCase_index in range(1, nCase +1):
    aData_case=list()
    #read data
    sFilename_configuration_in = realpath( sPath_parent +  '/examples/susquehanna/pyhexwatershed_susquehanna_square.json' )
    if os.path.isfile(sFilename_configuration_in):
        pass
    else:
        print('This configuration does not exist: ', sFilename_configuration_in )
    oPyhexwatershed = pyhexwatershed_read_model_configuration_file(sFilename_configuration_in,\
      iCase_index_in=iCase_index,\
          sDate_in= sDate)  
    sWorkspace_output_hexwatershed = oPyhexwatershed.sWorkspace_output_hexwatershed
    for iWatershed in range(1, 2):#there is only one watershed in this study
        pBasin_hexwatershed=   oPyhexwatershed.aBasin[iWatershed-1]
        sWatershed = "{:04d}".format(iWatershed) 
        pBasin_hexwatershed = oPyhexwatershed.aBasin[0]
        sWorkspace_output_basin = pBasin_hexwatershed.sWorkspace_output_basin       
        sFilename_json = pBasin_hexwatershed.sFilename_watershed_json
        with open(sFilename_json) as json_file:
            data = json.load(json_file)  
            ncell = len(data)
            lID =0 
            for i in range(ncell):
                pcell = data[i]
                lCellID = int(pcell['lCellID'])
                iSegment = int(pcell['iSegment'])
                dSlope_between=float(pcell['dSlope_between'])  
                if iSegment >=1:
                    aData_case.append(dSlope_between)
                else:
                    pass
    
    aData.append(aData_case)
sFilename_out = sPath_parent + '/' + 'figures' + '/' + 'channel_slope.png'
sLabel_x = 'Channel slope (percent)'
ridgeplot_data_density(case_dict, aData, sFilename_out, dMin_x_in =0, dMax_x_in= 0.01, sLabel_x_in = sLabel_x)
