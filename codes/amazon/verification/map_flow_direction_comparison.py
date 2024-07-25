import os, stat
from pathlib import Path
from os.path import realpath
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import matplotlib.cm as cm
from PIL import Image

from pyhexwatershed.configuration.read_configuration_file import pyhexwatershed_read_configuration_file
plt.rcParams["font.family"] = "Times New Roman"
class OOMFormatter(mpl.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1e", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        mpl.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$\mathdefault{%s}$' % self.format


# Open the images
sDate = '20240101'
nrow = 2
ncolumn = 2
iCase_start = 10
iCase_end = 13
iFlag_colorbar = 1
iFlag_scientific_notation_colorbar =0

#===========================
#setup workspace path
#===========================
sPath_parent = str(Path(__file__).parents[3]) # data is located two dir's up
sPath_data = realpath( sPath_parent +  '/data/amazon' )
sWorkspace_input =  str(Path(sPath_data)  /  'input')
sWorkspace_output = '/compyfs/liao313/04model/pyhexwatershed/amazon'

sFilename_configuration_in = realpath( sPath_parent +  '/data/amazon/input/pyhexwatershed_amazon_dggrid.json' )
if os.path.isfile(sFilename_configuration_in):
    pass
else:
    print('This configuration does not exist: ', sFilename_configuration_in )
    exit()

sFilename_png = 'flow_direction_w_observation_manaus.png'

aImage = list()
for iCase in range(iCase_start, iCase_end + 1):
    iCase_index = iCase
    oPyhexwatershed = pyhexwatershed_read_configuration_file(sFilename_configuration_in,
            iCase_index_in=iCase_index,
                sDate_in= sDate)
    pBasin_hexwatershed = oPyhexwatershed.aBasin[0]
    sWorkspace_output_basin = pBasin_hexwatershed.sWorkspace_output_basin
    print(sWorkspace_output_basin)
    sFilename = os.path.join(  sWorkspace_output_basin, sFilename_png )

    image_dummy = Image.open(sFilename)
    aImage.append(image_dummy)


# Create a figure and subplots


fig, axs = plt.subplots(nrow, ncolumn, figsize=(16, 12), gridspec_kw={'width_ratios': [5,5]},dpi=300)
plt.subplots_adjust(hspace=0.0, wspace=0.0, top=0.96)  # Adjust spacing here
# Plot each image on a subplot

for irow in range(1, nrow+1):
    for icolumn in range(1, ncolumn+1):

        iCase_index = (irow-1)*ncolumn + icolumn
        ax_dummy = axs[irow-1, icolumn-1]

        if iCase_index < 15:
            ax_dummy.imshow(aImage[iCase_index-1])

        ax_dummy.axis('off')


# Add a common title above the subplots
#anchored_text = AnchoredText("Surface elevation", loc='upper center', frameon=False, prop=dict(fontsize=16))

# Save the merged image with titles

fig.suptitle("Flow direction with observation", fontsize=16)
sFilename_out = '/qfs/people/liao313/workspace/python/liao_2023_scidata_dggs/figures/amazon/flow_direction_comparison_manaus.png'
#plt.show()
plt.savefig(sFilename_out,  bbox_inches='tight')

