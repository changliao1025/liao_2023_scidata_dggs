from pye3sm.mosart.general.unstructured.map.mosart_map_variable_unstructured import mosart_map_variable_unstructured
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_e3sm_configuration_file
from pye3sm.shared.pye3sm_read_configuration_file import pye3sm_read_case_configuration_file
from pye3sm.shared.case import pycase
from pye3sm.shared.e3sm import pye3sm
import os
import math
from pathlib import Path
from os.path import realpath
import importlib
import numpy as np
import netCDF4 as nc  # read netcdf
from osgeo import osr, gdal, ogr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import cartopy as cpl
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
from pyearth.visual.formatter import OOMFormatter

from pyearth.system.define_global_variables import *
from pyearth.visual.formatter import log_formatter, MathTextSciFormatter
from pyearth.toolbox.math.stat.remap import remap
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import convert_360_to_180
else:
    from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180


sRegion = 'amazon'
sMesh_type = 'mpas'

res = 'MOS_USRDAT'
res = 'MOS_USRDAT_MPAS'
compset = 'RMOSGPCC'
project = 'esmd'

iCase_index_e3sm = 1

dResolution_meter = 5000
sDate = '20230601'
# this one should be replace

sModel = 'e3sm'
sWorkspace_scratch = '/compyfs/liao313'




#
pDriver = ogr.GetDriverByName('GeoJSON')

sWorkspace_output_png = '/qfs/people/liao313/workspace/python/liao_2023_scidata_dggs/figures/amazon'
iFlag_mosart = 1
pSpatial_reference_gcs = osr.SpatialReference()
pSpatial_reference_gcs.ImportFromEPSG(4326)
iFlag_draw_mini_earth = 1

iFlag_draw_dem = 1
iFlag_draw_flow_direction = 1
iFlag_colorbar = 1
iFlag_legend = 1


sFont = "Times New Roman"
plt.rcParams["font.family"] = sFont
sColormap = 'terrain'
iDPI = 600
iSize_x = 14
iSize_y = 12
iFont_size = 12
cmap = mpl.cm.get_cmap(sColormap)
fig = plt.figure(dpi=iDPI)
fig.set_figwidth(iSize_x)
fig.set_figheight(iSize_y)

# define the box of amamzon
dLongitude_left = -82
dLongitude_right = -50
dLatitude_bot = -20
dLatitude_top = 7

pRing = ogr.Geometry(ogr.wkbLinearRing)
pRing.AddPoint(dLongitude_left, dLatitude_top)
pRing.AddPoint(dLongitude_right, dLatitude_top)
pRing.AddPoint(dLongitude_right, dLatitude_bot)
pRing.AddPoint(dLongitude_left, dLatitude_bot)
pRing.AddPoint(dLongitude_left, dLatitude_top)
pBoundary = ogr.Geometry(ogr.wkbPolygon)
pBoundary.AddGeometry(pRing)
pBoundary_wkt = pBoundary.ExportToWkt()
pBoundary = ogr.CreateGeometryFromWkt(pBoundary_wkt)
cmap = mpl.cm.get_cmap(sColormap)
pProjection_map = cpl.crs.Orthographic(central_longitude=0.50*(
    dLongitude_left+dLongitude_right),  central_latitude=0.50*(dLatitude_bot+dLatitude_top), globe=None)


ax_globe = fig.add_axes([0.8, 0.7, 0.1, 0.1], projection=pProjection_map)
# now draw the data

ax_cb = fig.add_axes([0.845, 0.11, 0.015, 0.57])    

# plot a mini earth
ax_globe.set_global()
if iFlag_draw_mini_earth == 1:    

    #draw a dggrid mesh

    sFilename_dggrid  ='/compyfs/liao313/04model/pyflowline/global/pyhexwatershed20231201001/dggrid.geojson'
    pDataset = pDriver.Open(sFilename_dggrid, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    aPolygon_dggrid = list()
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        # get attribute
        if sGeometry_type == 'POLYGON':            
            aCoords_gcs = get_geometry_coordinates(pGeometry_in)
            aPolygon_dggrid.append(aCoords_gcs[:, 0:2])
        else:
            if sGeometry_type == 'MULTIPOLYGON':
                nLine = pGeometry_in.GetGeometryCount()
                for i in range(nLine):
                    pBoundary_ogr = pGeometry_in.GetGeometryRef(i)    
                    aCoords_gcs = get_geometry_coordinates(pBoundary_ogr)
                    aPolygon_dggrid.append(aCoords_gcs[:, 0:2])
            
                
    
    aPatch = [Polygon(poly, closed=True) for poly in aPolygon_dggrid]
    pPC = PatchCollection(aPatch, cmap=cmap, alpha=1.0, edgecolor='blue',
                          facecolor='none', linewidths=0.1,
                          transform=cpl.crs.Geodetic())
    
    ax_globe.add_collection(pPC)

    # draw a boundary of amazon
    sFilename_boundary = '/qfs/people/liao313/data/hexwatershed/amazon/vector/hydrology/mesh_boundary_buffer.geojson'
    pDataset = pDriver.Open(sFilename_boundary, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        # get attribute
        if sGeometry_type == 'MULTIPOLYGON':
            nLine = pGeometry_in.GetGeometryCount()
            for i in range(nLine):
                pBoundary_ogr = pGeometry_in.GetGeometryRef(i)    
                aCoords_gcs = get_geometry_coordinates(pBoundary_ogr)    
                polygon = mpatches.Polygon(aCoords_gcs[:, 0:2], closed=True, linewidth=1.0,
                                           alpha=0.8, edgecolor='red', facecolor='none',
                                           transform=cpl.crs.Geodetic())
                ax_globe.add_patch(polygon)
                

    ax_globe.coastlines(color='black', linewidth=1)

    
    # draw the openstreet mep (optional)

    gl = ax_globe.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                         linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 5, 'color': 'k', 'rotation': 0, 'ha': 'right'}
    gl.ylabel_style = {'size': 5, 'color': 'k',
                   'rotation': 90, 'weight': 'normal'}    


dValue_max = 5500 
dValue_min = 0.0


dLat_min = 90
dLat_max = -90
dLon_min = 180
dLon_max = -180
if iFlag_draw_dem == 1:

    sFilename_in = '/compyfs/liao313/04model/pyhexwatershed/amazon/pyhexwatershed20230801013/hexwatershed/00000001/variable_polygon.geojson'
       

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)

    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

    aValue = list()
    sVariable = 'elevation'
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        dValue = float(pFeature.GetField(sVariable))
        aValue.append(dValue)
        aCoords_gcs = get_geometry_coordinates(pGeometry_in)
        dLon_max = np.max([dLon_max, np.max(aCoords_gcs[:, 0])])
        dLon_min = np.min([dLon_min, np.min(aCoords_gcs[:, 0])])
        dLat_max = np.max([dLat_max, np.max(aCoords_gcs[:, 1])])
        dLat_min = np.min([dLat_min, np.min(aCoords_gcs[:, 1])])


    aValue = np.array(aValue)
    dValue_max = 5500 #np.max(aValue)   
    dValue_min = 0.0

    dLongitude_left = dLon_min
    dLongitude_right = dLon_max
    dLatitude_bot = dLat_min
    dLatitude_top = dLat_max

    pProjection_map = cpl.crs.Orthographic(central_longitude=0.50*(
    dLongitude_left+dLongitude_right),  central_latitude=0.50*(dLatitude_bot+dLatitude_top), globe=None)
    ax_amazon = fig.add_axes([0.11, 0.1, 0.68, 0.7], projection=pProjection_map)
    ax_amazon.set_global()

    aPolygon_dem = list()
    aColor = list()
    aValue = list()
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        # get attribute
        dValue = float(pFeature.GetField(sVariable))
        
        aValue.append(dValue)
        if dValue != -9999:
            if dValue > dValue_max:
                dValue = dValue_max
            if dValue < dValue_min:
                dValue = dValue_min

            dValue0 = dValue
            if dValue0 <= 0.0:
                iColor_index = 0
            else:
                dValue = dValue0
                iColor_index = int((dValue - dValue_min) /
                                   (dValue_max - dValue_min) * 255)
                           

            if sGeometry_type == 'POLYGON':
                aCoords_gcs = get_geometry_coordinates(pGeometry_in)
                aColor.append(cmap(iColor_index))
                aPolygon_dem.append(aCoords_gcs[:, 0:2])

        else:
            pass

    aValue = np.array(aValue)
    print(np.max(aValue))
    aPatch = [Polygon(poly, closed=True) for poly in aPolygon_dem]
    pPC = PatchCollection(aPatch, cmap=cmap, alpha=0.5, edgecolor='none',
                          facecolor=aColor, linewidths=0.1,
                          transform=cpl.crs.Geodetic())
    ax_amazon.add_collection(pPC)

dField_min = 0.0
dField_max = 5.0E12
iThickness_min = 0.1
iThickness_max = 1.5
aThickness=list()
if iFlag_draw_flow_direction == 1:
    # draw the flow direction geojson

    sFilename_flow_direction = '/compyfs/liao313/04model/pyhexwatershed/amazon/pyhexwatershed20230801013/hexwatershed/00000001/flow_direction.geojson'
    pDataset = pDriver.Open(sFilename_flow_direction, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)

    aPolyline_flow_direction = list()

    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        dField = pFeature.GetField('drainage_area')
        if sGeometry_type == 'LINESTRING':
            aCoords_gcs = get_geometry_coordinates(pGeometry_in)
            nvertex = len(aCoords_gcs)
            codes = np.full(nvertex, mpl.path.Path.LINETO, dtype=int)
            codes[0] = mpl.path.Path.MOVETO
            path = mpl.path.Path(aCoords_gcs, codes)
            x, y = zip(*path.vertices)
            iThickness = remap( dField, dField_min, dField_max, iThickness_min, iThickness_max )
            aThickness.append(iThickness)
            # aPolyline_flow_direction.append(Line2D(x, y))  # Use Line2D
            aPolyline_flow_direction.append(list(zip(x, y)))

    pLC = LineCollection(aPolyline_flow_direction,  alpha=0.8, edgecolor='black',
                         facecolor='black', linewidths=aThickness, transform=cpl.crs.Geodetic())

    ax_amazon.add_collection(pLC)

marginx = 0.0 #(dLongitude_right - dLongitude_left) / 20
marginy = 0.0 #(dLatitude_top - dLatitude_bot) / 20
aExtent = [dLongitude_left - marginx, dLongitude_right + marginx,
           dLatitude_bot - marginy, dLatitude_top + marginy]
ax_amazon.set_extent(aExtent)
ax_amazon.coastlines(color='black', linewidth=1)
# ax_amazon.set_title(sTitle)

if iFlag_legend == 1:
    aLegend_in = ['Flow direction']
    nlegend = len(aLegend_in)
    for i in range(nlegend):
        sText = aLegend_in[i]
        dLocation = 0.06 + i * 0.04
        ax_amazon.text(0.03, dLocation, sText,
                       verticalalignment='top', horizontalalignment='left',
                       transform=ax_amazon.transAxes,
                       color='black', fontsize=iFont_size)
        pass

if iFlag_colorbar == 1:    

    formatter = OOMFormatter(fformat="%1.1e")
    # formatter = MathTextSciFormatter("%1.2e")
    cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                   cmap=cmap,
                                   norm=mpl.colors.Normalize(
                                       dValue_min, dValue_max),  # vmax and vmin
                                   extend='max', format=formatter)

    sUnit = 'Elevation: m'
    cb.ax.get_yaxis().set_ticks_position('right')
    cb.ax.get_yaxis().labelpad = 5
    cb.ax.set_ylabel(sUnit, rotation=90, fontsize=iFont_size-2)
    cb.ax.get_yaxis().set_label_position('left')
    cb.ax.tick_params(labelsize=iFont_size-2)

gl = ax_amazon.gridlines(crs=cpl.crs.PlateCarree(), draw_labels=True,
                         linewidth=1, color='gray', alpha=0.5, linestyle='--')
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation': 0, 'ha': 'right'}
gl.ylabel_style = {'size': 10, 'color': 'k',
                   'rotation': 90, 'weight': 'normal'}    

sFilename_output_in = os.path.join(sWorkspace_output_png, 'flow_direction.pdf')
sDirname = os.path.dirname(sFilename_output_in)
pDataset = pLayer = pFeature = None
if sFilename_output_in is None:
    plt.show()
else:
    sFilename = os.path.basename(sFilename_output_in)
    sFilename_out = os.path.join(sDirname, sFilename)
    sExtension = os.path.splitext(sFilename)[1]
    if sExtension == '.png':
        plt.savefig(sFilename_out, bbox_inches='tight')
    else:
        if sExtension == '.pdf':
            plt.savefig(sFilename_out, bbox_inches='tight')
        else:
            plt.savefig(sFilename_out, bbox_inches='tight', format='ps')
    plt.close('all')
    plt.clf()
