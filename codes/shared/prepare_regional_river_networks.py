from pyearth.toolbox.analysis.extract.clip_vector_by_polygon import clip_vector_by_polygon_file
from pyearth.toolbox.conversion.convert_vector_to_geojson import convert_vector_to_geojson
def prepare_regional_river_networks(sFilename_river_networks_in, 
                                    sFilename_watershed_boundary_in, 
                                    sFilename_river_networks_out):

    #generat a temporary file based on the output filename be replace the extension to .shp
    sFilename_temp = sFilename_river_networks_out.replace('.geojson', '.shp')


    clip_vector_by_polygon_file(sFilename_river_networks_in, sFilename_watershed_boundary_in, sFilename_temp)

    #now convert the shapefile to geojson using gdal python api
    convert_vector_to_geojson (sFilename_temp, sFilename_river_networks_out)



    return