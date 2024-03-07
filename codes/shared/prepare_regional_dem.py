from pyearth.toolbox.analysis.extract.clip_raster_by_polygon import clip_raster_by_geojson
def prepare_regional_dem(sFilename_dem_in, sFilename_watershed_boundary_in, sFilename_dem_out):

    clip_raster_by_geojson(sFilename_dem_in, sFilename_watershed_boundary_in, sFilename_dem_out)

    return