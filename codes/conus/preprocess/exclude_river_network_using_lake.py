from pyearth.system.define_global_variables import *
from pyearth.toolbox.analysis.extract.exclude_vector_by_polygon_files import exclude_vector_by_polygon_files

sFilename_river_network = '/qfs/people/liao313/workspace/python/liao_2023_scidata_dggs/data/conus/dggrid12/river_networks.geojson'
aLake_boundary = list()

aLake_boundary.append('/qfs/people/liao313/data/hexwatershed/greatlakes/vector/hydrology/lake_erie_new.geojson')
aLake_boundary.append('/qfs/people/liao313/data/hexwatershed/greatlakes/vector/hydrology/lake_huron_new.geojson')
aLake_boundary.append('/qfs/people/liao313/data/hexwatershed/greatlakes/vector/hydrology/lake_michigan_new.geojson')
aLake_boundary.append('/qfs/people/liao313/data/hexwatershed/greatlakes/vector/hydrology/lake_ontario_new.geojson')
aLake_boundary.append('/qfs/people/liao313/data/hexwatershed/greatlakes/vector/hydrology/lake_superior_new.geojson')

sFilename_river_network_out = '/qfs/people/liao313/workspace/python/liao_2023_scidata_dggs/data/conus/dggrid12/river_networks_wo_greatlakes.geojson'

exclude_vector_by_polygon_files(sFilename_river_network,
                                aLake_boundary,
                                sFilename_river_network_out)