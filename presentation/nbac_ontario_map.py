import ee
import geemap

# Authenticate and initialize Earth Engine
ee.Authenticate()
ee.Initialize(project='ee-jandrewgoldman')

# Load the NBAC raster
nbac_raster = ee.Image("projects/sat-io/open-datasets/CA_FOREST/NBAC/NBAC_MRB_1972_to_2023")

# Load Ontario boundary (assuming the image has a 'year' band with burn year values)
ontario = ee.FeatureCollection("FAO/GAUL/2015/level1").filter(ee.Filter.eq('ADM1_NAME', 'Ontario'))

# Clip the raster to Ontario
nbac_ontario = nbac_raster.clip(ontario)

# Define visualization parameters
vis_params = {
    'min': 1972,
    'max': 2023,
    'palette': ['blue', 'green', 'yellow', 'orange', 'red']  # Color scale for years
}

# Create a map
Map = geemap.Map()
Map.centerObject(ontario, 6)
Map.addLayer(nbac_ontario, vis_params, 'Fires in Ontario (1972-2023)')

# Add a colorbar legend
Map.add_colorbar(vis_params, label='Burn Year')

# Display the map
Map