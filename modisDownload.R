# packages
devtools::install_github("16EAGLE/getSpatialData")
library(sf)
library(getSpatialData)
library(gdalUtilities)
library(ggplot2)

# download data
# romania
romania <- sf::st_read("vectorFiles/roSHP/ro_judete_poligon.shp")

romania <- sf::st_transform(romania, crs = 4326)

username_USGS <- ""

password_USGS <- ""

username_earthdata <- ""

password_earthdata <- ""

login_USGS(username = username_USGS, password = password_USGS)

login_earthdata(username = username_earthdata, password = password_earthdata)

modis <- getSpatialData::get_products()

# NDVI and EVI Vegetation Indices at 16 days 250m = lpcs_modis_mod13q1
# source https://lpdaac.usgs.gov/products/mod13q1v006/
modis <- modis[modis %in% c("lpcs_modis_mod13q1")]

set_aoi(aoi = as(sf::st_geometry(romania), "Spatial"))

modis <- getSpatialData::get_records(time_range = c("2017-01-01", "2022-06-30"),
                                     products = modis, simplify_cols = TRUE)

# select a date
# Observation MODIS tile 194 covers most of Romania
# Modis tile 204 covers the rest (Black Sea region)
# source https://modis-land.gsfc.nasa.gov/MODLAND_grid.html
modis <- getSpatialData::check_availability(modis)

modis <- getSpatialData::get_data(modis, dir_out = "data/")

modisFiles <- list.files(pattern = ".hdf", recursive = TRUE)

lapply(modisFiles, function(x){

  ndvi <- gdalUtils::get_subdatasets(x)

  gdalUtils::gdal_translate(src_dataset = ndvi[1], dst_dataset = paste0(gsub("\\.", "", x), ".tif"))
})

# mosaic and tiff
modisFiles <- list.files(".", pattern = ".tif", recursive = TRUE)

gdalbuildvrt(modisFiles, "MODIS_NDVI.vrt")

gdal_translate("MODIS_NDVI.vrt", "MODIS_NDVI.tif")

modis <- stack("MODIS_NDVI.tif")




