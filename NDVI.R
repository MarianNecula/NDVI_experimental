################################ NDVI computation ##############################
########################### libraries and input files ##########################
# vector and raster
library(raster)
library(sf)
library(sen2r)
library(RStoolbox)

# utils
library(xml2)
library(dplyr)

# plotting
library(ggplot2)
library(dplyr)
library(grid)
library(patchwork)
library(tidytext)

# util functions 
# raster_percent
# raster percent for binary underlying values
# @param rasterFile = raster file
# @param shapeFile = shapeFile extent of raster
raster_percent <- function(rasterFile = NULL, shapeFile = NULL){
  
  resedintaArea <- st_area(shapeFile)
  
  # km^2
  resedintaArea <- round(resedintaArea/10^6, 2)
  
  # percentage dichotomous underlying raster variable
  v <- sum(rasterFile[] == 1L, na.rm = TRUE)
  
  nonv <- sum(rasterFile[] == 0L, na.rm = TRUE)
  
  vPercent <- round(v/(v + nonv), 2) 
  
  vKm2 <- round(resedintaArea * vPercent, 2)
  
  return(c(resedintaName= shapeFile$DENUMIRE,
              resedintaArea = resedintaArea,
              vPercent = vPercent,
              vKm2 = vKm2)
         )
}

# s2filter_cloud_nodata
# exclude cloud and nodata datasets
# @inputFile - vector containing string of SAFE file name 
# @cloud percent numeric between 0 and 100
# @nodataPercent numeric between 0 and 100
s2filer_cloud_nodata <- function(inputFile = NULL,
                                 cloudPercent = NULL, 
                                 nodataPercent = NULL){
  
  SAFEDiscard <- lapply(inputFile, function(SAFE){
  
    bandsPath <- paste0("data/", SAFE)
    
    myXML <- read_xml(paste0(bandsPath, "/MTD_MSIL2A.xml"))
    
    valueNODATA <- xml_find_all(myXML, ".//NODATA_PIXEL_PERCENTAGE") %>% 
    xml_text() %>% 
    as.numeric()
    
    valueCloud <- xml_find_all(myXML, ".//CLOUDY_PIXEL_OVER_LAND_PERCENTAGE") %>% 
    xml_text() %>% 
    as.numeric()
  
    if(valueNODATA > nodataPercent || valueCloud > cloudPercent){
      return(c(SAFE, "NOTOK", valueCloud, valueNODATA))
    } else {
      return(c(SAFE, 'OK', valueCloud, valueNODATA))}
    }
  )
  
  SAFEDiscard <- as.data.frame(do.call('rbind', SAFEDiscard), stringsAsFactors = FALSE)
  
  return(SAFEDiscard)
}


################################################################################
##################################### input files ##############################

# shapefiles
shapefilesPath <- "vectorFiles/resedinte/"

shapeFiles <- list.files(path = shapefilesPath, pattern = ".shp")

shapeFiles <- lapply(paste0(shapefilesPath, shapeFiles), st_read)

shapeFiles <- do.call('rbind', shapeFiles)

# add S2 tile intersection
shapeFiles$tiles <- NA

for(i in 1:nrow(shapeFiles)){
  
  shapeFiles[i, ]$tiles <- tiles_intersects(shapeFiles[i, ])
  
}

shapeFiles$tiles <- paste0("T", shapeFiles$tiles)

# bands path
bandsPath <- "data/"

# add S2 dataset names
# SAFE files
safeFiles <- list.files(pattern = ".SAFE$", path = "data/")

# extract TILE names from SAFE file name
tiles <-
  sapply(safeFiles, function(safeFiles)
    strsplit(safeFiles, "_")[[1]][6])

safeTiles <- data.frame(SAFE = safeFiles, tiles = tiles)

# create input dataset for sentinel
shapeFiles <-
  left_join(shapeFiles, safeTiles, by = c('tiles' = 'tiles'))


# filter S2 SAFE files
SAFEDiscard <- s2filer_cloud_nodata(shapeFiles$SAFE, cloudPercent = 5, nodataPercent = 10)

# SAFEDiscard <- SAFEDiscard[SAFEDiscard$V2 == 'OK',]

shapeFiles <- left_join(shapeFiles, SAFEDiscard, by = c('SAFE' = 'V1'))

shapeFiles$V3 <- as.numeric(shapeFiles$V3)

shapeFiles$V4 <- as.numeric(shapeFiles$V4)

shapeFiles$sumQ <- shapeFiles$V3 + shapeFiles$V4

# keep only one observation per resedinta
shapeFiles <- shapeFiles %>% 
  group_by(DENUMIRE) %>% filter(sumQ == min(sumQ)) %>%
  ungroup()

shapeFiles <- shapeFiles[!duplicated(shapeFiles$DENUMIRE),]

# MODIS input files
MODIS <- raster("data/MODIS_NDVI.tif")

# divide by 10^8 to get NDVI between -1 and 1
MODIS <- MODIS/10^8

# NDVI threshold
threshold <- read.csv("data/parametru.csv")

tholdMODIS <- mean(threshold$modis)

tholdS2 <- mean(threshold$s2)


################################################################################
######################### computation of NDVI and plots ########################
results <- lapply(unique(shapeFiles$DENUMIRE), function(DENUMIRE){
  
  resedinta <- shapeFiles[shapeFiles$DENUMIRE == DENUMIRE,]
  
  # resedinta area for output before CRS change
  rArea <- st_area(resedinta)
  
  rArea <- round(rArea/10^6, 2)
  
  resedintaMODIS <- st_transform(resedinta, st_crs(MODIS))
  
  ############################ MODIS ##########################################
  # MODIS NDVI
  modisNDVI <- crop(MODIS, resedintaMODIS)
  
  modisNDVI <- projectRaster(modisNDVI, crs = 32634)
  
  resedintaMODIS <- st_transform(resedinta, 32634)
  
  modisNDVI <- mask(modisNDVI, resedintaMODIS)
  
  # set threshold 
  modisNDVI[modisNDVI <= tholdMODIS] <- 0L
  
  modisNDVI[modisNDVI > tholdMODIS] <- 1L
  
  vegetationValuesMODIS <- raster_percent(rasterFile = modisNDVI, 
                                          shapeFile = resedintaMODIS)
  
  ################################### SENTINEL 2 ###############################
  # Sentinel 2 SAFE FILE 
  bandPATH <- paste0(bandsPath, resedinta$SAFE)
  
  lf <- list.files(bandPATH, pattern = "B(0[2348]_10m).jp2$", recursive = TRUE)
  
  s2 <- stack(paste0(bandPATH, "/", lf))
  
  resedintaS2 <- st_transform(resedinta, st_crs(s2))
  
  s2 <- crop(s2, resedintaS2)
  
  s2 <- mask(s2, resedintaS2)
  
  names(s2) <- c("B02", "B03", "B04", "B08")
  
  S2NDVI <- (s2$B08 - s2$B04) / (s2$B08 + s2$B04)
  
  S2NDVI[S2NDVI <= tholdS2] <- 0L
  
  S2NDVI[S2NDVI > tholdS2] <- 1L
  
  vegetationValuesS2 <- raster_percent(rasterFile = S2NDVI, 
                                       shapeFile = resedintaS2)
  
  ################################# Plots ######################################
  plotMODIS <- ggR(modisNDVI, geom_raster = TRUE, stretch = 'none', forceCat = TRUE) +
    scale_fill_manual(
      name = "NDVI",
      values = c("grey", "green"),
      labels = c('Other', 'Vegetation'),
      na.value = 'transparent',
      na.translate = FALSE) + 
     geom_sf(data = resedintaS2, alpha = 0.01) +
     xlab("Longitudine") +
    ylab("Latitudine") +
    ggtitle("MODIS NDVI") +
    theme_minimal() +
    coord_sf(xlim = st_bbox(resedintaMODIS)[c(1, 3)],
             ylim = st_bbox(resedintaMODIS)[c(2, 4)],
             expand = TRUE)


  plotS2 <-  ggR(S2NDVI, geom_raster = TRUE, 
                 stretch = 'none', 
                 forceCat = TRUE) +
    scale_fill_manual(
      name = "NDVI",
      values = c("grey", "green"),
      labels = c('Other', 'Vegetation'),
      na.value = 'transparent',
      na.translate = FALSE) +
      geom_sf(data = resedintaS2, alpha = 0.01) +
      xlab("Longitudine") +
      ylab("Latitudine") +
      ggtitle("SENTINEL 2 NDVI") +
      theme_minimal() + 
    coord_sf(xlim = st_bbox(resedintaS2)[c(1, 3)],
             ylim = st_bbox(resedintaS2)[c(2, 4)], 
             expand = TRUE)
  
  S2RGB <- ggRGB(
    s2,
    r = 3,
    g = 2,
    b = 1,
    stretch = "lin") +
    geom_sf(data = resedintaS2, alpha = 0.01) +
    ggtitle("SENTINEL 2 RGB") +
    xlab("Longitudine") +
    ylab("Latitudine") +
    theme_minimal() +
    coord_sf(xlim = st_bbox(resedintaS2)[c(1, 3)],
             ylim = st_bbox(resedintaS2)[c(2, 4)],
             expand = TRUE)
  
  t <- textGrob(paste0("Resedinta: ", vegetationValuesMODIS[1], "\n",
                       "Suprafata(Km^2):", vegetationValuesMODIS[2], "\n",
                       "MODIS Suprafata vegetatie(km^2): ", vegetationValuesMODIS[4], "\n",
                       "SENTINEL 2 Suprafata vegetatie(km^2): ", vegetationValuesS2[4], "\n",
                       "MODIS Vegetatie(%): ", as.numeric(vegetationValuesMODIS[3]) * 100, "\n",
                       "SENTINEL 2 Vegetatie(%): ", as.numeric(vegetationValuesS2[3]) * 100, "\n")
                )

  plots <- plotMODIS + plotS2 + S2RGB + t 
  
  ggsave(filename = paste0("plots/",DENUMIRE, ".png"), plot = plots, device = "png", width = 16, height = 9, units = "in")  
  
  return(c(Nume = resedinta$DENUMIRE, 
           Suprafata = rArea, 
           MODISSupVeg = vegetationValuesMODIS[3],
           SENTINEL2SupVeg = vegetationValuesS2[3],
           MODISVPercent = vegetationValuesMODIS[4],
           SENTINEL2VPercent = vegetationValuesS2[4]))
  })


tableVegetation <- do.call('rbind', results)

tableVegetation <- data.frame(tableVegetation)

colnames(tableVegetation) <- c("Nume", "Suprafata (Km2)", 
                               "Suprafata vegetatie (%) - MODIS", 
                               "Suprafata vegetatie (%) - Sentinel 2", 
                               "Suprafata vegetatie (km2) - MODIS",
                               "Suprafata vegetatie (km2) - Sentinel 2")

tableVegetation$`Suprafata vegetatie (%) - MODIS` <- as.numeric(tableVegetation$`Suprafata vegetatie (%) - MODIS`) * 100

tableVegetation$`Suprafata vegetatie (%) - Sentinel 2` <- as.numeric(tableVegetation$`Suprafata vegetatie (%) - Sentinel 2`) * 100

write.csv(tableVegetation, "tabelVegetatie.csv", row.names = FALSE)

tableVegetation <- reshape2::melt(tableVegetation, id.vars = c("Nume", "Suprafata (Km2)"), 
                                  measure.vars = c("Suprafata vegetatie (%) - MODIS", "Suprafata vegetatie (%) - Sentinel 2"))

tableVegetation$variable <- as.factor(tableVegetation$variable)

tableVegetation$value <- round(as.numeric(tableVegetation$value), 2)

rankingPlot <- ggplot(data = tableVegetation) + 
  geom_col(aes(x = forcats::fct_rev(Nume), y =  value), fill = 'darkgreen') + 
  facet_wrap(~variable, scales = 'free') +
  coord_flip() +
  theme_minimal() +
  xlab("Resedinta") +
  ylab("Procent") +
  theme(text = element_text(size = 16))

rankingPlot

ggsave(filename = "ranking.png", plot = rankingPlot, device = 'png', width = 16, height = 9, unit = 'in')

