# devtools::install_github('ranghetti/sen2r')

library(sen2r)

username <- ""

password <- ""

write_scihub_login(username = username, password = password)

shapefiles <- 'vectorFiles/resedinte/'

in.files <- list.files(shapefiles, pattern = ".shp")

for(i in in.files){

in.file <- paste0(shapefiles, i)

try({
  
sen2r(param_list = NULL, gui = FALSE, preprocess = FALSE, s2_levels = 'l2a',
      sel_sensor = c('s2a', 's2b'), online = TRUE, 
      overwrite_safe = FALSE, rm_safe = FALSE, step_atmcorr = 'auto', 
      max_cloud_safe = 1, timewindow = c(as.Date('2022-06-01'), as.Date('2022-07-30')), 
      extent = in.file, 
      timeperiod = 'full',
      clip_on_extent = TRUE, 
      extent_as_mask = TRUE,res_s2 = '10m', index_source = 'BOA', 
      resampling = 'bilinear', path_l2a = 'data/', 
      path_out = 'data/', parallel = TRUE)
}, silent = TRUE)
}


