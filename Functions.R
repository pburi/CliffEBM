################################################################################
# FUNCTIONS FILE
#
# 
# 2023/05/19
#
#
# Pascal Buri | High Mountain Glaciers and Hydrology | 
#  Swiss Federal Institute for Forest, Snow and Landscape Research, WSL |
#  Zürcherstrasse 111, 8903 Birmensdorf | pascal.buri@wsl.ch
#
#


################################################################################
# 'CLIFF_CLOSEHORIZON': CALCULATES CLOSE HORIZON FOR EACH POINT          
#
# v:                  vector of structure: 'position' 'elevation' 'X'  'Y'
# dem_r:              digital elevation model (raster layer)
# cliff_r:            raster representation of ice cliff
# cliff_poly:         spatial polygon of ice cliff
# resol:		          ORIGINAL resolution of DEM
# newgridsize:        desired gridsize where local horizon has to be considered
#                      (if highresDEM = TRUE) [m]
CLIFF_CLOSEHORIZON<-function(v,dem_r,cliff_r,cliff_poly,resol,newgridsize){
  
  library(raster,zoo)
  # endresult<-vector('list')
  
  # ###TESTING###
  # v<-horcalc[i,]
  # dem_r<-gl_geom_st$Elevation     
  # cliff_r<-cl_geom_r
  # cliff_poly<-cliffs_p
  # resol<-resol
  # newgridsize<-newgridsize
  # ID_cliffs<-ID_cliffs
  ###########
  DFR_x<-v[1,3]
  DFR_y<-v[1,4]
  DFR_ele<-v[1,2]
  
  newminx<-DFR_x-((newgridsize/resol)/2)*resol
  newmaxx<-DFR_x+((newgridsize/resol)/2)*resol
  newminy<-DFR_y-((newgridsize/resol)/2)*resol
  newmaxy<-DFR_y+((newgridsize/resol)/2)*resol
  dem_cropped<-crop(dem_r,raster::extent(newminx,newmaxx,
                                         newminy,newmaxy))
  # newgrid DEM
  xy<-c(DFR_x,DFR_y)
  dist_r<-distanceFromPoints(dem_cropped,xy)
  # subtract elevation from point to create base dem
  z_r<-dem_cropped-DFR_ele
  # calculate tangent from point to each cell
  tan_r<-z_r/dist_r
  
  # avoid -Inf for central cell & replace with 0
  tan_r[is.infinite(tan_r)]<-0
  
  # projection has to be defined for terrain analysis
  projection(dist_r)<-projection(dem_r)
  # classification of 360 zones, one per degree of view
  zone_r<-terrain(dist_r,opt='aspect',unit='degrees',neighbors=8)
  # create raster stack
  stacked_r<-stack(tan_r,dist_r,z_r)
  names(stacked_r)<-c('tan','d','z')
  # get maximum value (= horizon) for each view direction
  # (older R-version: 'stat=' instead of 'fun=')
  newgridsize_hor<-data.frame(zonal(stacked_r$tan,zone_r,max,
                                    na.rm=TRUE))
  
  # cliff DEM
  # query polygon which contains specific i-pixel coordinate
  xy_pt<-SpatialPoints(cbind(xy[1],xy[2]),
                       proj4string=CRS(projection(dem_r)))
  pol<-over(xy_pt,cliff_poly,returnList=FALSE) 
  # ->resulting variable should be a single integer 
  #   & is stored as a data frame
  
  if(typeof(pol) == 'integer'){
    pol<-which(names(ID_cliffs) == pol[1])# if stored as 'Named Int'
  }
  if(typeof(pol) == 'list'){
    pol<-which(ID_cliffs == pol$ID[1])# if stored as a data frame 
    # pol<-which(ID_cliffs == pol[1,1])# if stored as a data frame 
  }      
  pol<-cliff_poly[pol,]
  # make same extent
  st_r<-crop(stacked_r,extent(cliff_r[[1]]))
  # mask stacked raster with cliff polygon
  cl_r<-mask(st_r,pol)
  # classification of 360 zones, one per degree of view
  cl_zone_r<-terrain(st_r$d,opt='aspect',unit='degrees',neighbors=8)
  # get maximum value (= horizon) for each view direction
  # (older R-version: 'stat=' instead of 'fun=')
  cl_hor<-data.frame(zonal(cl_r$tan,cl_zone_r,max,na.rm=TRUE))
  # replace -Inf values with NA, otherwise conversion
  #  of -Inf to -90.0 further down while rad->deg converting
  is.na(cl_hor)<-do.call(cbind,lapply(cl_hor,is.infinite))
  
  # both
  # list both horizons from newgrid and from cliff only
  # ->creating df with:
  #  'zone'  'horangle_newgrid'  'horangle_cliff'
  #    (many NA's in 'horangle_cliff' as some azimuth
  #     directions are not defined due to the lack of pixels)
  hor<-merge(newgridsize_hor,cl_hor,by='zone',all=TRUE)
  colnames(hor)<-c('zone','all_hor','cliff_hor')
  
  # converting radians to degrees
  # (older R-version: 'hor$max' instead of 'hor$value')
  hor$all_hor<-(atan(hor$all_hor)*180)/pi
  hor$cliff_hor<-(atan(hor$cliff_hor)*180)/pi
  hor<-round(hor,3)
  
  # use angle from cliff horizon if cliff_hor > all_hor
  #  (might happen rarely due to different zone-arrangement)
  idx<-which(hor$all_hor < hor$cliff_hor)
  hor[idx,2]<-hor[idx,3]
  
  # find where zones leap azimuth directions & interpolate there
  full.circ<-0:360
  idx_pres<-full.circ[full.circ %in% hor$zone]  #present azimuths
  #idx_miss<-full.circ[!full.circ %in% hor$zone]  #missing azimuths
  full.circ<-data.frame(full.circ)
  names(full.circ)<-'zone'
  # merge ideal circle with present data
  t1<-merge(full.circ,hor,by='zone',all=TRUE)
  #to avoid possible error in newer R version?
  idx<-which(is.na(t1$zone))
  if(length(idx) > 0){t1<-t1[-idx,]}
  rm(idx)
  # interpolate missing 'hor_all' linearly between present values
  t2<-zoo::na.approx(t1[min(idx_pres):max(idx_pres)+1,2])
  # create data frame with interpolated data
  #  (evt. start & end still missing)
  part.circ<-full.circ[min(idx_pres):max(idx_pres)+1,1]
  t3<-data.frame(cbind(part.circ,t2))
  names(t3)<-c('zone','all_hor')
  # merge again ideal circle with partly interpolated data
  #  & add also cliff horizons
  t4<-merge(full.circ,t3,by='zone',all=TRUE)
  t4<-data.frame(cbind(t4,t1$cliff_hor))
  # convert raster package directions (0?=south) 
  #  to 'normal' directions (0?=north)
  t5<-rbind(t4[180:360,],t4[1:179,])
  t6<-data.frame(cbind(1:360,t5[,2:3]))
  colnames(t6)<-c('az','el_all','el_ice')
  # interpolate still missing 'el_all' linearly 
  #  between present values
  #  (around 180?/S for northerly facing cliffs) 
  t7<-zoo::na.approx(t6$el_all)
  t7<-data.frame(cbind(t6$az,t7,t6$el_ice))
  colnames(t7)<-c('az','el_all','el_ice')
  t8<-as.matrix(t7)
  # important: replace NA/Inf's with a number!
  t8[is.na(t8)]<- -9999
  
  ###TESTING###
  # plot(t7$el_all,type='n')
  # lines(t7$el_all,lwd=5,col='orange')
  # lines(t7$el_ice,lwd=5,col='red')
  # plot(dem_cropped)
  # plot(cliffs_p,add=T)
  # points(x=xy[1],y=xy[2],pch=3,lwd=3)
  # plot(atan(tan_r)*180/pi)
  # plot(cliffs_p,add=T)
  # points(x=xy[1],y=xy[2],pch=3,lwd=1)
  #############
  
  endresult<-t8[,2:3]
  return(endresult)
}

################################################################################
CLIFF_FARHORIZON<-function(v,dem_r,resol,resol2){
  
  library(raster,zoo)
  ###TESTING###
  # dfr<-horcalc
  # dem_r<-dem_ASTERclip     
  # resol<-res_ASTER
  # resol2<-resol
  # ID<-modelID
  # i<-296
  #############
  
  DFR_x<-v[1,3]
  DFR_y<-v[1,4]
  
  xy<-c(DFR_x,DFR_y)
  dist_r<-distanceFromPoints(dem_r,xy)
  # subtract elevation from point to create base dem
  # ?
  #####################
  #z_r<-dem_r-dfr[,2]
  z_r<-dem_r-v[1,2]
  #####################
  # ?
  # calculate tangent from point to each cell
  tan_r<-z_r/dist_r
  
  # avoid -Inf for central cell & replace with 0
  tan_r[is.infinite(tan_r)]<-0
  
  # projection has to be defined for terrain analysis
  projection(dist_r)<-projection(tan_r)
  # classification of 360 zones, one per degree of view
  zone_r<-terrain(dist_r,opt='aspect',
                  unit='degrees',neighbors=8)
  # create raster stack
  st_r<-stack(tan_r,dist_r,z_r)
  names(st_r)<-c('tan','d','z')
  # get maximum value (= horizon) for each view direction
  # (older R-version: 'stat=' instead of 'fun=')
  hor<-data.frame(zonal(st_r$tan,zone_r,max,na.rm=TRUE))
  colnames(hor)<-c('zone','all_hor')
  # converting radians to degrees
  # (older R-version: 'hor$max' instead of 'hor$value')
  hor$all_hor<-(atan(hor$all_hor)*180)/pi
  hor<-round(hor,3)  
  
  # find where zones leap azimuth directions & interpolate there
  full.circ<-0:360
  #present azimuths
  idx_pres<-full.circ[full.circ %in% hor$zone]
  #missing azimuths
  #idx_miss<-full.circ[!full.circ %in% hor$zone]
  full.circ<-data.frame(full.circ)
  names(full.circ)<-'zone'
  # merge ideal circle with present data
  t1<-merge(full.circ,hor,by='zone',all=TRUE)
  # interpolate missing 'hor_all' linearly 
  #  between present values
  t2<-zoo::na.approx(t1[min(idx_pres):max(idx_pres)+1,2])
  # create data frame with interpolated data
  #  (evt. start & end still missing)
  part.circ<-full.circ[min(idx_pres):max(idx_pres)+1,1]
  t3<-data.frame(cbind(part.circ,t2))
  names(t3)<-c('zone','all_hor')
  # merge again ideal circle with partly interpolated data
  t4<-merge(full.circ,t3,by='zone',all=TRUE)
  # convert raster package directions (0?=south) 
  #  to 'normal' directions (0?=north)
  t5<-rbind(t4[180:360,],t4[1:179,])
  t6<-data.frame(cbind(1:360,t5[,2]))
  colnames(t6)<-c('az','el_all')
  # interpolate still missing 'el_all' (around 180?/S) linearly
  #  between present values
  t7<-zoo::na.approx(t6$el_all)
  t7<-data.frame(cbind(t6$az,t7))
  colnames(t7)<-c('az','el_all')
  t8<-as.matrix(t7)
  
  ###TESTING###
  # plot(t7$el_all,type='n')
  # lines(t7$el_all,lwd=5,col='orange')
  # plot(dem_r)
  # #e<-drawExtent()
  # plot(crop(dem_r,e))
  # plot(cliffs_p,add=T)
  # plot(crop(atan(tan_r)*180/pi,e))
  # plot(cliffs_p,add=T)
  #############
  
  # keep only horizon value & remove azimuth
  endresult<-t8[,2]
  return(endresult)
}


################################################################################
# CALCULATE DIURNAL CYCLE (+SD) PER CLIFF FROM MATRIX ->TS x CLIFF PIXELS
#
# cliff_IDs:   Vector with cliff ID per pixel
# ts_v:        Vector with timesteps (numeric or POSIX-format) 
# Flux_m:      Matrix with rows as timesteps (equal to "ts_v") 
#               and columns for each cliff pixel (of different cliffs)

DC_PC<-function(cliff_IDs,ts_v,Flux_m){
  
  unique_ids<-unique(cliff_IDs)
  date_h<-as.POSIXlt(ts_v,origin='1970-01-01')$hour
  
  Flux_dc_ls<-list()
  Flux_sd_ls<-list()
  for(cl in 1:length(unique_ids)){
    # create index to select all pixels of a single cliff
    idx<-which(cliff_IDs == unique_ids[cl])
    
    # aggregate flux values to diurnal cycle per cliff pixel            
    Flux_dc<-aggregate(Flux_m[,idx],list(date_h),mean)
    
    # derive one diurnal cycle for entire cliff
    Flux_dc_ls[[cl]]<-apply(Flux_dc,1,mean)
    
    # calculate standard deviation between all pixels of a cliff 
    #  per hour of diurnal cycle
    Flux_sd_ls[[cl]]<-apply(Flux_dc,1,sd)
  }
  output_ls<-list(Flux_dc_ls,Flux_sd_ls)
  names(output_ls)<-c('DC','SD')
  return(output_ls)
}




################################################################################
# CALCULATE DIURNAL CYCLE (+SD) PER DATA-VECTOR ->TS x VARIABLES
#
# dataframe:   Dataframe with timesteps (numeric or POSIX-format) 
#               & values for each timestep (each column one variable)
#
# e.g.:
#  Date                  L_s   L_d   H      Q_m   melt
#  2013-05-19 00:00:00   220.3 127.5 6.793  48.44 0.0006
#  2013-05-19 01:00:00   216.2 126.2 3.954  40.19 0.0005
#  2013-05-19 02:00:00   210.8 122.5 2.172  29.27 0.0004
#  2013-05-19 03:00:00   208.3 119.5 2.453  24.14 0.0003
#  2013-05-19 04:00:00   210.6 120.2 3.912  28.48 0.0004
#  2013-05-19 05:00:00   212.3 121.6 4.608  32.30 0.0004

DC<-function(dataframe){
  
  # extract only hours per timestep
  date_h<-as.POSIXlt(dataframe$Date,origin='1970-01-01')$hour
  
  # remove 'Date'-column
  dataframe<-dataframe[!colnames(dataframe) %in% 'Date']
  
  # aggregate flux values to diurnal cycle (and compute sd)
  Flux_dc_ls<-list()
  Flux_sd_ls<-list()
  for(v in 1:ncol(dataframe)){          
    Flux_dc_ls[[v]]<-aggregate(dataframe[[v]],list(date_h),mean,na.rm=TRUE)[,2]
    Flux_sd_ls[[v]]<-aggregate(dataframe[[v]],list(date_h),sd,na.rm=TRUE)[,2]
  }
  
  # diurnal cycle                                       
  DC_df<-data.frame(1:24,do.call(cbind,Flux_dc_ls))
  colnames(DC_df)<-c('hours',colnames(dataframe)) 
  
  # standard deviation
  SD_df<-data.frame(1:24,do.call(cbind,Flux_sd_ls))
  colnames(SD_df)<-c('hours',colnames(dataframe)) 
  
  output_ls<-list(DC_df,SD_df)
  names(output_ls)<-c('DC','SD')
  return(output_ls)
}




################################################################################################
# recursive function to remove name from all levels of list
# 
# from: https://stackoverflow.com/questions/37853679/removing-elements-in-a-nested-r-list-by-name

stripname <- function(x, name) {
  thisdepth <- depth(x)
  if (thisdepth == 0) {
    return(x)
  } else if (length(nameIndex <- which(names(x) == name))) {
    x <- x[-nameIndex]
  }
  return(lapply(x, stripname, name))
}



################################################################################################
# function to find depth of a list element
# 
# see: http://stackoverflow.com/questions/13432863/determine-level-of-nesting-in-r

depth <- function(this, thisdepth=0){
  if (!is.list(this)) {
    return(thisdepth)
  } else{
    return(max(unlist(lapply(this,depth,thisdepth=thisdepth+1))))    
  }
}
