################################################################################
# ICE CLIFF ENERGY BALANCE MODEL
#
# Authors: Pascal Buri(1), Marin Kneib(1,2)
#          (1) Swiss Federal Institute for Forest, Snow and Landscape Research, WSL
#          (2) Université Grenoble Alpes, France
#   
#
# Model details: Buri et al., 2016a (https://doi.org/10.3189/2016AoG71A059)
#                Kneib et al., 2022 (https://doi.org/10.5194/tc-16-4701-2022)
#
#
#
# Contact: Dr. Pascal Buri
#          Swiss Federal Institute for Forest, Snow and Landscape Research, WSL
#          Zürcherstrasse 111, 8903 Birmensdorf, Switzerland
#          pascal.buri@wsl.ch
################################################################################

#===============================================================================
# INITIALIZE WORKSPACE
#===============================================================================
## clear entire workspace (excl. packages)
rm(list = ls())
gc()

## define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')


################################################################################
################################################################################
# PRIMARY DEFINITIONS (TO CHANGE REGULARLY)
################################################################################
################################################################################
## define cliff ID
CL<-2

## define simulation folder name (will be created) 
sim_folder<-'Lirung-C-'&CL

## define root path of all data
root<-'E:/CliffEBM'

## define the reduction of cores to be used (max. number - reduction) 
reduction<-4
# reduction<-0

## define parameter file name 
##  Parameter file defines for single or multiple (e.g. for sensitivity analysis, Monte Carlo analysis etc.) runs: 
##  albedos, emissivities, structural parameters, meteorological offsets
parfn<-'param_20220103.txt'

## define cliff-shapefile name
clifn<-'Cliff2_May2013.shp'

## define glacier shapefile name
glafn<-'glacier_mask_2014.shp'

## define primary dem1 file name (high resolution, covering cliff and close surrounding)
demfn1<-'UAV-DEM_1m_May2013_Cliff2.tif'

## define secondary dem1 file name (coarse resolution, covering glacier and relevant far surrounding)
demfn2<-'ALOS100m.tif'

## define meteodata file name
metfn<-'meteodata_AWS2013.txt'

## define target resolution for primary DEM [m]
resol_calc<-1   

## define if additional outputs should be generated from model (TRUE) or not (FALSE)
##  (e.g. mean & sd of fluxes per season, model parameters, horizons per cell, rasters)
additionalOutput<-TRUE
# additionalOutput<-FALSE

## set date-format to site-specific time (to avoid summertime)
targetTZ<-'Asia/Kathmandu'
localTZ<-'Europe/Zurich'

# define melt model loop (has to be inside of meteodata-timestamps)
#  [mm-dd HH:MM:SS]
simStart<-"2013-05-19 00:00:00" 
simStop<-"2013-10-22 23:00:00"             

## define interval of intermediate outputs
##  e.g. daily ('%j'), weekly ('%W'), monthly ('%m') or all timesteps ('ALL')
##  (call '?strptime' to see other formats)
# interv<-'%W'    #weekly
interv<-'%m'  #monthly
# interv<-'ALL' #run over all timesteps


################################################################################
################################################################################
# SECONDARY DEFINITIONS (TO CHANGE RARELY)
################################################################################
################################################################################
## define path to data folder (existing)
path_data<-root&'/Data'

## define path to horizon files folder (existing)
path_hor<-path_data&'/Horizons'

## define path to model output folder (existing)
path_out<-path_data&'/ModelOutputs'

## define path to code folder (existing)
path_code<-root&'/Code'

## define path to functions file (existing)
path_func<-path_code&'/Functions.r'

## define how many neighboring cells to use to compute slope/aspect for any cell:
## ?raster::terrain
## 4: 'Fleming and Hoffer (1979)', 'Ritter (1987)', better for smoother surfaces
## 8: 'Horn (1981)', best for rough surfaces
neighbCells<-8

# desired squared gridsize in which surrounding terrain has to be considered 
#  for debrisview-calculations [m], e.g.
#  200: uses 200m x 200m grid
newgridsize<-100

## define UTM projection of all geographical data (otherwise specified within code)
projec<-'+proj=utm +zone=45 +datum=WGS84 +units=m +no_defs'


################################################################################
################################################################################
# LOAD PACKAGES
################################################################################
################################################################################
library(doParallel)                                                                             
library(foreach)                                                                                 
library(grDevices)                                                                               
library(iterators)                                                                               
library(methods)                                                                                 
library(parallel)                                                                                
library(raster)                                                                                  
library(rgdal)                                                                                   
library(rgeos)                                                                                  
library(sf)                                                                                      
library(sp)                                                         
library(stats)                                                                                 
library(utils)
library(zoo)   #used in function file only
## install missing packages using the install.packages() function (e.g.: install.packages("zoo") )

## define no. of digits printed in console
options('scipen'=100, 'digits'=4)

## change local time from german to english (check in sessionInfo()) for plotting
Sys.setlocale('LC_TIME', 'C')
Sys.setenv(TZ=targetTZ)

## no. of cores to be used
n.cores<-parallel::detectCores(all.tests = FALSE, logical = TRUE) - reduction

## register cores for parallel functionality
registerDoParallel(cores=n.cores)   #put "stopImplicitCluster()" at the end


#===============================================================================
# LOAD GLACIER-MASK POLYGON
#===============================================================================
p<-st_read(path_data&'/Polygons/'&glafn,quiet=T)
p<-as(p,'Spatial')
projection(p)<-projec
buffer_gla_msk<-buffer(p,50) # add buffer around glacier [m]
rm(p)


#===============================================================================
# READ PARAMETER-FILE
#===============================================================================
paramfile<-read.csv(path_data&'/Parameters/'&parfn,header=TRUE,dec='.',sep='\t')
paramfile$Run<-as.character(paramfile$Run)


#===============================================================================
# LOAD CLIFF POLYGONS
#===============================================================================
cliff_p<-st_read(path_data&'/Polygons/'&clifn,quiet=T)
cliff_p<-as(cliff_p,'Spatial')
projection(cliff_p)<-projec


#===============================================================================
# READ PRIMARY DEM
#===============================================================================
dem1<-raster(paste(path_data,'/DEMs/',demfn1,sep=''))
f<-resol_calc/round(res(dem1)[1],digits=3) # conversion factor
##resample if needed
if(f > 1){dem1<-aggregate(dem1,fact=f,method='bilinear')}
if(f == 1){dem1<-dem1}
projection(dem1)<-projec
resol1<-res(dem1)[1]


#===============================================================================
# READ SECONDARY DEM
#===============================================================================
dem2<-raster(paste(path_data,'/DEMs/',demfn2,sep=''))
projection(dem2)<-projec
dem2<-mask(dem2,buffer_gla_msk,inverse=TRUE) #mask glacier from dem2
resol2<-res(dem2)[1]




################################################################################
################################################################################
# LOOP OVER SENSITIVITY RUNS 
# (N = NUMBER OF ROWS IN PARAMETERS-FILE)
################################################################################
################################################################################
for(SRN in 1:nrow(paramfile)){
  runID<-paramfile[SRN,'Run']
  id<-'C-'&CL&'___SR-'&runID
  print('Modelling '&id)
  
  
  #===============================================================================
  # EXTRACT PARAMETERS
  #===============================================================================
  # << ad >> 
  # ice albedo [-]
  alpha_i<-paramfile[SRN,'ai']
  # << ai >> 
  # debris albedo [-]
  alpha_d<-paramfile[SRN,'ad']
  # << ed >> 
  # debris emissivity [-]
  epsilon_d<-paramfile[SRN,'ed']
  # << ei >> 
  # ice emissivity [-]
  epsilon_i<-paramfile[SRN,'ei']
  # << ft >> 
  # value for filtering (smoothing) of dem1 in terms of aspect:
  #   "filter_val x filter_val" mean filter around pixel (e.g. '3')
  filter_val<-paramfile[SRN,'ft'] 
  # << lw >> 
  # longwave radiation deviation [W/m2]
  LWdev<-paramfile[SRN,'lw']
  # << rh >> 
  # relative humidity deviation [%]
  RHdev<-paramfile[SRN,'rh']
  # << sw >> 
  # shortwave radiation deviation [W/m2]
  SWdev<-paramfile[SRN,'sw']
  # << ta >> 
  # T air deviation [K]
  TAdev<-paramfile[SRN,'ta'] 
  # << ts >> 
  # T surface deviation [K]
  TSdev<-paramfile[SRN,'ts']
  # << ws >> 
  # wind speed multiplier [-]
  WSmult<-1  
  # << z0 >> 
  # surface roughness [m]
  z0<-paramfile[SRN,'z0'] 
  
  ###cliff polygon
  # simplify with Ramer-Douglas-Peucker algorithm
  #  https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm
  cliff_p<-rgeos::gSimplify(cliff_p,tol=0.1,topologyPreserve=TRUE)
  #plot(dem1)
  #plot(cliff_p,add=TRUE)
  
  ## Access polygon information (ID, area,...)
  area<-as.numeric(sapply(slot(cliff_p,'polygons'),function(x)slot(x,'area')))
  centr<-sapply(slot(cliff_p,'polygons'),function(x)slot(x,'labpt'))
  centrX<-centr[1,]
  centrY<-centr[2,]
  cliffs_df<-as.data.frame(cbind(area,centrX,centrY))
  cliffs_df$ID<-CL
  rm(area,centr,centrX,centrY)
  
  
  #===============================================================================
  # CREATE DIRECTORIES & TIMESTAMPS
  #===============================================================================
  # Model output directory    
  path_newdir<-path_out&'/'&sim_folder&'/Run-'&id
  dir.create(path_newdir,recursive=TRUE,showWarnings=FALSE)
  
  # specific output folders
  path_outrast<-path_newdir&'/Rasters'
  dir.create(path_outrast,showWarnings=FALSE)
  path_outpol<-path_newdir&'/Polygons'
  dir.create(path_outpol,showWarnings=FALSE)
  path_outtab<-path_newdir&'/Tables'
  dir.create(path_outtab,showWarnings=FALSE)
  path_outlist<-path_newdir&'/Lists'
  dir.create(path_outlist,showWarnings=FALSE)
  path_outplot<-path_newdir&'/Plots'
  dir.create(path_outplot,showWarnings=FALSE)
  
  #===============================================================================
  # CALCULATE GEOMETRY
  #===============================================================================
  # calculate terrain values of clipped initial dem1
  slope_r<-raster::terrain(dem1,opt='slope',unit='degrees',
                           neighbors=neighbCells)
  aspect_r<-raster::terrain(dem1,opt='aspect',unit='degrees',
                            neighbors=neighbCells)
  
  # stack geometry raster layers
  gl_geom_st<-raster::stack(dem1,slope_r,aspect_r)
  names(gl_geom_st)<-c('Elevation','Slope','Aspect')
  rm(slope_r,aspect_r)
  # plot(gl_geom_st$Elevation)
  # plot(cliff_p,add=TRUE)
  
  # rasterize cliff polygons (clipped dem1-raster needed only for extent)
  cliff_r<-rasterize(cliff_p,gl_geom_st)
  # plot(cliff_r)
  
  # mask clipped rasters with cliff raster
  cl_geom_r<-mask(gl_geom_st,cliff_r)
  
  
  #=================================================================
  # FILTER OUT ABNORMAL ASPECT VALUES
  #=================================================================
  # Aspect, component-wise:
  # convert degrees to radians & calculate sines
  #  (x-axis component, W-E, zonal direction)
  WEcomp_r<-sin(cl_geom_r$Aspect*pi/180)
  # plot(WEcomp_r)
  
  # convert degrees to radians & calculate cosines
  #  (y-axis component, S-N, meridional direction)
  SNcomp_r<-cos(cl_geom_r$Aspect*pi/180)
  # plot(SNcomp_r)
  
  # take the arctan of the averages, 
  #  rotate to positive values, 
  #  and then convert back to degrees
  aspect_r<-atan2(WEcomp_r,SNcomp_r)
  rm(WEcomp_r,SNcomp_r)
  # plot(aspect_r)              
  
  # Median filter window (e.g. 9x9)
  aspect_r<-focal(aspect_r,w=matrix(1,nrow=filter_val,ncol=filter_val),fun=median,
                  na.rm=TRUE)
  # plot(aspect_r) 
  
  aspect_v<-as.vector(aspect_r)
  idx<-which(aspect_v<0,arr.ind=FALSE)
  aspect_v[idx]<-aspect_v[idx]+(2*pi)
  aspect_v<-aspect_v*180/pi
  aspect_r<-setValues(aspect_r,aspect_v)
  
  # plot(crop(cl_geom_r$Aspect,extent(cliff_p)+10)) # old
  # plot(cliff_p,add=TRUE) 
  # plot(crop(aspect_r,extent(cliff_p)+10))# filtered
  # plot(cliff_p,add=TRUE)
  
  cl_geom_r$Aspect<-aspect_r
  rm(aspect_v,aspect_r)
  
  # mask clipped rasters with cliff raster
  cl_geom_r<-mask(cl_geom_r,cliff_r)
  
  
  #=======================================================================
  # GET GEOMETRY VALUES + CLIFF-ID FOR EVERY CLIFF PIXEL
  #=======================================================================
  # extract elevation/slope/aspect values from raster (incl. NAs)
  cliffs_geom_df<-data.frame(raster::extract(cl_geom_r,extent(cl_geom_r)))
  # select only cliff-pixels
  cliffs_geom_df<-cliffs_geom_df[complete.cases(cliffs_geom_df),]
  # use pixel numbers (relative to "new_geometry_r") as index
  idx<-as.numeric(rownames(cliffs_geom_df))
  # extract coordinates of cliff pixels
  cliff_coord_m<-xyFromCell(cl_geom_r,idx)
  rm(idx)
  
  # dataframe with elevation/slope/aspect/x/y
  cliffs_geom_df<-data.frame(cbind(cliffs_geom_df,cliff_coord_m))
  colnames(cliffs_geom_df)<-c('Elevation','Slope','Aspect','x','y')
  
  # add cliff ID
  # query cliff polygon which contains specific pixel coordinate
  xy_pts<-SpatialPoints(cbind(cliffs_geom_df[,'x'],cliffs_geom_df[,'y']),
                        proj4string=CRS(projec))
  projection(cliff_p)<-projec
  cliffs_geom_df$ID<-over(xy_pts,cliff_p)
  
  horcalc<-cliffs_geom_df[,c(1,4:6)] # only elevation, coordinates & ID
  horcalc<-cbind(1:nrow(cliffs_geom_df),horcalc)# add specific number per point
  horcalc[,] <- as.numeric(as.character(unlist(horcalc[,])))
  colnames(horcalc)<-c('cell','elevation','X','Y','ID')
  
  #plot(gl_geom_st$Elevation)
  #plot(cliff_p,add=TRUE)
  
  gc()
  
  
  #===============================================================================
  # GEOMETRY CALCULATIONS
  #===============================================================================
  ID_cliffs<-as.character(cliffs_df$ID)
  calculateInitHorizon<-TRUE
  
  cliff_p$ID<-ID_cliffs
  
  if (calculateInitHorizon == TRUE) {
    #################################
    # horizon angle from primary DEM  ->comput. expensive!
    #################################
    N<-nrow(horcalc)
    
    source(path_func)
    print('*** CLIFF_CLOSEHORIZON @t0   |||   '&id&'    |||   '&N&'px ***')
    horfile_dem1<-foreach(i = icount(N),.combine=cbind) %dopar% 
      CLIFF_CLOSEHORIZON(horcalc[i,],gl_geom_st$Elevation,cl_geom_r,cliff_p,
                         resol1,newgridsize)
    print('************************************************************************* DONE ***')
    horfile_dem1<-data.frame(horfile_dem1)
    horfile_dem1<-cbind(1:360,horfile_dem1)
    
    horfile_dem1[horfile_dem1== -9999]<-NA
    
    
    ###################################
    # horizon angle from secondary DEM  ->comput. expensive!
    ###################################
    print('**** CLIFF_FARHORIZON @t0   |||   '&id&'    |||   '&N&'px ****')
    horfile_dem2<-foreach(i = icount(N),.combine=cbind) %dopar% 
      CLIFF_FARHORIZON(horcalc[i,],dem2,resol2,resol1)
    print('************************************************************************* DONE ***')
    rm(N,horcalc)
    
    horfile_dem2<-data.frame(horfile_dem2)
    horfile_dem2<-cbind(as.numeric(1:360),horfile_dem2)
    elem<-ncol(horfile_dem2)-1
    colnames(horfile_dem2)<-c('Azimuth',rep('Px_'&1:elem,1))
    
    
    
    ###########################################
    # define 360° skyview angle for every point
    ###########################################
    # select only 'el_all' (highest elevation/horizon angle overall) from UAV
    horHRDEM_all<-horfile_dem1[,seq(2,ncol(horfile_dem1),by=2)]
    
    # select only 'el_ice' (highest elevation/horizon angle only on ice) from UAV
    horHRDEM_ice<-horfile_dem1[,seq(1,ncol(horfile_dem1),by=2)]
    rm(horfile_dem1)
    horHRDEM_ice<-horHRDEM_ice[,-1] # to get rid of azimuth column
    
    resulthor<-vector('list')
    for (cell in 1:nrow(cliffs_geom_df)){ # cells 1-n on cliff
      # horizon angle for shortwave radiation:
      # [hor_skyI]
      # combine horizon angles derived from ASTER & UAV
      #  & choose highest angle for every degree of azimuth
      hor1<-data.frame(cbind(horfile_dem2[,cell+1],horHRDEM_all[,cell]))
      colnames(hor1)<-c('horang_ASTER','horang_UAV')
      hor_skyI<-apply(hor1,1,max,na.rm=TRUE)
      hor_skyL<-hor1$horang_UAV
      
      # debris portion angle (for longwave radiation from debris)
      # [prt_deb]
      # combine horizon angles derived UAV
      #  & choose effective debris portion angle for every degree of azimuth
      hor2<-data.frame(cbind(horHRDEM_all[,cell],horHRDEM_ice[,cell]))
      prt_deb<-cbind(hor2,rep(NA,nrow(hor2)))
      colnames(prt_deb)<-c('horang_HRDEMall','horang_HRDEMice','debris_portion')
      #  [case 1] azimuth direction where horizon angle is defined by ice surface:
      idx1<-which(prt_deb$horang_HRDEMall == prt_deb$horang_HRDEMice)
      prt_deb[idx1,3]<-0
      #  [case 2] azimuth direction where horizon angle is higher than ice surface:
      idx2<-which(prt_deb$horang_HRDEMall > prt_deb$horang_HRDEMice)
      prt_deb[idx2,3]<-prt_deb[idx2,1]-prt_deb[idx2,2]
      #  [case 3] azimuth direction where horizon angle only for ice is not defined
      #   (debris can be assumed in this direction):
      idx3<-which(is.na(prt_deb$horang_HRDEMice))
      prt_deb[idx3,3]<-prt_deb[idx3,1]+90
      
      hor3<-data.frame(cbind(hor_skyI,hor_skyL,prt_deb[,3]))
      colnames(hor3)<-c('horizon_angle_I','horizon_angle_L','deb_prt_angle')
      resulthor[[cell]]<-hor3
      names(resulthor[[cell]])<-c('hor_angleI_'&cell,'hor_angleL_'&cell,'deb_prt_angle'&cell)
    }
    
    hor_all<-do.call('cbind',resulthor)
    rm(resulthor)
    hor_all<-data.frame(hor_all)
    
    # output-file of horizon/debris angles per cliff cell:
    # convert output to 3 digits after decimal point & avoid spacing between columns
    hor_all_output<-format(hor_all,digits=4,trim=TRUE)
    hor_all_output<-cbind(1:360,hor_all_output)
    colnames(hor_all_output)[1]<-'Azimuth'
    filename_horDEF<-path_hor&'/horDEF___C-'&ID_cliffs&'___'&
      round(resol1,digits=1)&'m.txt'
    write.table(hor_all_output,filename_horDEF,dec='.',sep=',',quote=FALSE,
                row.names=FALSE)
    rm(hor_all_output,filename_horDEF)
    
    hor_all<-cbind(1:360,hor_all)
    colnames(hor_all)[1]<-'Azimuth'
    
    horDEF<-hor_all
    rm(hor_all)
    
    
    ##################################################################
    # define skyviewfactor (I, L) & debrisviewfactor from every point
    ##################################################################
    horDEF<-horDEF[,-1] # get rid of azimuth column
    # horizon angle for sky SHORTWAVE (0?=horizontal, 90?=vertical upwards)
    obstr_s_I<-horDEF[,seq(1, ncol(horDEF), by = 3)]
    # horizon angle for sky LONGWAVE (0?=horizontal, 90?=vertical upwards)
    obstr_s_L<-horDEF[,seq(2, ncol(horDEF), by = 3)]
    # obstruction angle for debris (angle between lower and upper limit of debrisview)
    obstr_d<-horDEF[,seq(3, ncol(horDEF), by = 3)]
    
    skyview_factor_I<-vector('list')
    skyview_factor_L<-vector('list')
    debrisview_factor<-vector('list')
    
    for (cell in 1:ncol(obstr_s_I)){ # cells 1-n on cliff
      # skyview factor SHORTWAVE
      # ->convert horizon angles I appropriately as used in literature
      #   (=angle between vertical up (solar zenith) and horizon)
      sky_I<-abs(obstr_s_I[,cell]-90)
      # calculate portion of visible sky per azimuth direction
      sky_I<-sky_I/90
      # get mean of every portion of skyview around 360?
      skyview_factor_I[[cell]]<-mean(sky_I,na.rm=TRUE)
      
      # skyview factor LONGWAVE
      # ->convert horizon angles L appropriately as used in literature
      #   (=angle between vertical up (solar zenith) and horizon)
      sky_L<-abs(obstr_s_L[,cell]-90)
      # calculate portion of visible sky per azimuth direction
      sky_L<-sky_L/90  #
      # get mean of every portion of skyview around 360?
      skyview_factor_L[[cell]]<-mean(sky_L,na.rm=TRUE)
      
      
      # debrisview  factor
      # ->calculate portion of visible debris per azimuth direction
      # ->cliffs: lowest horizon is -90?, maximum 90?
      debris<-obstr_d[,cell]/180
      # get mean of every portion of skyview around 360?:
      debrisview_factor[[cell]]<-mean(debris,na.rm=TRUE)
    }
    
    # skyview factors I
    VsI<-do.call('rbind',skyview_factor_I)
    VsI<-round(VsI,3)
    # skyview factors L
    VsL<-do.call('rbind',skyview_factor_L)
    VsL<-round(VsL,3)
    # debrisview factors
    Vd<-do.call('rbind',debrisview_factor)
    Vd<-round(Vd,3)
    # combine with cliff geometry
    cliffs_geom_df<-cbind(cliffs_geom_df,VsI,VsL,Vd)
    rm(skyview_factor_I,skyview_factor_L,debrisview_factor,VsI,VsL,Vd)
    
    # add cliff ID
    cliffs_geom_df$ID<-as.character(cliff_p$ID)
    rm(xy_pts)
    
    # define format of output
    cliffs_geom_df<-within(cliffs_geom_df,
                           {Elevation<-formatC(Elevation,format='f',digits=2)
                           Slope<-formatC(Slope,format='f',digits=1)
                           Aspect<-formatC(Aspect,format='f',digits=1)
                           ID<-formatC(ID,format='f',digits=0)
                           x<-formatC(x,format='f',digits=3)
                           y<-formatC(y,format='f',digits=3)
                           VsI<-formatC(VsI,format='f',digits=3)
                           VsL<-formatC(VsL,format='f',digits=3)
                           Vd<-formatC(Vd,format='f',digits=3)
                           })
    filename_geom<-path_hor&'/geom___C-'&ID_cliffs&'___'&
      round(resol1,digits=1)&'m.txt'
    #
    write.table(cliffs_geom_df,filename_geom,dec='.',sep=',',
                quote=FALSE,row.names=FALSE)
    rm(filename_geom)
  }
  
  
  gc()
  
  #===============================================================================
  # READ GEOMETRY CALCULATIONS IF NOT CALCULATED BEFORE
  #===============================================================================
  t_load_horizons<-system.time(
    if (calculateInitHorizon == FALSE) {
      # read horizon files if not calculated above
      horDEF<-read.csv(file=path_data&'horizons/horDEF___C-'&ID_cliffs&'___'&
                         round(resol1,digits=1)&'m.txt',header=TRUE,dec='.',sep=',')
      horDEF<-horDEF[,-1] # get rid of azimuth column
    }
  )
  t_load_horizons
  
  # read geometry file anyway (numeric)
  cliffs_geom_df<-read.csv(file=path_hor&'/geom___C-'&ID_cliffs&'___'&
                             round(resol1,digits=1)&'m.txt',header=TRUE,dec='.',sep=',')
  
  # write file also into model output folder
  filename_geom<-path_newdir&'/geom_ts0___'&id&'_'&
    round(resol1,digits=1)&'m.txt'
  write.table(cliffs_geom_df,filename_geom,dec='.',sep=',',
              quote=FALSE,row.names=FALSE)
  rm(filename_geom)
  
  
  
  #===============================================================================
  # READ METEOROLOGICAL DATA
  #===============================================================================
  # read hourly data from AWS
  AWSdata<-read.csv(file=paste(path_data,'/Meteodata/'&metfn,sep=''),
                    header=TRUE,dec='.',sep='\t')
  
  # merge columns to date (POSIXct-format)
  AWSdata<-within(AWSdata,Date<-as.POSIXct(paste(year,doy,sprintf("%04d",time)),
                                           format='%Y %j %k'))
  
  # subset data according to defined start and end dates
  AWSdata<-subset(AWSdata,Date >= simStart)
  AWSdata<-subset(AWSdata,Date <= simStop)
  
  # plot(AWSdata$RH,type='l') # test
  
  # create vectors of AWSdata-columns
  Date<-AWSdata$Date
  
  # # check if no duplicated dates are in the dataset:
  # which(duplicated(Date))
  
  year<-as.numeric(strftime(Date,format='%Y'))
  doy<-as.numeric(strftime(Date,format='%j'))
  hour<-as.numeric(strftime(Date,format='%H'))
  ws<-AWSdata$WS
  rH<-AWSdata$RH
  SWin<-AWSdata$SWIN
  LWin<-AWSdata$LWIN
  LWout<-AWSdata$LWOUT
  Ta<-AWSdata$TA
  sigma <- 5.6704e-8;   # Stefan-Boltzmann constant (W m^-2 K^-4)
  # Tsurf<-((LWout - ((1-epsilon_d)*LWin))/(epsilon_d*sigma))^(1/4)-273.15
  Tsurf<-AWSdata$TS
  
  
  #===============================================================================
  # APPLY DEVIATIONS TO METEO-VARIABLES FROM UNCERTAINTY ANALYSIS
  #===============================================================================
  LWin<-LWin+LWdev
  rH<-rH+RHdev
  rH[rH<0]<-0
  rH[rH>100]<-100
  SWin<-SWin+SWdev
  Ta<-Ta+TAdev
  Ts<-Tsurf+TSdev
  ws<-ws*WSmult
  
  # avoid zero wind speed (numerical problems for turbulent fluxes)
  ws[ws < 0.01]<-0.01
  
  
  
  #===============================================================================
  # DEFINE REMAINING CONSTANTS
  #===============================================================================
  # define constants
  Env_LR<-0.0065    # Environmental temperature lapse rate [K m^-1] 
                    #  (used for pressure in turbulent flux)
  Lat<-28.21        # Latitude (pos. for N-hemisphere) [°]
  Lon<-85.56        # Longitude (pos. for W-hemisphere) [°]
  sigma<-5.67e-8    # Stefan-Boltzmann constant [W m^-2 K^-4]
  Ti<-0             # Ice cliff surface temperature [C]
  g<-9.81           # Gravitational acceleration [m s^-2]
  k<-0.41           # Von Karman's constant
  Rgas<-8.31447     # Gas constant [J mol^-1 K^-1]
  Mair<-0.0289644   # Molar mass of dry air [kg mol^-1]
  P_0<-101.3        # Standard sea level air pressure [kPa]
  T_0<-288.15       # Standard sea level air temperature [K]
  rho_0<-1.29       # Standard sea level air density [kg m^-3]
  L_e<-2514000        # Latent heat of evaporation of water [J kg^-1]
  L_f<-334000       # Latent heat of fusion of ice [J kg^-1]
  c_p<-1004         # Specific heat capacity of air at constant pressure
  rho_i<-900        # Ice density [kg m^-3]
  z<-2              # Height of meteorological measurements [m]
  e_s<-0.611        # Saturated vapour pressure at 0C [kPa]
  nu<-1.35e-5       # Viscosity of air
  
  
  #===============================================================================
  # PREPARE MODEL OUTPUT TIMESTEPS
  #===============================================================================
  # predefine timesteps in a vector, where intermediate results should be
  #  performed resp. the cliff geometry should be updated
  PERIOD<-Date[1:length(Date)]
  
  # intermediate output interval
  if(interv == 'ALL'){iOutput_ts<-length(PERIOD)} ##run over all timesteps
  if(interv != 'ALL'){
    iOutput_ts<-as.numeric(format(PERIOD,interv))
    iOutput_ts<-which(diff(iOutput_ts) != 0)
    iOutput_ts<-c(iOutput_ts,length(PERIOD))
  }
  
  
  #===========================================================================
  # PREPARE RASTERS ETC. FOR MODEL
  #===========================================================================
  
  ##
  # Raster brick with horizon angles
  # get cell numbers from cliff-cells for embedding back in rectangular raster
  cliffcell_nr<-cellFromXY(gl_geom_st$Elevation,cliff_coord_m)
  rm(cliff_coord_m)
  
  # create index selecting column with shortwave sky horizon for every point
  columnsSkyangle_I<-seq(1,ncol(horDEF),by=3)
  hors_I<-horDEF[,columnsSkyangle_I]
  rm(horDEF)
  
  hors_I<-data.frame(t(hors_I))
  #structure 'hors' with only 3 cliff-cells::
  #                   V1       V2      (...)   V360
  #         9         34.75    35.02   (...)   35.51
  #         13        34.66    34.94   (...)   35.37
  #         17        34.52    34.79   (...)   35.05
  
  # create data frame with NAs in every cell
  grid_df<-data.frame(replicate(360,sample(NA,ncell(gl_geom_st$Elevation),
                                           rep=TRUE)))
  #strucute 'grid_df':
  #                   X1        X2     (...)    X360
  #         1         NA        NA     (...)    NA
  #         2         NA        NA     (...)    NA
  #        (...)     (...)     (...)   (...)   (...)
  #         ncells    NA        NA     (...)    NA
  
  # merge all cell numbers
  grid_df[cliffcell_nr,]<-hors_I
  rm(hors_I)
  
  # create rasterbrick with 360 layers, each with the horizon angle per
  #   cliff cell in the correspondent azimuth direction
  horI_m<-data.matrix(grid_df)
  rm(grid_df)
  clip_br<-brick(gl_geom_st$Elevation)  #used later in model!
  # horI_br<-setValues(clip_br,horI_m)   #expensive!
  
  # horI_df<-as.data.frame(horI_br)
  horI_m<-horI_m[complete.cases(horI_m),] # only cliff-pixels
  
  ##
  # Raster with skyview/debrisview factor
  viewfactors<-data.frame(cbind(as.numeric(cliffs_geom_df$VsI),
                                as.numeric(cliffs_geom_df$VsL),
                                as.numeric(cliffs_geom_df$Vd)))
  
  colnames(viewfactors)<-c('VsI','VsL','Vd')
  #> head(viewfactors)
  #  VsI   VsL    Vd
  # 0.680 0.812 0.510
  # 0.660 0.777 0.514
  # 0.640 0.746 0.531
  
  # create data frame with NAs in every cell
  grid_df2<-data.frame(replicate(3,sample(NA,ncell(gl_geom_st$Elevation),
                                          rep=TRUE)))
  colnames(grid_df2)<-c('VsI','VsL','Vd')
  #> head(grid_df2)
  # VsI VsL Vd
  #  NA  NA NA
  #  NA  NA NA
  #  NA  NA NA
  
  # merge all cell numbers
  grid_df2[cliffcell_nr,]<-viewfactors
  rm(cliffcell_nr,viewfactors)
  
  # create rasterbrick with 3 layers, one with the skyview factor I,
  #  one with skyviewfactor L and one with the debrisview factor
  #  per cliff cell
  view_m<-data.matrix(grid_df2)
  view_br<-setValues(clip_br,view_m)
  rm(clip_br,grid_df2,view_m)
  
  # stack rasters & convert slope & aspect from deg to rad
  #  (names: elevation, slope, aspect, VsI, VsL, Vd)
  cl_geom_r<-stack(cl_geom_r,view_br)
  cl_geom_r$Slope<-cl_geom_r$Slope*pi/180
  cl_geom_r$Aspect<-cl_geom_r$Aspect*pi/180
  #plot(cl_geom_r)
  
  ##
  # prepare rasterstacks for calculation of
  #  diurnal cycle in melt model loop
  # make a vector of the unique hour values to see which are there
  hrs<-as.POSIXlt(Date)$hour
  hrs<-unique(hrs)
  
  # take an existing raster to create (24) layers with the correct extent
  empty_r<-cl_geom_r$Aspect
  empty_r<-setValues(empty_r,0)
  stack_dc<-stack(mget(rep('empty_r',length(hrs))))
  DC_names<-paste('Hour',hrs,sep='')
  names(stack_dc)<-DC_names
  
  # write initial cliff polygon to shapefile
  fn<-path_outpol&'/Cliffs_ts0_p.shp'
  st_write(st_as_sf(cliff_p),fn,driver="ESRI Shapefile",
           delete_dsn=TRUE,quiet=TRUE)
  rm(fn)
  
  # plot cliffs
  plot(cliff_p,col='orange',border='orange')
  
  # Histograms
  # Aspect:
  png(file=path_outplot&'/CliffAspect_ts0_hist.png',units='in',
      width=6,height=6,res=300)
  ti='Aspect  '&as.character(Date[1])&' (after filtering)'
  hist(cl_geom_r$Aspect*180/pi,main=ti,xlab='Aspect',
       breaks=seq(0,360,by=45),labels=c('NNE','ENE','ESE','SSE','SSW','WSW',
                                        'WNW','NNW'))
  dev.off()
  # Slope
  png(file=path_outplot&'/CliffSlope_ts0_hist.png',units='in',
      width=6,height=6,res=300)
  ti='Slope  '&as.character(Date[1])
  hist(cl_geom_r$Slope*180/pi,main=ti,xlab='Slope [°]',breaks=seq(0,90,by=5))
  dev.off()
  
  ##
  # Write model settings to file
  # switch to local time zone for current machine timestamp
  Sys.setenv(TZ=localTZ)
  generated<-as.character(Sys.time())
  print(generated)
  
  #switch back to target timezone
  Sys.setenv(TZ=targetTZ)
  print(Sys.time())
  
  simulation<-
    sink(path_newdir&'/MODELSETTINGS_'&sim_folder&'_'&id&'.txt',append=FALSE)
  sim_folder
  generated
  # modelversion
  print(ls.str(),max.level=0)
  sink()
  
  gc()
  
  ################################################################################
  ################################################################################
  #===============================================================================
  # DISTRIBUTED MELT MODEL
  #===============================================================================
  ################################################################################
  ################################################################################
  stop<-FALSE
  
  # LOOP FOR INTERMEDIATE GEOMETRY UPDATES
  for(loop in 1:length(iOutput_ts)){
    ## for testing only:
    # loop=1
    # loop=2
    
    ifelse(loop == 1,
           TS<-1:iOutput_ts[loop],
           TS<-(iOutput_ts[loop-1]+1):iOutput_ts[loop]
    )         
    
    # YEAR
    idx_yr<-1900+unique(as.POSIXlt(Date[TS])$year)
    
    #===============================================================================
    # CLIFF ALTITUDE
    #===============================================================================
    ### MEDIAN CLIFF ELEVATION ###
    # calculate median elevation per cliff
    MedianCliffEle<-aggregate(cliffs_geom_df[,'Elevation'],
                              list(cliffs_geom_df$ID),median)
    colnames(MedianCliffEle)<-c('ID','ele')
    # write median cliff elevations to file for documentation
    MedianCliffEle_tab<-within(MedianCliffEle,
                               {ele<-formatC(ele,format='f',digits=2)})
    filename_medelev<-path_outtab&'/MedianCliffEle_ts'&TS[1]-1&'_tab.txt'
    write.table(MedianCliffEle,filename_medelev,dec='.',sep=',',quote=FALSE,
                row.names=FALSE)
    rm(filename_medelev,MedianCliffEle_tab)
    
    
    #===============================================================================
    # SHORTWAVE RADIATION
    #===============================================================================
    
    library(cleaRskyQuantileRegression)
    
    I_E<-calc_PotRadiation_CosineResponsePower(doy = doy[TS],
                                               hour = hour[TS],
                                               latDeg = 28.21081,
                                               longDeg = 85.56169,
                                               timeZone = 0,
                                               isCorrectSolartime = FALSE,
                                               cosineResponsePower = 1.2 )
                                               
    
    # SWin can not be larger than I_E (in cases where this happens,
    #   numerical probelms occur, leading to model inaccuracies)
    SW_IN<-SWin[TS]
    I_in_obs<-pmin(SW_IN,I_E)
    # SWin can not be smaller than 0
    I_in_obs<-pmax(0,I_in_obs)
    
    # k_t:
    # Clearness index
    k_t<-I_in_obs / I_E
    # to avoid wrong diffuse fraction calculations
    k_t[is.na(k_t)]<-0
    
    ##
    # diffuse fraction calculations (Reindl-Model)
    ##
    # k_d:
    # Diffuse fraction
    
    # create empty k_d vector
    k_d<-k_t*0-9999
    
    # condition 1:
    idxCond1<-which(k_t <= 0.3)
    k_d[idxCond1]<-1.02-0.254*k_t[idxCond1]+0.0123*SIN_h[idxCond1]
    
    # condition 2:
    idxCond2<-which(k_t > 0.3 & k_t <= 0.78)
    k_d[idxCond2]<-1.4-1.749*k_t[idxCond2]+0.177*SIN_h[idxCond2]
    
    # condition 3:
    idxCond3<-which(k_t > 0.78)
    k_d[idxCond3]<-0.486*k_t[idxCond3]-0.182*SIN_h[idxCond3]
    
    ##
    # D_h:
    # Diffuse irradiance to a horizontal surface
    D_h<-k_d*I_in_obs
    # I_b:
    # Direct normal irradiance
    I_b<-(I_in_obs-D_h)/SIN_h
    
    # predefine rasters for radiation calculation (values in radians)
    cliffs_geom_df$COS_aspect<-cos(cliffs_geom_df$Aspect*pi/180)
    cliffs_geom_df$COS_slope<-cos(cliffs_geom_df$Slope*pi/180)
    cliffs_geom_df$SIN_slope<-sin(cliffs_geom_df$Slope*pi/180)
    cliffs_geom_df$SIN_aspect<-sin(cliffs_geom_df$Aspect*pi/180)
    
    # COSinci:
    # Cosine of solar incidence angle relative to
    #   vector normal to the (sloped) surface (Garnier & Ohmura, 1968)
    #  (as a factor to reduce/enhance radiation at a certain point in time
    
    A<- SIN_phi*COS_omega
    B<- -1*cliffs_geom_df$COS_aspect*cliffs_geom_df$SIN_slope
    C<- SIN_omega
    D<- cliffs_geom_df$SIN_aspect*cliffs_geom_df$SIN_slope
    E<- COS_phi*COS_omega
    F<- cliffs_geom_df$COS_slope
    G<- COS_delta
    H<- COS_phi
    I<- cliffs_geom_df$COS_aspect*cliffs_geom_df$SIN_slope
    K<- SIN_phi*cliffs_geom_df$COS_slope
    L<- SIN_delta
    
    AB<-sapply(B,function(x) A*x)
    CD<-sapply(D,function(x) C*x)               
    EF<-sapply(F,function(x) E*x)               
    
    Part1<-(AB-CD+EF)*G
    Part2<-sapply((H*I+K),function(x) x*L)
    
    COSinci_m<-Part1+Part2
    
    # I_b*COSinci:
    #  reduction/enhancement of the direct solar radiation depending on
    #  the deviation of the surface to the direction of the sun beam
    # I_s:
    # Direct solar irradiance from sky
    #  (reduced/enhanced depending on the deviation of the surface
    #  to the direction of the sun beam
    I_s_m<-I_b*COSinci_m
    
    # convert solar elevation angle from radians into degree for comparison with
    #  horizon angle
    hh<-h*180/pi
    
    # find at which timesteps sun elevation is lower than topography horizon
    #  ->hh [ timesteps(h) ] 
    #  ->solaraz [ timesteps(h) ] 
    #  ->horI_m [ pixels(n), aspect(deg) ]
    #  ->I_s_m [ timesteps(h), pixels(n) ]
    #IMPROVE/REPLACE LOOP (EXPENSIVE)
    for(t in 1:length(hh)){
      hSun<-hh[t]
      azSun<-round(360-solaraz[t])
      for(px in 1:nrow(horI_m)){
        hTopo<-horI_m[px,azSun]
        if(hTopo > hSun){I_s_m[t,px]<-0}
      }
    }
    
    # I_s radiation can not be negative
    I_s_m[I_s_m < 0]<-0
    
    ##TESTING
    # plot(360-round(solaraz[1:72]),type='n')
    # lines(360-round(solaraz[1:72]),col='orange',lwd=5)
    # plot(horI_m[102,],type='n')
    # lines(horI_m[102,],col='orange',lwd=5)
    # plot(horI_m[102,90:180],type='n')
    # lines(horI_m[102,90:180],col='orange',lwd=5)
    # I_s<-cbind(cliffs_geom_df[,c('x','y')],data.frame(t(I_s_m[7:17,])))
    # I_s<-rasterFromXYZ(I_s[,c('x','y','X2')])
    # projection(I_s)<-projec
    # plot(I_s)
    # cellStats(I_s,'mean')
    # pt<-cbind(cliffs_geom_df[,c('x','y')],rep(NA,nrow(cliffs_geom_df)))
    # pt[102,3]<-1000
    # pt<-rasterFromXYZ(pt)
    # projection(tt)<-projec
    # plot(I_s)
    # plot(pt,col='red',legend=FALSE,add=TRUE)
    
    # D_s:
    # Diffuse irradiance from sky
    D_s_m<-sapply(cliffs_geom_df$VsI,function(x) D_h*x)
    
    # D_t:
    # Diffuse irradiance from terrain
    D_t_m<-sapply(cliffs_geom_df$VsI,function(x) alpha_d*I_in_obs*(1-x))
    
    # I_n:
    # Net shortwave flux (normal to cliff face)
    I_n_m<-(I_s_m+D_s_m+D_t_m)*(1-alpha_i)
    
    # I_in:
    # Incoming shortwave flux (normal to cliff face)
    I_in_m<-I_s_m+D_s_m+D_t_m
    
    
    #===============================================================================
    # LONGWAVE RADIATION
    #===============================================================================
    # e_s_air:
    # Saturation vapour pressure of air (Teten's equation), uses ?C
    # Monteith and Unsworth (2008)
    e_s_air_df<-as.data.frame(610.8*exp(17.27*Ta[TS]/(237.3+Ta[TS]))) #[Pa]
    
    for(column in 1:length(unique(cliffs_geom_df$ID))){
      colnames(e_s_air_df)[column]<-'Cliff_'&sort(unique(cliffs_geom_df$ID))[column]
    }  
    
    # e_a:     
    # Actual vapour pressure of air (from rH and e_s_air)
    e_a_df<-rH[TS]*e_s_air_df/(100*1000) # [kPa]
    
    # L_d:
    # Longwave from surrounding terrain (debris)
    #L_t(x) = epsilon_d*sigma*(T_s(x)+273.15)^4*(1-V_s)
    TS_4<-as.data.frame((Ts[TS]+273.15)^4)
    # specify coluumn names (needed for cases where nCliffs=1)
    for(column in 1:length(unique(cliffs_geom_df$ID))){
      colnames(TS_4)[column]<-'Cliff_'&sort(unique(cliffs_geom_df$ID))[column]
    } 
    
    L_d_m_ls<-list()				 
    for(i in 1:ncol(TS_4)){
      cliffid<-sort(unique(cliffs_geom_df$ID))[i]
      idx<-which(cliffs_geom_df$ID == cliffid)
      Ld_tsteps<-list()
      for(tstep in 1:nrow(TS_4)){
        Ld_tsteps[[tstep]]<-cliffs_geom_df[idx,'Vd']*
          epsilon_d*sigma*TS_4[tstep,i]
      }
      L_d_m_ls[[i]]<-do.call('rbind',Ld_tsteps) 
    }
    # ->L_d_m [ timesteps(h), pixels(n) ]
    L_d_m<-as.matrix(do.call('cbind',L_d_m_ls))	
    
    ## TESTING
    # Ld<-cbind(cliffs_geom_df[,c('x','y')],data.frame(t(L_d_m[1:24,])))
    # Ld<-rasterFromXYZ(Ld[,c('x','y','X1')])
    # projection(Ld)<-projec
    # plot(Ld)
    
    # L_o:
    # Longwave outgoing from ice surface
    L_o<-epsilon_i*sigma*(Ti+273.15)^4
    
    # L_s:
    # Longwave from sky   
    L_s_m_ls<-list()				 
    for(i in 1:ncol(TS_4)){
      cliffid<-sort(unique(cliffs_geom_df$ID))[i]
      idx<-which(cliffs_geom_df$ID == cliffid)
      Ls_tsteps<-list()
      for(tstep in 1:length(TS)){
        Ls_tsteps[[tstep]]<-cliffs_geom_df[idx,'VsL']*LWin[TS][tstep]
      }
      L_s_m_ls[[i]]<-do.call('rbind',Ls_tsteps) 
    }
    # ->L_s_m [ timesteps(h), pixels(n) ]
    L_s_m<-as.matrix(do.call('cbind',L_s_m_ls))	
    dimnames(L_s_m)<-NULL
    
    # L_n:
    # Net longwave radiation
    L_n_m<-L_s_m+L_d_m-L_o
    
    
    #===============================================================================
    # TURBULENT FLUXES
    #===============================================================================
    # lambda_f(x) = 0.05*sin(beta)*(1+cos(A-W_D(x)));
    # z_0(x) = 0.5*h_c*lambda_f(x);	# according to Lettau (1969) modified
    # for high-relief terrain --> Lettau's approach not used, conceptual
    # error - see thesis JS for details
    
    # z_0:
    # Roughness length for momentum [m]
    # ->Literature value (0.001 m, e.g. Brock 2006 or Pelliciotti 2005)
    
    # ustar:
    # Friction velocity
    ustar<-ws[TS]*k/(log(z/z0))
    
    # Re:
    # Reynolds number
    Re<-ustar*z0/nu
    
    # z_0T:
    # Roughness length for temperature profiles [m]
    # z_0e:
    # Roughness length for a logarithmic profile of vapor pressure [m]
    #  (z_0T & z_0e are assumed to be equal and are derived
    #  using Andreas'(1987) model)
    LOG_Re<-log(Re)
    z0T<-z0*exp(0.317-0.565*LOG_Re-0.183*(LOG_Re)^2)
    
    z0e<-z0T
    
    #     z0(x) = 0.0012; # turned off by Reid
    #     z0T(x) = z0(x); # ignoring Andreas model # turned off by Reid
    #     z0e(x) = z0(x); # ignoring Andreas model # turned off by Reid
    
    # K_s:
    # Exchange coefficient for a neutral boundary layer with
    #  logarithmic profiles for wind speed and temperature
    LOG_zz0<-log(z/z0)
    LOG_zz0T<-log(z/z0T)
    LOG_zz0e<-log(z/z0e)
    K_s<-c_p*k^2*rho_0/(P_0*LOG_zz0*LOG_zz0T)
    
    # P:
    # Air pressure [kPa], based on altitude
    P<-P_0*((1-(Env_LR*MedianCliffEle[,'ele']/T_0))^(g*Mair/(Rgas*Env_LR)))
    names(P)<-'Cliff_'&sort(unique(cliffs_geom_df$ID))
    # K_L:
    # Exchange coefficient for latent heat
    K_L<-0.623*L_e*k^2*rho_0/(P_0*LOG_zz0*LOG_zz0e)
    
    # H:
    # Sensible heat flux
    H1<-((Ta[TS]+273.15)-(Ti+273.15))*ws[TS] 
    H2<-sapply(P['Cliff_'&cliffs_geom_df$ID],function(x) K_s*x)
    H_m<-data.matrix(H1*H2)
    dimnames(H_m)<-NULL
    
    # LE:
    # Latent heat flux
    LE_m<-data.matrix(K_L*(e_a_df[,'Cliff_'&cliffs_geom_df$ID]-e_s)*ws[TS])
    dimnames(LE_m)<-NULL
    
    
    #===============================================================================
    # MELT
    #===============================================================================
    # All fluxes are calculated perpendicular to the surface.
    
    # Q_m:
    # Heat for ice ablation perpendicular to cliff surface
    Q_m_m<-I_n_m+L_n_m+H_m+LE_m
    
    # melt:
    # Hourly melt rate [m w.e.] perpendicular to cliff surface
    melt_m<-(Q_m_m*3600/(rho_i*L_f))
    # avoid negative melt values (resulting from neg. Q_m-budget)
    melt_m[melt_m<0]<-0            
    
    # CumuMelt:
    # Cumulative vertical melt [m w.e.] over all timesteps with the same geometry 
    #  (e.g. 1 week/month etc.) until geometric correction
    CumuMelt_v<-as.vector(base::colSums(melt_m,na.rm=T))
    Q_m_m_v<-as.vector(base::colSums(Q_m_m,na.rm=T))
    I_n_m_v<-as.vector(base::colSums(I_n_m,na.rm=T))
    I_s_m_v<-as.vector(base::colSums(I_s_m,na.rm=T))
    D_s_m_v<-as.vector(base::colSums(D_s_m,na.rm=T))
    D_t_m_v<-as.vector(base::colSums(D_t_m,na.rm=T))
    L_n_m_v<-as.vector(base::colSums(L_n_m,na.rm=T))
    L_s_m_v<-as.vector(base::colSums(L_s_m,na.rm=T))
    L_d_m_v<-as.vector(base::colSums(L_d_m,na.rm=T))
    H_m_v<-as.vector(base::colSums(H_m,na.rm=T))
    LE_m_v<-as.vector(base::colSums(LE_m,na.rm=T))
    cliffs_geom_df$Cummelt<-CumuMelt_v
    cliffs_geom_df$Q_m_m<-Q_m_m_v
    cliffs_geom_df$I_n_m<-I_n_m_v
    cliffs_geom_df$I_s_m<-I_s_m_v
    cliffs_geom_df$D_s_m<-D_s_m_v
    cliffs_geom_df$D_t_m<-D_t_m_v
    cliffs_geom_df$L_n_m<-L_n_m_v
    cliffs_geom_df$L_s_m<-L_s_m_v
    cliffs_geom_df$L_d_m<-L_d_m_v
    cliffs_geom_df$H_m<-H_m_v
    cliffs_geom_df$LE_m<-LE_m_v
    
    # dt [days] between last update and actual date
    ddaysMelt<-as.numeric(Date[TS[length(TS)]]-Date[TS[1]])
    
    # Convert back from "m w.e" 
    #  to "m" distance in ice:
    CumuMelt_dist<-cliffs_geom_df$Cummelt/(rho_i/1000)
    #summary(CumuMelt_dist)
    
    
    #=========================================================================
    # CALCULATE VOLUME LOSS OVER WHOLE CLIFF 
    #  [VERTICAL MELT x PIXEL AREA (HORIZONTALLY-PROJECTED)]
    #=========================================================================
    # calculate projected area of a pixel [m2] (as seen in dem1)
    projAreaPPx<-resol1*resol1 #[m2]
    
    # water volume loss per day [m^3 w.e.]
    cliffs_geom_df$VolLossPDay<-(cliffs_geom_df$Cummelt*projAreaPPx)/ddaysMelt
    #summary(cliffs_geom_df$VolLossPDay) 
    
    
    #=============================================================================== 
    # MELT PER DAY PER PIXEL (PERPENDICLAR TO ICE CLIFF)
    #===============================================================================
    cliffs_geom_df$MPD<-cliffs_geom_df$Cummelt/ddaysMelt #[m w.e./d]
    #summary(cliffs_geom_df$MPD)
    
    
    #=============================================================================== 
    # INCLINED (TRUE) AREA OF ICE CLIFF PIXELS (USED FOR STATISTICS ONLY)
    #===============================================================================
    # calculate real inclined area for each pixel [m2]
    cliffs_geom_df$inclArea<-projAreaPPx/cos(cliffs_geom_df$Slope*pi/180)
    #summary(cliffs_geom_df$inclArea)
    
    
    #=============================================================================== 
    # dz: VERTICAL SHIFT [m] DUE TO ATMOSPHERIC MELT
    #===============================================================================
    # dz:
    # Raster of vertical shift [m, distance in ice] due to melt
    #  (slope in radians)
    cliffs_geom_df$dz<- -1*CumuMelt_dist*cos(cliffs_geom_df$Slope*pi/180)
    #summary(cliffs_geom_df$dz)
    
    
    #=============================================================================== 
    # dxy: HORIZONTAL SHIFT [m] DUE TO ATMOSPHERIC MELT
    #===============================================================================
    #  (slope in radians)
    cliffs_geom_df$dxy<-CumuMelt_dist*sin(cliffs_geom_df$Slope*pi/180)
    #summary(cliffs_geom_df$dxy)
    
    
    #============================================================================= 
    # z: ELEVATION [m] AFTER VERTICAL MELT
    #=============================================================================
    cliffs_geom_df$z<-cliffs_geom_df$Elevation+cliffs_geom_df$dz
    
    
    #============================================================================== 
    # bA_v: BACKAZIMUTH ANGLES (OPPOSITE OF CLIFF AZIMUTH) [rad]
    #==============================================================================
    # corr_v (in degrees):
    # Correction vector with values 180?/-180? where azimuth is smaller/bigger
    #  than 180?
    corr_v<-cliffs_geom_df$Aspect
    corr_v[cliffs_geom_df$Aspect >= 180 & cliffs_geom_df$Aspect <= 360]<- -180
    corr_v[cliffs_geom_df$Aspect >= 0   & cliffs_geom_df$Aspect <  180]<-  180
    
    # bA_v:
    bA_v<-(cliffs_geom_df$Aspect + corr_v)*pi/180
    rm(corr_v)
    
    
    #============================================================================== 
    # dx: DISPLACEMENT IN EAST DIRECTION [m]
    #==============================================================================
    cliffs_geom_df$dx<-cliffs_geom_df$dxy*sin(bA_v)
    #summary(cliffs_geom_df$dx)
    
    
    #============================================================================== 
    # dx: DISPLACEMENT IN NORTH DIRECTION [m]
    #==============================================================================
    cliffs_geom_df$dy<-cliffs_geom_df$dxy*cos(bA_v)
    #summary(cliffs_geom_df$dy)
    
    # create matrix with old cliff values
    xyz_m<-cbind(cliffs_geom_df$x,cliffs_geom_df$y,cliffs_geom_df$Elevation,
                 cliffs_geom_df$ID)
    colnames(xyz_m)<-c('x','y','elevation','ID')
    
    
    # compute new coordinates and add its elevations (for all raster cells)
    new_xyz_m<-cbind(cliffs_geom_df$x+cliffs_geom_df$dx,
                     cliffs_geom_df$y+cliffs_geom_df$dy,
                     cliffs_geom_df$z,cliffs_geom_df$ID)
    colnames(new_xyz_m)<-c('x','y','elevation','ID')
    
    # put coordinates into raster to drop points which are in the same pixel
    cliffID_r<-gl_geom_st[[1]]
    cliffID_r<-setValues(cliffID_r,NA)
    cliffID_r2<-cliffID_r
    cliffID_r<-rasterize(new_xyz_m[,1:2],cliffID_r,new_xyz_m[,'elevation'])
    cliffID_r2<-rasterize(new_xyz_m[,1:2],cliffID_r2,new_xyz_m[,'ID'])
    new_xyz_clean_m<-xyFromCell(cliffID_r,1:ncell(cliffID_r))
    elev<-as.vector(cliffID_r)
    ids<-as.vector(cliffID_r2)
    new_xyz_clean_m<-cbind(new_xyz_clean_m,elev,ids)
    new_xyz_clean_m<-new_xyz_clean_m[complete.cases(new_xyz_clean_m),]      
    colnames(new_xyz_clean_m)<-c('x','y','elevation','ID')
    
    # old cliff cells
    pts1<-SpatialPoints(cbind(cliffs_geom_df$x,cliffs_geom_df$y))
    projection(pts1)<-projec
    pts1_data<-as.data.frame(cliffs_geom_df$ID)
    colnames(pts1_data)<-'Cliff_ID'
    pts1_spdf<-SpatialPointsDataFrame(pts1,pts1_data)
    projection(pts1_spdf)<-projec
    
    
    #=========================================================================
    # CALCULATE STATISTICS PER CLIFF (PC)
    #=========================================================================
    ### AVERAGES/TOTALS PER CLIFF
    # mean melt per day [m w.e./d]
    MPD_pc<-aggregate(cliffs_geom_df$Cummelt,list(cliffs_geom_df$ID),
                      mean)
    MPD_pc$x<-MPD_pc$x/ddaysMelt
    
    # cumulative water volume loss [m3 w.e.]
    VolLoss_pc<-aggregate(cliffs_geom_df$Cummelt,list(cliffs_geom_df$ID),
                          sum)
    VolLoss_pc$x<-VolLoss_pc$x*projAreaPPx
    
    # total projected area [m2]
    projA_v<-table(cliffs_geom_df$ID)*projAreaPPx
    
    # total inclined area [m2]
    inclA_v<-projAreaPPx/cos(cliffs_geom_df$Slope*pi/180)
    inclA_pc<-aggregate(inclA_v,list(cliffs_geom_df$ID),sum)
    
    # mean aspect & slope [?] (weighted vectorial mean)
    ids<-unique(cliffs_geom_df$ID)
    Slope_v<-vector()
    Aspect_v<-vector()
    for(i in 1:length(ids)){
      idx<-which(cliffs_geom_df$ID == ids[i])
      # mean slope weighted with inclined area
      Slope_v[i]<-sum(cliffs_geom_df[idx,'Slope']*inclA_v[idx])/
        sum(inclA_v[idx])
      
      # vectorial mean aspect, weighted with inclined area
      rads<-cliffs_geom_df[idx,'Aspect']*pi/180    
      avsin<-sum(sin(rads)*inclA_v[idx])/sum(inclA_v[idx])
      avcos<-sum(cos(rads)*inclA_v[idx])/sum(inclA_v[idx])
      avaz<-atan2(avsin,avcos)
      if (avaz<0){avaz<-avaz+(2*pi)}
      Aspect_v[i]<-avaz*180/pi
    }
    rm(ids)    
    
    # mean elevation [m a.s.l.]
    Elevation_pc<-aggregate(cliffs_geom_df$Elevation,
                            list(cliffs_geom_df$ID),mean)
    
    # mean SW skyview factor [-]
    VsI_pc<-aggregate(cliffs_geom_df$VsI,list(cliffs_geom_df$ID),mean)
    
    # mean LW skyview factor [-]
    VsL_pc<-aggregate(cliffs_geom_df$VsL,list(cliffs_geom_df$ID),mean)
    
    # mean debrisview factor [-]
    Vd_pc<-aggregate(cliffs_geom_df$Vd,list(cliffs_geom_df$ID),mean)
    
    
    #=========================================================================
    # OUTPUT: MEAN GEOMETRY/MPD & TOTAL VOLUME LOSS PER CLIFF (PC)
    #=========================================================================
    Stats_df<-as.data.frame(cbind('Cliff_ID'=sort(unique(cliffs_geom_df$ID)),
                                  'Elevation'=Elevation_pc$x,
                                  'Slope'=Slope_v,
                                  'Aspect'=Aspect_v,
                                  'projA'=projA_v,
                                  'inclA'=inclA_pc$x,                                   
                                  'VsI'=VsI_pc$x,
                                  'VsL'=VsL_pc$x,
                                  'Vd'=Vd_pc$x,
                                  'MPD'=MPD_pc$x,
                                  'VolLoss'=VolLoss_pc$x))
    
    # define format of output
    Stats_df<-within(Stats_df,
                     {Cliff_ID<-formatC(Cliff_ID,format='f',digits=0)
                     Elevation<-formatC(Elevation,format='f',digits=2)
                     Slope<-formatC(Slope,format='f',digits=2)
                     Aspect<-formatC(Aspect,format='f',digits=2)
                     projA<-formatC(projA,format='f',digits=2)
                     inclA<-formatC(inclA,format='f',digits=2)
                     VsI<-formatC(VsI,format='f',digits=2)
                     VsL<-formatC(VsL,format='f',digits=2)
                     Vd<-formatC(Vd,format='f',digits=2)
                     MPD<-formatC(MPD,format='f',digits=4)
                     VolLoss<-formatC(VolLoss,format='f',digits=2)
                     })
    
    # units
    units_v<-c('[-]','[m a.s.l.]','[?]','[?]','[m2]','[m2]','[-]','[-]','[-]',
               '[m w.e./d]','[m3 w.e.]')
    Stats_df<-rbind(units_v,Stats_df)
    
    # write to text file
    filename_Stats<-path_newdir&'/PC_Stats_ts'&TS[1]&'to'&TS[length(TS)]&'_'&id&'_'&
      round(resol1,1)&'m.txt'
    write.table(Stats_df,filename_Stats,quote=FALSE,row.names=FALSE,
                dec='.',sep='\t')                   
    
    
    #=========================================================================
    # OUTPUT: MELT RASTER & FLUXES CONTRIBUTION RASTER
    #=========================================================================   
    # create raster from dataframe
    MPD_r<-gl_geom_st[[1]]
    MPD_r<-setValues(MPD_r,NA)
    names(MPD_r)<-'MPD'
    MPD_r<-rasterize(cliffs_geom_df[,c('x','y')],
                     MPD_r,cliffs_geom_df[,'MPD'])
    
    filename_newMPD<-path_outrast&'/MPD_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(MPD_r,filename=filename_newMPD,format='GTiff',overwrite=TRUE)
    
    Q_m_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'Q_m_m'])
    
    filename_newQ_m_m<-path_outrast&'/Q_m_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(Q_m_m_r,filename=filename_newQ_m_m,format='GTiff',overwrite=TRUE)
    
    I_n_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'I_n_m'])
    
    filename_newI_n_m<-path_outrast&'/I_n_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(I_n_m_r,filename=filename_newI_n_m,format='GTiff',overwrite=TRUE)
    
    I_s_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'I_s_m'])
    
    filename_newI_s_m<-path_outrast&'/I_s_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(I_s_m_r,filename=filename_newI_s_m,format='GTiff',overwrite=TRUE)
    
    D_s_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'D_s_m'])
    
    filename_newD_s_m<-path_outrast&'/D_s_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(D_s_m_r,filename=filename_newD_s_m,format='GTiff',overwrite=TRUE)
    
    D_t_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'D_t_m'])
    
    filename_newD_t_m<-path_outrast&'/D_t_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(D_t_m_r,filename=filename_newD_t_m,format='GTiff',overwrite=TRUE)
    
    L_n_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'L_n_m'])
    
    filename_newL_n_m<-path_outrast&'/L_n_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(L_n_m_r,filename=filename_newL_n_m,format='GTiff',overwrite=TRUE)
    
    L_s_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'L_s_m'])
    
    filename_newL_s_m<-path_outrast&'/L_s_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(L_s_m_r,filename=filename_newL_s_m,format='GTiff',overwrite=TRUE)
    
    L_d_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'L_d_m'])
    
    filename_newL_d_m<-path_outrast&'/L_d_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(L_d_m_r,filename=filename_newL_d_m,format='GTiff',overwrite=TRUE)
    
    H_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                     MPD_r,cliffs_geom_df[,'H_m'])
    
    filename_newH_m<-path_outrast&'/H_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(H_m_r,filename=filename_newH_m,format='GTiff',overwrite=TRUE)
    
    LE_m_r<-rasterize(cliffs_geom_df[,c('x','y')],
                      MPD_r,cliffs_geom_df[,'LE_m'])
    
    filename_newLE_m<-path_outrast&'/LE_m_ts'&TS[length(TS)]&'_r.tif'             
    writeRaster(LE_m_r,filename=filename_newLE_m,format='GTiff',overwrite=TRUE)
    
    
    if(loop == 1){
      #=========================================================================
      # OUTPUT: VIEW ANGLES RASTER
      #========================================================================= 
      VsI_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'VsI'])
      
      filename_newVsI<-path_outrast&'/VsI_ts'&TS[length(TS)]&'_r.tif'             
      writeRaster(VsI_r,filename=filename_newVsI,format='GTiff',overwrite=TRUE)
      
      VsL_r<-rasterize(cliffs_geom_df[,c('x','y')],
                       MPD_r,cliffs_geom_df[,'VsL'])
      
      filename_newVsL<-path_outrast&'/VsL_ts'&TS[length(TS)]&'_r.tif'             
      writeRaster(VsL_r,filename=filename_newVsL,format='GTiff',overwrite=TRUE)
      
      Vd_r<-rasterize(cliffs_geom_df[,c('x','y')],
                      MPD_r,cliffs_geom_df[,'Vd'])
      
      filename_newVd<-path_outrast&'/Vd_ts'&TS[length(TS)]&'_r.tif'             
      writeRaster(Vd_r,filename=filename_newVd,format='GTiff',overwrite=TRUE)
    }
    
    
    #===============================================================================
    # OUTPUT: HOURLY VALUES OF NON-GRIDDED MODEL VARIABLES
    #=============================================================================== 
    # 'modVariables':
    modVariables_df<-data.frame(cbind(Gamma,E_0,omega,delta,theta,solaraz,I_E,
                                      k_t,k_d,D_h,I_b,hh,ustar,Re,z0T,z0e,K_s,K_L))
    
    # solar elevation angle: can not be named 'h' in excel file
    #  as 'H' already used!  ->'SEA'                            
    colnames(modVariables_df)<-c('Gamma','E_0','omega','delta','theta','solaraz',
                                 'I_E','k_t','k_d','D_h','I_b','SEA','ustar','Re',
                                 'z0T','z0e','K_s','K_L')
    
    # define format of output
    modVariables_df<-within(modVariables_df,
                            {Gamma<-formatC(Gamma,format='f',digits=4)
                            E_0<-formatC(E_0,format='f',digits=4)
                            omega<-formatC(omega,format='f',digits=4)
                            delta<-formatC(delta,format='f',digits=4)
                            theta<-formatC(theta,format='f',digits=4)
                            solaraz<-formatC(solaraz,format='f',digits=4)
                            I_E<-formatC(I_E,format='f',digits=4)
                            k_t<-formatC(k_t,format='f',digits=4)
                            k_d<-formatC(k_d,format='f',digits=4)
                            D_h<-formatC(D_h,format='f',digits=4)
                            I_b<-formatC(I_b,format='f',digits=4)
                            SEA<-formatC(SEA,format='f',digits=4)
                            ustar<-formatC(ustar,format='f',digits=4)
                            Re<-formatC(Re,format='f',digits=4)
                            z0T<-formatC(z0T,format='f',digits=4)
                            z0e<-formatC(z0e,format='f',digits=4)
                            K_s<-formatC(K_s,format='f',digits=4)
                            K_L<-formatC(K_L,format='f',digits=4)
                            })
    
    modVariables_df<-cbind(as.character(Date[TS]),modVariables_df)
    colnames(modVariables_df)[1]<-'Date'
    
    filename_modVar<-path_outtab&'/Stats_NonGriddedData_ts'&TS[1]&'to'&
      TS[length(TS)]&'_'&id&'_'&round(resol1,digits=1)&'m.txt'
    write.table(modVariables_df,filename_modVar,dec='.',sep='\t',
                quote=FALSE,row.names=FALSE)
    
    
    #===============================================================================
    # OUTPUT: HOURLY VALUES OF DISTRIBUTED ENERGY FLUXES
    #         ->LIST OF MATRICES (TS x CLIFFPIXELS)
    #===============================================================================
    distribEnergyFluxes_ls<-list(I_s_m,D_s_m,D_t_m,I_n_m,I_in_m,L_s_m,L_d_m,
                                 H_m,LE_m,Q_m_m,melt_m)
    
    names(distribEnergyFluxes_ls)<-c('I_s','D_s','D_t','I_n','I_in',
                                     'L_s','L_d','H','LE','Q_m','melt')
    
    # calculate hourly mean per cliff
    meanFluxes_ls<-list()
    for(i in 1:length(distribEnergyFluxes_ls)){
      meanFluxes_ls[[i]]<-rowMeans(distribEnergyFluxes_ls[[i]])
    }
    #meanFluxes_df<-data.frame(distribEnergyFluxes_ls[[1]],do.call(cbind,meanFluxes_ls))
    meanFluxes_df<-data.frame(do.call(cbind,meanFluxes_ls))
    colnames(meanFluxes_df)<-names(distribEnergyFluxes_ls)
    rm(meanFluxes_ls,distribEnergyFluxes_ls)
    
    # check
    # plot(meanFluxes_df[1:100,2],type='o')
    
    # define format of output
    meanFluxes_df<-within(meanFluxes_df,
                          {I_s<-formatC(I_s,format='f',digits=4)
                          D_s<-formatC(D_s,format='f',digits=4)
                          D_t<-formatC(D_t,format='f',digits=4)
                          I_n<-formatC(I_n,format='f',digits=4)
                          I_in<-formatC(I_in,format='f',digits=4)
                          L_s<-formatC(L_s,format='f',digits=4)
                          L_d<-formatC(L_d,format='f',digits=4)
                          H<-formatC(H,format='f',digits=4)
                          LE<-formatC(LE,format='f',digits=4)
                          Q_m<-formatC(Q_m,format='f',digits=4)
                          melt<-formatC(melt,format='f',digits=4)
                          })
    
    meanFluxes_df<-cbind(as.character(Date[TS]),meanFluxes_df)
    colnames(meanFluxes_df)[1]<-'Date'
    
    filename_meanFluxes<-path_outtab&'/Stats_MeanFluxes_ts'&TS[1]&'to'&
      TS[length(TS)]&'_'&id&'_'&round(resol1,digits=1)&'m.txt'
    write.table(meanFluxes_df,filename_meanFluxes,dec='.',sep='\t',
                quote=FALSE,row.names=FALSE)
    
    
    #===============================================================================
    # OUTPUT: DIURNAL CYCLES OF ENERGY FLUXES (AVERAGED PER CLIFF & PER PIXEL)
    #===============================================================================
    source(path_func)
    # Direct SW (DC averaged for entire cliff)
    I_s_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],I_s_m)   
    I_s_dc<-do.call('cbind',I_s_dc_ls[['DC']])
    I_s_sd<-do.call('cbind',I_s_dc_ls[['SD']])
    # Direct SW (DC per cliff pixel)
    li<-apply(I_s_m,2, function(x) DC(cbind.data.frame(Date=Date[TS],x)))
    li<-stripname(li,'SD')
    li<-stripname(li,'hours')
    li<-do.call(Map,c(f=cbind,li))
    I_s_dc_pp<-do.call('cbind',li[['DC']])
    rm(li)
    
    # Diffuse SW from sky (DC averaged for entire cliff)
    D_s_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],D_s_m)
    D_s_dc<-do.call('cbind',D_s_dc_ls[['DC']])
    D_s_sd<-do.call('cbind',D_s_dc_ls[['SD']])
    # Diffuse SW (DC per cliff pixel)
    li<-apply(D_s_m,2, function(x) DC(cbind.data.frame(Date=Date[TS],x)))
    li<-stripname(li,'SD')
    li<-stripname(li,'hours')
    li<-do.call(Map,c(f=cbind,li))
    D_s_dc_pp<-do.call('cbind',li[['DC']])
    rm(li)
    
    # Diffuse SW reflected from terrain (DC averaged for entire cliff)
    D_t_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],D_t_m)
    D_t_dc<-do.call('cbind',D_t_dc_ls[['DC']])
    D_t_sd<-do.call('cbind',D_t_dc_ls[['SD']])
    # Diffuse SW reflected from terrain (DC per cliff pixel)
    li<-apply(D_t_m,2, function(x) DC(cbind.data.frame(Date=Date[TS],x)))
    li<-stripname(li,'SD')
    li<-stripname(li,'hours')
    li<-do.call(Map,c(f=cbind,li))
    D_t_dc_pp<-do.call('cbind',li[['DC']])
    rm(li)
    
    # Net SW (DC averaged for entire cliff)
    I_n_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],I_n_m)
    I_n_dc<-do.call('cbind',I_n_dc_ls[['DC']])
    I_n_sd<-do.call('cbind',I_n_dc_ls[['SD']])
    # Net SW (DC per cliff pixel)
    li<-apply(I_n_m,2, function(x) DC(cbind.data.frame(Date=Date[TS],x)))
    li<-stripname(li,'SD')
    li<-stripname(li,'hours')
    li<-do.call(Map,c(f=cbind,li))
    I_n_dc_pp<-do.call('cbind',li[['DC']])
    rm(li)
    
    # Total incoming SW (DC averaged for entire cliff)
    I_in_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],I_in_m)
    I_in_dc<-do.call('cbind',I_in_dc_ls[['DC']])
    I_in_sd<-do.call('cbind',I_in_dc_ls[['SD']])
    # Total incoming SW (DC per cliff pixel)
    li<-apply(I_in_m,2, function(x) DC(cbind.data.frame(Date=Date[TS],x)))
    li<-stripname(li,'SD')
    li<-stripname(li,'hours')
    li<-do.call(Map,c(f=cbind,li))
    I_in_dc_pp<-do.call('cbind',li[['DC']])
    rm(li)
    
    # LW from sky (DC averaged for entire cliff)
    L_s_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],L_s_m)
    L_s_dc<-do.call('cbind',L_s_dc_ls[['DC']])
    L_s_sd<-do.call('cbind',L_s_dc_ls[['SD']])
    # LW from sky (DC per cliff pixel)
    li<-apply(L_s_m,2, function(x) DC(cbind.data.frame(Date=Date[TS],x)))
    li<-stripname(li,'SD')
    li<-stripname(li,'hours')
    li<-do.call(Map,c(f=cbind,li))
    L_s_dc_pp<-do.call('cbind',li[['DC']])
    rm(li)
    
    # LW from debris (DC averaged for entire cliff)
    L_d_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],L_d_m)
    L_d_dc<-do.call('cbind',L_d_dc_ls[['DC']])
    L_d_sd<-do.call('cbind',L_d_dc_ls[['SD']])
    # LW from debris (DC per cliff pixel)
    li<-apply(L_d_m,2, function(x) DC(cbind.data.frame(Date=Date[TS],x)))
    li<-stripname(li,'SD')
    li<-stripname(li,'hours')
    li<-do.call(Map,c(f=cbind,li))
    L_d_dc_pp<-do.call('cbind',li[['DC']])
    rm(li)
    
    # Sensible heat (DC averaged for entire cliff)
    H_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],H_m)
    H_dc<-do.call('cbind',H_dc_ls[['DC']])
    H_sd<-do.call('cbind',H_dc_ls[['SD']])
    
    # Latent heat (DC averaged for entire cliff)
    LE_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],LE_m)
    LE_dc<-do.call('cbind',LE_dc_ls[['DC']])
    LE_sd<-do.call('cbind',LE_dc_ls[['SD']])
    
    # Energy available for melt (DC averaged for entire cliff)
    Q_m_dc_ls<-DC_PC(cliffs_geom_df$ID,Date[TS],Q_m_m)
    Q_m_dc<-do.call('cbind',Q_m_dc_ls[['DC']])
    Q_m_sd<-do.call('cbind',Q_m_dc_ls[['SD']])
    # Energy available for melt (DC per cliff pixel)
    li<-apply(Q_m_m,2, function(x) DC(cbind.data.frame(Date=Date[TS],x)))
    li<-stripname(li,'SD')
    li<-stripname(li,'hours')
    li<-do.call(Map,c(f=cbind,li))
    Q_m_dc_pp<-do.call('cbind',li[['DC']])
    rm(li)
    
    # combine all DC (DC averaged for entire cliff)
    DC_fluxes_ls<-list(I_s_dc,D_s_dc,D_t_dc,I_n_dc,I_in_dc,L_s_dc,L_d_dc,
                       H_dc,LE_dc,Q_m_dc)
    names(DC_fluxes_ls)<-c('I_s','D_s','D_t','I_n','I_in','L_s','L_d','H','LE','Q_m')
    # combine all DC (DC per cliff pixel)
    DC_fluxes_pp_ls<-list(I_s_dc_pp,D_s_dc_pp,D_t_dc_pp,I_n_dc_pp,I_in_dc_pp,L_s_dc_pp,
                          L_d_dc_pp,Q_m_dc_pp)
    names(DC_fluxes_pp_ls)<-c('I_s','D_s','D_t','I_n','I_in','L_s','L_d','Q_m')
    
    # combine all SD of DC
    SD_fluxes_ls<-list(I_s_sd,D_s_sd,D_t_sd,I_n_sd,I_in_sd,L_s_sd,L_d_sd,
                       H_sd,LE_sd,Q_m_sd)
    names(SD_fluxes_ls)<-names(DC_fluxes_ls)
    
    # (DC averaged for entire cliff)
    filename_DC_fluxes<-path_outlist&'/DC_Fluxes_ts'&TS[1]&'to'&TS[length(TS)]&
      '_'&id&'_'&round(resol1,digits=1)&'m_ls.rda'
    save(DC_fluxes_ls,file=filename_DC_fluxes)
    # (DC per cliff pixel)
    filename_ppDC_fluxes<-path_outlist&'/perPx_DC_Fluxes_ts'&TS[1]&'to'&TS[length(TS)]&
      '_'&id&'_'&round(resol1,digits=1)&'m_ls.rda'
    save(DC_fluxes_pp_ls,file=filename_ppDC_fluxes)
    
    # (SD of DC averaged for entire cliff)
    filename_SD_fluxes<-path_outlist&'/SD_Fluxes_ts'&TS[1]&'to'&TS[length(TS)]&
      '_'&id&'_'&round(resol1,digits=1)&'m_ls.rda'
    save(SD_fluxes_ls,file=filename_SD_fluxes)
    
    #===============================================================================
    # OUTPUT: MEAN MELT RATE PER CLIFF (PC)
    #===============================================================================
    # Time series of hourly melt [m w.e./h]
    ts_Melt_pc<-aggregate(t(melt_m),list(cliffs_geom_df$ID),
                          mean)
    ts_Melt_pc<-ts_Melt_pc[,-1] #remove first column ("row" later on)
    ts_Melt_pc<-as.data.frame(t(ts_Melt_pc))
    ts_Melt_pc<-round(ts_Melt_pc,digits=5)
    colnames(ts_Melt_pc)<-paste0('Cliff_',1:ncol(ts_Melt_pc))
    ts_Melt_pc$Date<-as.character(Date[TS])
    filename_tsMelt<-path_outlist&'/PC_Melt_ts'&TS[1]&'to'&TS[length(TS)]&
      '_'&id&'_'&round(resol1,1)&'m.txt'
    write.table(ts_Melt_pc,filename_tsMelt,quote=FALSE,row.names=FALSE,
                dec='.',sep='\t')          
    
    #plot(cumsum(ts_Melt_pc[,1]))
    if(stop){break}
  }
} ### CLOSING LOOP OVER "SENSITIVITY RUN" 

gc()

