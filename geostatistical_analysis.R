setwd("~/Desktop/1^LM ing mat/Applied Statistic/Applied Statistic Project 2")

library(data.table)

# Importing data
ITA_18_full_dataset = readRDS("newdata.rds")
setDT(ITA_18_full_dataset)

attach(ITA_18_full_dataset)

# Clearing data
ITA_18_full_dataset$event_id = as.factor(event_id)
ITA_18_full_dataset$ev_nation_code = as.factor(ev_nation_code)
ITA_18_full_dataset$sof = as.factor(sof)
ITA_18_full_dataset$station_id = as.factor(station_id)
ITA_18_full_dataset = ITA_18_full_dataset[!is.na(ITA_18_full_dataset$magnitude)]
ITA_18_full_dataset = ITA_18_full_dataset[!is.na(ITA_18_full_dataset$ev_nation_code)]
ITA_18_full_dataset = ITA_18_full_dataset[!is.na(ITA_18_full_dataset$SA_0)]
data_ev_nation_IT <- subset(ITA_18_full_dataset, ev_nation_code == "IT")

### mapping of italy 

library(maps)
world_map <- map_data("world")
italy <- as.data.frame(subset(world_map, world_map$region=="Italy"))

coordinates(italy) <- c('long','lat')

Sardinia <- as.data.frame(subset(italy, italy$subregion=="Sardinia"))
Sicily <- as.data.frame(subset(italy, italy$subregion=="Sicily"))

italy<- italy[-which(italy$subregion =='Isola di Pantelleria'),]
italy<- italy[-which(italy$subregion =='Asinara'),]
italy<- italy[-which(italy$subregion =='Forio'),]
italy<- italy[-which(italy$subregion =="Isola d'Elba"),]
italy<- italy[-which(italy$subregion =="Sant'Antonio"),]
italy<- italy[-which(italy$subregion =='Sardinia'),]
italy<- italy[-which(italy$subregion =='Sicily'),]

coordinates(Sardinia) <- c('long','lat')
coordinates(Sicily) <- c('long','lat')

italy.lst <- list(Polygons(list(Polygon(italy[,1:2])), "Italy"))
italy.sr <- SpatialPolygons(italy.lst)
sard.lst <- list(Polygons(list(Polygon(Sardinia[,1:2])), "Sardinia"))
sard.sr <- SpatialPolygons(sard.lst)
sic.lst <- list(Polygons(list(Polygon(Sicily[,1:2])), "Sicily"))
sic.sr <- SpatialPolygons(sic.lst)

full.ita.sr<- rbind(italy.sr, sic.sr, sard.sr)

# making the grid for prediction
grd<- makegrid(full.ita.sr, n=10000)
coordinates(grd)<-c('x1','x2')
grd_in<-grd[full.ita.sr,]
# but we have to remember that each location has a specific soil
# estimating vs30 for all the grid points
data_ev_nation_IT <- subset(ITA_18_full_dataset, ev_nation_code == "IT")
x.train<- data.frame(scale(data_ev_nation_IT[,c(5,6)]))
y.train<- data.frame(data_ev_nation_IT$vs30)
test<- data.frame(scale(coordinates(grd_in)))
vs30.knn<-FNN::knn.reg(train=x.train[,1:2],
                       test=test,
                       y=y.train,
                       k=15
                       )
grid<- data.frame(grd_in,
                  vs30=vs30.knn$pred)
coordinates(grid)<-c('x1','x2')
grid<- as(grid, 'SpatialPixelsDataFrame')

coordinates(data_ev_nation_IT)<-c('station_lon', 'station_lat')
plot(grid, main='reconstruction of vs30 via knn (15)')
########################
### geostat analysis ###
########################

library(sp)        
library(lattice)   
library(gstat)  

# We noticed that using the original data we have issues in kriging: this is due
# to the fact that we have obviously the same coordinates for the same stations 
# but different values for our variables. So we have to jitter the coordinates
# to avoid colinearity-like problems

new_lat <- jitter(data_ev_nation_IT$station_lat)
new_lon <- jitter(data_ev_nation_IT$station_lon)
geostat_ita_data<- data.frame(new_lat, new_lon, data_ev_nation_IT)
coordinates(geostat_ita_data)<-c( 'new_lon','new_lat')

stations<-data.frame(coordinates(data_ev_nation_IT))

## variogram estimation

v <- variogram( log(SA_0.01) ~1 , 
                data = geostat_ita_data)
v.fit<-fit.variogram(v, vgm(2, "Exp", 2, 3) )
plot(v, v.fit, pch = 19)

ita.gstat <- gstat(formula = log(SA_0.01) ~ magnitude+vs30 +
                     distance,
                   data = geostat_ita_data, nmax = 50, model=v.fit, 
                   set = list(gls=1))
v.gls<-variogram(ita.gstat)
v.gls.fit <- fit.variogram(v.gls, vgm(1, "Exp", 1, 0.7))
par(bg = "#F4F4F4", mfrow = c(3, 3))
plot(v.gls, v.gls.fit, pch = 19, col='#E15269')

ita.gstat <- gstat(formula = log(SA_0.01) ~ magnitude + vs30 + distance,
                   data = geostat_ita_data, nmax = 50, 
                   model=v.gls.fit, set = list(gls=1))

# we transform everything in points of class sf.
# in this way we could compute the distance in km starting from the spatial
# coordinates in WGS84

library(sf) 

# 1) computing geographic distance between epicenter and stations

# creation of the event as a spatial point 
event<- data.frame(x=14.33,y=41.74)
ev<-st_as_sf(event, coords = c('x','y'),crs = 4326)

# creation of the stations as spatial points
stat<- st_as_sf(stations, coords = c('station_lon','station_lat'),crs = 4326)

# the distance expressed in meters [m]
dist_from_event<-st_distance(ev, stat )

# 2) prediction for all the stations
# s0<- data.frame(coordinates(geostat_ita_data),
#                 magnitude=4,
#                 vs30=geostat_ita_data$vs30,
#                 distance=as.numeric(dist_from_event/1000))
# coordinates(s0)<- c('new_lon', 'new_lat')
# 
# lz.uk <- predict(ita.gstat, s0, BLUE = F)
# spplot(lz.uk)

#3) prediction over the grid 
grd_new<-st_as_sf(grid , coords = c('x1','x2'))
grd_new<-st_set_crs(grd_new, 4326)
grid_dist_from_ev<- st_distance(ev , grd_new)

z0<-data.frame(coordinates(grid),
               vs30=grid$vs30, 
               magnitude=6,
               distance=as.numeric(grid_dist_from_ev)/1000)
coordinates(z0)<-c('x1','x2')
z0<- as(z0, 'SpatialPixelsDataFrame')
lz.uk <- predict(ita.gstat, z0, BLUE = F)
coordinates(event)<-c('x','y')
spplot(lz.uk[,1], 
       sp.layout = list("sp.points", event, pch = 16, cex=1, col = "green"),
       main=paste("Earthquake in lon",coordinates(event)[,1], "and lat",coordinates(event)[,2], sep=" "))

##############################
## ground motion model ITA18 #
##############################

library(data.table)
library(car)
library(glmnet)
library(caret)
library(MASS)
###########################
### ITA 18 linear model ###
###########################

#source function
Mh = 5.7
source <- ITA_18_full_dataset$magnitude - Mh
selected_obs <- which( ITA_18_full_dataset$magnitude <= Mh )
dummy1 <- rep(0, dim( ITA_18_full_dataset )[1])
for (i in c(1:dim( ITA_18_full_dataset )[1])){
  for (j in c(1:length(selected_obs))){
    if (i == selected_obs[j]){
      dummy1[i] <- 1
    }
  }
}
selected_obs <- which( ITA_18_full_dataset$magnitude > Mh )
dummy2 <- rep(0, dim( ITA_18_full_dataset )[1])
for (i in c(1:dim( ITA_18_full_dataset )[1])){
  for (j in c(1:length(selected_obs))){
    if (i == selected_obs[j]){
      dummy2[i] <- 1
    }
  }
}
#path function
Mref = 4.5
path1 <- (ITA_18_full_dataset$magnitude - Mref) * log(ITA_18_full_dataset$distance, base=10)
path2 <- log(ITA_18_full_dataset$distance, base = 10) 
#site function
site <- rep(0,dim(ITA_18_full_dataset)[1])
site[which(ITA_18_full_dataset$vs30 < 1500)] <- log(ITA_18_full_dataset[which(ITA_18_full_dataset$vs30 < 1500)]$vs30/800)
site[which(ITA_18_full_dataset$vs30 >= 1500)] <- log(1500/800)
# - Final Design Matrix (remove the units with distance = 0)
Z <- cbind(source, dummy1, dummy2, path1, path2, site, ITA_18_full_dataset)
remove <- which(ITA_18_full_dataset$distance == 0)
Z <- Z[-remove,]

model <- lm(log(SA_0.01) ~ source:dummy1 + source:dummy2 + sof + distance + path1 + path2 + site, data = Z)
res<-residuals(model)  
pred<-model$fitted.values

############################
### ITA 18 geostat model ###
############################

# Since our believe is that there exist a spatial dependence which is not captured
# by the linear model we tried to apply geostatistics to ITA18.
# We thought that the spatial dependence was hidden in the residuals of the model
# and so we proceeded in estimating those latter via ordinary kriging and finally
# computing the sum between those residuals and the (mean) prediction given by
# the linear model over the whole italian grid

# fit variogram on the residuals
X<- data.frame(res,lat=jitter(Z$station_lat),lon=jitter(Z$station_lon))
coordinates(X)<-c('lon','lat')

v<-autofitVariogram(res~1 , input_data=X, model = c("Sph", "Exp", "Gau", "Ste","Per","Hol"))

par(bg = 'grey')
plot(v$exp_var, v$var_model,pch=19,col='#E15269')

plot(v)
ita.gstat <- gstat(formula=res~1, 
                   data = X, nmax = 50, model=v$var_model)
o.k<-predict(ita.gstat, grid)
spplot(o.k, main='prediction of residuals',
       sp.layout = list("sp.points", event, pch = 16, cex=1, col = "green"))

# computing the prediction over the grid

#source function
Mh = 5.7
source <- z0$magnitude - Mh
selected_obs <- which( z0$magnitude <= Mh )
dummy1 <- rep(0, dim( z0 )[1])
for (i in c(1:dim( z0 )[1])){
  for (j in c(1:length(selected_obs))){
    if (i == selected_obs[j]){
      dummy1[i] <- 1
    }
  }
}
selected_obs <- which( z0$magnitude > Mh )
dummy2 <- rep(0, dim( z0 )[1])
for (i in c(1:dim( z0 )[1])){
  for (j in c(1:length(selected_obs))){
    if (i == selected_obs[j]){
      dummy2[i] <- 1
    }
  }
}
#path function
Mref = 4.5
path1 <- (z0$magnitude - Mref) * log(z0$distance, base=10)
path2 <- log(z0$distance, base = 10) 
#site function
site <- rep(0,dim(z0)[1])
site[which(z0$vs30 < 1500)] <- log(ITA_18_full_dataset[which(z0$vs30 < 1500)]$vs30/800)
site[which(z0$vs30 >= 1500)] <- log(1500/800)
# - Final Design Matrix (remove the units with distance = 0)
Z_z0 <- data.frame(source, dummy1, dummy2, path1, path2, site, 
                   coordinates(grid),
                   vs30=grid$vs30, 
                   magnitude=6,
                   sof='NF',
                   distance=as.numeric(grid_dist_from_ev)/1000)
pred<-predict(model,Z_z0)

# summing the prediction of the linear model over the grid to the prediction 
# of the residuals given by the geaostat model

total<- o.k$var1.pred+pred

grid_total<-data.frame(coordinates(grid),total)
coordinates(grid_total)<-c('x1','x2')
grid_total<- as(grid_total, 'SpatialPixelsDataFrame')
par(mfrow=c(1,1), bg='white')
plot(grid_total)
points(event, col='green', pch=19)
title(main=paste("Earthquake in lon",coordinates(event)[,1], "and lat",coordinates(event)[,2], sep=" "))
