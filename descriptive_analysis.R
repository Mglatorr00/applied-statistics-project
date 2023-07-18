####################################################################
#                   THE ITALIC SEISMIC DATASET                     #
####################################################################

library(data.table)

# Importing data
ITA_18_full_dataset = readRDS("newdata.rds")
setDT(ITA_18_full_dataset)

attach(ITA_18_full_dataset)

# Clearing data
ITA_18_full_dataset$event_id = as.factor(event_id)
ITA_18_full_dataset$ev_nation_code = as.factor(ev_nation_code)
ITA_18_full_dataset$sof = as.factor(ITA_18_full_dataset$sof)
ITA_18_full_dataset$station_id = as.factor(station_id)

ITA_18_full_dataset = ITA_18_full_dataset[!is.na(ITA_18_full_dataset$magnitude)]
ITA_18_full_dataset = ITA_18_full_dataset[!is.na(ITA_18_full_dataset$ev_nation_code)]
ITA_18_full_dataset = ITA_18_full_dataset[!is.na(ITA_18_full_dataset$SA_0)]

data_ev_nation_IT <- subset(ITA_18_full_dataset, ev_nation_code == "IT")
data_rest_world <- subset(ITA_18_full_dataset, ev_nation_code != "IT")

SA = data_ev_nation_IT[, 11:47]

#-------------------------------------------------------------------
# 1. DATA VISUALIZATION AND EXPLORATION
#-------------------------------------------------------------------
library(GGally)
library(ggplot2)
library(maps)

# - Geo - visualization of the Italian data
world_map <- map_data("world")
italy <- as.data.frame(subset(world_map, world_map$region=="Italy"))

# - Plot of the Moment Magnitude on the Italian Map
x11()
ggplot() +
  geom_polygon(data=italy,aes(x=long, y=lat,group=group), color = 'black', fill = 'white') + 
  geom_point(data = data_ev_nation_IT, aes(x = event_lon, y = event_lat, color = magnitude), size = 4) +
  scale_color_gradient(low = "yellow", high = "red") +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude", color = "magnitude") + theme_bw()
ggsave("Magnitude_Italian_Map.jpeg",scale = 2)

# - Plot of type of faulting with respect to the Moment Magnitude on the Italian Map
ggplot() +
  geom_polygon(data=italy,aes(x=long, y=lat,group=group), color = 'black', fill = 'white') + 
  geom_point(data = data_ev_nation_IT, aes(x = event_lon, y = event_lat, size = magnitude, color = sof)) +
  scale_size(range = c(1, 5)) +
  scale_color_discrete() +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude", color = "Type of faulting", size="Moment Magnitude")
ggsave("Style_of_faulting_Italian_Map.jpeg",scale = 2)

# - Plot of the stations on the Italian Map
ggplot() +
  geom_polygon(data=italy,aes(x=long, y=lat,group=group), color = 'black', fill = 'white') + 
  geom_point(data = data_ev_nation_IT, aes(x = station_lon, y = station_lat)) +
  scale_color_gradient(low = "yellow", high = "red") +
  coord_equal() +
  labs(x = "Longitude", y = "Latitude", color = "vs30") + theme_bw()
ggsave("Stations_Italian_Map.jpeg",scale = 2)

# - Boxplot of Moment Magnitude by style of faulting in the world
ggplot(data = ITA_18_full_dataset, mapping = aes(x = sof, y = magnitude, fill = sof)) +
  geom_boxplot() + 
  scale_fill_manual(values = topo.colors(4)) +
  xlab("Style of faulting") + 
  ylab("Moment Magnitude") + 
  ggtitle("Boxplot of Moment Magnitude by style of faulting in the world")

# - Scatterplot 
flag_it = ITA_18_full_dataset$ev_nation_code == "IT"

ggplot(data = ITA_18_full_dataset) +
  geom_point(mapping = aes(x = distance, y = magnitude, color = flag_it )) +
  scale_fill_manual(values = topo.colors(11))

# - Plot of the ratio of SA with respect to periods
na = is.na(SA)
ratio_na = vector(mode = "numeric", length = 36)

periods = c(0, 0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3,
            0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.2, 1.4,
            1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10)

for (j in 1:37){
  ratio_na[j] = 1 - sum(na[, j])
}

plot(periods, ratio_na, pch = 16, type = "o", col = "green", xlab = "T", ylab="Ratio",
     main = "Ratio of non missing points in the Spectral Acceleration", ylim = c(0, 1), grid = TRUE)
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

# - Plot of the Spectral Acceleration as function of periods
color = hcl.colors(5522, palette="Sunset")
x11()
par(bg = "#F4F4F4")
plot(log10(periods), SA[1, ], type = 'l', lwd = 2, xlab = "log(T)", ylab="SA", 
     main="Plot of Spectral Acceleration", col = color[1], ylim = c(0, max(SA)))
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

for (j in 2:5522){
  points(log10(periods), SA[j, ], lwd = 2, type = 'l', col = color[j])
}

#-------------------------------------------------------------------
# 2. ANALYSIS OF MAGNITUDE BY STYLE OF FAULTING
#-------------------------------------------------------------------

# ANOVA testing for the mean of Moment Magnitude by style of faulting
# Test:
#   H0: The style of faulting don't have effect
#   H1: At least one style of faulting has effect

# - Boxplot of Moment Magnitude by style of faulting in Italy
x11()
ggplot(data = data_ev_nation_IT, mapping = aes(x = sof, y = magnitude, fill = sof)) +
  geom_boxplot() + 
  scale_fill_manual(values = topo.colors(4)) +
  xlab("Style of faulting") + 
  ylab("Moment Magnitude") + 
  ggtitle("Boxplot of Moment Magnitude by style of faulting in Italy")

# - Plot under different assumptions
x11()
par(mfrow=c(1,2))
barplot(rep(mean(data_ev_nation_IT$magnitude),3), names.arg=levels(data_ev_nation_IT$sof), ylim=c(0, max(data_ev_nation_IT$magnitude)),
        las=2, col='grey85', main='Model under H0')
barplot(tapply(data_ev_nation_IT$magnitude, data_ev_nation_IT$sof, mean), names.arg=levels(data_ev_nation_IT$sof), ylim=c(0, max(data_ev_nation_IT$magnitude)),
        las=2, col=topo.colors(4), main='Model under a special case of H1 with 3 populations')

# - Saving the problem dimensions
n = length(data_ev_nation_IT$sof)
ng = table(data_ev_nation_IT$sof)
treat = levels(data_ev_nation_IT$sof)
g = length(treat)     

# - Verify the assumptions for ANOVA:
#   1. Normality (univariate) in each group (6 tests) separately
Ps <- c(shapiro.test(data_ev_nation_IT$magnitude[data_ev_nation_IT$sof==treat[1]])$p,
        shapiro.test(data_ev_nation_IT$magnitude[data_ev_nation_IT$sof==treat[2]])$p,
        shapiro.test(data_ev_nation_IT$magnitude[data_ev_nation_IT$sof==treat[3]])$p)
Ps # -> Accept H0

# - Verify the covariance structure (= same sigma^2 - by using the hypothesis test)
Var <- c(var(data_ev_nation_IT$magnitude[data_ev_nation_IT$sof==treat[1]]),
         var(data_ev_nation_IT$magnitude[data_ev_nation_IT$sof==treat[2]]),
         var(data_ev_nation_IT$magnitude[data_ev_nation_IT$sof==treat[3]])) 
Var

# Test of homogeneity of variances for normal samples (Bartlett's test)
#     H0: sigma.1 = sigma.2 = sigma.3 
#     H1: there exist i,j s.t. sigma.i!=sigma.j

bartlett.test(data_ev_nation_IT$magnitude, data_ev_nation_IT$sof) # -> Accept H0

# One-way ANOVA 
fit <- aov(data_ev_nation_IT$magnitude ~ data_ev_nation_IT$sof)
summary(fit)

# We reject the test, so we have evidence to state that the style of faulting
# is responsible to differences in the magnitude (in mean) on at least one direction

# Now we can proceed in constructing confidence intervals for the means to spot
# in which directions there is a effect of the treatment

k <- g*(g-1)/2
alpha = 0.05
group_means  <- tapply(data_ev_nation_IT$magnitude, data_ev_nation_IT$sof, mean) # group-wise means
SSres <- sum(residuals(fit)^2)
S <- SSres/(n-g)
ICrange=NULL
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    print(paste(treat[i],"-",treat[j]))        
    print(as.numeric(c(group_means[i]-group_means[j] - qt(1-alpha/(2*k), n-g) * sqrt(S * ( 1/ng[i] + 1/ng[j])),
                       group_means[i]-group_means[j] + qt(1-alpha/(2*k), n-g) * sqrt(S * ( 1/ng[i] + 1/ng[j])))))
    ICrange=rbind(ICrange,as.numeric(c(group_means[i]-group_means[j] - qt(1-alpha/(2*k), n-g) * sqrt(S * (1/ng[i] + 1/ng[j])),
                                       group_means[i]-group_means[j] + qt(1-alpha/(2*k), n-g) * sqrt(S * (1/ng[i] + 1/ng[j])))))
  }}

x11()
plot(data_ev_nation_IT$sof, data_ev_nation_IT$magnitude, xlab='treat', ylab='magnitude', col = topo.colors(4), las=2)
h <- 1
x11()
plot(c(1,g*(g-1)/2),range(ICrange), pch='',xlab='pairs treat', ylab='Conf. Int. tau magnitude')
for(i in 1:(g-1)) {
  for(j in (i+1):g) {
    ind <- (i-1)*g-i*(i-1)/2+(j-i)
    lines (c(h,h), c(ICrange[ind,1],ICrange[ind,2]), col='grey55'); 
    points(h, group_means[i]-group_means[j], pch=16, col='grey55'); 
    points(h, ICrange[ind,1], col=topo.colors(4)[j], pch=16); 
    points(h, ICrange[ind,2], col=topo.colors(4)[i], pch=16); 
    h <- h+1
  }}
abline(h=0)

# We can conclude that the style of faulting has an effect on the mean magnitude

rm(n)
rm(ng)
rm(treat)
rm(g)
rm(Ps)
rm(Var)
rm(fit)

#-------------------------------------------------------------------
# 3. TEST FOR REPEATED MEASURES (full dataset)
#-------------------------------------------------------------------
# numero di eventi
n <- dim(ITA_18_full_dataset)[1]

# periodi
q <- dim(ITA_18_full_dataset[, 11:47])[2]

# sample mean
M <- sapply(ITA_18_full_dataset[, 11:47], mean) 

# covariance matrix
S <- cov(ITA_18_full_dataset[, 11:47]) 

#contrast matrix (q-1)xq
C<- rep(0,37)
C[1]=-1
C[2]=1
for (i in 2:36)
{
  row<- rep(0,37)
  row[i]=-1
  row[i+1]=1
  C=rbind(C, row)
}

# Test: H0: C %% mu == rep(0,36)  vs H1: C %% mu != rep(0,36)
alpha   <- .05
delta.0 <- rep(0,36)

Md <- C %*% M 
Sd <- C %% S %% t(C) 
Sdinv <- solve(Sd)

# Hotelling T2 statistics
T2 <- n * t(Md - delta.0) %% Sdinv %% (Md - delta.0)

# (q-1)*(n-1)/(n-(q-1)) times the 1-alpha Fisher quantile with q-1 and n-q+1 df 
cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha, (q-1), n-(q-1)) 

# Testing if we are in the rejection region
T2 < cfr.fisher # Testing if we are in the rejection region
cfr.fisher

# T2 is much higher than cfr.fisher => the p-value will be very small
P <- 1 - pf(T2*(n-(q-1))/((q-1)*(n-1)), (q-1), n-(q-1))

#CONFIDENCE INTERVALS
# Simultaneous T2 intervals in the direction of the contrasts h.* - h.0
IC.T2 <- cbind(Md - sqrt(cfr.fisher*diag(Sd)/n),
               Md,
               Md + sqrt(cfr.fisher*diag(Sd)/n))


# Bonferroni intervals 
k <- q-1   # number of increments 
cfr.t <- qt(1-alpha/(2*k), n-1)

IC.BF <- cbind(Md - cfr.t*sqrt(diag(Sd)/n),
               Md,
               Md + cfr.t*sqrt(diag(Sd)/n))

x11()
matplot(t(matrix(1:36, 36, 3)), t(IC.BF), type='b',ylim=c(-10,20),  pch='', xlab='',
        ylab='', main='Confidence intervals')

# Plotting the Bonferroni intervals
segments(matrix(1:36, 36, 1), IC.BF[,1], matrix(1:36, 36, 1), IC.BF[,3],
         col='green', lwd=2)
points(1:36, IC.BF[,2], col='green', pch=16)

# Plotting delta.0 under H0 (delta.0 == c(0, 0, 0))
points(1:36+.05, delta.0, col='black', pch=16)

# Plotting the simultaneous T2
segments(matrix(1:36+.3,36,1),IC.T2[,1],matrix(1:36+.3,36,1),IC.T2[,3], col='blue', lwd=2)
points(1:36+.3,IC.T2[,2], col='blue', pch=16)

legend('topright', c('Bonf. IC', 'Sim-T2 IC'), col=c('green', 'blue'), lty=1, lwd=2)


#changing the contrast matrix...
C1<- rep(1,36)
for (i in 1:36)
{
  col<- rep(0,36)
  col[i]=-1
  C1=cbind(C1, col)
}

Md <- C1 %*% M 
Sd <- C1 %% S %% t(C1) 
Sdinv <- solve(Sd)

# Hotelling T2 statistics
T2 <- n * t(Md - delta.0) %% Sdinv %% (Md - delta.0)

# (q-1)*(n-1)/(n-(q-1)) times the 1-alpha Fisher quantile with q-1 and n-q+1 df 
cfr.fisher <- ((q-1)*(n-1)/(n-(q-1)))*qf(1-alpha, (q-1), n-(q-1)) 

# Testing if we are in the rejection region
T2 < cfr.fisher # Testing if we are in the rejection region
cfr.fisher

# T2 is much higher than cfr.fisher => the p-value will be very small
P <- 1 - pf(T2*(n-(q-1))/((q-1)*(n-1)), (q-1), n-(q-1))

#CONFIDENCE INTERVALS
# Simultaneous T2 intervals in the direction of the contrasts h.* - h.0
IC.T2 <- cbind(Md - sqrt(cfr.fisher*diag(Sd)/n),
               Md,
               Md + sqrt(cfr.fisher*diag(Sd)/n))


# Bonferroni intervals 
k <- q-1   # number of increments 
cfr.t <- qt(1-alpha/(2*k), n-1)

IC.BF <- cbind(Md - cfr.t*sqrt(diag(Sd)/n),
               Md,
               Md + cfr.t*sqrt(diag(Sd)/n))


x11()
matplot(t(matrix(1:36, 36, 3)), t(IC.BF), type='b', pch='', ylim=c(-60,60), xlab='',
        ylab='', main='Confidence intervals')

# Plotting the Bonferroni intervals
segments(matrix(1:36, 36, 1), IC.BF[,1], matrix(1:36, 36, 1), IC.BF[,3],
         col='green', lwd=2)
points(1:36, IC.BF[,2], col='green', pch=16)

# Plotting delta.0 under H0 (delta.0 == c(0, 0, 0))
points(1:36+.05, delta.0, col='black', pch=16)

# Plotting the simultaneous T2
segments(matrix(1:36+.3,36,1),IC.T2[,1],matrix(1:36+.3,36,1),IC.T2[,3], col='blue', lwd=2)
points(1:36+.3,IC.T2[,2], col='blue', pch=16)

legend('topright', c('Bonf. IC', 'Sim-T2 IC'), col=c('green', 'blue'), lty=1, lwd=2)
