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

#-------------------------------------------------------------------
# 1. IMPLEMENTATION OF ITA18 GROUND MOTION MODEL
#-------------------------------------------------------------------

library(data.table)
library(car)
library(glmnet)
library(caret)
library(MASS)

numeric_values <- cbind(ITA_18_full_dataset[, c(7:10)])
SA_log = log(ITA_18_full_dataset[, 11:47])
initial_dataset <- data.frame(cbind(numeric_values, SA_log ) )

# - Construction of the design matrix:
#source function
Mh = 5.7
source <- initial_dataset$magnitude - Mh

selected_obs <- which( initial_dataset$magnitude <= Mh )
dummy1 <- rep(0, dim( initial_dataset )[1])
for (i in c(1:dim( initial_dataset )[1])){
  for (j in c(1:length(selected_obs))){
    if (i == selected_obs[j]){
      dummy1[i] <- 1
    }
  }
}

selected_obs <- which( initial_dataset$magnitude > Mh )
dummy2 <- rep(0, dim( initial_dataset )[1])
for (i in c(1:dim( initial_dataset )[1])){
  for (j in c(1:length(selected_obs))){
    if (i == selected_obs[j]){
      dummy2[i] <- 1
    }
  }
}

#path function
Mref = 4.5
#c1
path1 <- (initial_dataset$magnitude - Mref) * log(initial_dataset$distance, base=10)
#c2
path2 <- log(initial_dataset$distance, base = 10) 

#site function
site <- rep(0,dim(initial_dataset)[1])
site[which(initial_dataset$vs30 < 1500)] <- log(initial_dataset[which(initial_dataset$vs30 < 1500),4]/800)
site[which(initial_dataset$vs30 >= 1500)] <- log(1500/800)

# - Final Design Matrix (remove the units with distance = 0)
Z <- data.frame(source, dummy1, dummy2 , initial_dataset$sof, initial_dataset$distance, path1, path2, site, 
                initial_dataset[,5:41])
remove <- which(initial_dataset$distance == 0)
Z <- Z[-remove,]

# - Linear model of ITA18 GROUND MOTION MODEL
min_tot=c()
max_tot=c()
mean_tot=c()
res_std_err_tot=c()
coll <- c()
sd.ITA18.GMM = c()
for(i in c(9:45)){
  # one linear model for each period
  model <- lm(Z[,i] ~ source:dummy1 + source:dummy2 + initial_dataset.sof + initial_dataset.distance + 
                path1 + path2 + site, data = Z[,c(1,2,3,4,5,6,7,8,i)])
  
  #plot(model) #this is for the analysis of the residuals: they are all fine (normality of residuals is
  #verified, as well as absence of influential cases)
  sd.ITA18.GMM = c(sd.ITA18.GMM, sqrt(var(model$residuals)))
  cl <- vif(model)[,1]
  coll <- cbind(coll, cl)
  
  confint(model)
  min=confint(model)[,1]
  min_tot=rbind(min_tot,t(min))
  max=confint(model)[,2]
  max_tot=rbind(max_tot,t(max))
  mean=(max+min)/2
  mean_tot=rbind(mean_tot,t(mean))
  
  res_std_err=sqrt(deviance(model)/model$df.residual)
  res_std_err_tot=rbind(res_std_err_tot, t(res_std_err))
}

colnames(coll) <- names(ITA_18_full_dataset[, c(11:47)])
periods = log10(c(0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3,
                  0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.2, 1.4,
                  1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10))

names <- c("a", "f1", "f2", "c3", "c1", "c2", "k", "b1", "b2")
x11()
par(mfrow=c(3,3), bg = "white")
for(i in c(1:9)){
  plot(periods, mean_tot[2:37,i], ylim =c(min(min_tot[,i]) , max(max_tot[,i]) ),col = "#1B1B19", pch=19, ylab = '',xlab = '' ,
       main=paste(names[i]))
  points(periods, max_tot[2:37,i], col = "#1B1B19", pch = '_')
  points(periods, min_tot[2:37,i] , col = "#1B1B19", pch='_')
  for(j in c(2:37)){
    segments(x0=periods[j-1],y0=min_tot[j,i],x1=periods[j-1],y1=max_tot[j,i])
  }
  abline(h = 0)
}
# We can see very high values of VIF -> collinearity of the parameters
# so we try with a penalized regression -> RIDGE REGRESSION 

# - Ridge regressionof ITA18 GROUND MOTION MODEL
bestlam.ridge <- c()
coef.ridge <- c()

for(i in c(9:45)){
  xx <- model.matrix(Z[,i] ~ source:dummy1 + source:dummy2 + initial_dataset.sof + initial_dataset.distance + 
                       path1 + path2 + site, data = Z[,c(1,2,3,4,5,6,7,8)])
  yy <- Z[,i]
  lambda.grid <- 10^seq(5,-5,length=100)
  cv.ridge <- cv.glmnet(xx , yy , alpha = 0, lambda=lambda.grid) # alpha = 0 => ridge
  bestlam.ridge <- cbind(bestlam.ridge, cv.ridge$lambda.min)
  fit.ridge <- glmnet( xx , yy , alpha = 0, lambda = cv.ridge$lambda.min) # alpha = 0 => ridge
  # fit.ridge <- lm.ridge( Z[,i] ~ source:dummy1 + source:dummy2 + initial_dataset.sof + initial_dataset.distance + 
  # path1 + path2 + site, data = Z[,c(1,2,3,4,5,6,7,8)] ,lambda = cv.ridge$lambda.min) 
  coef.ridge <- cbind(coef.ridge, predict(fit.ridge, s=cv.ridge$lambda.min, type = 'coefficients')[1:10,])
}

coef.ridge <- coef.ridge[-2,]

names <- c("a", "f1", "f2", "c3", "c1", "c2", "k", "b1", "b2")
x11()
par(mfrow=c(3,3))
for(i in c(1:9)){
  plot(periods, coef.ridge[i,], ylim =c(min(coef.ridge[i,]) , max(coef.ridge[i,]) ),col = 'blue', ylab = '',xlab = '' ,
       main=paste(names[i],"ridge"), pch = 19)
  abline(h = 0)
}

#-------------------------------------------------------------------
# 2. CROSS VALIDATION BETWEEN ALL MODELS TO FINE THE BEST ONE !
#-------------------------------------------------------------------
# - Cross-validation on the new model with sof as variable  
library(caret)
numeric_values <- cbind(ITA_18_full_dataset[, c(7,8,10)], dummy_tf, dummy_ss)
SA_log = log(ITA_18_full_dataset[, 11:47])
data_cv <- cbind(numeric_values, SA_log )
attach(data_cv)


ctrl <- trainControl(method = "cv", number = 25)

SA_0_sof <- train(SA_0 ~ magnitude + distance + vs30 + dummy_tf + dummy_ss  +
                    dummy_tf:magnitude + dummy_tf:distance + dummy_tf:vs30 +  
                    dummy_ss:magnitude + dummy_ss:distance + dummy_ss:vs30, data = data_cv[, c(1,2,3,4,5,6)], method = "lm", trControl = ctrl) 
print(SA_0_sof)
detach(data_cv)

RMSE_sof <- as.data.frame(SA_0_sof$results$RMSE)

copie <- data_cv[,7:42]
colnames(copie) <- rep("SA", 36)

for (j in c(7:42)){
  data_j <- data_cv[ , 1:5]
  data_regression_sof <- cbind(data_j, copie[,1])
  attach(data_regression_sof)
  regression_sof <- train( SA ~ magnitude + distance + vs30 + dummy_tf + dummy_ss  +
                             dummy_tf:magnitude + dummy_tf:distance + dummy_tf:vs30 +  
                             dummy_ss:magnitude + dummy_ss:distance + dummy_ss:vs30, data = data_regression_sof, method = "lm", trControl = ctrl) 
  
  RMSE_sof <- cbind(RMSE_sof, regression_sof$results$RMSE)
  copie[,1] <- NULL
  detach(data_regression_sof)
}

colnames(RMSE_sof) <- names(ITA_18_full_dataset[, c(11:47)])
rm(data_j)
rm(data_regression_sof)
rm(regression_sof)
rm(copie)

# - Cross-validation on ITA18 GMM 
data_cv <- Z
attach(data_cv)

ita.control = trainControl(method = "cv", number = 25)

SA_0_ <- train(SA_0 ~ source:dummy1 + source:dummy2 + initial_dataset.sof + initial_dataset.distance +
                 path1 + path2 + site, data = data_cv[,c(1,2,3,4,5,6,7,8,9)], method = "lm", trControl = ita.control)

print(SA_0_)
detach(data_cv)

RMSE_ita <- as.data.frame(SA_0_$results$RMSE)
RMSE_ita


copie <- data_cv[,9:45]
colnames(copie) <- rep("SA", 36)


for (j in c(7:42)){
  data_j <- data_cv[ , 1:5]
  data_regression <- cbind(data_j, copie[,1])
  attach(data_regression)
  
  regression <- train(SA_0 ~ source:dummy1 + source:dummy2 + initial_dataset.sof + initial_dataset.distance + 
                        path1 + path2 + site, data = data_cv[,c(1,2,3,4,5,6,7,8,9)], method = "lm", trControl = ita.control)
  RMSE_ita <- cbind(RMSE_ita, regression$results$RMSE)
  copie[,1] <- NULL
  
  detach(data_regression)
}

colnames(RMSE_ita) <- names(data_cv[, c(9:45)])
rm(data_j)
rm(data_regression)
rm(regression)
rm(copie)



# - Cross-validation on Ridge-version of ITA18 GMM 
data_cv <- Z
attach(data_cv)

ridge.control = trainControl(method = "cv", number = 25)
ridge.grid = expand.grid(alpha = 0, lambda = cv.ridge$lambda.min) # alpha = 0 => ridge

SA_0_ <- train(SA_0 ~ source:dummy1 + source:dummy2 + initial_dataset.sof + initial_dataset.distance +
                 path1 + path2 + site, data = data_cv[,c(1,2,3,4,5,6,7,8,9)], method = "glmnet", trControl = ridge.control, tuneGrid = ridge.grid)

print(SA_0_)
detach(data_cv)

RMSE_ridge <- as.data.frame(SA_0_$results$RMSE)
RMSE_ridge


copie <- data_cv[,9:45]
colnames(copie) <- rep("SA", 36)


for (j in c(7:42)){
  data_j <- data_cv[ , 1:5]
  data_regression <- cbind(data_j, copie[,1])
  attach(data_regression)
  
  regression <- train(SA_0 ~ source:dummy1 + source:dummy2 + initial_dataset.sof + initial_dataset.distance + 
                        path1 + path2 + site, data = data_cv[,c(1,2,3,4,5,6,7,8,9)], method = "glmnet", trControl = ridge.control, tuneGrid = ridge.grid)
  RMSE_ridge <- cbind(RMSE_ridge, regression$results$RMSE)
  copie[,1] <- NULL
  
  detach(data_regression)
}

colnames(RMSE_ridge) <- names(data_cv[, c(9:45)])
rm(data_j)
rm(data_regression)
rm(regression)
rm(copie)


# - Plot comparing the RMSE of the models
periods = log10(c(0, 0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3,
                  0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.2, 1.4,
                  1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10))

x11()
par(bg = "white")
plot(periods[-1],RMSE_sof[,2:37], xlab = "periods", ylab = "RMSE", pch = 19, col = '#E15269', ylim = c(0.82,1.053), type = 'b')
points(periods[-1], RMSE_ridge[,2:37], col = "#66121F", pch = 19, type = 'b')
points(periods[-1], RMSE_ita[,2:37], col = "#1B1B19", pch = 19, type = 'b')
legend("topright", legend = c("RMSE with SOF","RMSE with ridge","RMSE of ITA18" ),
       col = c("#E15269", "#66121F", "#1B1B19"), pch = 19)

RMSE_sof <- as.numeric(RMSE_sof)
RMSE_ridge <- as.numeric(RMSE_ridge)
RMSE_ita <- as.numeric(RMSE_ita)
mean(RMSE_sof)
mean(RMSE_ridge)
mean(RMSE_ita)

#-------------------------------------------------------------------
# 3. FUNCTIONAL DATA ANALYSIS AND FUNCTIONAL REGRESSION
#-------------------------------------------------------------------
library(fda)
library(KernSmooth)
library(rgl)
library(fields)

periods = log10(c(0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3,
                  0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.2, 1.4,
                  1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10))
# - Select the best smoothing basis by looking at one seismic event
Xobs0 <- t(SA_log[1,])

NT <- length(periods) # number of locations of observations

basis <- create.bspline.basis(periods, norder=4)
lambda <- 10^seq(-15,0,by = 0.5)
gcv <- numeric(length(lambda))
for (i in 1:length(lambda)){
  functionalPar <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda[i])  
  gcv[i] <- smooth.basis(periods, Xobs0, functionalPar)$gcv
}
par(mfrow=c(1,1))
plot(log10(lambda),gcv)
opt.lambda <- lambda[which.min(gcv)]    
opt.lambda


basis <- create.bspline.basis(periods, norder=4)
functionalPar <- fdPar(fdobj=basis, Lfdobj=2, lambda= opt.lambda) 
Xss <- smooth.basis(periods, Xobs0, functionalPar)
Xss0 <- eval.fd(periods, Xss$fd, Lfd=0)
Xss1 <- eval.fd(periods, Xss$fd, Lfd=1)  #1st derivative
Xss2 <- eval.fd(periods, Xss$fd, Lfd=2)  #2nd derivative
rappincX1 <- (Xobs0[3:NT]-Xobs0[1:(NT-2)])/(periods[3:NT]-periods[1:(NT-2)])
rappincX2 <- ((Xobs0[3:NT]-Xobs0[2:(NT-1)])/(periods[3:NT]-periods[2:(NT-1)])-(Xobs0[2:(NT-1)]-Xobs0[1:(NT-2)])/(periods[2:(NT-1)]-periods[1:(NT-2)]))*2/(periods[3:(NT)]-periods[1:(NT-2)])

x11()
par(bg ="#F4F4F4", mfrow=c(1,3))
plot(periods,Xobs0,xlab="log(T)",main="observed data", ylab=" ", pch = 19, col = "#1B1B19", type="b")
points(periods,Xss0 ,type="l",col="#66121F",lwd=2)
plot(periods[2:(NT-1)],rappincX1,xlab="log(T)",main = "first differences x", type="b", ylab=" ",pch = 19, col = "#1B1B19")
points(periods,Xss1 ,type="l",col="#66121F",lwd=2)
plot(periods[2:(NT-1)],rappincX2,xlab="log(T)",main = "second differences x",type="b", ylab=" ",pch = 19, col = "#1B1B19")
points(periods,Xss2 ,type="l",col="#66121F",lwd=2)


# - Penalized smoothing of the data
# We are using cubic Bspline basis with coefficient of penalization 1e-05 on the second derivative

Y <- t(SA_log)
data.fd <- Data2fd(y = Y , argvals = periods, basisobj = basis, lambda = opt.lambda ) #create fda object

color = hcl.colors(6000, palette="Sunset")
x11()
par(bg = "#F4F4F4")
plot.fd(data.fd, main = "Functional Data after Smoothing", col = color)


# - Functional Regression
n <- dim(SA_log)[1]

xlist <- list("intercept" = rep(1,n),
              "source1"   = source*dummy1,
              "source2"   = source*dummy2,
              "sof"       = as.numeric(sof),
              "distance"  = R,
              "path1"     = path1,
              "path2"     = path2,
              "site"      = site)

basis  <- data.fd$basis
lambda <- opt.lambda # provate a cambiare e vedete che succede
fPar   <- fdPar(fdobj=basis, Lfdobj=2, lambda=lambda)
blist  <- list(fPar,fPar,fPar,fPar,fPar,fPar,fPar,fPar) # tante volte quanti sono le covariate + l'intercetta

mod <- fRegress(y        = data.fd,
                xfdlist  = xlist,
                betalist = blist)

beta0.hat.fd <- mod$betaestlist[[1]]
beta1.hat.fd <- mod$betaestlist[[2]]
beta2.hat.fd <- mod$betaestlist[[3]]
beta3.hat.fd <- mod$betaestlist[[4]]
beta4.hat.fd <- mod$betaestlist[[5]]
beta5.hat.fd <- mod$betaestlist[[6]]
beta6.hat.fd <- mod$betaestlist[[7]]
beta7.hat.fd <- mod$betaestlist[[8]]

y.hat.fd <- mod$yhatfdobj

x11()
plot(y.hat.fd)
title(main="Estimated functional response")

tt <- seq(-2,1, length.out=100)
beta0.vals <- eval.fd(tt, beta0.hat.fd$fd)
beta1.vals <- eval.fd(tt, beta1.hat.fd$fd)
beta2.vals <- eval.fd(tt, beta2.hat.fd$fd)
beta3.vals <- eval.fd(tt, beta3.hat.fd$fd)
beta4.vals <- eval.fd(tt, beta4.hat.fd$fd)
beta5.vals <- eval.fd(tt, beta5.hat.fd$fd)
beta6.vals <- eval.fd(tt, beta6.hat.fd$fd)
beta7.vals <- eval.fd(tt, beta7.hat.fd$fd)

x11()
par(mfrow=c(2,4))
plot(tt, beta0.vals, type='l', lwd=3,main="Estimated beta0")
plot(tt, beta1.vals, type='l', lwd=3,main="Estimated beta1")
plot(tt, beta2.vals, type='l', lwd=3,main="Estimated beta2")
plot(tt, beta3.vals, type='l', lwd=3,main="Estimated beta3")
plot(tt, beta4.vals, type='l', lwd=3,main="Estimated beta4")
plot(tt, beta5.vals, type='l', lwd=3,main="Estimated beta5")
plot(tt, beta6.vals, type='l', lwd=3,main="Estimated beta6")
plot(tt, beta7.vals, type='l', lwd=3,main="Estimated beta7")

# - Compare the scalar coefficients of ITA18 with the functional ones
# confronto dei beta_hat con ITA18 scalare
periods = log10(c(0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3,
                  0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.2, 1.4,
                  1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10))

x11()
par(mfrow=c(2,4))
plot(periods, mean_tot[2:37,1], ylim =c(min(min_tot[,1]) , max(max_tot[,1]) ),col = 'red', ylab = '',xlab = '' ,
     main=paste(names[1]), pch = 19)
abline(h = 0)
plot(periods, mean_tot[2:37,8], ylim =c(min(min_tot[,8]) , max(max_tot[,8]) ),col = 'black', ylab = '',xlab = '' ,
     main=paste(names[8]), pch = 19)
abline(h = 0)
plot(periods, mean_tot[2:37,9], ylim =c(min(min_tot[,9]) , max(max_tot[,9]) ),col = 'red', ylab = '',xlab = '' ,
     main=paste(names[9]), pch = 19)
abline(h = 0)
plot(periods, mean_tot[2:37,3], ylim =c(min(min_tot[,3]) , max(max_tot[,3]) ),col = 'red', ylab = '',xlab = '' ,
     main=paste(names[3]), pch = 19)
abline(h = 0)
for(i in c(4:7)){
  plot(periods, mean_tot[2:37,i], ylim =c(min(min_tot[,i]) , max(max_tot[,i]) ),col = 'red', ylab = '',xlab = '' ,
       main=paste(names[i]), pch = 19)
  abline(h = 0)
}

x11()
plot(tt, beta0.vals, type='l', lwd=3,main="Estimated beta0")
points(periods, mean_tot[2:37,1], ylim =c(min(min_tot[,1]) , max(max_tot[,1]) ),col = 'red', pch = 19)

# - Cross validation on the functional model

rands <- runif(500, min=1, max=n ) #extract randomly 500 observations
Y.cv <- t(SA_log[rands, ])
data.fd.cv <- Data2fd(y = Y.cv , argvals = periods, basisobj = basis, lambda = opt.lambda )
xlist <- as.data.frame(xlist)
xlist.cv <- as.list(xlist[rands, ])

mod.cv <- fRegress.CV(y        = data.fd.cv,
                      xfdlist  = xlist.cv,      # perform a LOO CV on the extracted samples
                      betalist = blist)
mean.MSE <- mean.fd(mod.cv$errfd.cv)   #mean of the errors


errors <- eval.fd(tt,mean.MSE)      #evaluate the errors in order to plot them
errors <- sqrt(abs(errors))     

periods_err <- eval.fd(periods,mean.MSE)     #MSE on the sampled periods
periods_err <- sqrt(abs(periods_err))

x11()
plot(tt, errors,type='l', lwd=3, xlab = "periods", ylab = "RMSE",
     main = 'LOO CV on functional ITA18 model')   #plot of the MSE throughout continuous periods
points(periods, periods_err, col = 'red', pch=19)
legend("bottomleft", legend = c("RMSE of functional ITA18 on sampled period"),
       col = c("red"), pch = 19)



x11()
par(bg ="white")
plot(periods, periods_err, ylim = c(0,1), col = "#1B1B19", pch=19)
points(periods, RMSE_ita[2:37], col="#66121F", pch=19)
legend("topright", legend = c("RMSE of functional ITA18", "RMSE of scalar ITA18"),
       col = c("#1B1B19", "#66121F"), pch = 19)

abline(h=0.81)

#-------------------------------------------------------------------
# 4. MIXED EFFECT MODEL
#-------------------------------------------------------------------

library(lme4)
library(insight)
library(lattice)

# - EVENT_ID
center_CI = c()
min_CI = c()
max_CI = c()
PVRE1 = c()
for(i in c(7:43)){
  sd.mod <- lmer(data[,i] ~ source:dummy1 + source:dummy2 + data$sof + data$distance + 
                   path1 + path2 + site + (1|event_id), data = data[,c(1,2,3,4,5,6,7,8,i)])
  center_CI = cbind(center_CI, sd.mod@beta)
  CI = confint(sd.mod, oldNames=TRUE)
  min_CI = cbind(min_CI, CI[, 1])
  max_CI = cbind(max_CI, CI[, 2])
  sigma2_eps <- as.numeric(get_variance_residual(sd.mod))
  sigma2_b <- as.numeric(get_variance_random(sd.mod))
  PVRE1 = c(PVRE1, (sigma2_b/(sigma2_b+sigma2_eps)))
}

# Addictional material
sd.mod1_T1 = lmer(data$SA_1 ~ source:dummy1 + source:dummy2 + data$sof + data$distance + 
                    path1 + path2 + site + (1|data$event_id))
x11()
dotplot(ranef(sd.mod1_T1, condVar=T))

x11()
boxplot(residuals(sd.mod1_T1)~data$event_id, main='boxplot of the residual wrt event_id T=1')

x11()
plot(sd.mod1_T1, main="Residual for T=1", ylab=" ", xlab="T")

x11()
qqnorm(resid(sd.mod1_T1))
qqline(resid(sd.mod1_T1), col='red', lwd=2)

dev.off()

# - STATION_ID
center_CI = c()
min_CI = c()
max_CI = c()
PVRE2 = c()
for(i in c(7:43)){
  sd.mod <- lmer(data[,i] ~ source:dummy1 + source:dummy2 + data$sof + data$distance + 
                   path1 + path2 + site + (1|station_id), data = data[,c(1,2,3,4,5,6,7,8,i)])
  center_CI = cbind(center_CI, sd.mod@beta)
  #CI = confint(sd.mod, oldNames=TRUE)
  #min_CI = cbind(min_CI, CI[, 1])
  #max_CI = cbind(max_CI, CI[, 2])
  sigma2_eps <- as.numeric(get_variance_residual(sd.mod))
  sigma2_b <- as.numeric(get_variance_random(sd.mod))
  PVRE2 = c(PVRE2, (sigma2_b/(sigma2_b+sigma2_eps)))
}

x11()
par(mfrow=c(3,3))
for(i in c(1:9)){
  plot(periods, center_CI[i, ], ylim = c(min(min_CI[,i]), max(max_CI[,i])), col = 'red', ylab = '', xlab = 'T',
       main=paste(names[i]), pch = 19, type='b')
  abline(h = 0)
}

# Addictional material
sd.mod2_T1 = lmer(data$SA_1 ~ source:dummy1 + source:dummy2 + data$sof + data$distance + 
                    path1 + path2 + site + (1|data$station_id))
x11()
dotplot(ranef(sd.mod2_T1, condVar=T))

x11()
plot(sd.mod2_T1, main="Residual for T=1", ylab=" ", xlab="T")

x11()
qqnorm(resid(sd.mod2_T1))
qqline(resid(sd.mod2_T1), col='red', lwd=2)

dev.off()

# - EVENT_ID + STATION_ID
center_CI = c()
min_CI = c()
max_CI = c()
sd.err = c()
PVRE3 = c()
for(i in c(8:43)){
  sd.mod <- lmer(data[,i] ~ source:dummy1 + source:dummy2 + data$sof + data$distance + 
                   path1 + path2 + site + (1|event_id) + (1|station_id), data = data[,c(1,2,3,4,5,6,7,8,i)])
  sd.err = c(sd.err, sqrt(var(residuals(sd.mod))))
  # CI = confint(sd.mod, oldNames=TRUE)
  # min_CI = cbind(min_CI, CI[, 1])
  # max_CI = cbind(max_CI, CI[, 2])
  # center_CI = cbind(center_CI, (CI[, 1] + CI[, 2])/2)
  sigma2_eps <- as.numeric(get_variance_residual(sd.mod))
  sigma2_b <- as.numeric(get_variance_random(sd.mod))
  PVRE3 = c(PVRE3, (sigma2_b/(sigma2_b+sigma2_eps)))
}

x11()
par(bg = "white")
plot(periods, sd.err,  col='#E15269', type='b', pch=19, ylim=c(min(sd.ITA18.GMM[-1], sd.err), max(sd.ITA18.GMM[-1], sd.err)))
points(periods, sd.ITA18.GMM[-1], col="#66121F", type='b', pch=19)
legend("topright", legend = c("LMM ITA18", "GMM ITA18 standard"), col = c("#E15269", "#66121F"), pch = 19)


x11()
par(bg = "white", mfrow = c(3, 3))
for(i in c(1:9)){
  plot(periods, center_CI[i, ], col='#E15269', ylab = '', xlab = 'T', main=paste(names[i]), ylim = c(min(min_CI[i, ]), max(max_CI[i, ])), pch = 19)
  points(periods, min_CI[i, ], col='#E15269', pch='_')
  points(periods, max_CI[i, ], col='#E15269', pch='_')
  for(j in c(1:36)){
    segments(x0=periods[j],y0=min_CI[i, j],x1=periods[j],y1=max_CI[i, j], col = '#E15269')
  }
  abline(h = 0, type = '-')
}

x11()
par(mfrow=c(3,3), bg = "#F4F4F4")
for(i in c(1:9)){
  plot(periods, mean_tot[2:37,i], ylim =c(min(min_tot[,i]) , max(max_tot[,i]) ),col = "#1B1B19", pch=19, ylab = '',xlab = '' ,
       main=paste(names[i]))
  points(periods, max_tot[2:37,i], col = "#1B1B19", pch = '_')
  points(periods, min_tot[2:37,i] , col = "#1B1B19", pch='_')
  for(j in c(2:37)){
    segments(x0=periods[j-1],y0=min_tot[j,i],x1=periods[j-1],y1=max_tot[j,i])
  }
  abline(h = 0)
}

# Addictional material
sd.mod3_T1 = lmer(data$SA_1 ~ source:dummy1 + source:dummy2 + data$sof + data$distance + 
                    path1 + path2 + site + (1|data$event_id) + (1|data$station_id))

x11()
par(bg = "#F4F4F4")
dotplot(ranef(sd.mod3_T1, condVar=T))

x11()
par(bg = "#F4F4F4")
boxplot(residuals(sd.mod3_T1)~data$event_id, main='boxplot of the residual wrt event_id and station_id T=1', col='#E15269')

x11()
par(bg = "#F4F4F4")
plot(sd.mod3_T1, main="Residual for T=1", ylab=" ", xlab="T", col='#E15269')

x11()
qqnorm(resid(sd.mod3_T1))
qqline(resid(sd.mod3_T1), col='#E15269', lwd=2)

dev.off()

x11()
par(bg = "white")
plot(periods, PVRE3[-1], main="Intraclass correlation (ICC)", xlab="log(T)", ylab="PVRE", col='#E15269', type = 'b', ylim=c(0, max(PVRE3)), pch=19)
points(periods, PVRE2[-1], col="#66121F", type="b", pch=19)
points(periods, PVRE1[-1], col="#1B1B19", type="b", pch=19)
legend("bottomleft", legend = c("event", "station", "both"), col = c("#1B1B19", "#66121F", "#E15269"), pch = 19)
