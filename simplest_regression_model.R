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
# 1. MULTIPLE LINEAR REGRESSION MODEL FOR EACH STYLE OF FAULTING 
#    CONSIDERING ONE-PERIOD-AT-THE-TIME
#-------------------------------------------------------------------
library(cowplot)

###  LINEAR REGRESSION FOR EACH PERIOD sof == "NF"  ###
# - Selection of variables and outputs
numeric_values = ITA_18_full_dataset[sof == "NF", c(7,8,10)]
SA_log = log(ITA_18_full_dataset[sof == "NF", 11:47])

world_SA_log_data <- cbind(numeric_values,SA_log)


# - Saving coefficients of linear regression model and residual std error for each period
SA_0_reg <- lm(SA_log$SA_0 ~ .,numeric_values)
summary(SA_0_reg)

# - Saving coefficients beta_hat of the linear regression and the residual std error
residual.std.err <- sqrt(deviance(SA_0_reg)/SA_0_reg$df.residual)
coeff = as.data.frame(c(SA_0_reg$coefficients, residual.std.err))
names_coeff <- c("intercept", "magnitude", "distance", "vs30","residual std err")

p <- SA_0_reg$rank
n <- dim(SA_log)[1]
alpha <- 0.05
t_alpha = qt( 1-alpha/2, n-p )

names_IC <- c("lwr bound","mean", "upr bound")

# - Saving IC for the coefficients beta_hat
beta_hat_magnitude = SA_0_reg$coefficients[2]
se_magnitude = summary(SA_0_reg)[[4]][2,2]
IC_magnitude <- as.data.frame(c( beta_hat_magnitude - t_alpha * se_magnitude, 
                                 beta_hat_magnitude,
                                 beta_hat_magnitude + t_alpha * se_magnitude ), row.names = names_IC)

beta_hat_distance = SA_0_reg$coefficients[3]
se_distance = summary(SA_0_reg)[[4]][3,2]
IC_distance <- as.data.frame(c( beta_hat_distance - t_alpha * se_distance,
                                beta_hat_distance,
                                beta_hat_distance + t_alpha * se_distance ), row.names = names_IC)

beta_hat_vs30 = SA_0_reg$coefficients[4]
se_vs30 = summary(SA_0_reg)[[4]][4,2]
IC_vs30 <- as.data.frame(c( beta_hat_vs30 - t_alpha * se_vs30,
                            beta_hat_vs30,
                            beta_hat_vs30 + t_alpha * se_vs30 ), row.names = names_IC)

copie = log(ITA_18_full_dataset[sof == "NF",c(12:47)])

for (i in c(12:47)){
  numeric = ITA_18_full_dataset[sof == "NF", c(7,8,10)]
  numeric$output = copie[,1]
  regression = lm(output ~ ., numeric)
  res_err <- sqrt(deviance(regression)/regression$df.residual)
  coeff_tmp = as.data.frame(c(regression$coefficients, res_err))
  coeff = cbind(coeff, coeff_tmp)
  
  beta_hat_magnitude= regression$coefficients[2]
  se_magnitude = summary( regression)[[4]][2,2]
  IC_m_temp <- as.data.frame(c( beta_hat_magnitude - t_alpha * se_magnitude,
                                beta_hat_magnitude,
                                beta_hat_magnitude + t_alpha * se_magnitude ), row.names = names_IC)
  IC_magnitude <- cbind(IC_magnitude, IC_m_temp)
  
  beta_hat_distance = regression$coefficients[3]
  se_distance = summary( regression )[[4]][3,2]
  IC_d_temp <- as.data.frame(c( beta_hat_distance - t_alpha * se_distance,
                                beta_hat_distance,
                                beta_hat_distance + t_alpha * se_distance ), row.names = names_IC)
  IC_distance <- cbind(IC_distance, IC_d_temp)
  
  beta_hat_vs30 = regression$coefficients[4]
  se_vs30 = summary( regression )[[4]][4,2]
  IC_v_temp <- as.data.frame(c( beta_hat_vs30 - t_alpha * se_vs30,
                                beta_hat_vs30,
                                beta_hat_vs30 + t_alpha * se_vs30 ), row.names = names_IC)
  IC_vs30 <- cbind(IC_vs30, IC_v_temp)
  
  copie[,1]= NULL
}

colnames(coeff) = names(ITA_18_full_dataset[, c(11:47)])
colnames(IC_magnitude) = names(ITA_18_full_dataset[, c(11:47)])
colnames(IC_distance) = names(ITA_18_full_dataset[, c(11:47)])
colnames(IC_vs30) = names(ITA_18_full_dataset[, c(11:47)])

rm(copie)
rm(coeff_tmp)
rm(IC_m_temp)
rm(IC_d_temp)
rm(IC_v_temp)
rm(beta_hat_magnitude)
rm(beta_hat_distance)
rm(beta_hat_vs30)


# - Creation of the 37 dataframe for the plot 

IC_0 <- as.data.frame(t(cbind(IC_magnitude$SA_0, IC_distance$SA_0, IC_vs30$SA_0)))
IC_0.01 <- as.data.frame(t(cbind(IC_magnitude$SA_0.01, IC_distance$SA_0.01, IC_vs30$SA_0.01)))
IC_0.025<- as.data.frame(t(cbind(IC_magnitude$SA_0.025, IC_distance$SA_0.025, IC_vs30$SA_0.025)))
IC_0.04<- as.data.frame(t(cbind(IC_magnitude$SA_0.04, IC_distance$SA_0.04, IC_vs30$SA_0.04)))
IC_0.05<- as.data.frame(t(cbind(IC_magnitude$SA_0.05, IC_distance$SA_0.05, IC_vs30$SA_0.05)))
IC_0.07<- as.data.frame(t(cbind(IC_magnitude$SA_0.07, IC_distance$SA_0.07, IC_vs30$SA_0.07)))
IC_0.1 <- as.data.frame(t(cbind(IC_magnitude$SA_0.1, IC_distance$SA_0.1, IC_vs30$SA_0.1)))
IC_0.15 <- as.data.frame(t(cbind(IC_magnitude$SA_0.15, IC_distance$SA_0.15, IC_vs30$SA_0.15)))
IC_0.2 <- as.data.frame(t(cbind(IC_magnitude$SA_0.2, IC_distance$SA_0.2, IC_vs30$SA_0.2)))
IC_0.25 <- as.data.frame(t(cbind(IC_magnitude$SA_0.25, IC_distance$SA_0.25, IC_vs30$SA_0.25)))
IC_0.3<- as.data.frame(t(cbind(IC_magnitude$SA_0.3, IC_distance$SA_0.3, IC_vs30$SA_0.3)))
IC_0.35 <- as.data.frame(t(cbind(IC_magnitude$SA_0.35, IC_distance$SA_0.35, IC_vs30$SA_0.35)))
IC_0.4 <- as.data.frame(t(cbind(IC_magnitude$SA_0.4, IC_distance$SA_0.4, IC_vs30$SA_0.4)))
IC_0.45 <- as.data.frame(t(cbind(IC_magnitude$SA_0.45, IC_distance$SA_0.45, IC_vs30$SA_0.45)))
IC_0.5 <- as.data.frame(t(cbind(IC_magnitude$SA_0.5, IC_distance$SA_0.5, IC_vs30$SA_0.5)))
IC_0.6<- as.data.frame(t(cbind(IC_magnitude$SA_0.6, IC_distance$SA_0.6, IC_vs30$SA_0.6)))
IC_0.7 <- as.data.frame(t(cbind(IC_magnitude$SA_0.7, IC_distance$SA_0.7, IC_vs30$SA_0.7)))
IC_0.75 <- as.data.frame(t(cbind(IC_magnitude$SA_0.75, IC_distance$SA_0.75, IC_vs30$SA_0.75)))
IC_0.8 <- as.data.frame(t(cbind(IC_magnitude$SA_0.8, IC_distance$SA_0.8, IC_vs30$SA_0.8)))
IC_0.9 <- as.data.frame(t(cbind(IC_magnitude$SA_0.9, IC_distance$SA_0.9, IC_vs30$SA_0.9)))
IC_1 <- as.data.frame(t(cbind(IC_magnitude$SA_1, IC_distance$SA_1, IC_vs30$SA_1)))
IC_1.2 <- as.data.frame(t(cbind(IC_magnitude$SA_1.2, IC_distance$SA_1.2, IC_vs30$SA_1.2)))
IC_1.4 <- as.data.frame(t(cbind(IC_magnitude$SA_1.4, IC_distance$SA_1.4, IC_vs30$SA_1.4)))
IC_1.6<- as.data.frame(t(cbind(IC_magnitude$SA_1.6, IC_distance$SA_1.6, IC_vs30$SA_1.6)))
IC_1.8 <- as.data.frame(t(cbind(IC_magnitude$SA_1.8, IC_distance$SA_1.8, IC_vs30$SA_1.8)))
IC_2 <- as.data.frame(t(cbind(IC_magnitude$SA_2, IC_distance$SA_2, IC_vs30$SA_2)))
IC_2.5 <- as.data.frame(t(cbind(IC_magnitude$SA_2.5, IC_distance$SA_2.5, IC_vs30$SA_2.5)))
IC_3 <- as.data.frame(t(cbind(IC_magnitude$SA_3, IC_distance$SA_3, IC_vs30$SA_3)))
IC_3.5 <- as.data.frame(t(cbind(IC_magnitude$SA_3.5, IC_distance$SA_3.5, IC_vs30$SA_3.5)))
IC_4 <- as.data.frame(t(cbind(IC_magnitude$SA_4, IC_distance$SA_4, IC_vs30$SA_4)))
IC_4.5 <- as.data.frame(t(cbind(IC_magnitude$SA_4.5, IC_distance$SA_4.5, IC_vs30$SA_4.5)))
IC_5<- as.data.frame(t(cbind(IC_magnitude$SA_5, IC_distance$SA_5, IC_vs30$SA_5)))
IC_6 <- as.data.frame(t(cbind(IC_magnitude$SA_6, IC_distance$SA_6, IC_vs30$SA_6)))
IC_7 <- as.data.frame(t(cbind(IC_magnitude$SA_7, IC_distance$SA_7, IC_vs30$SA_7)))
IC_8 <- as.data.frame(t(cbind(IC_magnitude$SA_8, IC_distance$SA_8, IC_vs30$SA_8)))
IC_9 <- as.data.frame(t(cbind(IC_magnitude$SA_9, IC_distance$SA_9, IC_vs30$SA_9)))
IC_10 <- as.data.frame(t(cbind(IC_magnitude$SA_10, IC_distance$SA_10, IC_vs30$SA_10)))



# - Plot of the confidence intervals for the thre variables + residuals

# magnitude
magnitude_IC = rbind(IC_0[1,], IC_0.01[1,], IC_0.025[1,], IC_0.04[1,], IC_0.05[1,],  IC_0.07[1,], IC_0.1[1,], IC_0.15[1,],IC_0.2[1,], IC_0.25[1,], IC_0.3[1,], IC_0.35[1,],IC_0.4[1,], IC_0.45[1,],  IC_0.5[1,],IC_0.6[1,], IC_0.7[1,], IC_0.75[1,], IC_0.8[1,],IC_0.9[1,], IC_1[1,])
magnitude_IC = rbind(magnitude_IC,IC_1.2[1,], IC_1.4[1,], IC_1.6[1,], IC_1.8[1,], IC_2[1,],  IC_2.5[1,], IC_3[1,], IC_3.5[1,],IC_4[1,], IC_4.5[1,], IC_5[1,], IC_6[1,], IC_7[1,],  IC_8[1,], IC_9[1,], IC_10[1,])
magnitude_IC$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

magnitude.plot = ggplot(magnitude_IC, aes(x = n, y = V2, ymin = V1, ymax = V3))+ geom_errorbar(width = 0.2, color = "#1B1B19") +
  geom_point(size = 1.5, color = "#1B1B19") + ylab("magnitude coefficient (NF)")  + xlab("log(T)") + 
  theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))


# distance
distance_IC = rbind(IC_0[2,], IC_0.01[2,], IC_0.025[2,], IC_0.04[2,], IC_0.05[2,],  IC_0.07[2,], IC_0.1[2,], IC_0.15[2,],IC_0.2[2,], IC_0.25[2,], IC_0.3[2,], IC_0.35[2,],IC_0.4[2,], IC_0.45[2,],  IC_0.5[2,],IC_0.6[2,], IC_0.7[2,], IC_0.75[2,], IC_0.8[2,],IC_0.9[2,], IC_1[2,])
distance_IC = rbind(distance_IC,IC_1.2[2,], IC_1.4[2,], IC_1.6[2,], IC_1.8[2,], IC_2[2,],  IC_2.5[2,], IC_3[2,], IC_3.5[2,],IC_4[2,], IC_4.5[2,], IC_5[2,], IC_6[2,], IC_7[2,],  IC_8[2,], IC_9[2,], IC_10[2,])
distance_IC$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

distance.plot = ggplot(distance_IC, aes(x = n, y = V2, ymin = V1, ymax = V3))+ geom_errorbar(width = 0.2, color = "#1B1B19") +
  geom_point(size = 1.5, color = "#1B1B19") + ylab("distance coefficient (NF)") + xlab("log(T)") + 
  theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))

# vs30
vs30_IC = rbind(IC_0[3,], IC_0.01[3,], IC_0.025[3,], IC_0.04[3,], IC_0.05[3,],  IC_0.07[3,], IC_0.1[3,], IC_0.15[3,],IC_0.2[3,], IC_0.25[3,], IC_0.3[3,], IC_0.35[3,],IC_0.4[3,], IC_0.45[3,],  IC_0.5[3,],IC_0.6[3,], IC_0.7[3,], IC_0.75[3,], IC_0.8[3,],IC_0.9[3,], IC_1[3,])
vs30_IC = rbind(vs30_IC,IC_1.2[3,], IC_1.4[3,], IC_1.6[3,], IC_1.8[3,], IC_2[3,],  IC_2.5[3,], IC_3[3,], IC_3.5[3,],IC_4[3,], IC_4.5[3,], IC_5[3,], IC_6[3,], IC_7[3,],  IC_8[3,], IC_9[3,], IC_10[3,])
vs30_IC$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

vs30.plot = ggplot(vs30_IC, aes(x = n, y = V2, ymin = V1, ymax = V3))+ geom_errorbar(width = 0.2, color = "#1B1B19") +
  geom_point(size = 1.5, color = "#1B1B19") + ylab("vs30 coefficient (NF)") + xlab("T")+ xlab("log(T)") + 
  theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))

x11()
par(bg = "#F4F4F4")
plot_grid(magnitude.plot, distance.plot, vs30.plot, nrow = 3)

# Residuals
Residuals = as.data.table(t(coeff[5,]))
Residuals$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

res.plot.NF = ggplot(Residuals, aes(x = n, y = V1))+ geom_point(size = 1.5, color = "#1B1B19") + ylab("residual std error (NF)") + xlab("log(T)") + 
  theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))


rm(coeff)
rm(vs30_IC)
rm(distance_IC)
rm(magnitude_IC)
rm(Residuals)

rm(IC_magnitude)
rm(IC_vs30)
rm(IC_distance)

###  LINEAR REGRESSION FOR EACH PERIOD sof == "TF"  ###
# - Selection of variables and outputs
numeric_values = ITA_18_full_dataset[sof == "TF", c(7,8,10)]
SA_log = log(ITA_18_full_dataset[sof == "TF", 11:47])

world_SA_log_data <- cbind(numeric_values,SA_log)


# - Saving coefficients of linear regression model and residual std error for each period
SA_0_reg <- lm(SA_log$SA_0 ~ .,numeric_values)
summary(SA_0_reg)

# - Saving coefficients beta_hat of the linear regression and the residual std error
residual.std.err <- sqrt(deviance(SA_0_reg)/SA_0_reg$df.residual)
coeff = as.data.frame(c(SA_0_reg$coefficients, residual.std.err))
names_coeff <- c("intercept", "magnitude", "distance", "vs30","residual std err")

p <- SA_0_reg$rank
n <- dim(SA_log)[1]
alpha <- 0.05
t_alpha = qt( 1-alpha/2, n-p )

names_IC <- c("lwr bound","mean", "upr bound")

# - Saving IC for the coefficients beta_hat
beta_hat_magnitude = SA_0_reg$coefficients[2]
se_magnitude = summary(SA_0_reg)[[4]][2,2]
IC_magnitude <- as.data.frame(c( beta_hat_magnitude - t_alpha * se_magnitude, 
                                 beta_hat_magnitude,
                                 beta_hat_magnitude + t_alpha * se_magnitude ), row.names = names_IC)

beta_hat_distance = SA_0_reg$coefficients[3]
se_distance = summary(SA_0_reg)[[4]][3,2]
IC_distance <- as.data.frame(c( beta_hat_distance - t_alpha * se_distance,
                                beta_hat_distance,
                                beta_hat_distance + t_alpha * se_distance ), row.names = names_IC)

beta_hat_vs30 = SA_0_reg$coefficients[4]
se_vs30 = summary(SA_0_reg)[[4]][4,2]
IC_vs30 <- as.data.frame(c( beta_hat_vs30 - t_alpha * se_vs30,
                            beta_hat_vs30,
                            beta_hat_vs30 + t_alpha * se_vs30 ), row.names = names_IC)

copie = log(ITA_18_full_dataset[sof == "TF",c(12:47)])

for (i in c(12:47)){
  numeric = ITA_18_full_dataset[sof == "TF", c(7,8,10)]
  numeric$output = copie[,1]
  regression = lm(output ~ ., numeric)
  res_err <- sqrt(deviance(regression)/regression$df.residual)
  coeff_tmp = as.data.frame(c(regression$coefficients, res_err))
  coeff = cbind(coeff, coeff_tmp)
  
  beta_hat_magnitude= regression$coefficients[2]
  se_magnitude = summary( regression)[[4]][2,2]
  IC_m_temp <- as.data.frame(c( beta_hat_magnitude - t_alpha * se_magnitude,
                                beta_hat_magnitude,
                                beta_hat_magnitude + t_alpha * se_magnitude ), row.names = names_IC)
  IC_magnitude <- cbind(IC_magnitude, IC_m_temp)
  
  beta_hat_distance = regression$coefficients[3]
  se_distance = summary( regression )[[4]][3,2]
  IC_d_temp <- as.data.frame(c( beta_hat_distance - t_alpha * se_distance,
                                beta_hat_distance,
                                beta_hat_distance + t_alpha * se_distance ), row.names = names_IC)
  IC_distance <- cbind(IC_distance, IC_d_temp)
  
  beta_hat_vs30 = regression$coefficients[4]
  se_vs30 = summary( regression )[[4]][4,2]
  IC_v_temp <- as.data.frame(c( beta_hat_vs30 - t_alpha * se_vs30,
                                beta_hat_vs30,
                                beta_hat_vs30 + t_alpha * se_vs30 ), row.names = names_IC)
  IC_vs30 <- cbind(IC_vs30, IC_v_temp)
  
  copie[,1]= NULL
}

colnames(coeff) = names(ITA_18_full_dataset[, c(11:47)])
colnames(IC_magnitude) = names(ITA_18_full_dataset[, c(11:47)])
colnames(IC_distance) = names(ITA_18_full_dataset[, c(11:47)])
colnames(IC_vs30) = names(ITA_18_full_dataset[, c(11:47)])

rm(copie)
rm(coeff_tmp)
rm(IC_m_temp)
rm(IC_d_temp)
rm(IC_v_temp)
rm(beta_hat_magnitude)
rm(beta_hat_distance)
rm(beta_hat_vs30)



# - Creation of the 37 dataframe for the plot 

IC_0 <- as.data.frame(t(cbind(IC_magnitude$SA_0, IC_distance$SA_0, IC_vs30$SA_0)))
IC_0.01 <- as.data.frame(t(cbind(IC_magnitude$SA_0.01, IC_distance$SA_0.01, IC_vs30$SA_0.01)))
IC_0.025<- as.data.frame(t(cbind(IC_magnitude$SA_0.025, IC_distance$SA_0.025, IC_vs30$SA_0.025)))
IC_0.04<- as.data.frame(t(cbind(IC_magnitude$SA_0.04, IC_distance$SA_0.04, IC_vs30$SA_0.04)))
IC_0.05<- as.data.frame(t(cbind(IC_magnitude$SA_0.05, IC_distance$SA_0.05, IC_vs30$SA_0.05)))
IC_0.07<- as.data.frame(t(cbind(IC_magnitude$SA_0.07, IC_distance$SA_0.07, IC_vs30$SA_0.07)))
IC_0.1 <- as.data.frame(t(cbind(IC_magnitude$SA_0.1, IC_distance$SA_0.1, IC_vs30$SA_0.1)))
IC_0.15 <- as.data.frame(t(cbind(IC_magnitude$SA_0.15, IC_distance$SA_0.15, IC_vs30$SA_0.15)))
IC_0.2 <- as.data.frame(t(cbind(IC_magnitude$SA_0.2, IC_distance$SA_0.2, IC_vs30$SA_0.2)))
IC_0.25 <- as.data.frame(t(cbind(IC_magnitude$SA_0.25, IC_distance$SA_0.25, IC_vs30$SA_0.25)))
IC_0.3<- as.data.frame(t(cbind(IC_magnitude$SA_0.3, IC_distance$SA_0.3, IC_vs30$SA_0.3)))
IC_0.35 <- as.data.frame(t(cbind(IC_magnitude$SA_0.35, IC_distance$SA_0.35, IC_vs30$SA_0.35)))
IC_0.4 <- as.data.frame(t(cbind(IC_magnitude$SA_0.4, IC_distance$SA_0.4, IC_vs30$SA_0.4)))
IC_0.45 <- as.data.frame(t(cbind(IC_magnitude$SA_0.45, IC_distance$SA_0.45, IC_vs30$SA_0.45)))
IC_0.5 <- as.data.frame(t(cbind(IC_magnitude$SA_0.5, IC_distance$SA_0.5, IC_vs30$SA_0.5)))
IC_0.6<- as.data.frame(t(cbind(IC_magnitude$SA_0.6, IC_distance$SA_0.6, IC_vs30$SA_0.6)))
IC_0.7 <- as.data.frame(t(cbind(IC_magnitude$SA_0.7, IC_distance$SA_0.7, IC_vs30$SA_0.7)))
IC_0.75 <- as.data.frame(t(cbind(IC_magnitude$SA_0.75, IC_distance$SA_0.75, IC_vs30$SA_0.75)))
IC_0.8 <- as.data.frame(t(cbind(IC_magnitude$SA_0.8, IC_distance$SA_0.8, IC_vs30$SA_0.8)))
IC_0.9 <- as.data.frame(t(cbind(IC_magnitude$SA_0.9, IC_distance$SA_0.9, IC_vs30$SA_0.9)))
IC_1 <- as.data.frame(t(cbind(IC_magnitude$SA_1, IC_distance$SA_1, IC_vs30$SA_1)))
IC_1.2 <- as.data.frame(t(cbind(IC_magnitude$SA_1.2, IC_distance$SA_1.2, IC_vs30$SA_1.2)))
IC_1.4 <- as.data.frame(t(cbind(IC_magnitude$SA_1.4, IC_distance$SA_1.4, IC_vs30$SA_1.4)))
IC_1.6<- as.data.frame(t(cbind(IC_magnitude$SA_1.6, IC_distance$SA_1.6, IC_vs30$SA_1.6)))
IC_1.8 <- as.data.frame(t(cbind(IC_magnitude$SA_1.8, IC_distance$SA_1.8, IC_vs30$SA_1.8)))
IC_2 <- as.data.frame(t(cbind(IC_magnitude$SA_2, IC_distance$SA_2, IC_vs30$SA_2)))
IC_2.5 <- as.data.frame(t(cbind(IC_magnitude$SA_2.5, IC_distance$SA_2.5, IC_vs30$SA_2.5)))
IC_3 <- as.data.frame(t(cbind(IC_magnitude$SA_3, IC_distance$SA_3, IC_vs30$SA_3)))
IC_3.5 <- as.data.frame(t(cbind(IC_magnitude$SA_3.5, IC_distance$SA_3.5, IC_vs30$SA_3.5)))
IC_4 <- as.data.frame(t(cbind(IC_magnitude$SA_4, IC_distance$SA_4, IC_vs30$SA_4)))
IC_4.5 <- as.data.frame(t(cbind(IC_magnitude$SA_4.5, IC_distance$SA_4.5, IC_vs30$SA_4.5)))
IC_5<- as.data.frame(t(cbind(IC_magnitude$SA_5, IC_distance$SA_5, IC_vs30$SA_5)))
IC_6 <- as.data.frame(t(cbind(IC_magnitude$SA_6, IC_distance$SA_6, IC_vs30$SA_6)))
IC_7 <- as.data.frame(t(cbind(IC_magnitude$SA_7, IC_distance$SA_7, IC_vs30$SA_7)))
IC_8 <- as.data.frame(t(cbind(IC_magnitude$SA_8, IC_distance$SA_8, IC_vs30$SA_8)))
IC_9 <- as.data.frame(t(cbind(IC_magnitude$SA_9, IC_distance$SA_9, IC_vs30$SA_9)))
IC_10 <- as.data.frame(t(cbind(IC_magnitude$SA_10, IC_distance$SA_10, IC_vs30$SA_10)))



# - Plot of the confidence intervals for the thre variables + residuals

# magnitude
magnitude_IC = rbind(IC_0[1,], IC_0.01[1,], IC_0.025[1,], IC_0.04[1,], IC_0.05[1,],  IC_0.07[1,], IC_0.1[1,], IC_0.15[1,],IC_0.2[1,], IC_0.25[1,], IC_0.3[1,], IC_0.35[1,],IC_0.4[1,], IC_0.45[1,],  IC_0.5[1,],IC_0.6[1,], IC_0.7[1,], IC_0.75[1,], IC_0.8[1,],IC_0.9[1,], IC_1[1,])
magnitude_IC = rbind(magnitude_IC,IC_1.2[1,], IC_1.4[1,], IC_1.6[1,], IC_1.8[1,], IC_2[1,],  IC_2.5[1,], IC_3[1,], IC_3.5[1,],IC_4[1,], IC_4.5[1,], IC_5[1,], IC_6[1,], IC_7[1,],  IC_8[1,], IC_9[1,], IC_10[1,])
magnitude_IC$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

magnitude.plot = ggplot(magnitude_IC, aes(x = n, y = V2, ymin = V1, ymax = V3))+ geom_errorbar(width = 0.2, color = "#66121F") +
  geom_point(size = 1.5, color = "#66121F") + ylab("magnitude coefficient (TF)")  + xlab("log(T)") + theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))

# distance
distance_IC = rbind(IC_0[2,], IC_0.01[2,], IC_0.025[2,], IC_0.04[2,], IC_0.05[2,],  IC_0.07[2,], IC_0.1[2,], IC_0.15[2,],IC_0.2[2,], IC_0.25[2,], IC_0.3[2,], IC_0.35[2,],IC_0.4[2,], IC_0.45[2,],  IC_0.5[2,],IC_0.6[2,], IC_0.7[2,], IC_0.75[2,], IC_0.8[2,],IC_0.9[2,], IC_1[2,])
distance_IC = rbind(distance_IC,IC_1.2[2,], IC_1.4[2,], IC_1.6[2,], IC_1.8[2,], IC_2[2,],  IC_2.5[2,], IC_3[2,], IC_3.5[2,],IC_4[2,], IC_4.5[2,], IC_5[2,], IC_6[2,], IC_7[2,],  IC_8[2,], IC_9[2,], IC_10[2,])
distance_IC$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

distance.plot = ggplot(distance_IC, aes(x = n, y = V2, ymin = V1, ymax = V3))+ geom_errorbar(width = 0.2, color = "#66121F") +
  geom_point(size = 1.5, color = "#66121F") + ylab("distance coefficient (TF)") + xlab("log(T)") + theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))

# vs30
vs30_IC = rbind(IC_0[3,], IC_0.01[3,], IC_0.025[3,], IC_0.04[3,], IC_0.05[3,],  IC_0.07[3,], IC_0.1[3,], IC_0.15[3,],IC_0.2[3,], IC_0.25[3,], IC_0.3[3,], IC_0.35[3,],IC_0.4[3,], IC_0.45[3,],  IC_0.5[3,],IC_0.6[3,], IC_0.7[3,], IC_0.75[3,], IC_0.8[3,],IC_0.9[3,], IC_1[3,])
vs30_IC = rbind(vs30_IC,IC_1.2[3,], IC_1.4[3,], IC_1.6[3,], IC_1.8[3,], IC_2[3,],  IC_2.5[3,], IC_3[3,], IC_3.5[3,],IC_4[3,], IC_4.5[3,], IC_5[3,], IC_6[3,], IC_7[3,],  IC_8[3,], IC_9[3,], IC_10[3,])
vs30_IC$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

vs30.plot = ggplot(vs30_IC, aes(x = n, y = V2, ymin = V1, ymax = V3))+ geom_errorbar(width = 0.2, color = "#66121F") +
  geom_point(size = 1.5, color = "#66121F") + ylab("vs30 coefficient (TF)")   + xlab("log(T)") +
  theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))

x11()
par(bg = "#F4F4F4")
plot_grid(magnitude.plot, distance.plot, vs30.plot, nrow = 3)

# Residuals
Residuals = as.data.table(t(coeff[5,]))
Residuals$n = c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10)

res.plot.TF = ggplot(Residuals, aes(x = n, y = V1))+ geom_point(size = 1.5, color = "#66121F") + ylab("residual std error (TF)") + xlab("log(T)") +
  theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))


rm(coeff)
rm(vs30_IC)
rm(distance_IC)
rm(magnitude_IC)
rm(Residuals)

rm(IC_magnitude)
rm(IC_vs30)
rm(IC_distance)

###  LINEAR REGRESSION FOR EACH PERIOD sof == "SS"  ###
# - Selection of variables and outputs
numeric_values = ITA_18_full_dataset[sof == "SS", c(7,8,10)]
SA_log = log(ITA_18_full_dataset[sof == "SS", 11:47])

world_SA_log_data <- cbind(numeric_values,SA_log)


# - Saving coefficients of linear regression model and residual std error for each period
SA_0_reg <- lm(SA_log$SA_0 ~ .,numeric_values)
summary(SA_0_reg)

# - Saving coefficients beta_hat of the linear regression and the residual std error
residual.std.err <- sqrt(deviance(SA_0_reg)/SA_0_reg$df.residual)
coeff = as.data.frame(c(SA_0_reg$coefficients, residual.std.err))
names_coeff <- c("intercept", "magnitude", "distance", "vs30","residual std err")

p <- SA_0_reg$rank
n <- dim(SA_log)[1]
alpha <- 0.05
t_alpha = qt( 1-alpha/2, n-p )

names_IC <- c("lwr bound","mean", "upr bound")

# - Saving IC for the coefficients beta_hat
beta_hat_magnitude = SA_0_reg$coefficients[2]
se_magnitude = summary(SA_0_reg)[[4]][2,2]
IC_magnitude <- as.data.frame(c( beta_hat_magnitude - t_alpha * se_magnitude, 
                                 beta_hat_magnitude,
                                 beta_hat_magnitude + t_alpha * se_magnitude ), row.names = names_IC)

beta_hat_distance = SA_0_reg$coefficients[3]
se_distance = summary(SA_0_reg)[[4]][3,2]
IC_distance <- as.data.frame(c( beta_hat_distance - t_alpha * se_distance,
                                beta_hat_distance,
                                beta_hat_distance + t_alpha * se_distance ), row.names = names_IC)

beta_hat_vs30 = SA_0_reg$coefficients[4]
se_vs30 = summary(SA_0_reg)[[4]][4,2]
IC_vs30 <- as.data.frame(c( beta_hat_vs30 - t_alpha * se_vs30,
                            beta_hat_vs30,
                            beta_hat_vs30 + t_alpha * se_vs30 ), row.names = names_IC)

copie = log(ITA_18_full_dataset[sof == "SS",c(12:47)])


for (i in c(12:47)){
  numeric = ITA_18_full_dataset[sof == "SS", c(7,8,10)]
  numeric$output = copie[,1]
  regression = lm(output ~ ., numeric)
  res_err <- sqrt(deviance(regression)/regression$df.residual)
  coeff_tmp = as.data.frame(c(regression$coefficients, res_err))
  coeff = cbind(coeff, coeff_tmp)
  
  beta_hat_magnitude= regression$coefficients[2]
  se_magnitude = summary( regression)[[4]][2,2]
  IC_m_temp <- as.data.frame(c( beta_hat_magnitude - t_alpha * se_magnitude,
                                beta_hat_magnitude,
                                beta_hat_magnitude + t_alpha * se_magnitude ), row.names = names_IC)
  IC_magnitude <- cbind(IC_magnitude, IC_m_temp)
  
  beta_hat_distance = regression$coefficients[3]
  se_distance = summary( regression )[[4]][3,2]
  IC_d_temp <- as.data.frame(c( beta_hat_distance - t_alpha * se_distance,
                                beta_hat_distance,
                                beta_hat_distance + t_alpha * se_distance ), row.names = names_IC)
  IC_distance <- cbind(IC_distance, IC_d_temp)
  
  beta_hat_vs30 = regression$coefficients[4]
  se_vs30 = summary( regression )[[4]][4,2]
  IC_v_temp <- as.data.frame(c( beta_hat_vs30 - t_alpha * se_vs30,
                                beta_hat_vs30,
                                beta_hat_vs30 + t_alpha * se_vs30 ), row.names = names_IC)
  IC_vs30 <- cbind(IC_vs30, IC_v_temp)
  
  copie[,1]= NULL
}

colnames(coeff) = names(ITA_18_full_dataset[, c(11:47)])
colnames(IC_magnitude) = names(ITA_18_full_dataset[, c(11:47)])
colnames(IC_distance) = names(ITA_18_full_dataset[, c(11:47)])
colnames(IC_vs30) = names(ITA_18_full_dataset[, c(11:47)])

rm(copie)
rm(coeff_tmp)
rm(IC_m_temp)
rm(IC_d_temp)
rm(IC_v_temp)
rm(beta_hat_magnitude)
rm(beta_hat_distance)
rm(beta_hat_vs30)




# - Creation of the 37 dataframe for the plot 

IC_0 <- as.data.frame(t(cbind(IC_magnitude$SA_0, IC_distance$SA_0, IC_vs30$SA_0)))
IC_0.01 <- as.data.frame(t(cbind(IC_magnitude$SA_0.01, IC_distance$SA_0.01, IC_vs30$SA_0.01)))
IC_0.025<- as.data.frame(t(cbind(IC_magnitude$SA_0.025, IC_distance$SA_0.025, IC_vs30$SA_0.025)))
IC_0.04<- as.data.frame(t(cbind(IC_magnitude$SA_0.04, IC_distance$SA_0.04, IC_vs30$SA_0.04)))
IC_0.05<- as.data.frame(t(cbind(IC_magnitude$SA_0.05, IC_distance$SA_0.05, IC_vs30$SA_0.05)))
IC_0.07<- as.data.frame(t(cbind(IC_magnitude$SA_0.07, IC_distance$SA_0.07, IC_vs30$SA_0.07)))
IC_0.1 <- as.data.frame(t(cbind(IC_magnitude$SA_0.1, IC_distance$SA_0.1, IC_vs30$SA_0.1)))
IC_0.15 <- as.data.frame(t(cbind(IC_magnitude$SA_0.15, IC_distance$SA_0.15, IC_vs30$SA_0.15)))
IC_0.2 <- as.data.frame(t(cbind(IC_magnitude$SA_0.2, IC_distance$SA_0.2, IC_vs30$SA_0.2)))
IC_0.25 <- as.data.frame(t(cbind(IC_magnitude$SA_0.25, IC_distance$SA_0.25, IC_vs30$SA_0.25)))
IC_0.3<- as.data.frame(t(cbind(IC_magnitude$SA_0.3, IC_distance$SA_0.3, IC_vs30$SA_0.3)))
IC_0.35 <- as.data.frame(t(cbind(IC_magnitude$SA_0.35, IC_distance$SA_0.35, IC_vs30$SA_0.35)))
IC_0.4 <- as.data.frame(t(cbind(IC_magnitude$SA_0.4, IC_distance$SA_0.4, IC_vs30$SA_0.4)))
IC_0.45 <- as.data.frame(t(cbind(IC_magnitude$SA_0.45, IC_distance$SA_0.45, IC_vs30$SA_0.45)))
IC_0.5 <- as.data.frame(t(cbind(IC_magnitude$SA_0.5, IC_distance$SA_0.5, IC_vs30$SA_0.5)))
IC_0.6<- as.data.frame(t(cbind(IC_magnitude$SA_0.6, IC_distance$SA_0.6, IC_vs30$SA_0.6)))
IC_0.7 <- as.data.frame(t(cbind(IC_magnitude$SA_0.7, IC_distance$SA_0.7, IC_vs30$SA_0.7)))
IC_0.75 <- as.data.frame(t(cbind(IC_magnitude$SA_0.75, IC_distance$SA_0.75, IC_vs30$SA_0.75)))
IC_0.8 <- as.data.frame(t(cbind(IC_magnitude$SA_0.8, IC_distance$SA_0.8, IC_vs30$SA_0.8)))
IC_0.9 <- as.data.frame(t(cbind(IC_magnitude$SA_0.9, IC_distance$SA_0.9, IC_vs30$SA_0.9)))
IC_1 <- as.data.frame(t(cbind(IC_magnitude$SA_1, IC_distance$SA_1, IC_vs30$SA_1)))
IC_1.2 <- as.data.frame(t(cbind(IC_magnitude$SA_1.2, IC_distance$SA_1.2, IC_vs30$SA_1.2)))
IC_1.4 <- as.data.frame(t(cbind(IC_magnitude$SA_1.4, IC_distance$SA_1.4, IC_vs30$SA_1.4)))
IC_1.6<- as.data.frame(t(cbind(IC_magnitude$SA_1.6, IC_distance$SA_1.6, IC_vs30$SA_1.6)))
IC_1.8 <- as.data.frame(t(cbind(IC_magnitude$SA_1.8, IC_distance$SA_1.8, IC_vs30$SA_1.8)))
IC_2 <- as.data.frame(t(cbind(IC_magnitude$SA_2, IC_distance$SA_2, IC_vs30$SA_2)))
IC_2.5 <- as.data.frame(t(cbind(IC_magnitude$SA_2.5, IC_distance$SA_2.5, IC_vs30$SA_2.5)))
IC_3 <- as.data.frame(t(cbind(IC_magnitude$SA_3, IC_distance$SA_3, IC_vs30$SA_3)))
IC_3.5 <- as.data.frame(t(cbind(IC_magnitude$SA_3.5, IC_distance$SA_3.5, IC_vs30$SA_3.5)))
IC_4 <- as.data.frame(t(cbind(IC_magnitude$SA_4, IC_distance$SA_4, IC_vs30$SA_4)))
IC_4.5 <- as.data.frame(t(cbind(IC_magnitude$SA_4.5, IC_distance$SA_4.5, IC_vs30$SA_4.5)))
IC_5<- as.data.frame(t(cbind(IC_magnitude$SA_5, IC_distance$SA_5, IC_vs30$SA_5)))
IC_6 <- as.data.frame(t(cbind(IC_magnitude$SA_6, IC_distance$SA_6, IC_vs30$SA_6)))
IC_7 <- as.data.frame(t(cbind(IC_magnitude$SA_7, IC_distance$SA_7, IC_vs30$SA_7)))
IC_8 <- as.data.frame(t(cbind(IC_magnitude$SA_8, IC_distance$SA_8, IC_vs30$SA_8)))
IC_9 <- as.data.frame(t(cbind(IC_magnitude$SA_9, IC_distance$SA_9, IC_vs30$SA_9)))
IC_10 <- as.data.frame(t(cbind(IC_magnitude$SA_10, IC_distance$SA_10, IC_vs30$SA_10)))


# - Plot of the confidence intervals for the thre variables + residuals

# magnitude
magnitude_IC = rbind(IC_0[1,], IC_0.01[1,], IC_0.025[1,], IC_0.04[1,], IC_0.05[1,],  IC_0.07[1,], IC_0.1[1,], IC_0.15[1,],IC_0.2[1,], IC_0.25[1,], IC_0.3[1,], IC_0.35[1,],IC_0.4[1,], IC_0.45[1,],  IC_0.5[1,],IC_0.6[1,], IC_0.7[1,], IC_0.75[1,], IC_0.8[1,],IC_0.9[1,], IC_1[1,])
magnitude_IC = rbind(magnitude_IC,IC_1.2[1,], IC_1.4[1,], IC_1.6[1,], IC_1.8[1,], IC_2[1,],  IC_2.5[1,], IC_3[1,], IC_3.5[1,],IC_4[1,], IC_4.5[1,], IC_5[1,], IC_6[1,], IC_7[1,],  IC_8[1,], IC_9[1,], IC_10[1,])
magnitude_IC$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

magnitude.plot = ggplot(magnitude_IC, aes(x = n, y = V2, ymin = V1, ymax = V3))+ geom_errorbar(width = 0.2, color = "#E15269") +
  geom_point(size = 1.5, color = "#E15269") + ylab("magnitude coefficient (SS)") + xlab("log(T)") + 
  theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))

# distance
distance_IC = rbind(IC_0[2,], IC_0.01[2,], IC_0.025[2,], IC_0.04[2,], IC_0.05[2,],  IC_0.07[2,], IC_0.1[2,], IC_0.15[2,],IC_0.2[2,], IC_0.25[2,], IC_0.3[2,], IC_0.35[2,],IC_0.4[2,], IC_0.45[2,],  IC_0.5[2,],IC_0.6[2,], IC_0.7[2,], IC_0.75[2,], IC_0.8[2,],IC_0.9[2,], IC_1[2,])
distance_IC = rbind(distance_IC,IC_1.2[2,], IC_1.4[2,], IC_1.6[2,], IC_1.8[2,], IC_2[2,],  IC_2.5[2,], IC_3[2,], IC_3.5[2,],IC_4[2,], IC_4.5[2,], IC_5[2,], IC_6[2,], IC_7[2,],  IC_8[2,], IC_9[2,], IC_10[2,])
distance_IC$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

distance.plot = ggplot(distance_IC, aes(x = n, y = V2, ymin = V1, ymax = V3))+ geom_errorbar(width = 0.2, color = "#E15269") +
  geom_point(size = 1.5, color = "#E15269") + ylab("distance coefficient (SS)") + xlab("log(T)") + 
  theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))

# vs30
vs30_IC = rbind(IC_0[3,], IC_0.01[3,], IC_0.025[3,], IC_0.04[3,], IC_0.05[3,],  IC_0.07[3,], IC_0.1[3,], IC_0.15[3,],IC_0.2[3,], IC_0.25[3,], IC_0.3[3,], IC_0.35[3,],IC_0.4[3,], IC_0.45[3,],  IC_0.5[3,],IC_0.6[3,], IC_0.7[3,], IC_0.75[3,], IC_0.8[3,],IC_0.9[3,], IC_1[3,])
vs30_IC = rbind(vs30_IC,IC_1.2[3,], IC_1.4[3,], IC_1.6[3,], IC_1.8[3,], IC_2[3,],  IC_2.5[3,], IC_3[3,], IC_3.5[3,],IC_4[3,], IC_4.5[3,], IC_5[3,], IC_6[3,], IC_7[3,],  IC_8[3,], IC_9[3,], IC_10[3,])
vs30_IC$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

vs30.plot = ggplot(vs30_IC, aes(x = n, y = V2, ymin = V1, ymax = V3))+ geom_errorbar(width = 0.2, color = "#E15269") +
  geom_point(size = 1.5, color = "#E15269") + ylab("vs30 coefficient (SS)") + xlab("log(T)") +
  theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))

x11()
par(bg = "#F4F4F4")
plot_grid(magnitude.plot, distance.plot, vs30.plot, nrow = 3)

# Residuals
Residuals = as.data.table(t(coeff[5,]))
Residuals$n = log10(c(0,0.01,0.025,0.04,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.75,0.8,0.9,1,1.2,1.4,1.6,1.8,2,2.5,3,3.5,4,4.5,5,6,7,8,9,10))

res.plot.SS = ggplot(Residuals, aes(x = n, y = V1))+ geom_point(size = 1.5, color = "#E15269") + ylab("residual std error (SS)") + xlab("log(T)") + theme_bw() + theme(plot.background = element_rect(fill = "#F4F4F4"))

x11()
par(bg = "#F4F4F4")
plot_grid(res.plot.NF, res.plot.TF, res.plot.SS, ncol = 3)


#-------------------------------------------------------------------
# 2. CROSS-VALIDATION FOR THE REGRESSION MODEL WITH STYLE OF 
#    FAULTING AS A FEATURE
#-------------------------------------------------------------------

# - Dummy variables for style of faulting 
tf_obs <- which(ITA_18_full_dataset$sof == "TF")
dummy_tf <- rep(0, dim(ITA_18_full_dataset)[1])
for (i in c(1:dim(ITA_18_full_dataset)[1])){
  for (j in c(1:length(tf_obs))){
    if (i == tf_obs[j]){
      dummy_tf[i] <- 1
    }
  }
}

ss_obs <- which(ITA_18_full_dataset$sof == "SS")
dummy_ss <- rep(0, dim(ITA_18_full_dataset)[1])
for (i in c(1:dim(ITA_18_full_dataset)[1])){
  for (j in c(1:length(ss_obs))){
    if (i == ss_obs[j]){
      dummy_ss[i] <- 1
    }
  }
}

detach(ITA_18_full_dataset)

# - Cross-validation on the new model with sof as variable  
library(caret)
numeric_values <- cbind(ITA_18_full_dataset[, c(7,8,10)], dummy_tf, dummy_ss)
SA_log = log(ITA_18_full_dataset[, 11:47])
data_cv <- cbind(numeric_values, SA_log )
attach(data_cv)


ctrl <- trainControl(method = "cv", number = 22)

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

periods = log10(c(0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3,
                  0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.2, 1.4,
                  1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10))

x11()
plot(periods,RMSE_sof[2:37], xlab = "periods", ylab = "RMS without sof", pch = 19)


# - Cross-validation on the model without sof as variable 
numeric_values <- ITA_18_full_dataset[, c(7,8,10)]
SA_log = log(ITA_18_full_dataset[, 11:47])
data_cv <- cbind(numeric_values, SA_log )
attach(data_cv)


ctrl <- trainControl(method = "cv", number = 22)

SA_0 <- train(SA_0 ~ magnitude + distance + vs30, data = data_cv[, c(1,2,3,4)], method = "lm", trControl = ctrl) 
print(SA_0)
detach(data_cv)

RMSE <- as.data.frame(SA_0$results$RMSE)

copie <- data_cv[,5:40]
colnames(copie) <- rep("SA", 36)

for (j in c(5:40)){
  data_j <- data_cv[ , 1:3]
  data_regression <- cbind(data_j, copie[,1])
  attach(data_regression)
  regression <- train( SA ~ magnitude + distance + vs30,data = data_regression, method = "lm", trControl = ctrl) 
  
  RMSE <- cbind(RMSE, regression$results$RMSE)
  copie[,1] <- NULL
  detach(data_regression)
}

colnames(RMSE) <- names(ITA_18_full_dataset[, c(11:47)])
plot(periods, RMSE[, 2:37], xlab = "periods", ylab = "RMSE without sof", pch = 19)

x11()
plot(log10(per), RMSE_sof[,1:26],xlab = "log(T)", ylab = "RMSE", pch = 19,col='#66121F', ylim = c(min(RMSE_sof[,1:26], RMSE[,1:26]),max(RMSE_sof[,1:26], RMSE[,1:26])), type = 'b')
points(log10(per), RMSE[,1:26], col = "#E15269", pch = 19, type='b')
legend("topright", legend = c("RMSE with SOF", "RMSE without SOF"), col = c("#66121F", "#E15269"), pch = 19)

mean(t(RMSE_sof[1,]))
mean(t(RMSE[1,]))
#so we can conclude that the error is greater for the model without style of faulting as a feature 
#By Cross-Validation we decide to consider the complete model.

#-------------------------------------------------------------------
# 3. MULTIPLE LINEAR REGRESSION CONSIDERING ONE-PERIOD-AT-THE-TIME
#-------------------------------------------------------------------

## LINEAR REGRESSION FOR EACH PERIOD WITH THE VARIABLE Style of Faulting AS A REGRESSOR ##
# log(SA_t) = b0 + b1*magnitude + b2*vs30 + b3*distance + b4*dummySS + 
#             + b5*dummyTF + (interactions between regressors) + Eps

numeric_values <- cbind(ITA_18_full_dataset[, c(7,8,10)], dummy_tf, dummy_ss)
SA_log = log(ITA_18_full_dataset[, 11:47])
data_cv <- data.frame(cbind(numeric_values, SA_log ) )
attach(data_cv)

# constructing one linear model for each period
min_tot=c()
max_tot=c()
mean_tot=c()
res_std_err_tot=c()
for(i in c(7:42)){
  # one linear model for each period
  model <- lm(data_cv[,i] ~ magnitude + distance+ vs30+
                dummy_ss + dummy_tf +
                vs30:dummy_ss + magnitude:dummy_ss + distance:dummy_ss +
                vs30:dummy_tf + magnitude:dummy_tf + distance:dummy_tf, data = data_cv[,c(1,2,3,4,5,i)])
  
  #plot(model) #this is for the analysis of the residuals: they are all fine (normality of residuals is
  #verified, as well as absence of influential cases)
  
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

# plot of the confidence intervals for the coefficients over the different periods

x11()
par(mfrow=c(3,4))
for(i in c(1:12)){
  plot(min_tot[,i], ylim =c(min(min_tot[,i]) , max(max_tot[,i]) ),col = 'red', ylab = '',xlab = '' ,
       main=paste("confidence interval for",labels(model$terms[i-1]), i) )
  points(max_tot[,i], col = 'red')
  points(mean_tot[,i] , col = 'blue')
  for(j in c(1:36)){
    segments(x0=j,y0=min_tot[j,i],x1=j,y1=max_tot[j,i])
  }
}

periods = log10(c(0, 0.01, 0.025, 0.04, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3,
                  0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.2, 1.4,
                  1.6, 1.8, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10))

x11()
par(bg = "white", mfrow=c(3, 4))
plot(log10(periods[-1]), min_tot[,1], ylim =c(min(min_tot[,1]) , max(max_tot[,1])), pch = '_', cex.lab = 1.2,col = "#1B1B19", ylab = 'CI intercept i = NF', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,1], col = '#1B1B19', pch='_')
points(log10(periods[-1]), mean_tot[,1] , col = '#1B1B19', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,1],x1=log10(periods[-1])[j],y1=max_tot[j,1])
}

plot(log10(periods[-1]), min_tot[,2], ylim =c(min(min_tot[,2]) , max(max_tot[,2])), pch = '_', cex.lab = 1.2, col = "#1B1B19", ylab = 'CI magnitude i = NF', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,2], col = '#1B1B19', pch='_')
points(log10(periods[-1]), mean_tot[,2] , col = '#1B1B19', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,2],x1=log10(periods[-1])[j],y1=max_tot[j,2])
}

plot(log10(periods[-1]), min_tot[,3], ylim =c(min(min_tot[,3]) , max(max_tot[,3])), pch = '_', cex.lab = 1.2, col = "#1B1B19", ylab = 'CI distance i = NF', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,3], col = '#1B1B19', pch='_')
points(log10(periods[-1]), mean_tot[,3] , col = '#1B1B19', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,3],x1=log10(periods[-1])[j],y1=max_tot[j,3])
}

plot(log10(periods[-1]), min_tot[,4], ylim =c(min(min_tot[,4]) , max(max_tot[,4])), pch = '_', cex.lab = 1.2, col = "#1B1B19", ylab = 'CI vs30 i = NF', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,4], col = '#1B1B19', pch='_')
points(log10(periods[-1]), mean_tot[,4] , col = '#1B1B19', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,4],x1=log10(periods[-1])[j],y1=max_tot[j,4])
}

par(bg = "#F4F4F4")
plot(log10(periods[-1]), min_tot[,6], ylim =c(min(min_tot[,6]) , max(max_tot[,6])), pch = '_', cex.lab = 1.2,col = "#66121F", ylab = 'CI intercept i = TF', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,6], col = '#66121F', pch='_')
points(log10(periods[-1]), mean_tot[,6] , col = '#66121F', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,6],x1=log10(periods[-1])[j],y1=max_tot[j,6])
}
par(bg = "white")
plot(log10(periods[-1]), min_tot[,11], ylim =c(min(min_tot[,11]) , max(max_tot[,11])), pch = '_', cex.lab = 1.2, col = "#66121F", ylab = 'CI magnitude i = TF', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,11], col = '#66121F', pch='_')
points(log10(periods[-1]), mean_tot[,11] , col = '#66121F', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,11],x1=log10(periods[-1])[j],y1=max_tot[j,11])
}

plot(log10(periods[-1]), min_tot[,12], ylim =c(min(min_tot[,12]) , max(max_tot[,12])), pch = '_', cex.lab = 1.2, col = "#66121F", ylab = 'CI distance i = TF', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,12], col = '#66121F', pch='_')
points(log10(periods[-1]), mean_tot[,12] , col = '#66121F', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,12],x1=log10(periods[-1])[j],y1=max_tot[j,12])
}

plot(log10(periods[-1]), min_tot[,10], ylim =c(min(min_tot[,10]) , max(max_tot[,10])), pch = '_', cex.lab = 1.2, col = "#66121F", ylab = 'CI vs30 i = NF', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,10], col = '#66121F', pch='_')
points(log10(periods[-1]), mean_tot[,10] , col = '#66121F', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,10],x1=log10(periods[-1])[j],y1=max_tot[j,10])
}

par(bg = "white")
plot(log10(periods[-1]), min_tot[,5], ylim =c(min(min_tot[,5]) , max(max_tot[,5])), pch = '_', cex.lab = 1.2, col = "#E15269", ylab = 'CI intercept i = SS', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,5], col = '#E15269', pch='_')
points(log10(periods[-1]), mean_tot[,5] , col = '#E15269', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,5],x1=log10(periods[-1])[j],y1=max_tot[j,5], col="#E15269")
}

par(bg = "#F4F4F4")
plot(log10(periods[-1]), min_tot[,8], ylim =c(min(min_tot[,8]) , max(max_tot[,8])), pch = '_', cex.lab = 1.2,col = "#E15269", ylab = 'CI magnitude i = SS', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,8], col = '#E15269', pch='_')
points(log10(periods[-1]), mean_tot[,8] , col = '#E15269', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,8],x1=log10(periods[-1])[j],y1=max_tot[j,8],col="#E15269")
}


plot(log10(periods[-1]), min_tot[,9], ylim =c(min(min_tot[,9]) , max(max_tot[,9])), pch = '_', cex.lab = 1.2, col = "#E15269", ylab = 'CI distance i = SS', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,9], col = '#E15269', pch='_')
points(log10(periods[-1]), mean_tot[,9] , col = '#E15269', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,9],x1=log10(periods[-1])[j],y1=max_tot[j,9],col="#E15269")
}

plot(log10(periods[-1]), min_tot[,7], ylim =c(min(min_tot[,7]) , max(max_tot[,7])), pch = '_', cex.lab = 1.2, col = "#E15269", ylab = 'CI vs30 i = SS', xlab = 'log(T)')
points(log10(periods[-1]), max_tot[,7], col = '#E15269', pch='_')
points(log10(periods[-1]), mean_tot[,7] , col = '#E15269', pch = 19)
for(j in c(1:36)){
  segments(x0=log10(periods[-1])[j],y0=min_tot[j,7],x1=log10(periods[-1])[j],y1=max_tot[j,7], col = '#E15269' )
}

graphics.off()


## LINEAR REGRESSION FOR EACH PERIOD WITHOUT THE VARIABLE SOF ##
# log(SA_t) = b0 + b1*magnitude + b2*vs30 + b3*distance + Eps
max_tot=c()
min_tot=c()
mean_tot=c()
res_std_err_tot_nosof=c()
for(i in c(6:36)){
  # one linear model for each period
  model <- lm(data_cv[,i] ~ magnitude + distance + vs30, data = data_cv[,c(1,2,3,4,5,i)])
  
  #plot(model) #this is for the analysis of the residuals: they are all fine (normality of residuals is
  #verified, as well as absence of influential cases)
  
  confint(model)
  min=confint(model)[,1]
  min_tot=rbind(min_tot,t(min))
  max=confint(model)[,2]
  max_tot=rbind(max_tot,t(max))
  mean=(max+min)/2
  mean_tot=rbind(mean_tot,t(mean))
  
  res_std_err_nosof=sqrt(deviance(model)/model$df.residual)
  res_std_err_tot_nosof=rbind(res_std_err_tot_nosof, t(res_std_err_nosof))
}

# plot of the confidence intervals for the coefficients over the different periods
x11()
par(mfrow=c(2,2))
for(i in c(2:4)){
  plot(min_tot[,i], ylim =c(min(min_tot[,i]) , max(max_tot[,i]) ),col = 'red', ylab = '',xlab = '' ,
       main=paste("confidence interval for",labels(model$terms[i-1])) )
  points(max_tot[,i], col = 'red')
  points(mean_tot[,i] , col = 'blue')
  for(j in c(1:31)){
    segments(x0=j,y0=min_tot[j,i],x1=j,y1=max_tot[j,i])
  }
}


x11()
par(bg = "#F4F4F4")
plot(res_std_err_tot, col = '#66121F', ylim = c(min(res_std_err_tot),max(res_std_err_tot_nosof)),
     ylab='RMSE', xlab='T', main= "residual standard error", type = 'b', pch=19)
points(res_std_err_tot_nosof, col='#E15269', type = 'b', pch=19)
legend(x="topright", legend = c("model with sof","model without sof"),
       fill= c('#66121F','#E15269') )


detach(data_cv)


