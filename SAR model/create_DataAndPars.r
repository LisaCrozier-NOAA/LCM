library(TMB)
library(mvtnorm)
library(matrixcalc)
library(corpcor)

load("pitTotal_Snake_condensed.Rdata") #Load the pit tag data.

#Subset the data based on the input variables
pitTotal$year <- 
  as.numeric(as.character(pitTotal$year))
pitTotal <- pitTotal[pitTotal$year<=lastYr,]

if(!is.null(myTrans)){
  pitTotal <-
    pitTotal[pitTotal$trans_dam==myTrans,]
}
if(!is.null(rear)){
  pitTotal <-
    pitTotal[pitTotal$rear_type==rear,]
}
k_name <- paste(pitTotal$trans_dam,pitTotal$rear_type)
pitTotal$k <- as.numeric(as.factor(k_name))-1
pitTotal$yr <- pitTotal$year - eStartYr

#Load the data, subset it and transform is based on the input in the wrapper_modelRuns.r
load("envData.rData")
tmpenv <- envdata[envdata$year>=eStartYr & envdata$year<=eLastYr,]
if(eLastYr>max(envdata$year)){
  projEnv <- as.data.frame(matrix(NA,eLastYr-max(envdata$year),ncol(tmpenv)))
  projEnv[,1] <- (max(envdata$year)+1):eLastYr
  names(projEnv) <- names(tmpenv)
  tmpenv <- rbind(tmpenv,projEnv)
}
names(tmpenv)[1] <- "Year"

# Jinky step to get the most recent coastal summer data
load(paste0("envData_7302018.Rdata"))
tmpenv$ersstWAcoast.sum[tmpenv$year%in%envdata$Year] <- envdata$ersstWAcoast.sum[tmpenv$year%in%envdata$Year]

#Get the scaling data
env_mu <- apply(tmpenv,2,function(x){return(mean(na.omit(x)))})
env_sc <- apply(tmpenv,2,function(x){return(sd(na.omit(x)))})
# env_sc <- rep(1,ncol(tmpenv))
env_mu_2000_2015 <- apply(tmpenv[tmpenv$Year>=2000 & tmpenv$Year<=2015,],2,function(x){mean(na.omit(x))})
env_sc_2000_2015 <- apply(tmpenv[tmpenv$Year>=2000 & tmpenv$Year<=2015,],2,function(x){sd(na.omit(x))})
env_mu_1980_2015 <- apply(tmpenv[tmpenv$Year>=1980 & tmpenv$Year<=2015,],2,function(x){mean(na.omit(x))})
env_sc_1980_2015 <- apply(tmpenv[tmpenv$Year>=1980 & tmpenv$Year<=2015,],2,function(x){sd(na.omit(x))})
Stationary_baseline_difference <- env_mu_1980_2015[myVars] - env_mu_2000_2015[myVars]

#Z-score variables and years
tmpenv[,names(tmpenv)!="Year"] <- t((t(tmpenv[,names(tmpenv)!="Year"]) - env_mu[names(tmpenv)!="Year"])/env_sc[names(tmpenv)!="Year"])

#Index years starting with zero
tmpenv$Year <- tmpenv$Year - eStartYr

#Get rid of NAs
subData <- tmpenv[,myVars]
subData[is.na(subData)] <- -10000


#Create increment limit
nte = nrow(subData)
nvar = ncol(subData)

k_name <- k_name[pitTotal$julian>=minJ & pitTotal$julian<=maxJ]
pitTotal <- pitTotal[pitTotal$julian>=minJ & pitTotal$julian<=maxJ,]

#marine variables for pit analysis using only 2000-2015 mean and standard deviation
xm <- t((t(pitTotal[,marVars])-env_mu_2000_2015[marVars])/
    env_sc_2000_2015[marVars])

load("gridData.rData")
#Timing across days and years
retroTiming <- apply(gridData["inriver",,,,],c(2,1),sum)
retroTiming <- t(retroTiming)/colSums(retroTiming)
retroTiming <- t(retroTiming[gridYears%in%eStartYr:(startYr-1),1:365%in%minJ:maxJ])

#This is the calibration timing for Lisa's calibration from ~1998:2015
calibration_timing <- NA
calibration_years <- NA
if(calibration_flag==1){
  calibration_timing <- read.table(calibration_file, header=TRUE)
  if(myTrans=="In-river"){
    calibration_timing <- calibration_timing[calibration_timing$MigrationType=="inriver",]
  }
  if(myTrans=="Transport"){
    calibration_timing <- calibration_timing[calibration_timing$MigrationType=="transportation",]
  }
  calibration_years <- calibration_timing$Year[calibration_timing$Year<=lastYr]-eStartYr
  
  calibration_timing <- calibration_timing[calibration_timing$Year<=lastYr,minJ:maxJ+3]
  calibration_timing <- calibration_timing/rowSums(calibration_timing)
}

#This is for projecting into the future.
obsTiming <- matrix(0,length(minJ:maxJ),length(startYr:lastYr))
for(j in minJ:maxJ){
  for(yy in startYr:lastYr){
    obsTiming[j-minJ+1,yy-startYr+1] <- sum(pitTotal$Total[pitTotal$year==yy & pitTotal$julian==j])
  }
}
obsTiming <- t(t(obsTiming)/colSums(obsTiming))
timing <- cbind(retroTiming,obsTiming)
projTiming <- matrix(rowSums(timing)/sum(timing),length(minJ:maxJ),length((lastYr+1):eLastYr))
timing <- cbind(timing,projTiming)

#You need and offest for years if you are doing the projection
if(retro)
  yShift <- 0
if(!retro)
  yShift <- startYr-eStartYr

#Data list
data <- list( yShift = yShift
             ,yr = pitTotal$yr
             ,j = pitTotal$julian-minJ
             ,s_n = pitTotal$Total
             ,s_k = pitTotal$resMatrix[,1] #pitTotal.txt
             ,k = pitTotal$k
             ,nvar = nvar
             ,env = as.matrix(subData)
             ,timing = timing 
             ,xm = as.matrix(xm)
             ,marVars = tmpMarVars-1
             ,sd=0.0001
             ,re_j = re_j
             ,re_t = re_t
             ,re_jt = re_jt
             ,cov_pars = cov_pars
             ,retro = retro
             ,env_mu = as.vector(env_mu[myVars])
             ,env_sc = as.vector(env_sc[myVars])
             ,env_mu_2000_2015 = as.vector(env_mu_2000_2015[myVars])
             ,env_sc_2000_2015 = as.vector(env_sc_2000_2015[myVars])
             ,calibration_flag = calibration_flag
             ,calibration_timing = as.matrix(calibration_timing)
             ,calibration_years = calibration_years
             ,Stationary_baseline_difference=Stationary_baseline_difference
)

# env_mu <- apply(tmpenv,2,function(x){return(mean(na.omit(x)))})
# env_sc <- apply(tmpenv,2,function(x){return(sd(na.omit(x)))})
# # env_sc <- rep(1,ncol(tmpenv))
# env_mu_2000_2015 <- apply(tmpenv[tmpenv$Year>=2000 & tmpenv$Year<=2015,],2,function(x){mean(na.omit(x))})
# env_sc_2000_2015 <- apply(tmpenv[tmpenv$Year>=2000 & tmpenv$Year<=2015,],2,function(x){sd(na.omit(x))})

nk_dim <- max(data$k)+1
#Parameter list
parameters <- list(mu_s  = rep(0,nk_dim)
                   ,frho_j = rep(0,nk_dim)
                   ,frho_t = rep(0,nk_dim)
                   ,frho1_jt = rep(0,nk_dim)
                   ,frho2_jt = rep(0,nk_dim)
                   ,fpsi_j = rep(0,nk_dim)
                   ,fpsi_t = rep(0,nk_dim)
                   ,fpsi_jt = rep(0,nk_dim)
                   ,eps_j = matrix(0,length(minJ:maxJ),nk_dim)
                   ,eps_t = matrix(0,length(eStartYr:eLastYr)-yShift,nk_dim)
                   ,eps_jt = array(0,c(nk_dim,length(minJ:maxJ),length(eStartYr:eLastYr)-yShift))
                   ,beta_mar = matrix(0,length(tmpMarVars),nk_dim)
                   ,frho_x = rep(0)
                   ,frho_Rx = rep(0,nvar*(nvar-1)/2)
                   ,fpsi_x = rep(0,nvar)
                   ,eps_x = matrix(0,nvar,length(eStartYr:eLastYr))
)
