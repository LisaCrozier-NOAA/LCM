#Call this file from wrapper_modelRuns_sameEnv_trans_inriver.r
#This is the wrapper for running climate change projections for 
# Crozier, Burke, Chasco, Widener and Zabel 2020
# Chinook salmon face climate change across multiple life stages putting iconic species in peril
# published in Communications Biology
####################################################
####################################################

library(mvtnorm)

set.seed(1000)

runSensitivity <- FALSE

load(files[1]) #It doesn't matter which file you load. They will all have the same environmental projection information
nvar <- length(myVars)

#These are the flow grid bins
flowSeq <- c(20,seq(50,350,50))
#Mapping to grid array
flowSeqID <- 1:length(flowSeq)
#Temperature bins
tempSeq <- 9:18
#Mapping to grid array
tempSeqID <- 1:length(tempSeq)

#Set up array for gridded data of flow and temp
gridData <- array(0,dim=c(2,length(1:80),365,length(flowSeq),length(tempSeq)),
                  dimnames=list(migrationType = c("inriver","transport"),
                                compassYears=1:80,
                                julian=1:365,
                                flow=flowSeq,
                                temp=tempSeq))

#Read in the gridded compass data
for(ff in flowSeq){
  for(tt in tempSeq){
    # ff <- flowSeq[1]
    # tt <- tempSeq[1]
if(NAA==TRUE)    cDat <- read.table(paste0(rootDir,"/DATA_COMPASS_FILES_NAA/ch1_BiOp-NAA-Flow-",ff,"_",tt,"deg_bonarrival_transport.out"), header = TRUE)
if(NAA==FALSE)   cDat <- read.table(paste0(rootDir,"/DATA_COMPASS_FILES_PA/ch1_Flow-",ff,"_",tt,"deg_bonarrival_transport.out"), header = TRUE)
    mcnt <- 1
    for(mm in c("inriver","transport")){
      # mm <- "inriver"
      tmpcDat <- cDat[cDat$MigrationType==mm,4:368]
      gridData[mcnt,,,flowSeqID[flowSeq==ff],tempSeqID[tempSeq==tt]] <- as.matrix(tmpcDat)
      mcnt <- mcnt + 1
    }
  }
}


#Marine and freshwater information 
rcps <- matrix(c(c("stationary","2.6","4.5","8.5"), #sstARC
                 c("stationary","2.6","4.5","8.5"), #sstWA
                 c("stationary","","",""), #upw
                 c("stationary","2.6","4.5","8.5"), #T20
                 c("stationary","2.6","4.5","8.5")), #F
               5,4,byrow=TRUE)
               
quants <- c("25","50","75")

#The RCP files require a lookup table to match the projections with the correct marine variable
RCP_lu <- read.csv(paste0(rootDir,"/RCP/RCP_lookup.csv"), header = TRUE)
rcp_varNames_pre <- RCP_lu[1,]
rcp_varNames_suf <- RCP_lu[2,]

#RCP data holder
RCP_rqvy <- array(0,dim=c(ncol(rcps),length(quants),nvar,ny),
             dimnames=list(RCP=rcps[1,],quants=quants,vars=dimnames(data$env)[[2]],y=1:ny))

#RCP directory
rcpDir <- paste0(rootDir,"/RCP/")

#Marine variable inputs from individual files
for(i in marVars){ #marine variables in the var map
  rcnt <- 2 #Need this since rcp is a character
  for(r in c("2.6","4.5","8.5")){
    qcnt <- 1 #Used to map the quantile (25,50,75) to the array index
    for(q in quants){
      tmp <- read.csv(paste0(rcpDir,
                  rcp_varNames_pre[1,i],
                  ".",
                  q,
                  ".delta",
                  r,
                  ".",
                  rcp_varNames_suf[1,i],
                  ".csv"),
           header=TRUE)
      
      tmp <- tmp[tmp$year>=2015 & tmp$year<=(2015+74),2]
      
      RCP_rqvy[r,q,i,] <- tmp
    }
  }
}


#Freshwater variable inputs. This is a different file and rcp structure compared to the marine variables
fw <- read.csv(paste0(rcpDir,"LGR_quantiles06242019.csv"),header=TRUE)

#For the marine variables, we read in the lookup table
#For the freshwater variables it's easier to just write it here.
rcp_FWvarNames_pre <- as.data.frame(matrix(c("T20.","F."),1,2))
rcp_FWvarNames_suf <- as.data.frame(matrix(c("",""),1,2))
names(rcp_FWvarNames_pre) <- c('aprmayjunetempLGR','aprmayjuneflowLGR')
names(rcp_FWvarNames_suf) <- c('aprmayjunetempLGR','aprmayjuneflowLGR')

for(i in c('aprmayjunetempLGR','aprmayjuneflowLGR')){
  for(r in c("2.6","4.5","8.5")){
    for(q in quants){
      myCol <- paste0(rcp_FWvarNames_pre[1,i],
                             q,
                             ".delta",
                             r,
                             rcp_FWvarNames_suf[1,i])
      tmp <- fw[fw$year>=2015 & fw$year<=(2015+74),names(fw)==myCol]
      #Get the RCP values
      RCP_rqvy[r,q,i,] <- tmp
    }
  }
}


#Environmental projections
#Replace the simulated values for the first year with the 2000:2015 mean values that are observed
#These environmental projections include all years.
envProj_rqvys <- array(0,dim=c(ncol(rcps),length(quants),nvar+3,ny,nsim),
                  dimnames = list(RCP=rcps[1,], quants=quants,vars=c(dimnames(data$env)[[2]],"gridYear","gridFlow","gridTemp"),years=1:ny,sim=1:nsim))
#unscaled env projections
sc_envProj_rqvys <- array(0,dim=c(ncol(rcps),length(quants),nvar+3,ny,nsim),
                       dimnames = list(RCP=rcps[1,], quants=quants,vars=c(dimnames(data$env)[[2]],"gridYear","gridFlow","gridTemp"),years=1:ny,sim=1:nsim))
#unscaled env projections with added RCP
sc_RCP_envProj_rqvys <- array(0,dim=c(ncol(rcps),length(quants),nvar+3,ny,nsim),
                       dimnames = list(RCP=rcps[1,], quants=quants,vars=c(dimnames(data$env)[[2]],"gridYear","gridFlow","gridTemp"),years=1:ny,sim=1:nsim))

#random day effects
rday_frqjys <- array(0,dim=c(length(files),ncol(rcps),length(quants),length(minJ:maxJ),ny,nsim),
                              dimnames = list(file=files,RCP=rcps[1,], quants=quants,day=minJ:maxJ,1:ny,sim=1:nsim))

#random dayXyear effects
rdayXyear_frqjys <- array(0,dim=c(length(files),ncol(rcps),length(quants),length(minJ:maxJ),ny,nsim),
               dimnames = list(file=files,RCP=rcps[1,], quants=quants,day=minJ:maxJ,1:ny,sim=1:nsim))

#random mu and betas
r_muBeta_frqjys <- array(0,dim=c(length(files),ncol(rcps),length(quants),1+length(tmpMarVars),ny,nsim),
                          dimnames = list(file=files,RCP=rcps[1,], quants=quants,fixed=c("mu",paste("beta",1:length(tmpMarVars))),1:ny,sim=1:nsim))

#covariance matrix
Sigma <- rep$`Sigma_x.cov()`#*(rep$psi_x%o%rep$psi_x)

#initEnv data from 2000 to 2015
initEnv <- colMeans(data$env[2000:2015-startYr+1,])

for(rr in simRCP){
  for(qq in simQuants){
    for(ss in 1:nsim){
      for(y in 1:ny){
        #Initial year
        if(y==1){
          envProj_rqvys[rr,qq,dimnames(data$env)[[2]],1,ss] <- 
            rmvnorm(1,rep(0,nvar),Sigma)
        }
        #Project the marine and freshwater env data. This is in the Z-score space
        if(y>1){
          envProj_rqvys[rr,qq,dimnames(data$env)[[2]],y,ss] <- 
            rmvnorm(1,rep$rho_x*envProj_rqvys[rr,qq,dimnames(data$env)[[2]],y-1,ss],Sigma*sqrt(1-rep$rho_x^2))   #SCALE by correlation coefficient
        }
  
        # #Now we add the climate trends to these simulations
        flow <- envProj_rqvys[rr,qq,"aprmayjuneflowLGR",y,ss] * #Projection
                  data$env_sc[colnames(data$env)=="aprmayjuneflowLGR"] + #Scale
                  data$env_mu[colnames(data$env)=="aprmayjuneflowLGR"] + #mean
                  RCP_rqvy[rr,qq,"aprmayjuneflowLGR",y] #rcp
        
        temp <- envProj_rqvys[rr,qq,"aprmayjunetempLGR",y,ss] * #Projection
                      data$env_sc[colnames(data$env)=="aprmayjunetempLGR"] + #Scale
                      data$env_mu[colnames(data$env)=="aprmayjunetempLGR"] + #mean
                      RCP_rqvy[rr,qq,"aprmayjunetempLGR",y] #rcp
        #             
        # 
        #Get the flow and temp that are closest to the compass files
        gridFlow <- c(20,seq(50,350,50))[abs(flow -  c(20,seq(50,350,50)))==min(abs(flow -  c(20,seq(50,350,50))))]
        gridTemp <- (9:18)[abs(temp -  9:18)==min(abs(temp -  9:18))]
        
        #This is a random draw of 80 years within the flow and temp compass simulation
        gridYear <- round(runif(1,1,80))

        #Save the flow, temp and grid year
        envProj_rqvys[rr,qq,"gridYear",y,ss] <- gridYear
        envProj_rqvys[rr,qq,"gridFlow",y,ss] <- gridFlow
        envProj_rqvys[rr,qq,"gridTemp",y,ss] <- gridTemp
      }
      
      # #Now unZ-score everything and add the RCP
      sc_envProj_rqvys[rr,qq,myVars,,ss] <-
        envProj_rqvys[rr,qq,myVars,,ss] * rep$psi_x * #Scale by the process standard deviation
        data$env_sc[colnames(data$env)==myVars] +
        data$env_mu[colnames(data$env)==myVars] +
        RCP_rqvy[rr,qq,myVars,]

      #The re-scale by the 2000-2015 data       
      sc_RCP_envProj_rqvys[rr,qq,myVars,,ss] <-
        ((sc_envProj_rqvys[rr,qq,myVars,,ss] + data$Stationary_baseline_difference * stationary_baseline_flag) - data$env_mu_2000_2015[colnames(data$env)==myVars])/
        data$env_sc_2000_2015[colnames(data$env)==myVars]
      
    }
  }
}

#Projected SAR output
sarProj_frqys <- array(0,
                  dim=c(length(files),ncol(rcps),length(quants),ny,nsim),
                  dimnames = list(file=files,rcps=rcps[1,],quants=quants,ny=1:ny,nsim=1:nsim))
sarBeta_frqvs <- array(NA,
                       dim=c(length(files),ncol(rcps),length(quants),length(marVars),nsim),
                       dimnames = list(file=files,rcps=rcps[1,],quants=quants,beta_mar=marVars,nsim=1:nsim))

# #Projected betas
for(f in files){
#   #Load the model
  load(f)
  #Number of env variables
  nvar <- ncol(data$env)
  for(rr in simRCP[1]){
    for(qq in simQuants[1]){
      #Get migration type
      ifelse(substr(f,25,25)=="T",
             migType <- "transport",
             migType <- "inriver")

      #project the environmental data
      for(ss in 1:nsim){
        print(paste(f))
        print(paste(rr,qq,ss))
        #If you are doing a sensitivity analysis, then you just use the same simulation for each model run
        sSens <- ss
        if(runSensitivity==TRUE){
          sSens <- SensitivitySimulation
        }

        #Get the env variables that are to be multiplied by the estimated beta parameters
        if(length(data$marVars)==1){
          tmp_env <- cbind(rep(1,ny),
                         sc_RCP_envProj_rqvys[rr,qq,marVars[tmpMarVars],,sSens])
        }
        if(length(data$marVars)>1){
          tmp_env <- cbind(rep(1,ny),
                           t(sc_RCP_envProj_rqvys[rr,qq,marVars[tmpMarVars],,sSens]))
        }


        #Names of the SD values
        sdNames <- names(SD$value)

        #Loop through all years and calculate the SAR
        for(y in 1:ny){
          #You need to find the flow, temperature and year that is most closely associated with the TMB simulation value
          temp <- envProj_rqvys[rr,qq,"gridTemp",y,sSens]
          temp <- tempSeqID[temp==tempSeq]
          flow <- envProj_rqvys[rr,qq,"gridFlow",y,sSens]
          flow <- flowSeqID[flow==flowSeq]
          year <- envProj_rqvys[rr,qq,"gridYear",y,sSens] #Grid year
          ifelse(migType=="inriver",mig <- 1, mig <- 2)

          #This is the timing based on the COMPASS model
          timing <- gridData[mig,y,,flow,temp]

          #random year
          ry <- round(runif(1,1,length(startYr:lastYr)))

          #Position of the covariates not including the dayXyear interaction
          sd_pos <- (1:length(names(SD$value)))[sdNames=="mu_s" |
                                                          sdNames=="beta_mar"]
          sd_pos_fixed <- (1:length(names(SD$value)))[sdNames=="mu_s" |
                                                  sdNames=="beta_mar"]
          if(data$re_j==1){
            sd_pos <- c(sd_pos,(1:length(names(SD$value)))[sdNames=="eps_j"])
            sd_pos_j <- (1:length(names(SD$value)))[sdNames=="eps_j"]
          }
          if(data$re_t==1){
            sd_pos <- c(sd_pos,(1:length(names(SD$value)))[sdNames=="eps_t"])
            sd_pos_t <- (1:length(names(SD$value)))[sdNames=="eps_t"]
          }
          if(data$re_jt==1){
            sd_pos <- c(sd_pos,(1:length(names(SD$value)))[sdNames=="eps_jt"][1:length(minJ:maxJ) + (ry-1)*length(minJ:maxJ)])
            sd_pos_jt <- (1:length(names(SD$value)))[sdNames=="eps_jt"][1:length(minJ:maxJ) + (ry-1)*length(minJ:maxJ)]
          }

          #Covariance matrix for all fixed and random effects
          COV <- SD$cov[sd_pos,
                        sd_pos]

          #Draw the random mu, beta, day, year, and dayXyear effect for a particular ry (random year)
          #Random draw
          # rd <- rmvnorm(1,
          #               mean=SD$value[sd_pos],
          #               sigma = COV)
          # rd <- as.data.frame(rd)

          #This scenario has no variability at all. Just the MLE estimate. This goes really, really fast
          rd_fix <- SD$value[sd_pos_fixed]
                        #rmvnorm(1,
                        #mean=SD$value[sd_pos_fixed],
                        #sigma = SD$cov[sd_pos_fixed,sd_pos_fixed])

          #This scenario assumes no covariance between parameters. This goes really fast
          # rd_fix <- rnorm(length(SD$value[sd_pos_fixed]),
          #     SD$value[sd_pos_fixed],
          #     diag(SD$cov[sd_pos_fixed,sd_pos_fixed]))

          #This scenario assumes full covariance between parameters. This goes really slow
          rd_fix <- rmvnorm(1,
                          SD$value[sd_pos_fixed],
                          SD$cov[sd_pos_fixed,sd_pos_fixed])


          rd_j <- rmvnorm(1,
                        mean=SD$value[sd_pos_j],
                        sigma = SD$cov[sd_pos_j,sd_pos_j]) #SD$value[sd_pos_j] #
          rd_jt <- rmvnorm(1,
                          mean=SD$value[sd_pos_jt],
                          sigma = SD$cov[sd_pos_jt,sd_pos_jt]) #SD$value[sd_pos_jt]#
          
          #Mean and marine covariates
          # rBeta <- rd_fix[names(rd_fix)=="beta_mar"]
          # rMu <- rd_fix[names(rd_fix)=="mu_s"]
          rBeta <- rd_fix[2:length(rd_fix)]
          rMu <- rd_fix[1]

          #Estimate the marine effect
          mar_ef <- tmp_env%*%as.matrix(unlist(c(rMu,rBeta)))

          #Save the random draws
          if(length(rBeta)==1){
            r_muBeta_frqjys[f,rr,qq,,y,ss] <- c(unlist(c(rMu,rBeta)),NA)
          }
          if(length(rBeta)==2){
            r_muBeta_frqjys[f,rr,qq,,y,ss] <- unlist(c(rMu,rBeta))
          }
          #*****************************
          #Random effects for day, year, and day X year
          #*****************************
          day_ef <- t(as.matrix(rep(0,length(minJ:maxJ))))
          if(data$re_j==1){
            # day_ef <- rd[names(rd)=="eps_j"] * exp(rep$fpsi_j)
            day_ef <- t(rd_j) * exp(rep$fpsi_j)
          }
          rday_frqjys[f,rr,qq,,y,ss] <- t(day_ef)

          #Random year effect
          year_ef <- t(as.matrix(rep(0,length(startYr:lastYr))))
          if(data$re_t==1){
            year_ef <- rd[names(rd)=="eps_t"]  * exp(rep$fpsi_t)
          }
          #Random day X year effect
          dayXyear_ef <- t(as.matrix(rep(0,length(minJ:maxJ))))
          if(data$re_jt==1){
            # dayXyear_ef <- rd[names(rd)=="eps_jt"]  * exp(rep$fpsi_jt)
            dayXyear_ef <- t(rd_jt)  * exp(rep$fpsi_jt)
          }
          rdayXyear_frqjys[f,rr,qq,,y,ss] <- t(dayXyear_ef)


          #The total effects of all fixed and random effects
          total_ef <- c(rep(mar_ef[y,1],length(1:(minJ-1))), #days before minJ
            (t(day_ef) + t(dayXyear_ef) + year_ef[1,ry] + mar_ef[y,1]), #minJ:maxJ
            rep(mar_ef[y,1],length((maxJ+1):365))) #after maxJ

          # total_ef <- c(rep(mar_ef[y,1],length(1:(minJ-1))),
          #               (t(day_ef) + t(dayXyear_ef) + year_ef[ry] + mar_ef[y,1]),
          #               rep(mar_ef[y,1],length((maxJ+1):365)))

          #SAR is the total effect
          sarProj_frqys[f,rr,qq,y,ss] <- sum(total_ef * timing/sum(timing))
        }
        #Because the number of betas will vary, you need to create a list of arrays
        sarBeta_frqvs[f,rr,qq,data$marVars+1,ss] <- as.matrix(rBeta)
      }
    }
  }
}


    filename<-paste0("sar_projections_",
      if(NAA==TRUE) "NAA_",
      if(NAA==FALSE) "PA_",
                     myTrans.tt[1],"_",myTrans.tt[2],"_","Nsim",nsim,"_Nmodel",length(files),"scenarios",length(simRCP)*length(simQuants),
                     ".rData")

save(sarProj_frqys
     ,sarBeta_frqvs
     ,envProj_rqvys
    ,sc_RCP_envProj_rqvys
     ,rday_frqjys
     ,rdayXyear_frqjys
     ,r_muBeta_frqjys
     ,gridData
     ,env_mu
     ,env_sc
     ,env_mu_2000_2015
     ,env_sc_2000_2015
     ,RCP_rqvy
     ,subData
     ,file=filename
  )
