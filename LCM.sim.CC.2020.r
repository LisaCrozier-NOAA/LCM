



####################################################

#This is the wrapper for running the LCM model for 
# Crozier, Burke, Chasco, Widener and Zabel 2020
# Iconic salmon populations face perilous challenges from climate change across their life cycle
# published in Communications Biology

####################################################
####################################################


set.seed(10000)
library(plyr)

# Model control parameters---------
#1.Select which model you want to run. Models are defined by which freshwater and marine covariates they include
      #e.g., Model 1:
          # FW: p2TsuFfall
          # In-river: ersstArc.win + ersstWAcoast.sum
          # Transport: ersstArc.win 

        model<-"Model1";filename<-"p2TsuFfall"; Fname="F.S3.ParrF"; FW.name="Ffall"
          # model<-"Model2"; filename<-"p2Fsu"; Fname="F.S2.ParrSu"; FW.name="Fsu"
          # model<-"Model3"; filename<-"p2TsuFfall"; Fname="F.S3.ParrF"; FW.name="Ffall"
                arrayname<-paste0(model,".array")
              load(file=paste("Input files/param",arrayname,sep="."),verbose=T) #param.array.1K

                 
#2. Decide how many simulations,populations and climate scenarios you want to run and build list of scenarios
        runs = 10
        years = 75 #2015:2089
        pop<-creeks<-"bearvalley"
#       pop<-creeks<-c("bearvalley","big","camas","loon","marsh","secesh","sulphur","valley")

            climate =c("stationary","4.5" , "8.5"   )
            climateshort =c("1","45" , "85"   )
            quant = c("25", "50", "75")

        nsim<-length(pop)*length(model)*length(climate)*length(quant)#;nsim
        scenario<-length(model)*length(climate)*length(quant);scenario
        scenario.name<-vector("character",length=scenario)
     cnt.scenario=0
        for(l in 1:length(model)){
            for(cc in 1:length(climate)){ 
            for(qq in 1:length(quant)){ 
                 cnt.scenario=cnt.scenario+1
               scenario.name[cnt.scenario]=paste(climate[cc],quant[qq],model[l],sep=".");scenario.name
            }}}
     

#Input files-------------------
    #FW
        source("LCM.fxn.r") #functions for model relationships
        popinit<-read.csv("Input files/Ninit.csv", header = T, row.names = 1)
        hydro<-read.csv(file="Input files/climategrid_survival_80yr_2020_PA.csv",header=TRUE,row.names = NULL);hydro$gridYear<-rep(1:80,length=nrow(hydro))
        FW<-read.csv(file="Input files/SALM.Qflow.csv",header=TRUE,row.names = NULL)#;head(FW)
            FW<-FW[FW$year>=2015,]
            scalefactors<-read.csv(file="Input files/scalefactors.corrected.May2019.csv",row.names=1)#;scalefactors
        su<-read.csv("Input files/AprilJuneAndSurv_LoopsTransport.csv",row.names = NULL,header=TRUE)#;head(su)
            su$LGRtemp<-round_any(su$LGrTempAJ,0.5)#;simtemp
            su$LGRflow<-round_any(su$LGrFlowAJ,20)#;simflow

        
    #SARs    
        if(model=="Model3") sarfile<-"Model3.SAR100.Rdata"  else sarfile<-"Model1.Model2.SAR100.Rdata"
            print(sarfile)
               load(paste0("Input files/",sarfile),verbose=TRUE) 
 
      #collect parameters used in each run
            param.marT<-as.data.frame(t(sarcoefT_rqvs["stationary","50",,]))
                  names(param.marT)<-paste("T",names(param.marT),sep=".")
            param.marI<-as.data.frame(t(sarcoefI_rqvs["stationary","50",,])) 
                  names(param.marI)<-paste("I",names(param.marI),sep=".")

# Set up temporary arrays to store data during model runs-----------
    N.new = N = array(0,5)
    surv.s1=surv.trib=surv.inriver=pt=pt.sar=pt.sup=surv.main=surv.s2=surv.s3=surv.sar=surv.sarI=surv.sarT=surv.sup=sup.spr.inriver=sup.spr.transport=sup.sum.inriver=sup.sum.transport=sup.spring=sup.summer=Spawners = Recruits = array(0,years) #indiv years within a run

    # Arrays to save output-----------
      spawner.array<-parrsmolt.array<-recruit.array<-array(0,dim=c(years,runs,length(pop),scenario),dimnames=list(2015:2089,paste0("sim",1:runs),pop,scenario.name))#;dim(spawner.array)
      response<-c("yrQET50","Mean2020","Mean2040","Mean2060","Mean2080","Meanallyears")
            allparam<-c( names(param.marI),names(param.marT),dimnames(param.array.1K)[[2]],response)# ;allparam
      param.array<-array(0,dim=c(runs,length(allparam),length(pop),scenario),dimnames=list(paste0("sim",1:runs),allparam,pop,scenario.name))#;dim(param.array)
    


#START LOOPING THROUGH SIMULATIONS ========================
    print(Sys.time())     
#k=1;l=1;m=1;cc=1;p=1;j=1;qq=1
cnt=0
for(k in 1:length(pop)){ # loop across populations
    cnt.scenario=0
  print(paste("Pop",k,pop[k],Sys.time()))
  ha= popinit[pop[k],"ha"]
  for(cc in 1:length(climate)){ 
  for(qq in 1:length(quant)){ #loop across quantiles of each climate scenario
     cnt.scenario=cnt.scenario+1;print(scenario.name[cnt.scenario])
  for (j in 1:runs){
   jj<-sample(length(param.array.1K[,1,1]),1)
   jj.sar<-sample(length(sarI_rqys[1,1,1,]),1)
        param<-as.list(param.array.1K[jj,,pop[k]])#;unlist(param)

        #collect parameters used in each run
            param.marT<-as.data.frame(t(sarcoefT_rqvs[climate[cc],quant[qq],,jj.sar]))  #;unlist(param.marT) #dimnames(sarcoefT_rqvs)
            param.marI<-as.data.frame(t(sarcoefI_rqvs[climate[cc],quant[qq],,jj.sar])) #dimnames(sarBeta_frqvs)[[1]][2]
            param.array[j,paste("T",names(param.marT),sep="."),k,scenario.name[cnt.scenario]]<-unlist(param.marT)
            param.array[j,paste("I",names(param.marI),sep="."),k,scenario.name[cnt.scenario]]<-unlist(param.marI)
            param.array[j,names(param),k,scenario.name[cnt.scenario]]<-unlist(param)

        #Load raw T and F projections from TMB model
              T1<- sc_envProj_rqvys[climate[cc],quant[qq],"T.S2.ParrSu",,jj.sar]#*env_sc["T.S2.ParrSu"]+env_mu["T.S2.ParrSu"]
              F1<- sc_envProj_rqvys[climate[cc],quant[qq],Fname,,jj.sar]#*env_sc[Fname]+env_mu[Fname]
          #Add trends
               if(cc>1) T1<-T1+FW[,paste0("Tsu.Q",quant[qq],".RCP",climate[cc])]
               if(cc>1) F1<-F1+FW[,paste0(FW.name,quant[qq],".RCP",climateshort[cc])]
          #Scale for gompertz model coefficients
               T1<-(T1-scalefactors["mean","T.S2.ParrSu"])/scalefactors["sd","T.S2.ParrSu"]
               F1<-(F1-scalefactors["mean",Fname])/scalefactors["sd",Fname]
               
               
#Initialize all age classes 
     N[1:5] = c(0,0,0,0,0)	
        for (i in 1:5){ 
          xx<-sample(1:dim(envProj_rqvys)[5],1)
          N.new[1] = popinit[pop[k],"stationaryMeanSp"]*exp(param$p1+param$c1*log(popinit[pop[k],"stationaryMeanSp"]/ha))	
          if(filename=="p2Fsu")        surv.trib[i] = exp(param$p2+param$c2*log(max(N[1]/ha,1))+param$BF*F1[i])
          if(length(grep("p2Tsu",filename))>0)   {surv.trib[i] = exp(param$p2+param$c2*log(max(N[1]/ha,1))+param$BT*T1[i]+param$BF*F1[i])}
          if(length(grep("c2",filename))>0)   {surv.trib[i] = exp(param$p2+param$c2*log(max(N[1]/ha,1))+param$BT*T1[i]*log(max(N[1]/ha,1))+param$BF*F1[i]*log(max(N[1]/ha,1)))}
          surv.inriver[i]=  hydro[hydro$gridYear==envProj_rqvys[climate[cc],quant[qq],"gridYear",i,xx] & 
                    hydro$flow==envProj_rqvys[climate[cc],quant[qq],"gridFlow",i,xx]
                  & hydro$temp==envProj_rqvys[climate[cc],quant[qq],"gridTemp",i,xx],"surv"]
          pt[i]=hydro[hydro$gridYear==envProj_rqvys[climate[cc],quant[qq],"gridYear",i,xx] & 
                  hydro$flow==envProj_rqvys[climate[cc],quant[qq],"gridFlow",i,xx]
                & hydro$temp==envProj_rqvys[climate[cc],quant[qq],"gridTemp",i,xx],"proptrans"]
          surv.main[i] =S2.mainstem(pt=pt[i], s.inriver=surv.inriver[i])
          pt.sar[i]<-pt[i]*0.98/((1-pt[i])*surv.inriver[i]+pt[i]*0.98)
          N.new[2] = N[1]*surv.trib[i]*surv.main[i]
          N.new[3] = N[2]*s3.fxn(pt=pt.sar[i],sarI=plogis(sarI_rqys[climate[cc],quant[qq],i,xx]),sarT=plogis(sarT_rqys[climate[cc],quant[qq],i,xx]),b3=param$b3,b4=param$b4,So=param$s0)$s3
          N.new[4] = N[3]*(1-param$b3)*param$s0
          N.new[5] = N[4]*(1-param$b4)*param$s0
          N = N.new
        };N

             for (i in 1:years){
          # Calculate effective spawners for recruit calculations
          #In the first year, we do not have pt.sup yet, so we are just using the average survival of upstream migrants
              if(i==1){
                      myvar="aprmayjunetempLGR";lgrtemp<-sc_envProj_rqvys[climate[cc],quant[qq],myvar,i,jj.sar]#*env_sc[myvar]+env_mu[myvar]
                      myvar="aprmayjuneflowLGR";lgrflow<-sc_envProj_rqvys[climate[cc],quant[qq],myvar,i,jj.sar]#*env_sc[myvar]+env_mu[myvar]
                      simtemp<-round_any(lgrtemp,0.5)#;simtemp
                      simflow<-round_any(lgrflow,20)#;simflow
                      su.Temp.matches<-su[which(abs(simtemp-(su$LGRtemp))==min(abs(simtemp-(su$LGRtemp)))),]#;su.Temp.matches
                      su.Flow.matches<-su.Temp.matches[which(abs(simflow-(su.Temp.matches$LGRflow))==min(abs(simflow-(su.Temp.matches$LGRflow)))),]#;su.Flow.matches
                      surv.sup[i]<-sample(su.Flow.matches$SpringSurv_BoToLG,1)                                         
                      if(k %in% c(6,8)) surv.sup[i]<-sample(su.Flow.matches$SummerSurv_BoToLG,1)
            
                      spawners = (param$b4*N[4] + param$f5*N[5])*surv.sup[i]*0.9
              }
               
          #In subsequent years, we use output from SARs in the previous year to determine the precent of upstream migrants that were transported as juveniles 
            if(i>1){    spawners = (param$b4*N[4] + param$f5*N[5])*surv.sup[i-1]*0.9 }
               
    #Survival in each age class           
          surv.s1[i]=exp(param$p1+param$c1*log(spawners/ha))
                 if(filename=="p2Fsu")        surv.trib[i] = exp(param$p2+param$c2*log(max(N[1]/ha,1))+param$BF*F1[i])
                 if(length(grep("p2Tsu",filename))>0)   {surv.trib[i] = exp(param$p2+param$c2*log(max(N[1]/ha,1))+param$BT*T1[i]+param$BF*F1[i])}
                 if(length(grep("c2",filename))>0)   {surv.trib[i] = exp(param$p2+param$c2*log(max(N[1]/ha,1))+param$BT*T1[i]*log(max(N[1]/ha,1))+param$BF*F1[i]*log(max(N[1]/ha,1)))}
          surv.inriver[i]=  hydro[hydro$gridYear==envProj_rqvys[climate[cc],quant[qq],"gridYear",i,jj.sar] &
                                  hydro$flow==envProj_rqvys[climate[cc],quant[qq],"gridFlow",i,jj.sar] &
                                  hydro$temp==envProj_rqvys[climate[cc],quant[qq],"gridTemp",i,jj.sar],"surv"]
          pt[i]=  hydro[hydro$gridYear==envProj_rqvys[climate[cc],quant[qq],"gridYear",i,jj.sar] & 
                        hydro$flow==envProj_rqvys[climate[cc],quant[qq],"gridFlow",i,jj.sar] &
                        hydro$temp==envProj_rqvys[climate[cc],quant[qq],"gridTemp",i,jj.sar],"proptrans"]
          surv.main[i] =S2.mainstem(pt=pt[i], s.inriver=surv.inriver[i])
          pt.sar[i]<-pt[i]*0.98/((1-pt[i])*surv.inriver[i]+pt[i]*0.98)
          surv.s2[i] =surv.trib[i]*surv.main[i] 
                s3 = s3.fxn(pt=pt.sar[i],sarI=plogis(sarI_rqys[climate[cc],quant[qq],i,jj.sar]),sarT=plogis(sarT_rqys[climate[cc],quant[qq],i,jj.sar]),b3=param$b3,b4=param$b4,So=param$s0)
          surv.s3[i]   = s3$s3
          surv.sarI[i] = s3$sarI
          surv.sarT[i] = s3$sarT
          surv.sar[i]  = s3$sar
      #tracking in-river and transported fish separately for upstream migration survival
          pt.sup[i]<-pt.sar[i]*surv.sarT[i]/surv.sar[i]

              myvar="aprmayjunetempLGR";lgrtemp<-sc_envProj_rqvys[climate[cc],quant[qq],myvar,i,jj.sar]#*env_sc[myvar]+env_mu[myvar]
              myvar="aprmayjuneflowLGR";lgrflow<-sc_envProj_rqvys[climate[cc],quant[qq],myvar,i,jj.sar]#*env_sc[myvar]+env_mu[myvar]
              simtemp<-round_any(lgrtemp,0.5)#;simtemp
              simflow<-round_any(lgrflow,20)#;simflow
              su.Temp.matches<-su[which(abs(simtemp-(su$LGRtemp))==min(abs(simtemp-(su$LGRtemp)))),]#;su.Temp.matches
              su.Flow.matches<-su.Temp.matches[which(abs(simflow-(su.Temp.matches$LGRflow))==min(abs(simflow-(su.Temp.matches$LGRflow)))),]#;su.Flow.matches
                          xx<-sample(1:nrow(su.Flow.matches),size=1);xx
                  sup.spr.inriver[i]<-su.Flow.matches$Spr_Pred_BOLG_inriver[xx]                                         
                  sup.spr.transport[i]<-su.Flow.matches$Spr_Pred_BOLG_transport[xx]                                         
                  sup.sum.inriver[i]<-su.Flow.matches$Sum_Pred_BOLG_inriver[xx]                                         
                  sup.sum.transport[i]<-su.Flow.matches$Sum_Pred_BOLG_transport[xx]                                         
                  
                  sup.spring[i]<-pt.sup[i]* sup.spr.transport[i]+(1-pt.sup[i])* sup.spr.inriver[i]                 
                  sup.summer[i]<-pt.sup[i]* sup.sum.transport[i]+(1-pt.sup[i])* sup.sum.inriver[i]                                        
          surv.sup[i]<-sup.spring[i]                                        
          if(k %in% c(6,8)) surv.sup[i]<-sup.summer[i]
          #Note that surv.sup[i] will be used in next year's calculation of spawners, not this year's N[1] 
          
          N.new[1] = spawners*surv.s1[i]	
          N.new[2] = N[1]*surv.s2[i]
          N.new[3] = N[2]*surv.s3[i]
          N.new[4] = N[3]*(1-param$b3)*param$s0
          N.new[5] = N[4]*(1-param$b4)*param$s0
          
          N = N.new
          
      # calculate spawners
          Spawners[i] = (param$b4*N[4] + N[5])*surv.sup[i]*0.9
          # calculate recruits referenced to brood year
          if (i > 4 && i <= (years-1))
               Recruits[i-4] = param$b4*N[4]*surv.sup[i]*0.9
          if (i > 5)
               Recruits[i-5] = Recruits[i-5] + N[5]*surv.sup[i]*0.9
        } #end years
        Spawners2<-Spawners[5:years]
        
        spawner.array[,j,k,cnt.scenario]<-Spawners
        parrsmolt.array[,j,k,cnt.scenario]<-surv.trib
        recruit.array[,j,k,cnt.scenario]<-Recruits
  }  #end runs
     
     
#Store spawner and recruit matrices and name them by scenario
      save(spawner.array,file=paste("output/spawner",filename,arrayname,sep="."))
      save(param.array,file=paste("output/param",filename,arrayname,sep="."))
      save(parrsmolt.array,file=paste("output/parrsmolt",filename,arrayname,sep="."))
      save(recruit.array,file=paste("output/recruit",filename,arrayname,sep="."))


      cnt=cnt+1
      
      
    }   #end quantiles qq
    }   #end climate cc
}    #end populations

#Add response metrics to param.array-------------
            years<-2015:2089 #dimnames(spawner.array)[[1]]
            runs<-dim(spawner.array)[[2]]
            pop<-dimnames(spawner.array)[[3]]
            scenarios<-dimnames(spawner.array)[[4]]
            spQET<-spawner.array
                spQET[1:3,,,]<-NA
                for(k in 1:length(pop)){
                    for(ss in 1:length(scenarios)){
                     for ( j in 1:runs){
                        spawners<-spawner.array[,j,k,scenarios[ss]];spawners
                        for (i in 4:75){ifelse (sum(spawners[(i-3):i]) < 200,spQET[i,j,k,ss]<-1,spQET[i,j,k,ss]<-0) }
                        param.array[j,"yrQET50",k,scenarios[ss]]<-ifelse(max(spQET[,j,k,scenarios[ss]],na.rm=TRUE)>0.5  ,   
                                    years[min(which(spQET[,j,k,scenarios[ss]]>0.5))]  ,  2115)
                        param.array[j,"Mean2020",k,scenarios[ss]]<-exp(mean(log(spawners[6:15])))
                        param.array[j,"Mean2040",k,scenarios[ss]]<-exp(mean(log(spawners[26:35])))
                        param.array[j,"Mean2060",k,scenarios[ss]]<-exp(mean(log(spawners[46:55])))
                        param.array[j,"Mean2080",k,scenarios[ss]]<-exp(mean(log(spawners[66:75])))
                        param.array[j,"Meanallyears",k,scenarios[ss]]<-exp(mean(log(spawners)))
                        } #end runs
                      }#end scenarios
                    }#end pop

      save(param.array,file=paste("output/param",filename,arrayname,sep="."))


 
