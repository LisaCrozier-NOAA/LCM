
rootDir <- path #root directory
retroDir <- paste0(rootDir,"/OUTPUT_LISA_SAR_RETRO/")

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
    for(mm in c("inriver","transport")){
      if(mm=="inriver"){
        cDat <- read.table(paste0(rootDir,"/DATA_COMPASS_FILES_PA/ch1_Flow-",ff,"_",tt,"deg_bonarrival.out"), header = TRUE)
        tmpcDat <- cDat[,3:367]
        gridData[mm,,,flowSeqID[flowSeq==ff],tempSeqID[tempSeq==tt]] <- as.matrix(tmpcDat)
      }
      if(mm=="transport"){
        cDat <- read.table(paste0(rootDir,"/DATA_COMPASS_FILES_PA/ch1_Flow-",ff,"_",tt,"deg_bonarrival_transport.out"), header = TRUE)
        tmpcDat <- cDat[cDat$MigrationType==mm,4:368]
        gridData[mm,,,flowSeqID[flowSeq==ff],tempSeqID[tempSeq==tt]] <- as.matrix(tmpcDat)
      }
    }
  }
}
gridYears <- unique(cDat[,1])
save(gridData,gridYears,file=paste0(rootDir,"gridData.rData"))
