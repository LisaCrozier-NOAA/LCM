myMap <- list()

if(retro)
  nYears <- length(eStartYr:eLastYr)
if(!retro)
  nYears <- length(startYr:eLastYr)

if(re_j == 0){
  map <- list(
    eps_j = as.factor(matrix(NA,nrow(parameters$eps_j),nk_dim))
    ,frho_j = as.factor(rep(NA,nk_dim))
    ,fpsi_j = as.factor(rep(NA,nk_dim))
  )
  myMap <- append(myMap,
                  map)
}

if(re_t == 0){
  map <- list(
    eps_t = as.factor(matrix(NA,nYears,nk_dim))
    ,frho_t = as.factor(rep(NA,nk_dim))
    ,fpsi_t = as.factor(rep(NA,nk_dim))
  )
  myMap <- append(myMap,
                  map)
  
}

if(re_jt == 0){
  map <- list(
    eps_jt = as.factor(array(NA,c(nk_dim,nrow(parameters$eps_j),nYears)))
    ,frho1_jt = as.factor(rep(NA,nk_dim))
    ,frho2_jt = as.factor(rep(NA,nk_dim))
    ,fpsi_jt = as.factor(rep(NA,nk_dim))
  )
  myMap <- append(myMap,
                  map)
}


if(fixed_mar == 0){
  map <- list(
    beta_mar = as.factor(matrix(NA,length(tmpMarVars),max(data$k)+1))
  )
  myMap <- append(myMap,
                  map)
}

if(mean_s == 0){
  map <- list(
    mu_s = as.factor(rep(NA,max(data$k)+1))
  )
  myMap <- append(myMap,
                  map)
}

if(cov_pars == 0){
  map <- list(
    frho_x = as.factor(NA)
    ,frho_Rx = as.factor(rep(NA,length(myVars)*(length(myVars)-1)/2))
    ,fpsi_x = as.factor(rep(NA,length(myVars)))
    ,eps_x = as.factor(matrix(NA,dim(parameters$eps_x)[1],dim(parameters$eps_x)[2]))
  )
  myMap <- append(myMap,
                  map)
}

obj <- MakeADFun(data = data,
                 parameters = parameters,
                 map = myMap,
                 random=c("eps_j"
                          ,"eps_t"
                          ,"eps_jt"
                          ,"eps_x"
                          ,"frho_Rx"
                 ),
                 silent = FALSE,
                 bias.correct=TRUE,
                 DLL = "integrated2")

# out <- tryCatch(TMBhelper::fit_tmb(obj = obj
#                                    ,getsd=TRUE
#                                    ,getJointPrecision = TRUE
#                                    ,newtonsteps = 3),
#                 error=function(e) NULL)

out <- nlminb(obj$par,obj$fn,obj$gr)
rep <- obj$report()
SD <- sdreport(obj, getJointPrecision = TRUE)

