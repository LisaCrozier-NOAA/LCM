S2.mainstem = function( pt, s.inriver,s2.fac=1){
	s2 = (pt*0.98 + (1.0-pt)*s.inriver)
	s2 = s2.fac*s2 # potential hydro improvements
	s2	
}


s3.fxn<-function(pt,sarI,sarT,b3=b3,b4=b4,So=So){
  a1=b3
  a2=(1-b3)*b4*So
  a3=(1-b3)*(1-b4)*So*So
  SAR=pt*sarT + (1-pt)*sarI
  s3=SAR/(a1+a2+a3)
  x<-list(s3=s3,sarI=sarI,sarT=sarT,sar=SAR)
  return(x)
}