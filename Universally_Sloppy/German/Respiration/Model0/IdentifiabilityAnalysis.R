# Load required packages
library(SoilR)
library(FME)
library(beeswarm)

#### Implementation of microbial model with radiocarbon
lambda=0.0001209681/(365*24) #radiocarbon decay constant on day^-1
pars=c(tau_M=0.0005,Vm=0.0126,Km=360,e=0.5) #Parameters from German et al. (2012, GCB 18:1468)
X0=c(SO=100,MIC=2,SO14=1*100,MIC14=1*2) #Initial condition for radiocarbon are in F*C


#The Michaelis-Menten model
MM=function(time, pars, X0){
  MMsys=function(time, pars, X0){
    with(as.list(c(pars,X0)),{
      dSO=tau_M*MIC-MIC*Vm*SO/(Km+SO)
      dMIC=(e*MIC*Vm*SO/(Km+SO))-tau_M*MIC
      
      dSO14=tau_M*MIC14-(MIC14*Vm*SO14/(Km+SO14))-lambda*SO14
      dMIC14=(e*MIC14*Vm*SO14/(Km+SO14))-(tau_M*MIC14)-lambda*MIC14
  
      return(list(c(dSO,dMIC,dSO14,dMIC14)))
    })
  }
  
  odeout=ode(y=X0,parms=pars,times=time,func=MMsys)
  R=(1-pars["e"])*odeout[,"MIC"]*pars["Vm"]*odeout[,"SO"]/(pars["Km"]+odeout[,"SO"])
  SOD14=((odeout[,"SO14"]/odeout[,"SO"])-1)*1000
  MICD14=((odeout[,"MIC14"]/odeout[,"MIC"])-1)*1000
  return(as.data.frame(cbind(odeout[,1:3],SOD14,MICD14,R))) #Output in Delta14C
}


# Create artificial bbservations
t_end=24*30*12 #One year with 12 months with 30 days each
Obst=append(seq(17.31462926,t_end,17.31462926), t_end)
ObsM=MM(Obst,pars,X0)

plot(ObsM$time, ObsM$R)

# Respiration Data
DataR <- cbind(time = ObsM$time,
               R = ObsM$R*(1 + 0.001*rnorm(sd = 1, n = length(ObsM$R))),
               sd = 0.001*max(ObsM$R))

DataR14 <- cbind(time = ObsM[seq(1, 499, 41.5833),c(1)],
                 MICD14 = ObsM[seq(1, 499, 41.5833),c(5)]*(1 + 0.001*rnorm(sd = 1, n = length(ObsM[seq(1, 499, 41.5833),c(5)]))),
                 sd = 0.001*max(abs(ObsM[seq(1, 499, 41.5833),c(5)])))

DataS14 <- cbind(time = ObsM[seq(0, 499, 499),c(1)],
                 SOD14 = ObsM[seq(0, 499, 499),c(4)]*(1 + 0.001*rnorm(sd = 1, n = length(ObsM[seq(0, 499, 499),c(4)]))),
                 sd = 0.001*max(abs(ObsM[seq(0, 499, 499),c(4)])))


#Cost function
costMM=function(parameters){
  out <- MM(time=seq(0,t_end),pars=parameters,X0=X0)
  cost1 <- modCost(model=out, obs=DataR, err="sd")
  cost2 <- modCost(model=out,obs=DataR14, err="sd", cost=cost1)
  return(modCost(model=out,obs=DataS14,err="sd", cost=cost2))
}


sensMM=sensFun(func=costMM,parms=pars)

gR=collin(sensMM,which="R")
gR14=collin(sensMM,which=c("R","MICD14"))
gS14=collin(sensMM,which=c("R","SOD14"))
gRS14=collin(sensMM)


#The reduced Michaelis-Menten model
MM1=function(time, pars, X0){
  MMsys=function(time, pars, X0){
    with(as.list(c(pars,X0)),{
      dSO=tau_M*MIC-MIC*Phi1*SO
      dMIC=(e*MIC*Phi1*SO)-tau_M*MIC
      
      dSO14=tau_M*MIC14-(MIC14*Phi1*SO14)-lambda*SO14
      dMIC14=(e*MIC14*Phi1*SO14)-(tau_M*MIC14)-lambda*MIC14
      
      return(list(c(dSO,dMIC,dSO14,dMIC14)))
    })
  }
  
  odeout=ode(y=X0,parms=pars,times=time,func=MMsys)
  R=(1-pars["e"])*odeout[,"MIC"]*pars["Phi1"]*odeout[,"SO"]
  SOD14=((odeout[,"SO14"]/odeout[,"SO"])-1)*1000
  MICD14=((odeout[,"MIC14"]/odeout[,"MIC"])-1)*1000
  return(as.data.frame(cbind(odeout[,1:3],SOD14,MICD14,R))) #Output in Delta14C
}  


#Cost function
costMM1=function(parameters){
  out <- MM1(time=seq(0,t_end),pars=parameters,X0=X0)
  cost1 <- modCost(model=out, obs=DataR, err="sd")
  cost2 <- modCost(model=out,obs=DataR14, err="sd", cost=cost1)
  return(modCost(model=out,obs=DataS14,err="sd", cost=cost2))
}

pars1=c(tau_M=0.000248,Phi1=2.268e-5,e=0.53)
sensMM1=sensFun(func=costMM1,parms=pars1)

gR1=collin(sensMM1,which="R")
gR141=collin(sensMM1,which=c("R","MICD14"))
gS141=collin(sensMM1,which=c("R","SOD14"))
gRS141=collin(sensMM1)


Measurement=c(rep("R",11), rep("R_red",4), rep("R-R14",11), rep("R-R14_red",4), rep("R-S14",11), rep("R-S14_red",4), rep("R-R14-S14",11), rep("R-R14-S14_red",4))
gR1$Km <- 0
gR1 <- gR1[, c("tau_M", "Phi1", "Km", "e","N","collinearity")]
names(gR1)[names(gR1) == 'Phi1'] <- 'Vm'
gR141$Km <- 0
gR141 <- gR141[, c("tau_M", "Phi1", "Km", "e","N","collinearity")]
names(gR141)[names(gR141) == 'Phi1'] <- 'Vm'
gS141$Km <- 0
gS141 <- gS141[, c("tau_M", "Phi1", "Km", "e","N","collinearity")]
names(gS141)[names(gS141) == 'Phi1'] <- 'Vm'
gRS141$Km <- 0
gRS141 <- gRS141[, c("tau_M", "Phi1", "Km", "e","N","collinearity")]
names(gRS141)[names(gRS141) == 'Phi1'] <- 'Vm'

gammaMM=rbind(gR, gR1, gR14, gR141,gS14, gS141,gRS14,gRS141)
gammaMM=cbind(gammaMM,Measurement)

#Figure 3
tmpcol = append(rep(1,11),rep(1,4))
tmpcol1 = append(tmpcol, rep(1,11))
tmpcol2 = append(tmpcol1, rep(1,4))
tmpcol3 = append(tmpcol2, rep(1, 11))
tmpcol4 = append(tmpcol3, rep(1,4))
tmpcol5 = append(tmpcol4, rep(1,11))
tmpcol6 = append(tmpcol5, rep(1,4))
tmpch = append(rep(2,11),rep(17,4))
tmpch1 = append(tmpch, rep(0,11))
tmpch2 = append(tmpch1, rep(15,4))
tmpch3 = append(tmpch2, rep(1, 11))
tmpch4 = append(tmpch3, rep(19, 4))
tmpch5 = append(tmpch4, rep(23, 11))
tmpch6 = append(tmpch5, rep(18, 4))

png("Plot3.png", width = 6.5, height = 5, units = 'in', res = 300)

beeswarm(collinearity~N,data=gammaMM,pwcol=tmpcol6,pwpch=tmpch6,ylab="",xlab="")
abline(h=16.5,lty=2); abline(v=1:3,lty=3,col="gray")
mtext("Collinearity index",side=2,line=3)
mtext("Number of parameters in a combination",side=1,line=3)
legend("left",c("R","R+R14","R+S14","R+R14+S14"),pch=c(2,0,1,23),col=1,bty="n")
dev.off()

# costs_R = c(242.464, 258.108, 1108.53)
# costs_R14 = c(246.291, 1437.225)
# costs_RS14 = c(243.743, 273.72)
# costs_all = c(246.313, 1929.24)
# N = c(4,3)
# N1 = c(4,3,2)
# 
# par(new = T)
# plot(N1, costs_R, type="o", axes=F, xlab=NA, ylab=NA, cex=1.2, pch=0, xlim=c(2,4),ylim=c(0,1400), col='red')
# par(new = T)
# plot(N, costs_R14, type="o", axes=F, xlab=NA, ylab=NA, cex=1.2, pch=1, xlim=c(2,4),ylim=c(0,1400), col='red')
# par(new = T)
# plot(N, costs_RS14, type="o", axes=F, xlab=NA, ylab=NA, cex=1.2, pch=23, xlim=c(2,4),ylim=c(0,1400), col='red')
# axis(side = 4)
# mtext(side = 4, line = 3, 'Cost')
