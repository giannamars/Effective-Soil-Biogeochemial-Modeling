# Load required packages
library(SoilR)
library(FME)

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

# Respiration Data
DataR <- cbind(time = ObsM$time,
               R = ObsM$R*(1 + 0.001*rnorm(sd = 1, n = length(ObsM$R))),
               sd = 0.001*max(ObsM$R))
               
               
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

ObsM1=MM1(Obst,pars1,X0)

#Cost function
costMM=function(parameters){
  out <- MM1(time=seq(0,t_end),pars=parameters,X0=X0)
  cost1 <- modCost(model=out, obs=DataR, err="sd")
  return(cost1)
}


pars1=c(tau_M=0.000248,Phi1=2.268e-5,e=0.53)
costMM(pars1)$model

sensMM=sensFun(func=costMM,parms=pars1)
ident <- collin(sensMM)

plot(ident, log = "y")
