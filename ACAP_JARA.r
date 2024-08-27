##########################################################################
#
# ACAP PRIORITY SPECIES POPULATION TRENDS
#
##########################################################################
# requested by Richard Phillips on 15 Aug 2024
# based on https://github.com/Henning-Winker/JARA
# description: https://www.biorxiv.org/content/10.1101/672899v3
# written by Steffen Oppel, August 2024

# updated on 27 Aug 2024 to reduce observation error after Richard Phillips mentioned that trends were too positive

# library(devtools)
# install_github("henning-winker/JARA")
library(JARA)
library(tidyverse)
library(runjags)
library(lubridate)
library(data.table)
library(readxl)
filter<-dplyr::filter
select<-dplyr::select




###################################################################################
##   1. READ IN DATA FROM DATABASE AND FILTER DATA FOR SPECIES OF INTEREST ####
###################################################################################

try(setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Marine\\ACAP_trend_assessment"), silent=T)
try(setwd("C:\\Users\\steffenoppel\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Marine\\ACAP_trend_assessment"), silent=T)

ACAP<-fread("PriorityPopsCountsForRichard13Aug2024.csv") %>%
  mutate(Year=as.numeric(substr(Year,1,4))) %>%
  separate(Count, c("min", "max"), "-") %>%
  mutate(min=as.numeric(min),max=as.numeric(max), Reliability=as.numeric(Reliability)) %>%
  mutate(Count=if_else(is.na(max),min,(min+max)/2)) %>%
  select(-min,-max) %>%
  filter(!is.na(Count)) %>%
  bind_rows(fread("TRAL_additions.csv")) %>%
  rename(ID=`bs id`, ACC=`Survey accuracy`,colony=`Part site name`) %>%
  select(ID,Species,site,colony,Year,Count,ACC) %>%
  filter(!(Species=="Diomedea dabbenena" & colony=="")) %>%
  mutate(ACC=if_else(ACC %in% c("High","H","high"),0.1,
                     if_else(ACC %in% c("Medio","Medium"),0.2,0.3)))   ### updated uncertainty to reduce obs error
  # mutate(ACC=if_else(ACC %in% c("High","H","high"),0.25,
  #                    if_else(ACC %in% c("Medio","Medium"),0.4,0.6)))
head(ACAP)  



#############################################################################
##   2. PREPARE THE POPULATION META DATA (for which an assessment is needed) ################
#############################################################################
## number of species/assessments
SPECIES<-unique(ACAP$Species)

ACAP %>% group_by(Species,colony) %>%
  summarize(nyears=length(unique(Year)), span=max(Year)-min(Year)) %>%
  ungroup() %>%
  group_by(Species) %>%
  summarize(nsites=length(unique(colony)), max_years=max(nyears), max_span=max(span), min_years=min(nyears), minspan=min(span))


table(ACAP$ACC)
## remove data from sites with <2 counts!



#########################################################################
##   3. RUN RED LIST ASSESSMENT (SSM MODEL) IN JARA FOR EACH SPECIES
#########################################################################

JARAout<-tibble()
SPECIES<-SPECIES[-(1:(length(list.dirs())-1))]   ### if loop is interrupted due to a plotting error retsart without re-running all completed species
for (sp in SPECIES) {
  obserr<-ACAP %>% filter(Species==sp) %>% summarise(err=mean(ACC))
  X<-ACAP %>% filter(Species==sp) %>%
    mutate(site=paste(site,colony,sep="_")) %>%
    group_by(site,Year) %>%
    summarise(N=mean(Count,na.rm=T)) %>%
    mutate(N=if_else(N==0,0.0001,N))  ## log(0) will result in an error, so replace 0 counts with tiny value
  bullshit<-X %>% ungroup() %>% group_by(site) %>%
    summarise(ny=length(unique(Year))) %>%
    filter(ny<2)
  X<-X %>% 
    filter(!(site %in% bullshit$site)) %>%
    spread(key=site, value=N)
  
  ## create complete timeseries
  timeseries<-seq(min(X$Year),max(X$Year,1))
  
  ## add missing rows
  xadd<-timeseries[which(!(timeseries %in% X$Year))]
  for(my in xadd){
    Xadd<-X[1,]
    Xadd[1,]<-NA
    Xadd$Year<-my
    X<-bind_rows(X,Xadd)
  }
  X<-X %>% arrange(Year)
  
  #--------------------------------------------------
  # Build JARA model
  #--------------------------------------------------
  dir.create(sprintf("./%s",sp),showWarnings = F)
  setwd(sprintf("%s/%s",getwd(),sp))
  jara.input = build_jara(I=X,
                          #se=ts$SE,
                          model.type = "census",
                          assessment = sp,
                          scenario = "noK",
                          GL=20,  # generation length in years - set to 20 so we get 20 year trend
                          fixed.obsE=as.numeric(obserr),  # numeric, average across all the counts
                          pk.prior = c(0.0001,0.1)
  )
  # Check input
  jrplot_indices(jara.input,as.png = T,output.dir = getwd())
  
  # Run JARA with csv file output
  fit = fit_jara(jarainput =jara.input,save.csvs = T,output.dir = getwd())
  
  # plot all routine output
  try(jara_plots(fit,output.dir = getwd()), silent=T)
  
  ## RETURN OUTPUT
  out<- X %>% gather(key=site, value=count,-Year) %>%
    group_by(site) %>%
    summarize(nyears=length(unique(Year)), span=max(Year)-min(Year)) %>%
    ungroup() %>%
    summarize(nsites=length(unique(site)), max_years=max(nyears), max_span=max(span)) %>%
    mutate(obs.error=as.numeric(obserr)) %>%
    mutate(Species=sp,
           full.trend=quantile(fit$posteriors$r.All.yrs,0.5),
           lcl.full.trend=quantile(fit$posteriors$r.All.yrs,0.025),
           ucl.full.trend=quantile(fit$posteriors$r.All.yrs,0.975)
           )
  if(dim(fit$posteriors)[2]>2){
    out<- out %>%
      mutate(Species=sp,
             y20.trend=quantile(fit$posteriors$r.1GL,0.5),
             lcl.y20.trend=quantile(fit$posteriors$r.1GL,0.025),
             ucl.y20.trend=quantile(fit$posteriors$r.1GL,0.975)
      )
  }
  
  try(setwd("C:\\STEFFEN\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Marine\\ACAP_trend_assessment"), silent=T)
  try(setwd("C:\\Users\\steffenoppel\\OneDrive - THE ROYAL SOCIETY FOR THE PROTECTION OF BIRDS\\STEFFEN\\RSPB\\Marine\\ACAP_trend_assessment"), silent=T)
  JARAout<-bind_rows(JARAout,out)
  fwrite(JARAout,"ACAP_priority_species_trend_assessment.csv")

  
}





#########################################################################
##   discarded scratchpad code chunks - necessary because JARA package required debugging
#########################################################################

jrplot_state



### troubleshoot function
jara=fit
type=NULL
ref.yr=NULL
extinction=0.01
credibility=0.95
output.dir=getwd()
as.png=FALSE
width=5
height=4.5
ylab = "Density"
xlab="Relative state"
xlim=NULL
plot.cex=1
legend.cex=0.9
legend.pos="right"
add=FALSE
  
cat(paste0("\n","><> jrplot_state() - %change relative to reference year <><","\n"))

Par = list(mfrow=c(1,1),mar = c(4, 4, 1, 1), mgp =c(2.5,1,0),mai = c(0.6, 0.6, 0.1, 0.1),mex=0.8, tck = -0.02,cex=plot.cex)
if(as.png==TRUE){png(file = paste0(output.dir,"/State_",jara$assessment,"_",jara$scenario,".png"), width = width, height = height,
                     res = 200, units = "in")}
if(add==FALSE) par(Par)

pdyn = jara$pop.posterior
yrs = 1:ncol(pdyn)
nyrs = length(yrs)
yr = as.numeric(names(pdyn))
if(is.null(ref.yr)) ref.yr = yr[1:3]
end.yr = max(jara$yr) 
prj.yr = max(jara$pyr)
pop.ref = apply(pdyn[,which(yr%in%ref.yr)],1,mean)
states =  cbind(pdyn[,which(yr%in%end.yr)]/pop.ref,pdyn[,which(yr%in%prj.yr)]/pop.ref)
states[is.nan(states)] <- 0  ### this is necessary if a population goes extinct to prevent an error in the following loop (missing values)
if(is.null(type)){
  type=ifelse(prj.yr-end.yr<3,"current","both") 
}


lymax=rymax = lxrange = lxmax =NULL # maximum and range for plotting
for(i in 1:2){
  if(i == 1 & type =="current" | type== "both" |i == 2 & type =="projected"){
    den = stats::density(states[,i],adjust=2)
    assign(paste0("xl",i),den$x)
    assign(paste0("yl",i),den$y)
    lymax=c(lymax,max(den$y))
    lxmax = c(lxmax,quantile(states[,i],0.99))
  }}

lxrange = ifelse(lxrange<0,0,lxrange)
if(is.null(xlim)) xlim = c(0,max(lxmax,1.1))
jcol = c(grey(0.4,0.6),rgb(1,0,0,0.6))
plot(0,0,type="n",ylab=ylab,xlab=xlab,xaxt="n",yaxt="n",cex.main=0.9,ylim=c(0,1.22*max(lymax)),xlim=xlim,xaxs="i",yaxs="i",frame=FALSE) 
for(i in 2:1){
  if(i == 1 & type =="current" | type== "both" |i == 2 & type =="projected"){
    x = get(paste0("xl",i))
    y = get(paste0("yl",i))
    xp = x[x>xlim[1] & x<xlim[2]]
    yp = y[x>xlim[1] & x<xlim[2]]
    polygon(c(xp,rev(xp)),c(yp,rep(0,length(yp))),col=jcol[i],border=NA)
    mu = round(median(states[,i]),10)
    lines(rep(mu,2),c(0,max(lymax*c(1.05,1.0)[i])),col=c(1,2)[i],lwd=1,lty=c(1))
    text(max(mu,0.05),max(lymax*c(1.11,1.05)[i]),c(end.yr,prj.yr)[i],cex=0.9)
  }}
axis(1,at=seq(0,ceiling(max(states)),0.2),cex.axis=0.9)
axis(2,cex.axis=0.9)

lines(rep(1,2),c(0,max(lymax*c(1.13,1.0)[1])),col=1,lwd=1,lty=2)
text(1,max(lymax*c(1.18)),min((ref.yr)),cex=0.9)

cnam = c(paste0("Cur = ",round(median(states[,1]),2)),paste0("Proj = ",round(median(states[,2]),2)))
if(type =="current") type.id = 1 
if(type =="projected") type.id = 2 
if(type =="both") type.id = 1:2 
legend(legend.pos,cnam[type.id],pch=15,col=c(jcol),box.col = "white",cex=legend.cex,y.intersp = 0.8,x.intersp = 0.8)
mu =apply(states,2,quantile,c(0.5))
quants = rbind(mu,HDInterval::hdi(states,credMass=credibility))
box()
state = NULL
state$state = data.frame(State=cnam,year=c(end.yr,prj.yr),median=quants[1,],lci=quants[2,],uci=quants[3,])
if(type=="current"){state$prob.pextinct = "Requires projection horizon"} else {
  state$prob.extinct = sum(ifelse(states[,2]<extinction,1,0))/nrow(states)
}