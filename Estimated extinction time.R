##################################################################################################################
## EXTINCTION TIME
##################################################################################################################
##Data obtained from bootstraps from "Population growth" script needed:
#pla.bootlam.all
#ant.bootlam.all
#tri.bootlam.all

#lambdas_cis

## data.plot_pla
## data.plot_ant
## data.plot_tri
## data.plot_cam

##################################################################################################################
## Line 27-67 "Predicted extinction time" - Line 68-213 "Transform bootstrapped population growth to extinction time"
## Line 214- 269 "BCPIs extinction time" - Line 270-309 "Pairwise comparisons extinction time"
## Line 310-393 "Trendlines for Fig.3" - Line 394-548 "Plot extinction time"
##################################################################################################################

## For explanation of sites and treatments, see script "Population growth"

## ! ## Note where scripts differ slightly between species:

## lines 354-368 - A.vulneraria and T.badium - limited trendline slope for plotting

##################################################################################################################

#################################################################################
##           1. Predicted extinction time
#################################################################################

#################################################################################
## Run one species at the time until Plotting (line.393)
#################################################################################

##subset for bootstrapped lambda dataset
dataset=pla.bootlam.all

##subset for species
species="pla"

##subset for initial population size

pop.size=49000 ##example P.alpina

### n0 = estimated population size (see Supporting information) 
## P.alpina (pla) = n0 = 49000
## A. alpestris (ant) = n0 = 147000 
## T. badium (tri) = n0 = 50100
## C. scheucheri = n0 = 228000


Nt=1
lambdas_cis_p = subset(lambdas_cis, subset = Species == species)
n0=pop.size
extinction_time<-log(Nt/n0)/log(lambdas_cis_p$lambda)

extinction_original<- extinction_time

name <- paste(species, "extinction_original",sep=".")
assign(name, extinction_original)

#### Extinction time for each species saved:

##pla.extinction_original
##ant.extinction_original
##tri.extinction_original
##cam.extinction_original

#################################################################################
##           2. TRANSFORM BOOTSTRAPPED LAMBDAS TO EXTINCTION TIME  
#################################################################################

## Sort boot data and add the corresponding treatment, site, elevation and species

dataset<-as.data.frame(dataset)

p_cs<-as.data.frame(dataset$p.cs)
colnames(p_cs)[1] <- "lambda"
p_cs$Treatment<- "Current species"
p_cs$Site<-"P"
p_cs$Species<-species
p_cs$elev<-"1400"

m_cs<-as.data.frame(dataset$m.cs)
colnames(m_cs)[1] <- "lambda"
m_cs$Treatment<- "Current species"
m_cs$Site<-"M"
m_cs$Species<-species
m_cs$elev<-"1750"

c_cs<-as.data.frame(dataset$c.cs)
colnames(c_cs)[1] <- "lambda"
c_cs$Treatment<- "Current species"
c_cs$Site<-"C"
c_cs$Species<-species
c_cs$elev<-"1950"

r_cs<-as.data.frame(dataset$r.cs)
colnames(r_cs)[1] <- "lambda"
r_cs$Treatment<- "Current species"
r_cs$Site<-"R"
r_cs$Species<-species
r_cs$elev<-"2200"

p_ns<-as.data.frame(dataset$p.ns)
colnames(p_ns)[1] <- "lambda"
p_ns$Treatment<- "Novel species"
p_ns$Site<-"P"
p_ns$Species<-species
p_ns$elev<-"1400"

m_ns<-as.data.frame(dataset$m.ns)
colnames(m_ns)[1] <- "lambda"
m_ns$Treatment<- "Novel species"
m_ns$Site<-"M"
m_ns$Species<-species
m_ns$elev<-"1750"

c_ns<-as.data.frame(dataset$c.ns)
colnames(c_ns)[1] <- "lambda"
c_ns$Treatment<- "Novel species"
c_ns$Site<-"C"
c_ns$Species<-species
c_ns$elev<-"1950"

dataset<-rbind(p_cs,m_cs, c_cs,r_cs,p_ns, m_ns,c_ns)
dataset<-as.data.frame(dataset)

## Transform bootstrapped population growth to extinction time

t= time=NULL

for(i in 1:50) {
  random = sample(i, 1)
  random = i
  p.cs = subset(dataset, Treatment == "Current species" & Site == "P" & Species == "pla"); dim(p.cs)
  m.cs = subset(dataset, Treatment == "Current species" & Site == "M"& Species == "pla"); dim(m.cs)
  c.cs = subset(dataset, Treatment == "Current species" & Site == "C"& Species == "pla"); dim(c.cs)
  r.cs = subset(dataset, Treatment == "Current species" & Site == "R"& Species == "pla"); dim(r.cs)
  p.ns = subset(dataset, Treatment == "Novel species" & Site == "P"& Species == "pla"); dim(p.ns)
  m.ns = subset(dataset, Treatment == "Novel species" & Site == "M"& Species == "pla"); dim(m.ns)
  c.ns = subset(dataset, Treatment == "Novel species" & Site == "C"& Species == "pla"); dim(c.ns)
  
  new.data = rbind(p.cs[random,],m.cs[random,],c.cs[random,],r.cs[random,],p.ns[random,],m.ns[random,],c.ns[random,])
  
  n0=pop.size
  Nt=1
  t<-log(Nt/n0)/log(new.data$lambda)
  t
  time = c(time,t)
  print(i) 
  
  
}

########## Organize the translated extinction time and exclude all 0s ########## 

time_p.cs =time[seq(1, length(time), 7)]
time_m.cs =time[seq(2, length(time), 7)]
time_c.cs =time[seq(3, length(time), 7)]
time_r.cs =time[seq(4, length(time), 7)]
time_p.ns =time[seq(5, length(time), 7)]
time_m.ns =time[seq(6, length(time), 7)]
time_c.ns =time[seq(7, length(time), 7)]

time_p.cs<-as.data.frame(time_p.cs)
time_p.cs$Treatment<-c("Current species")
time_p.cs$Site<-c("P")
time_p.cs$Species<-species

time_m.cs<-as.data.frame(time_m.cs)
time_m.cs$Treatment<-c("Current species")
time_m.cs$Site<-c("M")
time_m.cs$Species<-species

time_c.cs<-as.data.frame(time_c.cs)
time_c.cs$Treatment<-c("Current species")
time_c.cs$Site<-c("C")
time_c.cs$Species<-species

time_r.cs<-as.data.frame(time_r.cs)
time_r.cs$Treatment<-c("Current species")
time_r.cs$Site<-c("R")
time_r.cs$Species<-species

time_p.ns<-as.data.frame(time_p.ns)
time_p.ns$Treatment<-c("Novel species")
time_p.ns$Site<-c("P")
time_p.ns$Species<-species

time_m.ns<-as.data.frame(time_m.ns)
time_m.ns$Treatment<-c("Novel species")
time_m.ns$Site<-c("M")
time_m.ns$Species<-species

time_c.ns<-as.data.frame(time_c.ns)
time_c.ns$Treatment<-c("Novel species")
time_c.ns$Site<-c("C")
time_c.ns$Species<-species

time_p.cs_no_0 =subset(time_p.cs ,subset=time_p.cs> 0)
time_m.cs_no_0 =subset(time_m.cs ,subset=time_m.cs> 0)
time_c.cs_no_0 =subset(time_c.cs ,subset=time_c.cs> 0)
time_r.cs_no_0 =subset(time_r.cs ,subset=time_r.cs> 0)
time_p.ns_no_0 =subset(time_p.ns ,subset=time_p.ns> 0)
time_m.ns_no_0 =subset(time_m.ns ,subset=time_m.ns> 0)
time_c.ns_no_0 =subset(time_c.ns ,subset=time_c.ns> 0)

###delete all population growth above 0.99 (see Material and methods)

n0=pop.size
Nt=1
above_0.99<-log(Nt/n0)/log(0.99)

#################################################################################
##           3.            CALCULATE BCPIs   
#################################################################################

bcpi <- function(real_extinction, t, alpha) {
  B <- length(t)
  z0 <- qnorm(mean(t < real_extinction))
  a1 <- pnorm(2 * z0 + qnorm(alpha/2))
  a2 <- pnorm(2 * z0 + qnorm(1 - alpha/2))
  c1 <- quantile(t, a1)
  c2 <- quantile(t, a2)
  return(as.numeric(c(c1, c2)))
}

real_extinction <-extinction_original[1,1]
time_p.cs_del_out =subset(time_p.cs_no_0 ,subset=time_p.cs_no_0$time_p.cs< above_0.99)
bcpi_extinction_p.cs<-bcpi(real_extinction, time_p.cs_del_out$time_p.cs, 0.05)

real_extinction <-extinction_original[3,1]
time_m.cs_del_out =subset(time_m.cs_no_0 ,subset=time_m.cs_no_0$time_m.cs< above_0.99)
bcpi_extinction_m.cs<-bcpi(real_extinction, time_m.cs_del_out$time_m.cs, 0.05)

real_extinction <-extinction_original[5,1]
time_c.cs_del_out =subset(time_c.cs_no_0 ,subset=time_c.cs_no_0$time_c.cs< above_0.99)
bcpi_extinction_c.cs<-bcpi(real_extinction, time_c.cs_del_out$time_c.cs, 0.05)

real_extinction <-extinction_original[7,1]
time_r.cs_del_out =subset(time_r.cs_no_0 ,subset=time_r.cs_no_0$time_r.cs< above_0.99)
bcpi_extinction_r.cs<-bcpi(real_extinction, time_r.cs_del_out$time_r.cs, 0.05)

real_extinction <-extinction_original[2,1]
time_p.ns_del_out =subset(time_p.ns_no_0 ,subset=time_p.ns_no_0$time_p.ns< above_0.99)
bcpi_extinction_p.ns<-bcpi(real_extinction, time_p.ns_del_out$time_p.ns, 0.05)

real_extinction <-extinction_original[4,1]
time_m.ns_del_out =subset(time_m.ns_no_0 ,subset=time_m.ns_no_0$time_m.ns< above_0.99)
bcpi_extinction_m.ns<-bcpi(real_extinction, time_m.ns_del_out$time_m.ns, 0.05)

real_extinction <-extinction_original[6,1]
time_c.ns_del_out =subset(time_c.ns_no_0 ,subset=time_c.ns_no_0$time_c.ns< above_0.99)
bcpi_extinction_c.ns<-bcpi(real_extinction, time_c.ns_del_out$time_c.ns, 0.05)

extinction_bcpi<-rbind(bcpi_extinction_p.cs,bcpi_extinction_p.ns,
                       bcpi_extinction_m.cs,bcpi_extinction_m.ns,
                       bcpi_extinction_c.cs,bcpi_extinction_c.ns,
                       bcpi_extinction_r.cs)

name <- paste(species, "extinction_bcpi",sep=".")
assign(name, extinction_bcpi)

# Extinction time BCPIs saved:
##pla.extinction_bcpi
##ant.extinction_bcpi
##tri.extinction_bcpi
##cam.extinction_bcpi

#################################################################################
##           4.            PAIRWISE comparisons extinction time
#################################################################################

## compare generated extinction times between sites
## 
comp_p.cs_m.cs <-as.data.frame(time_p.cs_del_out$time_p.cs-time_m.cs_del_out$time_m.cs)
comp_p.cs_c.cs <-as.data.frame(time_p.cs_del_out$time_p.cs-time_c.cs_del_out$time_c.cs)
comp_p.cs_r.cs <-as.data.frame(time_p.cs_del_out$time_p.cs- time_r.cs_del_out$time_r.cs)
comp_m.cs_c.cs <-as.data.frame(time_m.cs_del_out$time_m.cs-time_c.cs_del_out$time_c.cs)
comp_m.cs_r.cs <-as.data.frame(time_m.cs_del_out$time_m.cs-time_r.cs_del_out$time_r.cs)
comp_c.cs_r.cs <-as.data.frame(time_c.cs_del_out$time_c.cs-time_r.cs_del_out$time_r.cs)

## compare generated extinction times between communities
comp_p.cs_ns <-as.data.frame(time_p.cs_del_out$time_p.cs- time_p.ns_del_out$time_p.ns)
comp_m.cs_ns <-as.data.frame(time_m.cs_del_out$time_m.cs- time_m.ns_del_out$time_m.ns)
comp_c.cs_ns <-as.data.frame(time_c.cs_del_out$time_c.cs- time_c.ns_del_out$time_c.ns)

dif_extinction<-cbind(comp_p.cs_m.cs,comp_p.cs_c.cs,comp_p.cs_r.cs,comp_m.cs_c.cs,comp_m.cs_r.cs,comp_c.cs_r.cs,
                      comp_p.cs_ns,comp_m.cs_ns,comp_c.cs_ns)


colnames(dif_extinction) <- c("comp_p.cs_m.cs","comp_p.cs_c.cs","comp_p.cs_r.cs","comp_m.cs_c.cs","comp_m.cs_r.cs","comp_c.cs_r.cs",
                              "comp_p.cs_nc","comp_m.cs_nc","comp_c.cs_nc")

dif_extinction<-as.data.frame(dif_extinction)

#### counting number of cases with differences

pairwise_comparison_extinction<-as.data.frame(sapply(dif_extinction, function(x) sum(x >0, na.rm = TRUE)/length(x)))

name <- paste(species, "pairwise_comparison_extinction",sep=".")
assign(name, pairwise_comparison_extinction)

# Pairwise comparisons exticntion time saved:
#pla.pairwise_comparison_extinction
#ant.pairwise_comparison_extinction
#tri.pairwise_comparison_extinction
#cam.pairwise_comparison_extinction

######################################################################################################
## 5.     Translate slopes from lambdas to change in extinction time per each 0.1 C
##        to be plotted in Fig.3  
######################################################################################################

##use this dataset where we have temperatures saved and subset:
data.slopes =data.plot_pla
#data.slopes =data.plot_ant
#data.slopes =data.plot_tri
##data.slopes =data.plot_cam

data.slopes =subset(data.slopes, subset = Site == "P"|Site == "M"|Site == "C")

## Now we use a model to estimate how population growth change with warming (i.e. exclude site 2200 m ("R"; see line above))

model<-lm(lambda~temp*Treatment, data=data.slopes)
sum<-as.data.frame(coef(model))
sum$Species <- species
colnames(sum)<-c("coef", "Species")

slope<-cbind(sum,(summary(model)$coef[,"Pr(>|t|)",drop=F]))

intercept = slope[1,1]
temperature = slope[2,1]
NS = slope[3,1]
int_NS = slope[4,1]

## Produce a temperature range from 1.6-5 and calculate extinction time for each 0.1 degree

range<-seq(1.6, 5, by = 0.1)
range<-as.data.frame(range)

temperature_range_cs<-data.frame(matrix(ncol=5, nrow=35))
colnames(temperature_range_cs)<- c("Species", "Treatment","temperature","lambda", "extinction_time")
temperature_range_cs[,3]<-range
temperature_range_cs[,1]=species
temperature_range_cs[,2]="Current species"
temperature_range_cs[,4]=(intercept +(temperature*(range)))
temperature_range_cs[,5]=log(Nt/n0)/log(temperature_range_cs$lambda)

######################################################################################################
## ! ## Special A. alpestris Current species,
## for plotting, see below line 387-391

#temperature_range_cs[10,3]<- "2.588"
#temperature_range_cs[10,5]<- "1173.367"

######################################################################################################
## ! ## Special T. badium Current species,
## for plotting, see below line 387-391

#temperature_range_cs[10,3]<- "2.0886"
#temperature_range_cs[10,5]<- "399.7184"
#temperature_range_cs= subset(temperature_range_cs, subset = extinction_time >500)
######################################################################################################

temperature_range_ns<-data.frame(matrix(ncol=5, nrow=35))
colnames(temperature_range_ns)<- c("Species", "Treatment","temperature","lambda", "extinction_time")
temperature_range_ns[,3]<-range
temperature_range_ns[,1]=species
temperature_range_ns[,2]="Novel species"
temperature_range_ns[,4]=((intercept+NS) +((temperature+int_NS)*(range)))
temperature_range_ns[,5]=log(Nt/n0)/log(temperature_range_ns$lambda)

temperature_range<-rbind(temperature_range_cs,temperature_range_ns)
## exclude "negative" extinction time i.e. lambdas >1
temperature_range= subset(temperature_range, subset = extinction_time >0)

name <- paste(species, "temperature_range",sep=".")
assign(name, temperature_range)

## Predicted time to extinction across temperature range saved:
#pla.temperature_range
#ant.temperature_range
#tri.temperature_range
#cam.temperature_range

#####################################################################################################################
## We limited the warming levels shown by the trend line for A. alpestris (max. 2.588 C) and T. badium (max. 2.0886 C) 
## for aesthetic reasons, as populations for these cases are not predicted to go extinct and create very steep slopes.
## see figure text Fig.3
#####################################################################################################################

#####################################################################################################################
## 6. Plot extinction time 
#####################################################################################################################
library(scales)

pla.extinction_original<-as.data.frame(pla.extinction_original)
pla.extinction_bcpi<-as.data.frame(pla.extinction_bcpi)
row.names(pla.extinction_original)<-NULL
row.names(pla.extinction_bcpi)<-NULL
pla_ext_comb<-cbind(pla.extinction_original, pla.extinction_bcpi)
colnames(pla_ext_comb)<-c("extinction_time", "lower_extinction", "upper_extinction")
pla_ext_comb$Species<- "pla"

ant.extinction_original<-as.data.frame(ant.extinction_original)
ant.extinction_bcpi<-as.data.frame(ant.extinction_bcpi)
row.names(ant.extinction_original)<-NULL
row.names(ant.extinction_bcpi)<-NULL
ant_ext_comb<-cbind(ant.extinction_original, ant.extinction_bcpi)
colnames(ant_ext_comb)<-c("extinction time", "lower_extinction", "upper_extinction")
ant_ext_comb$Species<- "ant"

tri.extinction_original<-as.data.frame(tri.extinction_original)
tri.extinction_bcpi<-as.data.frame(tri.extinction_bcpi)
row.names(tri.extinction_original)<-NULL
row.names(tri.extinction_bcpi)<-NULL
tri_ext_comb<-cbind(tri.extinction_original, tri.extinction_bcpi)
colnames(tri_ext_comb)<-c("extinction time", "lower_extinction", "upper_extinction")
tri_ext_comb$Species<- "tri"

cam.extinction_original<-as.data.frame(cam.extinction_original)
cam.extinction_bcpi<-as.data.frame(cam.extinction_bcpi)
row.names(cam.extinction_original)<-NULL
row.names(cam.extinction_bcpi)<-NULL
cam_ext_comb<-cbind(cam.extinction_original, cam.extinction_bcpi)
colnames(cam_ext_comb)<-c("extinction time", "lower_extinction", "upper_extinction")
cam_ext_comb$Species<- "cam"

my_theme<- theme_classic()+
  theme(
    legend.background=element_blank(),
    plot.title = element_text(hjust = 0.5,vjust=-8, size = 28, color="black"),   
    plot.caption = element_text(hjust = 0, face = "italic"), 
    legend.title=element_text(size=20),
    legend.text=element_text(size=20),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x =element_text(size=28,margin=margin(5,8,10,5,"pt")),
    axis.text.y =element_text(size=28, margin=margin(5,8,10,5,"pt")),
    axis.ticks.length = unit(-2, "mm"),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.ticks = element_line(colour = "black", size = 1.1),
    legend.position = c(0.8, 0.6),
    strip.text.x = element_text(
      size = 14,face = "italic"),
    legend.box.background = element_rect(colour = "black"))

pla_ext_comb$Site<-c("P","P","M","M","C","C","R")
pla_ext_comb$Treatment<-c("Current species","Novel species","Current species","Novel species","Current species","Novel species","Current species")
pla_ext_comb[data.plot_pla$Site == "R", "temperature"] = 0
pla_ext_comb[data.plot_pla$Site == "P", "temperature"] =  4.85 
pla_ext_comb[data.plot_pla$Site == "C", "temperature"] = 1.65 
pla_ext_comb[data.plot_pla$Site == "M", "temperature"] =  3.05 

points.pla = pla_ext_comb
points.pla =subset(points.pla, subset = Site == "P"|Site == "M"|Site == "C")

x<-expression(paste(italic("Plantago alpina")))

p1<-ggplot(data = pla.temperature_range, aes(x = temperature, y = extinction_time, color=Treatment)) + scale_color_manual(values=c("gold", "steelblue4"), labels=c("Current species", "Novel species"))+
  scale_fill_manual(values=c("gold", "steelblue4"),  labels=c("Current species", "Novel species"))+
  geom_line(data = pla.temperature_range,alpha=1, size=1.6)

p1<-p1+ geom_errorbar(data=points.pla,aes(x = temperature,ymin = lower_extinction, ymax = upper_extinction, group=Treatment), colour="black", position = position_dodge(0.4),  width = 0.2)+
  geom_point(data = points.pla, aes(x = temperature, y = extinction_time, fill=Treatment, group=Treatment), position = position_dodge(width=0.4),colour="black", shape=21, size = 6.4)+
  annotate("text", label = "*", x = 5.2, y = 10, size = 10, colour = "black")+
  theme(legend.title = element_text(size = 10, face="bold"))+ theme(legend.text = element_text(size = 10))+
  my_theme+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 1100),breaks = scales::pretty_breaks(n = 6))+guides(fill=guide_legend(title="Interacting with:"))+
  scale_y_continuous(trans = "log",breaks = log_breaks(n = 6))

pla<-p1

ant_ext_comb$Site<-c("P","P","M","M","C","C","R")
ant_ext_comb$Treatment<-c("Current species","Novel species","Current species","Novel species","Current species","Novel species","Current species")
ant_ext_comb[data.plot_ant$Site == "R", "temperature"] = 0
ant_ext_comb[data.plot_ant$Site == "P", "temperature"] =  4.85 
ant_ext_comb[data.plot_ant$Site == "C", "temperature"] = 1.65 
ant_ext_comb[data.plot_ant$Site == "M", "temperature"] =  3.05 

points.ant = ant_ext_comb
points.ant =subset(points.ant, subset = Site == "P"|Site == "M"|Site == "C")

x<-expression(paste(italic("Anthyllis alpestris")))

a1<-ggplot(data = ant.temperature_range, aes(x = temperature, y = extinction_time, color=Treatment)) + scale_color_manual(values=c("gold", "steelblue4"), labels=c("Current species", "Novel species"))+
  scale_fill_manual(values=c("gold", "steelblue4"),  labels=c("Current species", "Novel species"))+
  geom_line(data = ant.temperature_range,alpha=1, size=1.6)

a1<-p1+ geom_errorbar(data=points.ant,aes(x = temperature,ymin = lower_extinction, ymax = upper_extinction, group=Treatment), colour="black", position = position_dodge(0.4),  width = 0.2)+
  geom_point(data = points.ant, aes(x = temperature, y = extinction_time, fill=Treatment, group=Treatment), position = position_dodge(width=0.4),colour="black", shape=21, size = 6.4)+
  annotate("text", label = "*", x = 5.2, y = 10, size = 10, colour = "black")+
  theme(legend.title = element_text(size = 10, face="bold"))+ theme(legend.text = element_text(size = 10))+
  my_theme+
  #scale_y_continuous(expand = c(0, 0),limits = c(0, 1100),breaks = scales::pretty_breaks(n = 6))+guides(fill=guide_legend(title="Interacting with:"))+
  scale_y_continuous(trans = "log",breaks = log_breaks(n = 6))

ant<-p1

a1<-ggplot(data = ant.temperature_range, aes(x = temperature, y = extinction, color=Treatment)) + scale_color_manual(values=c("gold", "steelblue4"), labels=c("Current species", "Novel species"))+
  scale_fill_manual(values=c("gold", "steelblue4"),  labels=c("Current species", "Novel species"))+
  geom_line(data = data.ant_short_2,alpha=1, size=1.6)+
  geom_errorbar(data=points.ant,aes(ymin = ci_min_bcpi_del_ext, ymax = ci_max_bcpi_del_ext, group=Treatment), colour="black", position = position_dodge(0.4),  width = 0.2)+
  geom_point(data = points.ant, aes(x = temperature, y = extinction, fill=Treatment, group=Treatment), position = position_dodge(width=0.4),colour="black", shape=21, size = 6.4)+
  annotate("text", label = "*", x = 5.2, y = 7, size = 10, colour = "black")
theme(legend.title = element_text(size = 10, face="bold"))+ theme(legend.text = element_text(size = 10))+
  my_theme+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1400),breaks = scales::pretty_breaks(n = 6))+guides(fill=guide_legend(title="Interacting with:"))+
  scale_y_continuous(trans = "log", breaks = log_breaks(n = 6))
ant<-a1


points.tri = subset(combined_extinction, subset = species == "tri")
x<-expression(paste(italic("Trifolium badium")))
tri.temperature_range

t1<-ggplot(data = tri.temperature_range, aes(x = temperature, y = extinction, color=Treatment)) + scale_color_manual(values=c("gold", "steelblue4"), labels=c("Current species", "Novel species"))+
  scale_fill_manual(values=c("gold", "steelblue4"),  labels=c("Current species", "Novel species"))+
  geom_line(data = data.tri_short_2,alpha=1, size=1.6)+
  geom_errorbar(data=points.tri,aes(ymin = ci_min_bcpi_del_ext, ymax = ci_max_bcpi_del_ext, group=Treatment), colour="black", position = position_dodge(0.4),  width = 0.2)+
  geom_point(data = points.tri, aes(x = temperature, y = extinction, fill=Treatment, group=Treatment), position = position_dodge(width=0.4),colour="black", shape=21, size = 6.4)+
  annotate("text", label = "*", x = 5.2, y = 4, size = 10, colour = "black")
theme(legend.title = element_text(size = 10, face="bold"))+ theme(legend.text = element_text(size = 10))+
  my_theme+
  scale_y_continuous(expand = c(0, 0),limits = c(0, 1400),breaks = scales::pretty_breaks(n = 6))+guides(fill=guide_legend(title="Interacting with:"))+
  scale_y_continuous(trans = "log", breaks = log_breaks(n = 6))
tri<-t1

points.cam = subset(combined_extinction, subset = species == "cam")
x<-expression(paste(italic("Campanula scheuchzeri")))
cam.temperature_range

c1<-ggplot(data = cam.temperature_range, aes(x = temperature, y = extinction, color=Treatment)) + scale_color_manual(values=c("gold", "steelblue4"), labels=c("Current species", "Novel species"))+
  scale_fill_manual(values=c("gold", "steelblue4"),  labels=c("Current species", "Novel species"))+
  geom_line(data = data.cam_short_2,alpha=1, size=1.6)+
  geom_errorbar(data=points.cam,aes(ymin = ci_min_bcpi_del_ext, ymax = ci_max_bcpi_del_ext, group=Treatment), colour="black", position = position_dodge(0.4),  width = 0.2)+
  geom_point(data = points.cam, aes(x = temperature, y = extinction, fill=Treatment, group=Treatment), position = position_dodge(width=0.4),colour="black", shape=21, size = 6.4)+
  theme(legend.title = element_text(size = 10, face="bold"))+ theme(legend.text = element_text(size = 10))+
  my_theme+
  scale_y_continuous(trans = "log", breaks = log_breaks(n = 6),limits = c(2, 260))

cam<-c1

extinction_plot <-grid.arrange(arrangeGrob(pla, ant, tri, cam,
                                           nrow = 2,
                                           left = textGrob("Extinction time", rot = 90,gp = gpar(cex = 2.8)),
                                           bottom = textGrob("Climate warming (°C)",gp = gpar(cex = 2.8))))