##################################################################################################################
## LTREs
##################################################################################################################
## we need data obtained from "Population growth" script
## parameters_p
## max_min_values_all_p

##################################################################################################################
## For explanation of sites and treatments, see script "Population growth"
##################################################################################################################
## Line 15-299 "Estimate relative LTRE contribution" - Line 300-416 "Plot relative LTRE contribution"
##################################################################################################################
## 1. Estimate relative LTRE contribution
##################################################################################################################
## subset for species

#species = "pla"

#species = "ant"

#species = "tri"

#species = "cam"

parameters_p = subset(parameters, subset = Species == species)

max_min_values_all_p<-subset(max_min_values_all, subset = Species == species)
min.p<-as.data.frame(max_min_values_all_p[,3])
min.p<-t(min.p)
min.p<-as.data.frame(min.p)
min.p$parameter<-"min.size"
min.p$Species<-species
row.names(min.p)<-NULL
colnames(min.p)<-colnames(parameters_p)

max.p<-as.data.frame(max_min_values_all_p[,4])
max.p<-t(max.p)
max.p<-as.data.frame(max.p)
max.p$parameter<-"max.size"
max.p$Species<-species
row.names(max.p)<-NULL
colnames(max.p)<-colnames(parameters_p)

parameters_p<-rbind(parameters_p,min.p,max.p)

## High elevation with current species serves as "control"  
control.mpar = parameters_p$R.CS

## subset for treatment lines xx -xx for each treatment to obtain LTREs for each treatment:
## P.CS --> 1400 m Current species
## P.NS --> 1400 m Novel species
## M.CS --> 1750 m Current species
## M.NS --> 1750 m Novel species
## C.CS --> 1950 m Current species
## C.NS --> 1950 m Novel species

mid_treatment="C.NS"
treatment.mpar=parameters_p$C.NS

## calculate mid-parameters
midparams.Em<-(control.mpar[1:11]+treatment.mpar[1:11])/2 

grow.int=midparams.Em[1]
grow.z=midparams.Em[2]
grow.sd =midparams.Em[3]
surv.int=midparams.Em[4]
surv.z=midparams.Em[5]
flow.int=midparams.Em[6]
flow.z =midparams.Em[7]
seed.mean  = midparams.Em[8]
rcsz.int =midparams.Em[9]
rcsz.sd = midparams.Em[10]
p.r = midparams.Em[11]

names(m.par) <-c("grow.int","grow.z","grow.sd","surv.int","surv.z","flow.int","flow.z","seed.mean","rcsz.int","rcsz.sd","p.r")

midparams.Em = c(grow.int,
                 grow.z,
                 grow.sd,
                 surv.int,
                 surv.z,
                 flow.int,
                 flow.z,
                 seed.mean,
                 rcsz.int, 
                 rcsz.sd,
                 p.r )

names(midparams.Em) <- names(m.par)

## run IPM
G_z1z <- function (l1, l, midparams.Em) {
  dnorm(l1, mean = midparams.Em["grow.int"]+ midparams.Em["grow.z"]*l, sd = midparams.Em["grow.sd"])
  
}

s_z <- function(l, midparams.Em)
{
  linear.p <- midparams.Em["surv.int"] + midparams.Em["surv.z"] * l 
  p <- 1/(1+exp(-linear.p))                            
  return(p)
}

P_z1z <- function (l1, l, midparams.Em) {
  
  return((s_z(l, midparams.Em) * G_z1z(l1, l, midparams.Em)))
}

F_z1z <- function (l1, l, midparams.Em) {
  return(
    (1/(1+1/exp(midparams.Em["flow.int"] + midparams.Em["flow.z"]*l))) *
      midparams.Em["seed.mean"]*
      midparams.Em["p.r"] * 
      dnorm(l1, mean = midparams.Em["rcsz.int"], sd = midparams.Em["rcsz.sd"]))
}


nBigMatrix <- 250

##select the min and max size based on the control and treatment min/max sizes

min.size <- min(treatment.mpar[12], control.mpar[12], na.rm = TRUE)
max.size<- max(treatment.mpar[13], control.mpar[13], na.rm = TRUE)

U <- max.size *1.1
L<-min.size

h <- (U-L)/nBigMatrix 
meshpts <- L + (1:nBigMatrix)*h - h/2

mk_K <- function(nBigMatrix, midparams.Em, L, U) {
  h <- (U - L)/nBigMatrix
  meshpts <- L + ((1:nBigMatrix) - 1/2) * h
  
  
  P <- h * (outer(meshpts, meshpts, P_z1z, midparams.Em = midparams.Em))
  F. <- h * (outer(meshpts, meshpts, F_z1z, midparams.Em = midparams.Em))
  K <- P + F.
  
  
  return(list(K = K, meshpts = meshpts, P = P, F. = F.))
}

LTRE.K.mid <- mk_K(nBigMatrix, midparams.Em, L,U)
midpoint.lam<-Re(eigen(LTRE.K.mid$K)$values[1])


## calculate the sensitivities of th mid-kernel for each vital rate with a pertubation of 0.001

perturbation<-0.001
LTRE.sens<-c()

for(j in 1:11){ 
  LTRE.mid.mpar<-midparams.Em 
  LTRE.mid.mpar[j]<-midparams.Em[j]+perturbation 
  LTRE.K.control<-mk_K(nBigMatrix, LTRE.mid.mpar, L, U)
  LTRE.lam<-Re(eigen(LTRE.K.control$K)$values[1]) 
  LTRE.sens[j]<-(LTRE.lam-midpoint.lam)/perturbation
}

## calculate the LTREs for each vitalrate --> following the order as vital rates are arranged in the "parameters_p" file

growth.ltre<-sum(na.omit(LTRE.sens[1:3]*(treatment.mpar[1:3]-control.mpar[1:3])))
surv.ltre<-sum(na.omit(LTRE.sens[4:5]*(treatment.mpar[4:5]-control.mpar[4:5])))
flow.ltre<-sum(na.omit(LTRE.sens[6:7]*(treatment.mpar[6:7]-control.mpar[6:7])))
seed.ltre<-sum(na.omit(LTRE.sens[8]*(treatment.mpar[8]-control.mpar[8])))
rcsz.ltre<-sum(na.omit(LTRE.sens[9:10]*(treatment.mpar[9:10]-control.mpar[9:10])))
recruit.ltre<-sum(na.omit(LTRE.sens[11]*(treatment.mpar[11]-control.mpar[11])))

ltre <- rbind( "Survival" = surv.ltre, "Flowering" = flow.ltre, "Growth"= growth.ltre, "Recr.sz" =rcsz.ltre, "Seeds" =seed.ltre, "Recruitment"=recruit.ltre)
ltre<-as.data.frame(ltre)
colnames(ltre)<-c("value")

## calculate relative LTRE contribution

ltre$value_pos<-abs(ltre$value)
ltre_sum_positive<-sum(ltre[1,2],ltre[2,2],ltre[3,2],ltre[4,2],ltre[5,2],ltre[6,2])
ltre$sum_positive<-ltre_sum_positive
ltre$normalize<- (ltre$value_pos/ltre$sum_positive)
ltre$true_cont<-ifelse(ltre$value< 0, (ltre$normalize*-1), ltre$normalize)
ltre$Species = species
name <- paste(species, "ltre", sep=".")
name <- paste(name, mid_treatment, sep=".")

assign(name, ltre)

## When produced this for all treatments for a species you will have the files:
#pla.C.CS.ltre
#pla.M.CS.ltre
#pla.P.CS.ltre
#pla.C.NS.ltre
#pla.M.NS.ltre
#pla.P.NS.ltre

#ant.C.CS.ltre
#ant.M.CS.ltre
#ant.P.CS.ltre
#ant.C.NS.ltre
#ant.M.NS.ltre
#ant.P.NS.ltre

#tri.C.CS.ltre
#tri.M.CS.ltre
#tri.P.CS.ltre
#tri.C.NS.ltre
#tri.M.NS.ltre
#tri.P.NS.ltre

#cam.C.CS.ltre
#cam.M.CS.ltre
#cam.P.CS.ltre
#cam.C.NS.ltre
#cam.M.NS.ltre
#cam.P.NS.ltre

#### combine datasets

## P. alpina

ltre_pla_CS<-rbind(pla.ltre.P.CS,pla.ltre.M.CS,pla.ltre.C.CS)
ltre_pla_CS$Treatment<-"Current species"
ltre_pla_CS$Site<-c("P","P","P","P","P","P","M","M","M","M","M","M","C","C","C","C","C","C")
ltre_pla_CS$var<-c( "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment")

row.names(ltre_pla_CS)<-NULL

ltre_pla_NS<-rbind(pla.ltre.P.NS,pla.ltre.M.NS,pla.ltre.C.NS)
ltre_pla_NS$Treatment<-"Novel species"
ltre_pla_NS$Site<-c("P","P","P","P","P","P","M","M","M","M","M","M","C","C","C","C","C","C")
ltre_pla_NS$var<-c( "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment")

row.names(ltre_pla_NS)<-NULL

ltre_pla<-rbind(ltre_pla_CS,ltre_pla_NS)

## A. alpestris

ltre_ant_CS<-rbind(ant.ltre.P.CS,ant.ltre.M.CS,ant.ltre.C.CS)
ltre_ant_CS$Treatment<-"Current species"
ltre_ant_CS$Site<-c("P","P","P","P","P","P","M","M","M","M","M","M","C","C","C","C","C","C")
ltre_ant_CS$var<-c( "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment")
row.names(ltre_ant_CS)<-NULL

ltre_ant_NS<-rbind(ant.ltre.P.NS,ant.ltre.M.NS,ant.ltre.C.NS)
ltre_ant_NS$Treatment<-"Novel species"
ltre_ant_NS$Site<-c("P","P","P","P","P","P","M","M","M","M","M","M","C","C","C","C","C","C")
ltre_ant_NS$var<-c( "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment")
row.names(ltre_ant_NS)<-NULL

ltre_ant<-rbind(ltre_ant_CS,ltre_ant_NS)

## T. badium

ltre_tri_CS<-rbind(tri.ltre.P.CS,tri.ltre.M.CS,tri.ltre.C.CS)
ltre_tri_CS$Treatment<-"Current species"
ltre_tri_CS$Site<-c("P","P","P","P","P","P","M","M","M","M","M","M","C","C","C","C","C","C")
ltre_tri_CS$var<-c( "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment")
row.names(ltre_tri_CS)<-NULL

ltre_tri_NS<-rbind(tri.ltre.P.NS,tri.ltre.M.NS,tri.ltre.C.NS)
ltre_tri_NS$Treatment<-"Novel species"
ltre_tri_NS$Site<-c("P","P","P","P","P","P","M","M","M","M","M","M","C","C","C","C","C","C")
ltre_tri_NS$var<-c( "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment")
row.names(ltre_tri_NS)<-NULL

ltre_tri<-rbind(ltre_tri_CS,ltre_tri_NS)

## C. scheucheri

ltre_cam_CS<-rbind(cam.ltre.P.CS,cam.ltre.M.CS,cam.ltre.C.CS)
ltre_cam_CS$Treatment<-"Current species"
ltre_cam_CS$Site<-c("P","P","P","P","P","P","M","M","M","M","M","M","C","C","C","C","C","C")
ltre_cam_CS$var<-c( "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment")
row.names(ltre_cam_CS)<-NULL

ltre_cam_NS<-rbind(cam.ltre.P.NS,cam.ltre.M.NS,cam.ltre.C.NS)
ltre_cam_NS$Treatment<-"Novel species"
ltre_cam_NS$Site<-c("P","P","P","P","P","P","M","M","M","M","M","M","C","C","C","C","C","C")
ltre_cam_NS$var<-c( "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment",
                    "Survival","Flowering", "Growth", "Recr.sz" , "Seeds", "Recruitment")
row.names(ltre_cam_NS)<-NULL

ltre_cam<-rbind(ltre_cam_CS,ltre_cam_NS)

##################################################################################################################
## 2. Plot relative LTRE contribution
##################################################################################################################

##packages required
library(Rmisc)
library(gridExtra)
library(grid)

my_theme<- theme_classic()+
  theme(
    legend.background=element_blank(),
    plot.title = element_text(hjust = 0.5,vjust=-1, size = 22, color="black"),
    plot.caption = element_text(hjust = 0, face = "italic"), 
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x =element_text(size=28, angle = 90,margin=margin(5,8,10,5,"pt")),
    axis.text.y =element_text(size=32, margin=margin(5,8,10,5,"pt")),
    axis.ticks.length = unit(-2, "mm"),
    axis.line = element_line(colour = 'black', size = 1.1),
    axis.ticks = element_line(colour = "black", size = 1.1),
    legend.position = c(0.22, 0.22),
    strip.text.x = element_text(
      size = 14,face = "italic"))

## P. alpina

x<-expression(paste(italic("Plantago alpina")))

ltre_pla$var = factor(ltre_pla$var, levels = c("Recruitment", "Recr.sz", "Growth", "Flowering", "Seeds","Survival"))
ltre_pla$Site = factor(ltre_pla$Site, levels = c("P", "M", "C"))
cdata<- summarySE(ltre_pla, measurevar="true_cont", groupvars=c("Treatment", "var"),na.rm=TRUE)

dogde<-0.7

p1 <- ggplot(data=cdata, aes(x=var, y=true_cont, fill=Treatment)) +
  geom_hline(yintercept=0, size=1)+
  geom_bar(stat="identity", color="black", width=0.7, position=position_dodge(0.7))+
  scale_fill_manual(values=c("gold", "steelblue4"))+theme_light()+theme(axis.text.x=element_blank())+
  geom_point(data=ltre_pla,size=2.4, aes(y=true_cont, x=var,  fill=Treatment, group=Treatment, shape=Site), position = position_dodge(0.7))+
  scale_fill_manual(values=c("gold", "steelblue4"))+scale_shape_manual(values=c(22, 21, 24))+
  geom_errorbar(aes(ymin = true_cont-se, ymax = true_cont+se), width=0.4,colour = "black", position = position_dodge(0.7))+
  scale_y_continuous(limits = c(-1, 0.8) ,breaks = scales::pretty_breaks(n = 6))+theme(strip.text.x = element_text(
    size = 10,face = "italic"), panel.grid.major = element_blank())+ ggtitle(x)+
  my_theme+guides(fill=FALSE, color=FALSE, shape=FALSE)

pla<-p1

## A. alpestris

x<-expression(paste(italic("Anthyllis alpestris")))

ltre_ant$var = factor(ltre_ant$var, levels = c("Recruitment", "Recr.sz", "Growth", "Flowering", "Seeds","Survival"))
ltre_ant$Site = factor(ltre_ant$Site, levels = c("P", "M", "C"))
cdata<- summarySE(ltre_ant, measurevar="true_cont", groupvars=c("Treatment", "var"),na.rm=TRUE)

a1 <- ggplot(data=cdata, aes(x=var, y=true_cont, fill=Treatment)) +
  geom_hline(yintercept=0, size=1)+
  geom_bar(stat="identity", color="black", width=0.7, position=position_dodge(0.7))+
  scale_fill_manual(values=c("gold", "steelblue4"))+theme_light()+theme(axis.text.x=element_blank())+
  geom_point(data=ltre_ant,size=2.4, aes(y=true_cont, x=var,  fill=Treatment, group=Treatment, shape=Site), position = position_dodge(0.7))+
  scale_fill_manual(values=c("gold", "steelblue4"))+scale_shape_manual(values=c(22, 21, 24))+
  geom_errorbar(aes(ymin = true_cont-se, ymax = true_cont+se), width=0.4,colour = "black", position = position_dodge(0.7))+
  scale_y_continuous(limits = c(-1, 0.8) ,breaks = scales::pretty_breaks(n = 6))+theme(strip.text.x = element_text(
    size = 10,face = "italic"), panel.grid.major = element_blank())+ ggtitle(x)+
  my_theme+guides(fill=FALSE, color=FALSE, shape=FALSE)

ant<-a1

## T. badium

x<-expression(paste(italic("Trifolium badium")))

ltre_tri$var = factor(ltre_tri$var, levels = c("Recruitment", "Recr.sz", "Growth", "Flowering", "Seeds","Survival"))
ltre_tri$Site = factor(ltre_tri$Site, levels = c("P", "M", "C"))
cdata<- summarySE(ltre_tri, measurevar="true_cont", groupvars=c("Treatment", "var"),na.rm=TRUE)

t1 <- ggplot(data=cdata, aes(x=var, y=true_cont, fill=Treatment)) +
  geom_hline(yintercept=0, size=1)+
  geom_bar(stat="identity", color="black", width=0.7, position=position_dodge(0.7))+
  scale_fill_manual(values=c("gold", "steelblue4"))+theme_light()+theme(axis.text.x=element_blank())+
  geom_point(data=ltre_tri,size=2.4, aes(y=true_cont, x=var,  fill=Treatment, group=Treatment, shape=Site), position = position_dodge(0.7))+
  scale_fill_manual(values=c("gold", "steelblue4"))+scale_shape_manual(values=c(22, 21, 24))+
  geom_errorbar(aes(ymin = true_cont-se, ymax = true_cont+se), width=0.4,colour = "black", position = position_dodge(0.7))+
  scale_y_continuous(limits = c(-1, 0.8) ,breaks = scales::pretty_breaks(n = 6))+theme(strip.text.x = element_text(
    size = 10,face = "italic"), panel.grid.major = element_blank())+ ggtitle(x)+
  my_theme+guides(fill=FALSE, color=FALSE, shape=FALSE)

tri<-t1

## C. scheuchzeri

x<-expression(paste(italic("Campanula scheuchzeri")))

ltre_cam$var = factor(ltre_cam$var, levels = c("Recruitment", "Recr.sz", "Growth", "Flowering", "Seeds","Survival"))
ltre_cam$Site = factor(ltre_cam$Site, levels = c("P", "M", "C"))
cdata<- summarySE(ltre_cam, measurevar="true_cont", groupvars=c("Treatment", "var"),na.rm=TRUE)

dogde<-0.7

c1 <- ggplot(data=cdata, aes(x=var, y=true_cont, fill=Treatment)) +
  geom_hline(yintercept=0, size=1)+
  geom_bar(stat="identity", color="black", width=0.7, position=position_dodge(0.7))+
  scale_fill_manual(values=c("gold", "steelblue4"))+theme_light()+theme(axis.text.x=element_blank())+
  geom_point(data=ltre_cam,size=2.4, aes(y=true_cont, x=var,  fill=Treatment, group=Treatment, shape=Site), position = position_dodge(0.7))+
  scale_fill_manual(values=c("gold", "steelblue4"))+scale_shape_manual(values=c(22, 21, 24))+
  geom_errorbar(aes(ymin = true_cont-se, ymax = true_cont+se), width=0.4,colour = "black", position = position_dodge(0.7))+
  scale_y_continuous(limits = c(-1, 0.8) ,breaks = scales::pretty_breaks(n = 6))+theme(strip.text.x = element_text(
    size = 10,face = "italic"), panel.grid.major = element_blank())+ ggtitle(x)+
  my_theme+guides(fill=FALSE, color=FALSE, shape=FALSE)

cam<-c1

all <-grid.arrange(arrangeGrob(pla, ant, tri, cam,
                               nrow = 2, #
                               left = textGrob("Relative LTRE contribution", rot = 90,gp = gpar(cex = 2.8))))
