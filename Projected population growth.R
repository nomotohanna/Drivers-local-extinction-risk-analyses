##################################################################################################################
## POPULATION GROWTH
##################################################################################################################
## Packages required
library(MASS)
library(data.table)

## Read in data
data = fread("data_fit_models.csv", data.table = F)

##################################################################################################################
#### Explanation of sites and treatments throughout scripts
##################################################################################################################

################################# Sites ##################################

# Site P = 1400 m
# Site M = 1750 m
# Site C = 1950 m
# Site R = 2200 m

######################## Treatment abbreviations ########################

# P.CS = 1400 m site, Current species
# P.NS = 1400 m site, Novel species
# M.CS = 1750 m site, Current species
# M.NS = 1750 m site, Novel species
# C.CS = 1950 m site, Current species
# C.NS = 1950 m site, Novel species
# R.CS = 2200 m site, Current species

######################## Species abbreviations ########################

# "pla" = Plantago alpina
# "ant" = Anthyllis alpestris
# "tri" = Trifolium badium
# "cam" = Campanula scheuchzeri

##################################################################################################################
## ! ## Note where scripts differ slightly between species:

## lines 294-322 - P.alpina and T.badium have no recruit size for site "P"
## lines 357 - 401 C. scheuchzeri recruit data fitted as constants, 
## lines 451-544 C. scheuchzeri minimum size based on recruit data from other species
## Bootstrap:
## lines 970-1002 - P.alpina and T.badium have no recruit size for site "P", 
## lines 11038-1039 - P.alpina and T.badium have no recruit size for site "P"
## lines 1529-1534 - C. scheuchzeri recruit data fitted as constants

## For each section of the script;
## Line 89-625 "Fitting models" - Line ,628- 773 "Estimating population growth"
## Line 775-1665 "Bootstrap population growth" -Line 775-1665 "Pairwise comparisons"
## Line 1728-1872 "Plot population growth"
## Sections should one species at the time, thus the subsetting for species in beginning of each section
##################################################################################################################

##################################################################################################################
## Organize data
##################################################################################################################

data$Site = factor(data$Site, levels = c("P", "M", "C", "R"))
data$Treatment = factor(data$Treatment, levels = c( "Current species",  "Novel species","Alpine soil","Low elevation soil"))
data$Species = factor(data$Species, levels = c( "pla","ant","tri","cam"))
data$Site <-as.factor(data$Site)
data$Treatment <-as.factor(data$Treatment)
data$size<- as.numeric(data$size)
data$Repr<- as.numeric(data$Repr)
data$Surv<- as.numeric(data$Surv)
data$size_1<- as.numeric(data$size_1)

## Log-transform size
data$size<- log(data$size)
data$size_1<- log(data$size_1)

##Subset for only community treatments 
data = subset(data, subset = Treatment == "Current species"|Treatment == "Novel species")

##### do one species at the time

species = "pla"
#species = "tri"
#species = "ant"
#species = "cam"

data.p = subset(data, subset = Species == species)
head(data.p)

##################################################################################################################
##   1.       Fitting models to obtain parameters to implement in IPM 
##################################################################################################################

#################################################################################
## 1. GROWTH
#################################################################################

growth_model<-lm(size_1~size*Site*Treatment, data=data.p)

#################################################################################
## Extract model coefficients
#################################################################################

coef_1<-growth_model$coefficients

growth_int<- data.frame(
  P.CS = c(coef_1["(Intercept)"]),
  P.NS = c(sum(coef_1["(Intercept)"]+coef_1["TreatmentNovel species"])),
  M.CS = c(sum(coef_1["(Intercept)"]+coef_1["SiteM"])),
  M.NS = c(sum(coef_1["(Intercept)"]+coef_1["TreatmentNovel species"]
               +coef_1["SiteM"]+coef_1["SiteM:TreatmentNovel species"])),
  C.CS = c(sum(coef_1["(Intercept)"]+coef_1["SiteC"])), 
  C.NS = c(sum(coef_1["(Intercept)"]+coef_1["TreatmentNovel species"]
                   +coef_1["SiteC"]+coef_1["SiteC:TreatmentNovel species"])),
  R.CS = c(sum(coef_1["(Intercept)"]+coef_1["SiteR"])))

growth_int$parameter <-c("grow.int")

growth_slope<- data.frame(
  P.CS = c(coef_1["size"]),
  P.NS =  c(sum(coef_1["size"]+coef_1["size:TreatmentNovel species"])),
  M.CS = c(sum(coef_1["size"]+coef_1["size:SiteM"])),
  M.NS = c(sum(coef_1["size"]+coef_1["size:TreatmentNovel species"],
               coef_1["size:SiteM"]+coef_1["size:SiteM:TreatmentNovel species"])),
  C.CS = c(sum(coef_1["size"]+coef_1["size:SiteC"])),
  C.NS = c(sum(coef_1["size"]+coef_1["size:TreatmentNovel species"],
               coef_1["size:SiteC"]+coef_1["size:SiteC:TreatmentNovel species"])),
  R.CS = c(sum(coef_1["size"]+coef_1["size:SiteR"])))

growth_slope$parameter <-c("grow.z")

## SD around fitted growth slope

p.cs<-lm(size_1~size, data=data.p,subset = Site == "P" & Treatment=="Current species")
P.CS<-summary(p.cs)$sigma

m.cs<-lm(size_1~size, data=data.p,subset = Site == "M" & Treatment=="Current species")
M.CS<-summary(m.cs)$sigma

c.cs<-lm(size_1~size, data=data.p,subset = Site == "C" & Treatment=="Current species")
C.CS<-summary(c.cs)$sigma

r.cs<-lm(size_1~size, data=data.p,subset = Site == "R" & Treatment=="Current species")
R.CS<-summary(r.cs)$sigma

p.ns<-lm(size_1~size, data=data.p,subset = Site == "P" & Treatment=="Novel species")
P.NS<-summary(p.ns)$sigma

m.ns<-lm(size_1~size, data=data.p,subset = Site == "M" & Treatment=="Novel species")
M.NS<-summary(m.ns)$sigma

c.ns<-lm(size_1~size, data=data.p,subset = Site == "C" & Treatment=="Novel species")
C.NS<-summary(c.ns)$sigma

growth_sd<- data.frame(P.CS,P.NS,M.CS,M.NS,C.CS,C.NS,R.CS)

growth_sd$parameter <-c("grow.sd")

#################################################################################
## 2. Survival
#################################################################################

survival_model<-glm(Surv~size*Site*Treatment, data=data.p, family="binomial")

#################################################################################
## Extract model coefficients
#################################################################################

coef_1<-survival_model$coefficients

survival_int<- data.frame(
  P.CS = c(coef_1["(Intercept)"]),
  P.NS = c(sum(coef_1["(Intercept)"]+coef_1["TreatmentNovel species"])),
  M.CS = c(sum(coef_1["(Intercept)"]+coef_1["SiteM"])),
  M.NS = c(sum(coef_1["(Intercept)"]+coef_1["TreatmentNovel species"]
               +coef_1["SiteM"]+coef_1["SiteM:TreatmentNovel species"])),
  C.CS = c(sum(coef_1["(Intercept)"]+coef_1["SiteC"])), 
  C.NS = c(sum(coef_1["(Intercept)"]+coef_1["TreatmentNovel species"]
               +coef_1["SiteC"]+coef_1["SiteC:TreatmentNovel species"])),
  R.CS = c(sum(coef_1["(Intercept)"]+coef_1["SiteR"])))

survival_int$parameter <-c("surv.int")

survival_slope<- data.frame(
  P.CS = c(coef_1["size"]),
  P.NS =  c(sum(coef_1["size"]+coef_1["size:TreatmentNovel species"])),
  M.CS = c(sum(coef_1["size"]+coef_1["size:SiteM"])),
  M.NS = c(sum(coef_1["size"]+coef_1["size:TreatmentNovel species"],
               coef_1["size:SiteM"]+coef_1["size:SiteM:TreatmentNovel species"])),
  C.CS = c(sum(coef_1["size"]+coef_1["size:SiteC"])),
  C.NS = c(sum(coef_1["size"]+coef_1["size:TreatmentNovel species"],
               coef_1["size:SiteC"]+coef_1["size:SiteC:TreatmentNovel species"])),
  R.CS = c(sum(coef_1["size"]+coef_1["size:SiteR"])))

survival_slope$parameter <-c("surv.z")

#################################################################################
## 3. Flowering
#################################################################################

flowering_model<-glm(Repr~size*Site*Treatment, data=data.p, family="binomial")
summary(flowering_model)

#################################################################################
## Extract model coefficients
#################################################################################

coef_1<-flowering_model$coefficients

flowering_int<- data.frame(
  P.CS = c(coef_1["(Intercept)"]),
  P.NS = c(sum(coef_1["(Intercept)"]+coef_1["TreatmentNovel species"])),
  M.CS = c(sum(coef_1["(Intercept)"]+coef_1["SiteM"])),
  M.NS = c(sum(coef_1["(Intercept)"]+coef_1["TreatmentNovel species"]
               +coef_1["SiteM"]+coef_1["SiteM:TreatmentNovel species"])),
  C.CS = c(sum(coef_1["(Intercept)"]+coef_1["SiteC"])), 
  C.NS = c(sum(coef_1["(Intercept)"]+coef_1["TreatmentNovel species"]
               +coef_1["SiteC"]+coef_1["SiteC:TreatmentNovel species"])),
  R.CS = c(sum(coef_1["(Intercept)"]+coef_1["SiteR"])))

flowering_int$parameter <-c("flow.int")

flowering_slope<- data.frame(
  P.CS = c(coef_1["size"]),
  P.NS =  c(sum(coef_1["size"]+coef_1["size:TreatmentNovel species"])),
  M.CS = c(sum(coef_1["size"]+coef_1["size:SiteM"])),
  M.NS = c(sum(coef_1["size"]+coef_1["size:TreatmentNovel species"],
               coef_1["size:SiteM"]+coef_1["size:SiteM:TreatmentNovel species"])),
  C.CS = c(sum(coef_1["size"]+coef_1["size:SiteC"])),
  C.NS = c(sum(coef_1["size"]+coef_1["size:TreatmentNovel species"],
               coef_1["size:SiteC"]+coef_1["size:SiteC:TreatmentNovel species"])),
  R.CS = c(sum(coef_1["size"]+coef_1["size:SiteR"])))

flowering_slope$parameter <-c("flow.z")

#################################################################################
## 4. Seed production
#################################################################################

data.seeds = subset(data.p, subset = Seeds < 1000)

seed_model<-lm(Seeds ~ Site*Treatment, data=data.seeds)

seed_production <- data.frame(P.CS = coef(seed_model)[1],
                              P.NS = (coef(seed_model)[1]+coef(seed_model)[5]),
                              M.CS = (coef(seed_model)[1]+coef(seed_model)[2]),
                              M.NS = (coef(seed_model)[1]+coef(seed_model)[2]+coef(seed_model)[5]+coef(seed_model)[6]),
                              C.CS = (coef(seed_model)[1]+coef(seed_model)[3]),  
                              C.NS = (coef(seed_model)[1]+coef(seed_model)[3]+coef(seed_model)[5]+coef(seed_model)[7]),
                              R.CS = (coef(seed_model)[1]+coef(seed_model)[4]))
                             
seed_production$parameter <-c("seed.mean")

#################################################################################
## 4. Recruitment
#################################################################################

data.recruitment = fread("data_recruitment.csv", data.table = F)
data.recruitment$Site = factor(data.recruitment$Site, levels = c("P", "M", "C", "R"))
data.recruitment$Treatment = factor(data.recruitment$Treatment, levels = c("Current species", "Novel species","Alpine soil","Low elevation soil"))
data.recruitment = subset(data.recruitment, subset = Treatment == "Current species"|Treatment == "Novel species")
data.recruitment = subset(data.recruitment, subset = Species == species)

recruitment_model<-glm(recruitment ~ Site*Treatment, data=data.recruitment,family="binomial",na.action=na.exclude)

recruitment <- data.frame(P.CS = coef(recruitment_model)[1],
                          P.NS = (coef(recruitment_model)[1]+coef(recruitment_model)[5]),
                  M.CS = (coef(recruitment_model)[1]+coef(recruitment_model)[2]),
                  M.NS = (coef(recruitment_model)[1]+coef(recruitment_model)[2]+coef(recruitment_model)[5]+coef(recruitment_model)[6]),
                  C.CS = (coef(recruitment_model)[1]+coef(recruitment_model)[3]), 
                  C.NS = (coef(recruitment_model)[1]+coef(recruitment_model)[3]+coef(recruitment_model)[5]+coef(recruitment_model)[7]),
                  R.CS = (coef(recruitment_model)[1]+coef(recruitment_model)[4]))

recruitment <- 1/(1+exp(-recruitment))
recruitment$parameter <-c("p.r")

#################################################################################
## 5. Recruit size
#################################################################################

data.recruit_size = fread("data_recruit_size.csv", data.table = F)
data.recruit_size$Site = factor(data.recruit_size$Site, levels = c("P", "M", "C", "R"))
data.recruit_size$Treatment = factor(data.recruit_size$Treatment, levels = c("Current species", "Novel species","Alpine soil","Low elevation soil"))
data.recruit_size = subset(data.recruit_size, subset = Species == species)

##log transform size
data.recruit_size$size<- log(data.recruit_size$size)

###################################################################################################################
## ! ## Note that recruit size is fitted per site across treatment (see material method + supporting information)
###################################################################################################################

recruit_size_model<-lm(size ~ Site, data=data.recruit_size,na.action=na.exclude)

###############################################################################################################
## ! ##  Special case for P. alpina ("pla") and T.badium ("tri") where no size was obtained for 1400 site 
## ("P"), thus this site is excluded from the model
###############################################################################################################

recruit_size <- data.frame( P.CS =0, ##as not included in model =only NAs
                   P.NS =0, ##as not included in model = only NAs
                   M.CS =  coef(recruit_size_model)[1],
                   M.NS =  coef(recruit_size_model)[1],
                   C.CS = (coef(recruit_size_model)[1]+coef(recruit_size_model)[2]),
                   C.NS = (coef(recruit_size_model)[1]+coef(recruit_size_model)[2]),
                   R.CS = (coef(recruit_size_model)[1]+coef(recruit_size_model)[3]))                

recruit_size$parameter <-c("rcsz.int")

############################################### SD  #######################################################

P.CS<-0 ##as not included in model =only NAs
P.NS<-0 ##as not included in model =only NAs
M.CS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "M"))$sigma
M.NS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "M"))$sigma
C.CS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "C"))$sigma
C.NS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "C"))$sigma
R.CS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "R"))$sigma 

recruit_size_sd<- data.frame(P.CS,P.NS,M.CS, M.NS,C.CS,C.NS,R.CS)

recruit_size_sd$parameter <-c("rcsz.sd")

###############################################################################################################
## for A. alpestris ("ant"), we can obtain data from  1400 site ("P"), thus model looks like this instead:
###############################################################################################################

recruit_size <-data.frame(P.CS = coef(recruit_size_model)[1],
                          P.NS = coef(recruit_size_model)[1],
                          M.CS = (coef(recruit_size_model)[1]+coef(recruit_size_model)[2]),
                          M.NS = (coef(recruit_size_model)[1]+coef(recruit_size_model)[2]),
                          C.CS = (coef(recruit_size_model)[1]+coef(recruit_size_model)[3]),     
                          C.NS = (coef(recruit_size_model)[1]+coef(recruit_size_model)[3]),     
                          R.CS = (coef(recruit_size_model)[1]+coef(recruit_size_model)[4])) 

recruit_size$parameter <-c("rcsz.int")

P.CS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "P"))$sigma
P.NS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "P"))$sigma
M.CS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "M"))$sigma
M.NS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "M"))$sigma
C.CS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "C"))$sigma
C.NS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "C"))$sigma
R.CS<-summary(lm(size ~ 1, data=data.recruit_size,na.action=na.exclude, subset = Site == "R"))$sigma 

recruit_size_sd<- data.frame(P.CS,P.NS,M.CS,M.NS,C.CS,C.NS,R.CS)

recruit_size_sd$parameter <-c("rcsz.sd")

##################################################################################################################################################################

parameters<-rbind(growth_int, growth_slope, growth_sd,survival_int, survival_slope,flowering_int,flowering_slope,seed_production,recruit_size,recruit_size_sd,recruitment)
rownames(parameters)<-NULL

name <- paste(species, "parameters",sep=".")
assign(name, parameters)

###############################################################################################################
## ! ## RUN Campanula scheuchzeri ("cam") last when you have the datasets "pla.parameters","ant.parameters" and 
##"tri.parameters" as we have no data for recruitment or recruit size, where fixed values instead were implemented 
## (see Material and Method and Supporting information). For recruit size we use the minimum recruit size and SD
## based on the other species, thus we extract the parameters from the other species first 
###############################################################################################################

## run until recruitment above line xx", then
recruitment <- data.frame( P.CS = 0.001,
                           M.CS = 0.001,
                           C.CS = 0.001,
                           R.CS = 0.001,
                           
                           P.NS = 0.001,
                           M.NS = 0.001,
                           C.NS = 0.001)

recruitment$parameter <-c("p.r")

#####################################################################################################################################
## Implement recruit size and sd for C.scheuzheri based on the minimum values obtained for the other species (note that site 1400 m,
## "P", is deleted from pla.parameters and tri.parameters as here we had 0 size
#####################################################################################################################################


recruit_size <-data.frame( P.CS = min(pla.parameters[9,c(3:7)],ant.parameters[9,1:7], tri.parameters[9,c(3:7)]),
                           P.NS = min(pla.parameters[9,c(3:7)],ant.parameters[9,1:7], tri.parameters[9,c(3:7)]),
                           M.CS = min(pla.parameters[9,c(3:7)],ant.parameters[9,1:7], tri.parameters[9,c(3:7)]),
                           M.NS = min(pla.parameters[9,c(3:7)],ant.parameters[9,1:7], tri.parameters[9,c(3:7)]),
                           C.CS = min(pla.parameters[9,c(3:7)],ant.parameters[9,1:7], tri.parameters[9,c(3:7)]),
                           C.NS = min(pla.parameters[9,c(3:7)],ant.parameters[9,1:7], tri.parameters[9,c(3:7)]),
                           R.CS = min(pla.parameters[9,c(3:7)],ant.parameters[9,1:7], tri.parameters[9,c(3:7)]))


recruit_size$parameter <-c("rcsz.int")

recruit_size_sd <-data.frame(P.CS = min(pla.parameters[10,c(3:7)],ant.parameters[10,1:7], tri.parameters[10,c(3:7)]),
                             P.NS = min(pla.parameters[10,c(3:7)],ant.parameters[10,1:7], tri.parameters[10,c(3:7)]),
                             M.CS = min(pla.parameters[10,c(3:7)],ant.parameters[10,1:7], tri.parameters[10,c(3:7)]),
                             M.NS = min(pla.parameters[10,c(3:7)],ant.parameters[10,1:7], tri.parameters[10,c(3:7)]),
                             C.CS = min(pla.parameters[10,c(3:7)],ant.parameters[10,1:7], tri.parameters[10,c(3:7)]),
                             C.NS = min(pla.parameters[10,c(3:7)],ant.parameters[10,1:7], tri.parameters[10,c(3:7)]),
                             R.CS = min(pla.parameters[10,c(3:7)],ant.parameters[10,1:7], tri.parameters[10,c(3:7)]))

recruit_size_sd$parameter <-c("rcsz.sd")

parameters<-rbind(growth_int, growth_slope, growth_sd,survival_int, survival_slope,flowering_int,flowering_slope,seed_production,recruit_size,recruit_size_sd,recruitment)
rownames(parameters)<-NULL

name <- paste(species, "parameters",sep=".")
assign(name, parameters)

## replace "NA" seed production in the models with 0s (and also the extremely low values?)
cam.parameters[is.na(cam.parameters)] <- 0

#### Parameters for each species saved:

#pla.parameters --> Plantago alpina parameters
#tri.parameters --> Trifolium badium parameters
#ant.parameters --> Anthyllis alpestris parameters
#cam.parameters --> Campanula scheuchzeri parameters

pla.parameters$Species<-"pla"
ant.parameters$Species<-"ant"
tri.parameters$Species<-"tri"
cam.parameters$Species<-"cam"

parameters<-rbind(pla.parameters,ant.parameters,tri.parameters,cam.parameters)

########################################################################################################################
## Calculate minimum and maximum size
########################################################################################################################

########################################################################################################################
## Minimum size
########################################################################################################################

## Subset for one of the species
species = "pla"
#species = "tri"
#species = "ant"

##do cam last due to special case : see row xx
#species = "cam"

data.p = subset(data, subset = Species == species)

data.recruit_size = fread("data_recruit_size.csv", data.table = F)
data.recruit_size$Site = factor(data.recruit_size$Site, levels = c("P", "M", "C", "R"))
data.recruit_size$Treatment = factor(data.recruit_size$Treatment, levels = c("Current species", "Novel species","Alpine soil","Low elevation soil"))

data.recruit_size = subset(data.recruit_size, subset = Species == species)
data.recruit_size$size<- log(data.recruit_size$size)

########################################################################################################################
### ! ## special case C.scheuchzeri as we have no recruit data we use the minimum recruit size across the other species 
### and dont subset per site nor species, so the line 448 "data.recruit_size = subset(data.recruit_size, subset = Species 
### ==species)"should be ignored for C.scheuchzeri ("cam") and min.size_3 = "with(data.recruit_size, min(size,na.rm = TRUE))"
########################################################################################################################

min_value<-data.frame(matrix(ncol=3, nrow=7))
x<- c("Species", "C.S", "min.size")
colnames(min_value)<- x
min_value$Species=species

data.p_x= subset(data.p, subset = Treatment == "Current species")
data.p_x= subset(data.p_x, subset = Site == "P")
data.recruit_size.px<-subset(data.recruit_size, subset = Site == "P")

min.size_1<- with(data.p_x , min(size, na.rm = TRUE)) 
min.size_2 <- with(data.p_x , min(size_1, na.rm = TRUE)) 
min.size_3 <- with(data.recruit_size.px, min(size,na.rm = TRUE))
## For C.scheuchzeri run line below for min.size_3 instead of line above
#min.size_3 <- with(data.recruit_size, min(size,na.rm = TRUE)) 
min_value[1,3]<- min(min.size_1,min.size_2, min.size_3)
min_value[1,2]<-"P.CS"

data.p_x= subset(data.p, subset = Treatment == "Novel species")
data.p_x= subset(data.p_x, subset = Site == "P")
data.recruit_size.px<-subset(data.recruit_size, subset = Site == "P")

min.size_1<- with(data.p_x , min(size, na.rm = TRUE)) 
min.size_2 <- with(data.p_x , min(size_1, na.rm = TRUE)) 
min.size_3 <- with(data.recruit_size.px, min(size,na.rm = TRUE))
## For C.scheuchzeri run line below for min.size_3 instead of line above
#min.size_3 <- with(data.recruit_size, min(size,na.rm = TRUE))
min_value[2,3]<- min(min.size_1,min.size_2, min.size_3)
min_value[2,2]<-"P.NS"

data.p_x= subset(data.p, subset = Treatment == "Current species")
data.p_x= subset(data.p_x, subset = Site == "M")
data.recruit_size.px<-subset(data.recruit_size, subset = Site == "M")

min.size_1<- with(data.p_x , min(size, na.rm = TRUE)) 
min.size_2 <- with(data.p_x , min(size_1, na.rm = TRUE)) 
min.size_3 <- with(data.recruit_size.px, min(size,na.rm = TRUE))
## For C.scheuchzeri run line below for min.size_3 instead of line above
#min.size_3 <- with(data.recruit_size, min(size,na.rm = TRUE))
min_value[3,3]<- min(min.size_1,min.size_2, min.size_3)
min_value[3,2]<-"M.CS"

data.p_x= subset(data.p, subset = Treatment == "Novel species")
data.p_x= subset(data.p_x, subset = Site == "M")
data.recruit_size.px<-subset(data.recruit_size, subset = Site == "M")

min.size_1<- with(data.p_x , min(size, na.rm = TRUE)) 
min.size_2 <- with(data.p_x , min(size_1, na.rm = TRUE)) 
min.size_3 <- with(data.recruit_size.px, min(size,na.rm = TRUE))
## For C.scheuchzeri run line below for min.size_3 instead of line above
#min.size_3 <- with(data.recruit_size, min(size,na.rm = TRUE))
min_value[4,3]<- min(min.size_1,min.size_2, min.size_3)
min_value[4,2]<-"M.NS"

data.p_x= subset(data.p, subset = Treatment == "Current species")
data.p_x= subset(data.p_x, subset = Site == "C")
data.recruit_size.px<-subset(data.recruit_size, subset = Site == "C")

min.size_1<- with(data.p_x , min(size, na.rm = TRUE)) 
min.size_2 <- with(data.p_x , min(size_1, na.rm = TRUE)) 
min.size_3 <- with(data.recruit_size.px, min(size,na.rm = TRUE))
## For C.scheuchzeri run line below for min.size_3 instead of line above
#min.size_3 <- with(data.recruit_size, min(size,na.rm = TRUE))
min_value[5,3]<- min(min.size_1,min.size_2, min.size_3)
min_value[5,2]<-"C.CS"

data.p_x= subset(data.p, subset = Treatment == "Novel species")
data.p_x= subset(data.p_x, subset = Site == "C")
data.recruit_size.px<-subset(data.recruit_size, subset = Site == "C")

min.size_1<- with(data.p_x , min(size, na.rm = TRUE)) 
min.size_2 <- with(data.p_x , min(size_1, na.rm = TRUE)) 
min.size_3 <- with(data.recruit_size.px, min(size,na.rm = TRUE))
## For C.scheuchzeri run line below for min.size_3 instead of line above
#min.size_3 <- with(data.recruit_size, min(size,na.rm = TRUE))
min_value[6,3]<- min(min.size_1,min.size_2, min.size_3)
min_value[6,2]<-"C.NS"

data.p_x= subset(data.p, subset = Treatment == "Current species")
data.p_x= subset(data.p_x, subset = Site == "R")
data.recruit_size.px<-subset(data.recruit_size, subset = Site == "R")

min.size_1<- with(data.p_x , min(size, na.rm = TRUE)) 
min.size_2 <- with(data.p_x , min(size_1, na.rm = TRUE)) 
min.size_3 <- with(data.recruit_size.px, min(size,na.rm = TRUE))
## For C.scheuchzeri run line below for min.size_3 instead of line above
#min.size_3 <- with(data.recruit_size, min(size,na.rm = TRUE))
min_value[7,3]<- min(min.size_1,min.size_2, min.size_3)
min_value[7,2]<-"R.CS"

########################################################################################################################
## Maximum size
########################################################################################################################

max_value<-data.frame(matrix(ncol=3, nrow=7))
x<- c("Species", "C.S", "max.size")
colnames(max_value)<- x
max_value$Species=species

data.p_x= subset(data.p, subset = Treatment == "Current species")
data.p_x= subset(data.p_x, subset = Site == "P")

max.size_1<- with(data.p_x , max(size, na.rm = TRUE)) 
max.size_2 <- with(data.p_x , max(size_1, na.rm = TRUE)) 
max_value[1,3]<- max(max.size_1,max.size_2)
max_value[1,2]<-"P.CS"

data.p_x= subset(data.p, subset = Treatment == "Novel species")
data.p_x= subset(data.p_x, subset = Site == "P")

max.size_1<- with(data.p_x , max(size, na.rm = TRUE)) 
max.size_2 <- with(data.p_x , max(size_1, na.rm = TRUE)) 
max_value[2,3]<- max(max.size_1,max.size_2)
max_value[2,2]<-"P.NS"

data.p_x= subset(data.p, subset = Treatment == "Current species")
data.p_x= subset(data.p_x, subset = Site == "M")

max.size_1<- with(data.p_x , max(size, na.rm = TRUE)) 
max.size_2 <- with(data.p_x , max(size_1, na.rm = TRUE)) 
max_value[3,3]<- max(max.size_1,max.size_2)
max_value[3,2]<-"M.CS"

data.p_x= subset(data.p, subset = Treatment == "Novel species")
data.p_x= subset(data.p_x, subset = Site == "M")

max.size_1<- with(data.p_x , max(size, na.rm = TRUE)) 
max.size_2 <- with(data.p_x , max(size_1, na.rm = TRUE)) 
max_value[4,3]<- max(max.size_1,max.size_2)
max_value[4,2]<-"M.NS"

data.p_x= subset(data.p, subset = Treatment == "Current species")
data.p_x= subset(data.p_x, subset = Site == "C")

max.size_1<- with(data.p_x , max(size, na.rm = TRUE)) 
max.size_2 <- with(data.p_x , max(size_1, na.rm = TRUE)) 
max_value[5,3]<- max(max.size_1,max.size_2)
max_value[5,2]<-"C.CS"

data.p_x= subset(data.p, subset = Treatment == "Novel species")
data.p_x= subset(data.p_x, subset = Site == "C")

max.size_1<- with(data.p_x , max(size, na.rm = TRUE)) 
max.size_2 <- with(data.p_x , max(size_1, na.rm = TRUE)) 
max_value[6,3]<- max(max.size_1,max.size_2)
max_value[6,2]<-"C.NS"
max_value[c("max.size")]

data.p_x= subset(data.p, subset = Treatment == "Current species")
data.p_x= subset(data.p_x, subset = Site == "R")

max.size_1<- with(data.p_x , max(size, na.rm = TRUE)) 
max.size_2 <- with(data.p_x , max(size_1, na.rm = TRUE)) 
max_value[7,3]<- max(max.size_1,max.size_2)
max_value[7,2]<-"R.CS"

max_min_values<-cbind(min_value, max_value[c("max.size")])

name <- paste(species, "max_min_values",sep=".")

assign(name, max_min_values)

#### Minimum and maximum sizes for each species saved:
#pla.max_min_values
#ant.max_min_values
#tri.max_min_values
#cam.max_min_values

max_min_values_all<-rbind(pla.max_min_values,ant.max_min_values,tri.max_min_values,cam.max_min_values)

##################################################################################################################
##   2.       Estimate population growth 
##################################################################################################################

## subset for species

species = "pla"

#species = "ant"

#species = "tri"

#species = "cam"

parameters_p = pla.parameters
  ##subset(parameters, subset = Species == species)

## Combine the parameters with the minimum and maximum size

max_min_values_all_p<-read.csv("max_min_values_all_trial_resubmission_2.csv")
max_min_values_all_p<-subset(max_min_values_all_p, subset = Species == species)
min.p<-as.data.frame(max_min_values_all_p[,4])
min.p<-t(min.p)
min.p<-as.data.frame(min.p)
min.p$parameter<-"min.size"
min.p$Species<-species
row.names(min.p)<-NULL
colnames(min.p)<-colnames(parameters_p)

max.p<-as.data.frame(max_min_values_all_p[,5])
max.p<-t(max.p)
max.p<-as.data.frame(max.p)
max.p$parameter<-"max.size"
max.p$Species<-species
row.names(max.p)<-NULL
colnames(max.p)<-colnames(parameters_p)

parameters_p<-rbind(parameters_p,min.p,max.p)

 

dat<-parameters_p

counter=lambda=NULL

## Run the IPM for each treatment per species
for (i in 1:7) {
  # i-th element of `u1` squared into `i`-th position of `usq
  grow.int=dat[1,i]
  grow.z=dat[2,i]
  surv.int=dat[4,i]
  surv.z=dat[5,i]
  flow.int=dat[6,i]
  flow.z =dat[7,i]
  grow.sd =dat[3,i]
  seed.mean  = dat[8,i]
  rcsz.int =dat[9,i] 
  rcsz.sd = dat[10,i]
  p.r = dat[11,i]
  
  m.par.est = c(grow.int,
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
  
  names(m.par.est) <-c("grow.int","grow.z","grow.sd","surv.int","surv.z","flow.int","flow.z","seed.mean","rcsz.int","rcsz.sd","p.r")
  
  ## Growth kernel
  G_z1z <- function (l1, l, m.par.est) {
    dnorm(l1, mean = m.par.est["grow.int"]+ m.par.est["grow.z"]*l, sd = m.par.est["grow.sd"])
    
  }
  
  ## Survival kernel
  s_z <- function(l, m.par.est)
  {
    linear.p <- m.par.est["surv.int"] + m.par.est["surv.z"] * l 
    p <- 1/(1+exp(-linear.p))                            
    return(p)
  }
  
  P_z1z <- function (l1, l, m.par.est) {
    
    return((s_z(l, m.par.est) * G_z1z(l1, l, m.par.est)))
  }
  
  ## Reproductive kernel
  
  F_z1z <- function (l1, l, m.par.est) {
    return(
      (1/(1+1/exp(m.par.est["flow.int"] + m.par.est["flow.z"]*l))) *
        m.par.est["seed.mean"]*
        m.par.est["p.r"] * 
        dnorm(l1, mean = m.par.est["rcsz.int"], sd = m.par.est["rcsz.sd"]))
  }
  
  nBigMatrix <- 250
  min.size <-dat[12,i]
  max.size <-dat[13,i]
  U <- max.size*1.1
  L<-min.size
  h <- (U-L)/nBigMatrix 
  meshpts <- L + (1:nBigMatrix)*h - h/2
  
  mk_K <- function(nBigMatrix, m.par.est, L, U) {
    h <- (U - L)/nBigMatrix
    meshpts <- L + ((1:nBigMatrix) - 1/2) * h
    
    
    P <- h * (outer(meshpts, meshpts, P_z1z, m.par.est = m.par.est))
    F. <- h * (outer(meshpts, meshpts, F_z1z, m.par.est = m.par.est))
    K <- P + F.
    
    
    return(list(K = K, meshpts = meshpts, P = P, F. = F.))
  }
  
  IPM.est <- mk_K(nBigMatrix, m.par.est, L,U)
  lambda <- Re(eigen(IPM.est$K, only.values = TRUE)$values[1])
  
  
  counter = rbind(counter,lambda)
  
  print(i) 
  
  
}

lambdas<-as.data.frame(counter,row.names = FALSE)
colnames(lambdas)<- c("lambda")

name <- paste(species, "lambdas",sep=".")
assign(name, lambdas)

#### Population growth for each species saved:

#cam.lambdas
#tri.lambdas
#ant.lambdas
#pla.lambdas

lambdas_all<-rbind(pla.lambdas, ant.lambdas, tri.lambdas, cam.lambdas)
lambdas_max_min_all<-cbind(max_min_values_all,lambdas_all)

##################################################################################################################
##   3.       Bootstrap lambdas 
##################################################################################################################

## Packages required
library(boot)
library(overlap)
library(data.table)
library(gridExtra)
library(ggplot2)

## Read in data
data = fread("data_fit_models.csv", data.table = F)


#################################################################################
## Organize data
#################################################################################

data$Site = factor(data$Site, levels = c("P", "M", "C", "R"))
data$Treatment = factor(data$Treatment, levels = c( "Current species",  "Novel species","Alpine soil","Low elevation soil"))
data$Species = factor(data$Species, levels = c( "pla","ant","tri","cam"))

data$Site <-as.factor(data$Site)
data$Treatment <-as.factor(data$Treatment)
data$size<- as.numeric(data$size)
data$Repr<- as.numeric(data$Repr)
data$Surv<- as.numeric(data$Surv)
data$size_1<- as.numeric(data$size_1)

### Log-transform size

data$size<- log(data$size)
data$size_1<- log(data$size_1)

##Subset for only community treatments 
data = subset(data, subset = Treatment == "Current species"|Treatment == "Novel species")

#################################################################################
## Organize data recruitment
#################################################################################

data.recruitment = fread("data_recruitment.csv", data.table = F)

data.recruitment$Site = factor(data.recruitment$Site, levels = c("P", "M", "C", "R"))
data.recruitment$Treatment = factor(data.recruitment$Treatment, levels = c("Current species", "Novel species","Alpine soil","Low elevation soil"))

#################################################################################
## Organize data recruit size
#################################################################################

data.rs = fread("data_recruit_size.csv", data.table = F)
head(data.rs)
data.rs$Site = factor(data.rs$Site, levels = c("P", "M", "C", "R"))
data.rs$Treatment = factor(data.rs$Treatment, levels = c( "Current species",  "Novel species","Alpine soil","Low elevation soil"))

data.rs$size<- log(data.rs$size)

#################################################################################
## Subset for species
#################################################################################

## number of bootstraps
reps= 5000

##subset for species and run one species at the time

#species = "ant"

#species = "pla"

#species = "tri"

#species = "cam"

data.p = subset(data, subset = Species == species)

data.p= subset(data.p, subset = Treatment == "Current species" |Treatment == "Novel species" )
data.seeds = subset(data.p, subset = Seeds < 1000)
data.recruitment <- subset(data.recruitment, data.recruitment$Species == species)
data.recruitment <- subset(data.recruitment, subset = Treatment == "Current species" |Treatment == "Novel species" )

data.rs <- subset(data.rs, data.rs$Species == species)

#################################################################################
## Parametric bootstrap for each vital rate:
#################################################################################

#################################################################################
################################ Seed production ###############################
#################################################################################

seed_model<-lm(Seeds ~ Site*Treatment, data=data.seeds,na.action=na.exclude)

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  
  seed_model<-lm(Seeds ~ Site*Treatment,na.action=na.exclude,out)
  
  m.par.est <- c( seed_p.cs = coef(seed_model)[1],
                  seed_m.cs = (coef(seed_model)[1]+coef(seed_model)[2]),
                  seed_c.cs = (coef(seed_model)[1]+coef(seed_model)[3]),     
                  seed_r.cs = (coef(seed_model)[1]+coef(seed_model)[4]),                
                  
                  seed_p.ns = (coef(seed_model)[1]+coef(seed_model)[5]),
                  seed_m.ns = (coef(seed_model)[1]+coef(seed_model)[2]+coef(seed_model)[5]+coef(seed_model)[6]),
                  seed_c.ns = (coef(seed_model)[1]+coef(seed_model)[3]+coef(seed_model)[5]+coef(seed_model)[7]))
  
}

rgen<-function(data.seeds,mle){
  out<-data.seeds
  out$Seeds <-unlist(simulate(mle))
  
  return(out)
}

b2<-boot(data.seeds,foo,R=reps,sim="parametric",ran.gen=rgen,mle=seed_model)

data.b<-as.data.frame(b2$t)
data.p.out <- data.b
names(data.p.out) <- c("seeds.p.cs", "seeds.m.cs", "seeds.c.cs", "seeds.r.cs", "seeds.p.ns", "seeds.m.ns", "seeds.c.ns" )
##for campanula where no seeds produced
data.p.out[is.na(data.p.out)] <- 0

#################################################################################
################################### Recruitment #################################
#################################################################################

recruit_model<-glm(recruitment ~ Site*Treatment, data=data.recruitment,family="binomial",na.action=na.exclude)
p <- 1/(1+exp(-coef(m)))  

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  
  recruit_model<-glm(recruitment ~ Site*Treatment, data=out,family="binomial",na.action=na.exclude)
  
  m.par.est <- c( recr_p.cs = coef(recruit_model)[1],
                  recr_m.cs = (coef(recruit_model)[1]+coef(recruit_model)[2]),
                  recr_c.cs = (coef(recruit_model)[1]+coef(recruit_model)[3]),     
                  recr_r.cs = (coef(recruit_model)[1]+coef(recruit_model)[4]),                
                  
                  recr_p.ns = (coef(recruit_model)[1]+coef(recruit_model)[5]),
                  recr_m.ns = (coef(recruit_model)[1]+coef(recruit_model)[2]+coef(recruit_model)[5]+coef(recruit_model)[6]),
                  recr_c.ns = (coef(recruit_model)[1]+coef(recruit_model)[3]+coef(recruit_model)[5]+coef(recruit_model)[7]))
  
  m.par.est <- 1/(1+exp(-m.par.est))
  
}

rgen<-function(data.recruitment,mle){
  out<-data.recruitment
  out$recruitment <-unlist(simulate(mle))
  
  return(out)
}

b2<-boot(data.recruitment,foo,R=reps,sim="parametric",ran.gen=rgen,mle=recruit_model)

data.b<-as.data.frame(b2$t)

data.recr.out <- data.b
names(data.recr.out) <- c("recr.p.cs", "recr.m.cs", "recr.c.cs", "recr.r.cs", "recr.p.nc", "recr.m.nc", "recr.c.nc" )

#################################################################################
################################### Recruit size ################################
#################################################################################

recruit_size_model<-lm(size ~ Site, data=data.rs,na.action=na.exclude)

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  
  recruit_size_model<-lm(size ~ Site, data=out,na.action=na.exclude)
  
  m.par.est <- c( recr_p = coef(recruit_size_model)[1],
                  recr_m = (coef(recruit_size_model)[1]+coef(recruit_size_model)[2]),
                  recr_c = (coef(recruit_size_model)[1]+coef(recruit_size_model)[3]),     
                  recr_r = (coef(recruit_size_model)[1]+coef(recruit_size_model)[4]))                
}

rgen<-function(data.recruitment,mle){
  out<-data.rs
  out$size <-unlist(simulate(mle))
  
  return(out)
}

b2<-boot(data.rs,foo,R=reps,sim="parametric",ran.gen=rgen,mle=recruit_size_model)

data.b<-as.data.frame(b2$t)

data.rs.out <- data.b
names(data.rs.out) <- c("rs.p", "rs.m", "rs.c", "rs.r" )

################################################################################################################################
## ! ## Exception P.alpina (pla) and T.badium (tri), instead of section above, run below 
## to obtrain recruit size as recruit size and SD = 0
################################################################################################################################

recruit_size_model<-lm(size ~ Site, data=data.rs,na.action=na.exclude)


foo<-function(out,sample.index){
  out <- out[sample.index, ]
  
  recruit_size_model<-lm(size ~ Site, data=out,na.action=na.exclude)
  
  m.par.est <- c( recr_p =0, ##as not included in model =only NAs
                  recr_m =  coef(recruit_size_model)[1],
                  recr_c = (coef(recruit_size_model)[1]+coef(recruit_size_model)[2]),
                  recr_r = (coef(recruit_size_model)[1]+coef(recruit_size_model)[3]))                
}

rgen<-function(data.germ,mle){
  out<-data.rs
  out$size <-unlist(simulate(mle))
  
  return(out)
}

b2<-boot(data.rs,foo,R=reps,sim="parametric",ran.gen=rgen,mle=recruit_size_model)

data.b<-as.data.frame(b2$t)

data.rs.out <- data.b
names(data.rs.out) <- c("rs.p", "rs.m", "rs.c", "rs.r" )

##################### Recruit size - SD ###################

data.rp= subset(data.rs, subset = Site == "P")
data.rm= subset(data.rs, subset = Site == "M")
data.rc= subset(data.rs, subset = Site == "C")
data.rr= subset(data.rs, subset = Site == "R")

modrs_sd.p<-lm(size ~ 1, data=data.rs,na.action=na.exclude, subset = Site == "P") 
modrs_sd.m<-lm(size ~ 1, data=data.rs,na.action=na.exclude, subset = Site == "M") 
modrs_sd.c<-lm(size ~ 1, data=data.rs,na.action=na.exclude, subset = Site == "C") 
modrs_sd.r<-lm(size ~ 1, data=data.rs,na.action=na.exclude, subset = Site == "R") 

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  modrs_sd<-lm(size ~ 1, data=out, na.action=na.exclude)
  m.par.est <- c( 
    rs.sd = summary(modrs_sd)$sigma)
}

rgen<-function(data.recsize,mle){
  out<-data.recsize
  out$size <-unlist(simulate(mle))
  return(out)
}

################################################################################################################################
##!## Exception P.alpina (pla) and T.badium (tri) as recruit size and SD = 0
## replace row 1037 and 1040 with this below:
#sd2<-data.frame()[1:5000,1 ]
#sd2$t[1:5000]<-c(0)
#p.rs_sd<-as.data.frame(sd2$t) 
################################################################################################################################

data.recsize <- data.rp
sd2<-boot(data.recsize,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modrs_sd.p)
#sd2<-data.frame()[1:5000,1 ]
#sd2$t[1:5000]<-c(0)
p.rs_sd<-as.data.frame(sd2$t)

data.recsize <- data.rm
sd2<-boot(data.recsize,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modrs_sd.m)
m.rs_sd<-as.data.frame(sd2$t)

data.recsize <- data.rc
sd2<-boot(data.recsize,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modrs_sd.c)
c.rs_sd<-as.data.frame(sd2$t)

data.recsize <- data.rr
sd2<-boot(data.recsize,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modrs_sd.r)
r.rs_sd<-as.data.frame(sd2$t)

data.rs_sd.out <- cbind(p.rs_sd, m.rs_sd, c.rs_sd, r.rs_sd)
names(data.rs_sd.out) <- c("rs_sd.p", "rs_sd.m", "rs_sd.c", "rs_sd.r" )

#################################################################################
##################################### Growth ####################################
#################################################################################

growth_model<-lm(size_1~size*Site*Treatment, data=data.p,na.action=na.exclude)

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  
  growth_model<-lm(size_1~size*Site*Treatment,na.action=na.exclude,out)
  
  m.par.est <- c( grow.int_p.cc = coef(growth_model)[1],
                  grow.z_p.cc = coef(growth_model)[2],
                  
                  grow.int_p.nc = (coef(growth_model)[1]+coef(growth_model)[6]),
                  grow.z_p.nc = (coef(growth_model)[2]+coef(growth_model)[10]),
                  
                  grow.int_m.cc = (coef(growth_model)[1]+coef(growth_model)[3]),
                  grow.z_m.cc = (coef(growth_model)[2]+coef(growth_model)[7]),
                  
                  grow.int_m.nc = (coef(growth_model)[1]+coef(growth_model)[6]+coef(growth_model)[3]+coef(growth_model)[11]),
                  grow.z_m.nc = (coef(growth_model)[2]+coef(growth_model)[10]+coef(growth_model)[7]+coef(growth_model)[14]),
                  
                  grow.int_c.cc = (coef(growth_model)[1]+coef(growth_model)[4]),
                  grow.z_c.cc = (coef(growth_model)[2]+coef(growth_model)[8]),
                  
                  grow.int_c.nc = (coef(growth_model)[1]+coef(growth_model)[6]+coef(growth_model)[4]+coef(growth_model)[12]),
                  grow.z_c.nc = (coef(growth_model)[2]+coef(growth_model)[10]+coef(growth_model)[8]+coef(growth_model)[15]),
                  
                  grow.int_r.cc = (coef(growth_model)[1]+coef(growth_model)[5]),
                  grow.z_r.cc = (coef(growth_model)[2]+coef(growth_model)[9]))
  
}

rgen<-function(data.p,mle){
  out<-data.p
  out$size_1<-unlist(simulate(mle))
  
  return(out)
}

b2<-boot(data.p,foo,R=reps, sim="parametric",ran.gen=rgen,mle=growth_model)

data.b<-as.data.frame(b2$t)

p.cs_gr<-as.data.frame(data.b[1:2])

colnames(p.cs_gr) <- c("grow.int", "grow.z")

p.ns_gr<-as.data.frame(data.b[3:4])
head(p.ns_gr)
head(data.b)
colnames(p.ns_gr) <- c("grow.int", "grow.z")

m.cs_gr<-as.data.frame(data.b[5:6])
colnames(m.cs_gr) <- c("grow.int", "grow.z")
head(m.cs_gr)
head(data.b)

m.ns_gr<-as.data.frame(data.b[7:8])
colnames(m.ns_gr) <- c("grow.int", "grow.z")
head(m.ns_gr)
head(data.b)

c.cs_gr<-as.data.frame(data.b[9:10])
colnames(c.cs_gr) <- c("grow.int", "grow.z")
head(c.cs_gr)
head(data.b)

c.ns_gr<-as.data.frame(data.b[11:12])
colnames(c.ns_gr) <- c("grow.int", "grow.z")
head(c.ns_gr)
head(data.b)

r.cs_gr<-as.data.frame(data.b[13:14])
colnames(r.cs_gr) <- c("grow.int", "grow.z")
head(r.cs_gr)
head(data.b)

#################################################################################
################################### Survival ####################################
#################################################################################

data.s<-data.p[!is.na(data.p$Surv),]
survival_model<-glm(Surv~size*Site*Treatment, family="binomial", data=data.s, maxit=100)


foo<-function(out,sample.index){
  out <- out[sample.index, ]
  
  survival_model<-glm(Surv~size*Site*Treatment,family="binomial",out, maxit=100)
  
  m.par.est <- c( surv.int_p.cc = coef(survival_model)[1],
                  surv.z_p.cc = coef(survival_model)[2],
                  
                  surv.int_p.nc = (coef(survival_model)[1]+coef(survival_model)[6]),
                  surv.z_p.nc = (coef(survival_model)[2]+coef(survival_model)[10]),
                  
                  surv.int_m.cc = (coef(survival_model)[1]+coef(survival_model)[3]),
                  surv.z_m.cc = (coef(survival_model)[2]+coef(survival_model)[7]),
                  
                  surv.int_m.nc = (coef(survival_model)[1]+coef(survival_model)[6]+coef(survival_model)[3]+coef(survival_model)[11]),
                  surv.z_m.nc = (coef(survival_model)[2]+coef(survival_model)[10]+coef(survival_model)[7]+coef(survival_model)[14]),
                  
                  surv.int_c.cc = (coef(survival_model)[1]+coef(survival_model)[4]),
                  surv.z_c.cc = (coef(survival_model)[2]+coef(survival_model)[8]),
                  
                  surv.int_c.nc = (coef(survival_model)[1]+coef(survival_model)[6]+coef(survival_model)[4]+coef(survival_model)[12]),
                  surv.z_c.nc = (coef(survival_model)[2]+coef(survival_model)[10]+coef(survival_model)[8]+coef(survival_model)[15]),
                  
                  surv.int_r.cc = (coef(survival_model)[1]+coef(survival_model)[5]),
                  surv.z_r.cc = (coef(survival_model)[2]+coef(survival_model)[9]))
  
}

rgen<-function(data.s,mle){
  out<-data.s
  out$Surv<-unlist(simulate(mle))
  
  return(out)
}

s2<-boot(data.s,foo,R=reps,sim="parametric",ran.gen=rgen,mle=survival_model)

data.s<-as.data.frame(s2$t)
head(data.s)

p.cs_surv<-as.data.frame(data.s[1:2])
colnames(p.cs_surv) <- c("surv.int", "surv.z")

p.ns_surv<-as.data.frame(data.s[3:4])
colnames(p.ns_surv) <- c("surv.int", "surv.z")

m.cs_surv<-as.data.frame(data.s[5:6])
colnames(m.cs_surv) <- c("surv.int", "surv.z")

m.ns_surv<-as.data.frame(data.s[7:8])
colnames(m.ns_surv) <- c("surv.int", "surv.z")

c.cs_surv<-as.data.frame(data.s[9:10])
colnames(c.cs_surv) <- c("surv.int", "surv.z")

c.ns_surv<-as.data.frame(data.s[11:12])
colnames(c.ns_surv) <- c("surv.int", "surv.z")

r.cs_surv<-as.data.frame(data.s[13:14])
colnames(r.cs_surv) <- c("surv.int", "surv.z")


#################################################################################
################################## Flowering ###################################
#################################################################################

data.f<-data.p[!is.na(data.p$Repr),]
flowering_model<-glm(Repr~size*Site*Treatment, family="binomial", data=data.f, maxit=100)

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  
  flowering_model<-glm(Repr~size*Site*Treatment,family="binomial",out, maxit=100)
  
  m.par.est <- c( Repr.int_p.cs = coef(flowering_model)[1],
                  Repr.z_p.cs = coef(flowering_model)[2],
                  
                  Repr.int_p.ns = (coef(flowering_model)[1]+coef(flowering_model)[6]),
                  Repr.z_p.ns = (coef(flowering_model)[2]+coef(flowering_model)[10]),
                  
                  Repr.int_m.cs = (coef(flowering_model)[1]+coef(flowering_model)[3]),
                  Repr.z_m.cs = (coef(flowering_model)[2]+coef(flowering_model)[7]),
                  
                  Repr.int_m.ns = (coef(flowering_model)[1]+coef(flowering_model)[6]+coef(flowering_model)[3]+coef(flowering_model)[11]),
                  Repr.z_m.ns = (coef(flowering_model)[2]+coef(flowering_model)[10]+coef(flowering_model)[7]+coef(flowering_model)[14]),
                  
                  Repr.int_c.cs = (coef(flowering_model)[1]+coef(flowering_model)[4]),
                  Repr.z_c.cs = (coef(flowering_model)[2]+coef(flowering_model)[8]),
                  
                  Repr.int_c.ns = (coef(flowering_model)[1]+coef(flowering_model)[6]+coef(flowering_model)[4]+coef(flowering_model)[12]),
                  Repr.z_c.ns = (coef(flowering_model)[2]+coef(flowering_model)[10]+coef(flowering_model)[8]+coef(flowering_model)[15]),
                  
                  Repr.int_r.cs = (coef(flowering_model)[1]+coef(flowering_model)[5]),
                  Repr.z_r.cs = (coef(flowering_model)[2]+coef(flowering_model)[9]))
  
}

rgen<-function(data.f,mle){
  out<-data.f
  out$Repr<-unlist(simulate(mle))
  
  return(out)
}

f2<-boot(data.f,foo,R=reps,sim="parametric",ran.gen=rgen,mle=flowering_model)

data.f<-as.data.frame(f2$t)

p.cs_flow<-as.data.frame(data.f[1:2])
colnames(p.cs_flow) <- c("flow.int", "flow.z")

p.ns_flow<-as.data.frame(data.f[3:4])
head(p.ns_flow)

m.cs_flow<-as.data.frame(data.f[5:6])
colnames(m.cs_flow) <- c("flow.int", "flow.z")

m.ns_flow<-as.data.frame(data.f[7:8])
colnames(m.ns_flow) <- c("flow.int", "flow.z")

c.cs_flow<-as.data.frame(data.f[9:10])
colnames(c.cs_flow) <- c("flow.int", "flow.z")


c.ns_flow<-as.data.frame(data.f[11:12])
colnames(c.ns_flow) <- c("flow.int", "flow.z")


r.cs_flow<-as.data.frame(data.f[13:14])
colnames(r.cs_flow) <- c("flow.int", "flow.z")

#################################################################################
################################## Growth SD ###################################
#################################################################################

##################### P.CS ##################### 

data.p.cs= subset(data.p, subset = Treatment == "Current species")
data.p.cs= subset(data.p.cs, subset = Site == "P")

modgr_sd<-lm(size_1~size, data=data.p.cs,na.action=na.exclude) 
grow.sd = summary(modgr_sd)$sigma

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  modgr_sd<-lm(size_1~size, data=out, na.action=na.exclude)
  m.par.est <- c( 
    grow.sd = summary(modgr_sd)$sigma)
  
}

rgen<-function(data.p.cs,mle){
  out<-data.p.cs
  out$size_1<-unlist(simulate(mle))
  
  return(out)
}

sd2<-boot(data.p.cs,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modgr_sd)

p.cs_gr_sd<-as.data.frame(sd2$t)
colnames(p.cs_gr_sd) <- c("grow.sd")


##################### P.NS ##################### 

data.p.ns= subset(data.p, subset = Treatment == "Novel species")
data.p.ns= subset(data.p.ns, subset = Site == "P")

modgr_sd<-lm(size_1~size, data=data.p.ns, na.action=na.exclude) 
grow.sd = summary(modgr_sd)$sigma

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  modgr_sd<-lm(size_1~size, data=out, na.action=na.exclude)
  m.par.est <- c( 
    grow.sd = summary(modgr_sd)$sigma)
  
}

rgen<-function(data.p.ns,mle){
  out<-data.p.ns
  out$size_1<-unlist(simulate(mle))
  
  return(out)
}

sd2<-boot(data.p.ns,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modgr_sd)

p.ns_gr_sd<-as.data.frame(sd2$t)
colnames(p.ns_gr_sd) <- c("grow.sd")

##################### M.CS ##################### 

data.m.cs= subset(data.p, subset = Treatment == "Current species")
data.m.cs= subset(data.m.cs, subset = Site == "M")

modgr_sd<-lm(size_1~size, data=data.m.cs, na.action=na.exclude) 
grow.sd = summary(modgr_sd)$sigma

coef(modgr_sd)

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  modgr_sd<-lm(size_1~size, data=out, na.action=na.exclude)
  m.par.est <- c( 
    grow.sd = summary(modgr_sd)$sigma)
  
}

rgen<-function(data.m.cs,mle){
  out<-data.m.cs
  out$size_1<-unlist(simulate(mle))
  
  return(out)
}

sd2<-boot(data.m.cs,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modgr_sd)

m.cs_gr_sd<-as.data.frame(sd2$t)

colnames(m.cs_gr_sd) <- c("grow.sd")

##################### m.ns ##################### 

data.m.ns= subset(data.p, subset = Treatment == "Novel species")
data.m.ns= subset(data.m.ns, subset = Site == "M")

modgr_sd<-lm(size_1~size, data=data.m.ns, na.action=na.exclude) 
grow.sd = summary(modgr_sd)$sigma


foo<-function(out,sample.index){
  out <- out[sample.index, ]
  modgr_sd<-lm(size_1~size, data=out, na.action=na.exclude)
  m.par.est <- c( 
    grow.sd = summary(modgr_sd)$sigma)

}

rgen<-function(data.m.ns,mle){
  out<-data.m.ns
  out$size_1<-unlist(simulate(mle))
  
  return(out)
}

sd2<-boot(data.m.ns,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modgr_sd)

m.ns_gr_sd<-as.data.frame(sd2$t)

colnames(m.ns_gr_sd) <- c("grow.sd")

##################### C.CS ##################### 

data.c.cs= subset(data.p, subset = Treatment == "Current species")
data.c.cs= subset(data.c.cs, subset = Site == "C")

modgr_sd<-lm(size_1~size, data=data.c.cs, na.action=na.exclude) 
grow.sd = summary(modgr_sd)$sigma

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  modgr_sd<-lm(size_1~size, data=out, na.action=na.exclude)
  m.par.est <- c( 
    grow.sd = summary(modgr_sd)$sigma)
  
}

rgen<-function(data.c.cs,mle){
  out<-data.c.cs
  out$size_1<-unlist(simulate(mle))
  
  return(out)
}

sd2<-boot(data.c.cs,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modgr_sd)

c.cs_gr_sd<-as.data.frame(sd2$t)

colnames(c.cs_gr_sd) <- c("grow.sd")

##################### C.NS ##################### 

data.c.ns= subset(data.p, subset = Treatment == "Novel species")
data.c.ns= subset(data.c.ns, subset = Site == "C")

modgr_sd<-lm(size_1~size, data=data.c.ns, na.action=na.exclude) 
grow.sd = summary(modgr_sd)$sigma

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  modgr_sd<-lm(size_1~size, data=out, na.action=na.exclude)
  m.par.est <- c( 
    grow.sd = summary(modgr_sd)$sigma)
  
}

rgen<-function(data.c.ns,mle){
  out<-data.c.ns
  out$size_1<-unlist(simulate(mle))
  
  return(out)
}

sd2<-boot(data.c.ns,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modgr_sd)

c.ns_gr_sd<-as.data.frame(sd2$t)

colnames(c.ns_gr_sd) <- c("grow.sd")

##################### R.CS ##################### 

data.r.cs= subset(data.p, subset = Treatment == "Current species")
data.r.cs= subset(data.r.cs, subset = Site == "R")

modgr_sd<-lm(size_1~size, data=data.r.cs, na.action=na.exclude) 
grow.sd = summary(modgr_sd)$sigma

foo<-function(out,sample.index){
  out <- out[sample.index, ]
  modgr_sd<-lm(size_1~size, data=out, na.action=na.exclude)
  m.par.est <- c( 
    grow.sd = summary(modgr_sd)$sigma)
}

rgen<-function(data.r.cs,mle){
  out<-data.r.cs
  out$size_1<-unlist(simulate(mle))
  return(out)
}

sd2<-boot(data.r.cs,foo,R=reps,sim="parametric",ran.gen=rgen,mle=modgr_sd)
sd2$t

r.cs_gr_sd<-as.data.frame(sd2$t)

colnames(r.cs_gr_sd) <- c("grow.sd")

################### combine all ################### 

p.cs_sim_par<-cbind(p.cs_gr,p.cs_surv, p.cs_flow, p.cs_gr_sd, data.p.out[,1], data.recr.out[,1], data.rs.out[,1], data.rs_sd.out[,1])

p.ns_sim_par<-cbind(p.ns_gr,p.ns_surv, p.ns_flow, p.ns_gr_sd, data.p.out[,5], data.recr.out[,5], data.rs.out[,1], data.rs_sd.out[,1])

m.cs_sim_par<-cbind(m.cs_gr,m.cs_surv, m.cs_flow, m.cs_gr_sd, data.p.out[,2], data.recr.out[,2], data.rs.out[,2], data.rs_sd.out[,2])

m.ns_sim_par<-cbind(m.ns_gr,m.ns_surv, m.ns_flow, m.ns_gr_sd, data.p.out[,6], data.recr.out[,6], data.rs.out[,2], data.rs_sd.out[,2])

c.cs_sim_par<-cbind(c.cs_gr,c.cs_surv, c.cs_flow, c.cs_gr_sd, data.p.out[,3], data.recr.out[,3], data.rs.out[,3], data.rs_sd.out[,3])

c.ns_sim_par<-cbind(c.ns_gr,c.ns_surv, c.ns_flow, c.ns_gr_sd, data.p.out[,7], data.recr.out[,7], data.rs.out[,3], data.rs_sd.out[,3])

r.cs_sim_par<-cbind(r.cs_gr,r.cs_surv, r.cs_flow, r.cs_gr_sd, data.p.out[,4], data.recr.out[,4], data.rs.out[,4], data.rs_sd.out[,4])

##########################################################################################
## Generate population growth values from bootstrapped data
##########################################################################################

lambdas<-lambdas_max_min_all

lam.sp <- lambdas[lambdas$Species== species,]
lam.extract <- lambdas$lambda[lambdas$Species== species]

dat.comb <- list(p.cs_sim_par, p.ns_sim_par, m.cs_sim_par, m.ns_sim_par, c.cs_sim_par, c.ns_sim_par, r.cs_sim_par)
name.dat <- c("p.cs", "p.ns","m.cs", "m.ns","c.cs", "c.ns","r.cs")

head(dat.comb)
for(j in 1:7) {
  
  dat <- dat.comb[[j]]
  nam <- name.dat[j]

  counter=lambda=NULL
  
  for (i in 1:nrow(dat)) {
    grow.int=dat[i,1]
    grow.z=dat[i,2]
    surv.int=dat[i,3]
    surv.z=dat[i,4]
    flow.int=dat[i,5]
    flow.z =dat[i,6]
    grow.sd =dat[i,7]
    seed.mean  = dat[i,8]
    
    ####################################################################################
    ##! ## Exception for C.scheuchzeri as recruit size obtained from other species
    ## see Supporting information. Replace row 1536-1538 with: 
    # rcsz.int = 1.033246
    # rcsz.sd = 0.3071283
    # p.r = 0.001
    ####################################################################################
    rcsz.int =dat[i,10] 
    rcsz.sd = dat[i,11]
    p.r = dat[i,9]
    
    m.par.est = c(grow.int,
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
    
    names(m.par.est) <- names(m.par)
    
    G_z1z <- function (l1, l, m.par.est) {
      dnorm(l1, mean = m.par.est["grow.int"]+ m.par.est["grow.z"]*l, sd = m.par.est["grow.sd"])
      
    }
    
    s_z <- function(l, m.par.est)
    {
      linear.p <- m.par.est["surv.int"] + m.par.est["surv.z"] * l 
      p <- 1/(1+exp(-linear.p))                            
      return(p)
    }
    
    P_z1z <- function (l1, l, m.par.est) {
      
      return((s_z(l, m.par.est) * G_z1z(l1, l, m.par.est)))
    }
    
    F_z1z <- function (l1, l, m.par.est) {
      return(
        (1/(1+1/exp(m.par.est["flow.int"] + m.par.est["flow.z"]*l))) *
          m.par.est["seed.mean"]*
          m.par.est["p.r"] * 
          dnorm(l1, mean = m.par.est["rcsz.int"], sd = m.par.est["rcsz.sd"]))
    }
    
    nBigMatrix <- 250
    min.size <-lam.sp$min.size[j]
    max.size <-lam.sp$max.size[j]
    U <- max.size
    U <- max.size*1.1
    L<-min.size
    h <- (U-L)/nBigMatrix 
    meshpts <- L + (1:nBigMatrix)*h - h/2
    
    mk_K <- function(nBigMatrix, m.par.est, L, U) {
      h <- (U - L)/nBigMatrix
      meshpts <- L + ((1:nBigMatrix) - 1/2) * h
      
      
      P <- h * (outer(meshpts, meshpts, P_z1z, m.par.est = m.par.est))
      F. <- h * (outer(meshpts, meshpts, F_z1z, m.par.est = m.par.est))
      K <- P + F.
      
      
      return(list(K = K, meshpts = meshpts, P = P, F. = F.))
    }
    
    IPM.est <- mk_K(nBigMatrix, m.par.est, L,U)
    lambda <- Re(eigen(IPM.est$K, only.values = TRUE)$values[1])
    
    counter = rbind(counter,lambda)
    
    print(i) 
    
    
  }
  assign(nam, counter)
  
}

bootlam.all <- cbind(p.cs,p.ns,m.cs, m.ns,c.cs, c.ns,r.cs)
colnames(bootlam.all) <- c("p.cs", "p.ns","m.cs", "m.ns","c.cs", "c.ns","r.cs")

#### calculate BCPI confidence intervals

bcpi <- function(real_lambda, t, alpha) {
  B <- length(t)
  z0 <- qnorm(mean(t < real_lambda))
  a1 <- pnorm(2 * z0 + qnorm(alpha/2))
  a2 <- pnorm(2 * z0 + qnorm(1 - alpha/2))
  c1 <- quantile(t, a1)
  c2 <- quantile(t, a2)
  return(as.numeric(c(c1, c2)))
}

ci1 <-bcpi(lam.extract[1], bootlam.all[,1], 0.05)
ci2 <-bcpi(lam.extract[2], bootlam.all[,2], 0.05)
ci3 <-bcpi(lam.extract[3], bootlam.all[,3], 0.05)
ci4 <-bcpi(lam.extract[4], bootlam.all[,4], 0.05)
ci5 <-bcpi(lam.extract[5], bootlam.all[,5], 0.05)
ci6 <-bcpi(lam.extract[6], bootlam.all[,6], 0.05)
ci7 <-bcpi(lam.extract[7], bootlam.all[,7], 0.05)
cis <- rbind(ci1,ci2,ci3,ci4,ci5,ci6,ci7)

boot.lambdas <- data.frame(
  lambda = lam.extract,
  lower = cis[,1],
  upper = cis[,2])
rownames(boot.lambdas) <- colnames(bootlam.all)

name <- paste(species, "boot.lambdas",sep=".")
assign(name, boot.lambdas)

name <- paste(species, "bootlam.all",sep=".")
assign(name, bootlam.all)

#### Bootstrapped population growth for each species saved:

#pla.bootlam.all
#ant.bootlam.all
#tri.bootlam.all
#cam.bootlam.all

#### BCPIs for each species saved:
#pla.boot.lambdas
#ant.boot.lambdas
#tri.boot.lambdas
#cam.boot.lambdas

cis<-rbind(pla.boot.lambdas,ant.boot.lambdas,tri.boot.lambdas,cam.boot.lambdas)
lambdas_cis<-cbind(lambdas_max_min_all[,c(1,2,5)],cis[,2:3])
row.names(lambdas_cis)<-NULL

################################################################################################################################
## 4.  Pairwise comparisons 
################################################################################################################################

##subset for dataset
#dataset=pla.bootlam.all
#dataset=ant.bootlam.all
#dataset=tri.bootlam.all
#dataset=cam.bootlam.all

row.names(dataset)<-NULL
##subset for species

#species="pla"
#species="ant"
#species="tri"
#species="cam"

data_p_cs= as.data.frame(dataset[,1])
data_p_ns= as.data.frame(dataset[,2])
data.m_cs= as.data.frame(dataset[,3])
data.m_ns= as.data.frame(dataset[,4])
data.c_cs= as.data.frame(dataset[,5])
data.c_ns= as.data.frame(dataset[,6])
data.r_cs= as.data.frame(dataset[,7])

## compare bootstraps between sites

comp_p.cs_m.cs <-data_p_cs$`dataset[, 1]`-data.m_cs$`dataset[, 3]`
comp_p.cs_c.cs <-data_p_cs$`dataset[, 1]`-data.c_cs$`dataset[, 5]`
comp_p.cs_r.cs <-data_p_cs$`dataset[, 1]`-data.r_cs$`dataset[, 7]`
comp_m.cs_c.cs <-data.m_cs$`dataset[, 3]`-data.c_cs$`dataset[, 5]`
comp_m.cs_r.cs <-data.m_cs$`dataset[, 3]`-data.r_cs$`dataset[, 7]`
comp_c.cs_r.cs <-data.c_cs$`dataset[, 5]`-data.r_cs$`dataset[, 7]`

### compare bootstraps between communities

comp_p.cs_ns <-data_p_cs$`dataset[, 1]`-data_p_ns$`dataset[, 2]`
comp_m.cs_ns <-data.m_cs$`dataset[, 3]`-data.m_ns$`dataset[, 4]`
comp_c.cs_ns <-data.c_cs$`dataset[, 5]`-data.c_ns$`dataset[, 6]`

dif<-cbind(comp_p.cs_m.cs,comp_p.cs_c.cs,comp_p.cs_r.cs,comp_m.cs_c.cs,comp_m.cs_r.cs,comp_c.cs_r.cs,
           comp_p.cs_ns,comp_m.cs_ns,comp_c.cs_ns)

dif<-as.data.frame(dif)

## counting number of cases with differences - this shows how large proportion of differences are above 0.
## if more than 0.95 it means that it has positive effect, if less that 0.05 is has positive effect
## if in between 0.05 - 0.95 -->  non significant effect

pairwise_comparison<-as.data.frame(sapply(dif, function(x)sum(x >0, na.rm = TRUE)/length(x)))
name <- paste(species, "pairwise_comparison",sep=".")
assign(name, pairwise_comparison)

#### Population growth comparisons for each species saved:
## pla.pairwise_comparison
## ant.pairwise_comparison
## tri.pairwise_comparison
## cam.pairwise_comparison

#####################################################################################################################
## 5. Plot population growth 
#####################################################################################################################
data.plot=lambdas_cis

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
    legend.position = c(0.30, 0.22),
    strip.text.x = element_text(
      size = 14,face = "italic"),
    legend.box.background = element_rect(colour = "black"))

## P.alpina

data.plot_pla = subset(data.plot, subset = Species == "pla")
data.plot_pla$Site<-c("P","P","M","M","C","C","R")
data.plot_pla$Treatment<-c("Current species","Novel species","Current species","Novel species","Current species","Novel species","Current species")
data.plot_pla[data.plot_pla$Site == "R", "temp"] = 0
data.plot_pla[data.plot_pla$Site == "P", "temp"] =  4.85 
data.plot_pla[data.plot_pla$Site == "C", "temp"] = 1.65 
data.plot_pla[data.plot_pla$Site == "M", "temp"] =  3.05 

x<-expression(paste(italic("Plantago alpina")))

p1<-ggplot(data = data.plot_pla, aes(x = temp, y = lambda, group=Treatment)) + 
  geom_hline(yintercept=1, size=0.6,linetype="dashed")+
  scale_fill_manual(values=c("gold", "steelblue4", "grey"),labels = c("Current species", "Novel species"))+theme_light()+
  guides(fill=guide_legend(title="Interacting with:"))+
  xlab(expression("Temperature"))+
  geom_errorbar(data=data.plot_pla,aes(ymin = lower, ymax = upper, group=Treatment), position = position_dodge(0.4),  width = 0.2)+
  geom_point(data=data.plot_pla, colour="black", shape=21, size = 7, position = position_dodge(0.4),
             aes(fill = Treatment)) + 
  annotate("text", label = "*", x = 4.95, y = 0.62, size = 10, colour = "black")+
  theme(strip.text.x = element_text(
    size = 10,face = "italic"), panel.grid.major = element_blank())+ ggtitle(x)+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0, 1.2))+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  my_theme

pla<-p1
pla

## A. alpestris

data.plot_ant = subset(data.plot, subset = Species == "ant")
data.plot_ant$Site<-c("P","P","M","M","C","C","R")
data.plot_ant$Treatment<-c("Current species","Novel species","Current species","Novel species","Current species","Novel species","Current species")
data.plot_ant[data.plot_ant$Site == "R", "temp"] = 0
data.plot_ant[data.plot_ant$Site == "P", "temp"] =  4.85 
data.plot_ant[data.plot_ant$Site == "C", "temp"] = 1.65 
data.plot_ant[data.plot_ant$Site == "M", "temp"] =  3.05 

x<-expression(paste(italic("A.alpestris")))
a1<-ggplot(data = data.plot_pla, aes(x = temp, y = lambda, group=Treatment)) + 
  geom_hline(yintercept=1, size=0.6,linetype="dashed")+
  scale_fill_manual(values=c("gold", "steelblue4", "grey"),labels = c("Current species", "Novel species"))+theme_light()+
  guides(fill=guide_legend(title="Interacting with:"))+
  xlab(expression("Temperature"))+
  geom_errorbar(data=data.plot_pla,aes(ymin = lower, ymax = upper, group=Treatment), position = position_dodge(0.4),  width = 0.2)+
  geom_point(data=data.plot_pla, colour="black", shape=21, size = 7, position = position_dodge(0.4),
             aes(fill = Treatment)) + 
  annotate("text", label = "*", x = 1.75, y = 0.72, size = 10, colour = "black")+
  annotate("text", label = "*", x = 4.95, y = 0.52, size = 10, colour = "black")+
  theme(strip.text.x = element_text(
    size = 10,face = "italic"), panel.grid.major = element_blank())+ ggtitle(x)+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0, 1.2))+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  my_theme

ant<-a1

## T. badium

data.plot_tri = subset(data.plot, subset = Species == "tri")
data.plot_tri$Site<-c("P","P","M","M","C","C","R")
data.plot_tri$Treatment<-c("Current species","Novel species","Current species","Novel species","Current species","Novel species","Current species")
data.plot_tri[data.plot_tri$Site == "R", "temp"] = 0
data.plot_tri[data.plot_tri$Site == "P", "temp"] =  4.85 
data.plot_tri[data.plot_tri$Site == "C", "temp"] = 1.65 
data.plot_tri[data.plot_tri$Site == "M", "temp"] =  3.05 

x<-expression(paste(italic("T.badium")))
t1<-ggplot(data = data.plot_pla, aes(x = temp, y = lambda, group=Treatment)) + 
  geom_hline(yintercept=1, size=0.6,linetype="dashed")+
  scale_fill_manual(values=c("gold", "steelblue4", "grey"),labels = c("Current species", "Novel species"))+theme_light()+
  guides(fill=guide_legend(title="Interacting with:"))+
  xlab(expression("Temperature"))+
  geom_errorbar(data=data.plot_pla,aes(ymin = lower, ymax = upper, group=Treatment), position = position_dodge(0.4),  width = 0.2)+
  geom_point(data=data.plot_pla, colour="black", shape=21, size = 7, position = position_dodge(0.4),
             aes(fill = Treatment)) + 
  annotate("text", label = "*", x = 1.78, y = 1, size = 10, colour = "black")+
  annotate("text", label = "*", x = 5, y = 0.38, size = 10, colour = "black")+
  theme(strip.text.x = element_text(
    size = 10,face = "italic"), panel.grid.major = element_blank())+ ggtitle(x)+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0, 1.2))+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  my_theme

tri<-t1

## C. scheuzheri

data.plot_cam = subset(data.plot, subset = Species == "cam")
data.plot_cam$Site<-c("P","P","M","M","C","C","R")
data.plot_cam$Treatment<-c("Current species","Novel species","Current species","Novel species","Current species","Novel species","Current species")
data.plot_cam[data.plot_cam$Site == "R", "temp"] = 0
data.plot_cam[data.plot_cam$Site == "P", "temp"] =  4.85 
data.plot_cam[data.plot_cam$Site == "C", "temp"] = 1.65 
data.plot_cam[data.plot_cam$Site == "M", "temp"] =  3.05 

data.p = subset(data, subset = species == "cam")
x<-expression(paste(italic("C.scheuchzeri")))
c1<-ggplot(data = data.plot_pla, aes(x = temp, y = lambda, group=Treatment)) + 
  geom_hline(yintercept=1, size=0.6,linetype="dashed")+
  scale_fill_manual(values=c("gold", "steelblue4", "grey"),labels = c("Current species", "Novel species"))+theme_light()+
  guides(fill=guide_legend(title="Interacting with:"))+
  xlab(expression("Temperature"))+
  geom_errorbar(data=data.plot_pla,aes(ymin = lower, ymax = upper, group=Treatment), position = position_dodge(0.4),  width = 0.2)+
  geom_point(data=data.plot_pla, colour="black", shape=21, size = 7, position = position_dodge(0.4),
             aes(fill = Treatment)) + 
  theme(strip.text.x = element_text(
    size = 10,face = "italic"), panel.grid.major = element_blank())+ ggtitle(x)+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), limits = c(0, 1.2))+ 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  my_theme
cam<-c1

x<-ylab(expression("Population growth rate ("*lambda*")"))

all <-grid.arrange(arrangeGrob(pla,ant, tri,cam,
                               nrow = 2,
                               left = textGrob("Population growth", rot = 90,gp = gpar(cex = 2.8)),
                               bottom = textGrob("Climate warming (°C)",gp = gpar(cex = 2.8))))

