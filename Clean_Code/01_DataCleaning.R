#####################################################--
#Differences in wild turkey nest initiation and population growth rates between a hunted and nonhunted site in South Carolina, USA   
# This script loads the raw trapping and nesting data and cleans it for analysis. 
# It also graphs several bits of raw data 
################################--
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to source file location

library(tidyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(ggpubr)

#### Trapping Data ####
Trapping_Webb <- read.csv('./RawData/Turkey_Trapping_Webb.csv')
Trapping_Webb$Date <- as.Date(Trapping_Webb$Date, format = '%m/%d/%y')
Trapping_Webb$Site <- 'Webb Complex'
Trapping_SRS <- read.csv('./RawData/Turkey_Trapping_SRS.csv')
Trapping_SRS$Date <- as.Date(Trapping_SRS$Date, format = '%m/%d/%y')
Trapping_SRS$Note <- NA
Trapping_SRS$Site <- 'SRS'
Trapping <- rbind(Trapping_SRS, Trapping_Webb)
Trapping <- subset(Trapping, !is.na(Trapping$Transmitter) & (Trapping$CaptureMyopathy != 1 | is.na(Trapping$CaptureMyopathy)) & !is.na(Trapping$Sex))

## Nests ####
Nests_SRS <- read.csv('./RawData/Nests_SRS_2020_2025.csv')
Nests_Webb <- read.csv('./RawData/Nests_Webb.csv')
Nests_Webb$Sex <- NA
Nests_Webb$Date_Incubate <- yday(as.Date(Nests_Webb$Date_Incubate, format = "%m/%d/%y"))
Nests_Webb$Site <- 'Webb Complex'
Nests_SRS$Date_Incubate <- yday(as.Date(Nests_SRS$Date_Laying, format = "%m/%d/%y"))+12
Nests <- rbind(Nests_SRS[,c(2:6, 8:18)], Nests_Webb[,-c(7)])

Nests2 <- Nests %>% group_by(Site, Year, Age, Date_Incubate, ID, Total_Offspring, Total.Eggs, NestInit, H_inc, Hatch, H_brood, Brood_lived, Attempt, Est_Poults) %>% summarize(n.male.eggs = sum(Sex), n.tested = n()) 

## We want to know the sex ratio for all our eggs (when tested) per nest attempt
#Sex = 1 is male
Nests <- Nests %>% group_by(Site, Year, Age, Date_Incubate, ID, Total_Offspring, Total.Eggs, NestInit, H_inc, Hatch, H_brood, Brood_lived, Attempt, Est_Poults) %>% summarize(n.male.eggs = sum(Sex), n.tested = n()) 
Nests$n.tested[is.na(Nests$n.male.eggs)] <- 0
Nests$n.male.eggs[is.na(Nests$n.male.eggs)] <- 0
Nests$Age <- ifelse(is.na(Nests$Age), 1, Nests$Age)

## We need the per-bird probability of nesting as well:
nesting_birds <- Nests %>% group_by(ID, Year, Age, Attempt, H_inc, H_brood, Site) %>% summarize(Nested = sum(NestInit) >= 1)
n.birds <- length(unique(paste0(nesting_birds$ID, nesting_birds$Year)))

#Precipitation information for nests
precip.dat <- data.frame(year = c(rep(2014:2018, each = 2), rep(2021:2025, each = 2)),
                         precip = c(5.06, 3.16, 5.34, 1.28, 1.87, 8.92, 3.16, 7.37,3.87,6.09,
                                    1.66, 2.96, 4.58,  3.94, 6.74 , 4.78, 2.1, 3.11, 2.31, 8.61),
                         month = rep(4:5, 10)
)

precip2 <- precip.dat %>% group_by(year) %>% summarize(totP = sum(precip),
                                                       meanP = mean(precip))

nested <- nested2 <- survived_nest1 <- birds.age <- site <- precip <- array(0, n.birds)
nesting_birds$idyear <- paste0(nesting_birds$ID, nesting_birds$Year)
for(j in 1:n.birds){
  me <- nesting_birds[nesting_birds$idyear == unique(nesting_birds$idyear)[j],]
  nested[j] <- me$Nested[1]*1
  birds.age[j] <- me$Age[1]
  site[j] <- ifelse(me$Site[1] == 'SRS', 0, 1) #site 1 = Webb, site 0 = SRS
  precip[j] <- as.numeric(precip2[precip2$year == me$Year[1], 'meanP'])
}
nested[is.na(nested)] <- 0

## If nests  didn't hatch our egg count is very low and inaccurate. Remove Total.Eggs for anywhere Hatch = 0
Nests$Total.Eggs[Nests$Hatch == 0] <- NA

#Scale doy covariate:
mean(Nests$Date_Incubate, na.rm = T) #119.1574
sd(Nests$Date_Incubate, na.rm = T) # 20.90836
Nests$Date_Incubate_s  <- (Nests$Date_Incubate  - mean(Nests$Date_Incubate, na.rm = T))/sd(Nests$Date_Incubate, na.rm = T)
Nests$Date_Incubate_s <- ifelse(is.na(Nests$Date_Incubate_s), 0, Nests$Date_Incubate_s ) #the NAs are not actually in model, since nest isn't initiated, which means the 0's are also not in model, but this avoids an annoying warning

Trapping <- subset(Trapping, Trapping$BandID != '')
Trapping$CaptureMyopathy <- ifelse(is.na(Trapping$CaptureMyopathy), 0, Trapping$CaptureMyopathy)
Hen_dat <- subset(Trapping,Trapping$Sex == 0 & Trapping$CaptureMyopathy == 0)
Male_dat <- subset(Trapping,Trapping$Sex == 1 & Trapping$CaptureMyopathy == 0)

#### *Clean male data first ####
Male_dat$eff_id <- paste0(Male_dat$BandID, Male_dat$Year)
Male_dat$eff_id <- as.factor(Male_dat$eff_id)

## Turn dates into chunks of time:
#season starts = March 20 for most of the years
## Jan 1 - March 19 prehunt; 1 to 80
#March 20 to May 5 hunt; 81 to 125
#May 6 - July 31st posthunt; 126  to 212
#Aug1 - Jan #rest of year; 213 to 365

male.age <- first.male <- last.male <- malesite <- array(NA, length(levels(Male_dat$eff_id)))
z.male <- array(NA, c(length(levels(Male_dat$eff_id)), 5))
z.male[,1] <- 1 #have to be alive to be captured
harvested <- diedhunting <- array(0, length(levels(Male_dat$eff_id)))
remove <- NA
for(i in 1:length(levels(Male_dat$eff_id))){
  me <- Male_dat[as.numeric(Male_dat$eff_id) == i,]
  malesite[i] <- ifelse(me$Site[1] == 'SRS', 0, 1)
  #check for transmitter failure:
  if(any(me$TransmitterFail == 1) & !is.na(any(me$TransmitterFail == 1))){
    stophere <- which(me$TransmitterFail == 1)
    me <- me[1:stophere,]
  }
  if(nrow(me) == 1){
    remove <- c(remove, i)
    next
  }
  
  male.age[i] <- me$Age[1]
  me$days <- yday(as.Date(me$Date, format = '%Y-%m-%d'))
  
  
  me$days2 <- as.numeric(cut(me$days, breaks = c(0, 80, 125, 212, 365), 
                             labels = c("1", "2", "3", "4")))
  for(j in 1:nrow(me)){
    z.male[i,me$days2[j]+1] <- me$Alive.Dead[j]
  }
  
  first.male[i] <- 1 #pretty sure we never caught any transmitter birds except in the first time period
  last.alive <- max(which(z.male[i,] == 1))
  z.male[i,1:last.alive] <- 1
  if(any(z.male[i,] ==0) & !is.na(any(z.male[i,] ==0))){
    last.dead <- max(which(z.male[i,] == 0))
    z.male[i,1:(last.dead-1)] <- 1
    if(last.dead == 3){
      diedhunting[i] <- 1
      if(any(me$Harvested == 1) & !is.na(any(me$Harvested == 1))){
        harvested[i] <- 1
      }
    }
  }
  last.male[i] <- max(which(!is.na(z.male[i,])))
}
male.age1 <- male.age[-remove[-1]]
first.male1 <- first.male[-remove[-1]]
last.male1 <- last.male[-remove[-1]]
z.male1 <- z.male[-remove[-1],]
diedhunting1 <- diedhunting[-remove[-1]]
harvested1 <- harvested[-remove[-1]] #were they harvested
harvested2 <- array(0, c(length(harvested1), 5)) 
harvested2[,2] <- harvested1 #when in the year (obviously, we hope, during harvest season)
length(last.male1) == nrow(z.male1) #must be true
malesite1 <- malesite[-remove[-1]]
z.adults <- z.male1[male.age1 == 1,]

#### *females next ####
Hen_dat$eff_id <- paste0(Hen_dat$BandID, Hen_dat$Year)
Hen_dat$eff_id <- as.factor(Hen_dat$eff_id)

## Turn dates into chunks of time:
#Jan 1- March 19; 1-78 (pre-nest)
#March 20 to July 30; 79- 211 (nesting times)
#July 31 - Jan 1; 212 - 365

female.age <- first.female <- last.female <- nestingattempts <- female.year <-  fem.site <- array(NA, length(levels(Hen_dat$eff_id)))
broodattempts_nodeath <- neststobrood <- broodingattempts <- broodstolive  <- nestingattempts_nodeath  <- array(0, length(levels(Hen_dat$eff_id)))
z.female <- array(NA, c(length(levels(Hen_dat$eff_id)), 4))
z.female[,1] <- 1 #have to be alive to be captured

remove.f <- NA
Hen_dat$DiedBrood <- ifelse(is.na(Hen_dat$DiedBrood), NA, ifelse(Hen_dat$DiedBrood == 0, NA, 1))
Hen_dat$DiedNesting <- ifelse(is.na(Hen_dat$DiedNesting), NA, ifelse(Hen_dat$DiedNesting == 0, NA, 1))
for(i in 1:length(levels(Hen_dat$eff_id))){
  me <- Hen_dat[as.numeric(Hen_dat$eff_id) == i,]
  fem.site[i] <- ifelse(me$Site[1] == 'SRS', 0 ,1)
  #check for transmitter failure:
  if(any(me$TransmitterFail == 1) & !is.na(any(me$TransmitterFail == 1))){
    stophere <- which(me$TransmitterFail == 1)
    me <- me[1:stophere,]
  }
  if(nrow(me) == 1){
    remove.f <- c(remove.f, i)
    next
  }
  
  female.age[i] <- me$Age[1]
  female.year[i] <- me$Year[1] - 2013
  if(any(!is.na(me$NestAttempts))){
    nestingattempts[i] <- max(me$NestAttempts, na.rm =T) #could still be NA, that's okay
    nestingattempts_nodeath[i] <- nestingattempts[i]
    if(any(!is.na(me$DiedNesting))){
      nestingattempts_nodeath[i] <- nestingattempts[i] - 1
    }
    if(any(!is.na(me$BroodStage))){
      neststobrood[i] <- max(me$BroodStage, na.rm = T)
    }
  }
  
  if(any(!is.na(me$BroodStage))){
    broodingattempts[i] <- max(me$BroodStage, na.rm =T) #could still be NA, that's okay
    if(any(!is.na(me$DiedBrood))){
      broodattempts_nodeath[i] <- broodingattempts[i] - 1
    } else{
      broodattempts_nodeath[i]  <- broodingattempts[i] 
    }
  }
  
  me$days <- yday(as.Date(me$Date, format = '%Y-%m-%d')) 
  me$day <- ifelse(me$days > 365, 365, me$days)
  me$days2 <- as.numeric(cut(me$day, breaks = c(0, 79, 212, 365), 
                             labels = c("1", "2", "3")))
  if(max(me$days2) == 1){nestingattempts[i] <- 0}
  if(any(me$days >= 100) & is.na(nestingattempts[i])){
    nestingattempts[i] <- 0
  }
  for(j in 1:nrow(me)){
    z.female[i,me$days2[j]+1] <- me$Alive.Dead[j]
  }
  
  first.female[i] <- 1 #pretty sure we never caught any transmitter birds except in the first time period
  last.alive <- max(which(z.female[i,] == 1))
  z.female[i,1:last.alive] <- 1
  if(any(z.female[i,] ==0) & !is.na(any(z.female[i,] ==0))){
    last.dead <- max(which(z.female[i,] == 0))
    z.female[i,1:(last.dead-1)] <- 1
  }
  last.female[i] <- max(which(!is.na(z.female[i,])))
}

## we need to manually add in which broods lived, since there aren't that many...
Hen_dat$YearTransmitter <- paste0(Hen_dat$Transmitter, Hen_dat$Year)
GoodBroods <- subset(Nests, !is.na(Nests$Brood_lived) & Nests$Brood_lived == 1)
GoodBroods$eff_id <- Hen_dat[match(paste0(GoodBroods$ID, GoodBroods$Year), Hen_dat$YearTransmitter), 'eff_id']
GoodBroods$eff_id <- as.numeric(factor(GoodBroods$eff_id, levels = levels(Hen_dat$eff_id)))
broodstolive[GoodBroods$eff_id] <- 1

female.age1 <- female.age[-remove.f[-1]]
female.age1[is.na(female.age1)] <- 1
first.female1 <- first.female[-remove.f[-1]]
last.female1 <- last.female[-remove.f[-1]]
z.female1 <- z.female[-remove.f[-1],]
nestingattempts1 <- nestingattempts[-remove.f[-1]]
broodingattempts1 <- broodingattempts[-remove.f[-1]]
broodattempts_nodeath1 <- broodattempts_nodeath[-remove.f[-1]]
broodattempts_nodeath1[broodattempts_nodeath1 < 0] <- 0 #due to a spurious "0" in brooding count
neststobrood1 <- neststobrood[-remove.f[-1]] 
fem.site1 <- fem.site[-remove.f[-1]]

nestingattempts_nodeath1 <- nestingattempts_nodeath[-remove.f[-1]]
nestingattempts_nodeath1a <- nestingattempts_nodeath1 - 1 #this gives us 2nd attempts or later
nestingattempts_nodeath1a[nestingattempts_nodeath1a < 0] <- 0

##then for just first attempts:
nestingattempts_nodeath1b <- nestingattempts_nodeath1
nestingattempts_nodeath1b[nestingattempts_nodeath1b >0] <- 1

neststobroodfirst <- neststobrood1
neststobroodfirst[] <- 0
neststobroodfirst[nestingattempts_nodeath1b == 1 & neststobrood1 == 1] <- 1
neststobroodsecond <- neststobrood1
neststobroodsecond[neststobroodfirst == 1] <- 0

broodstolive1 <- broodstolive[-remove.f[-1]]
broodstolive1[128] <- 0 #error
female.year1 <- female.year[-remove.f[-1]]
length(last.female1) == nrow(z.female1) #must be true
zeronest <- (nestingattempts1 < 1)*1 #did this bird make a nest? 1 = 0 nests, 0 = at least one nest

Nests$NestInit[is.na(Nests$NestInit)] <- 0
### Build Nimble Objects ####
nim.dat <- list(nest.init = Nests$NestInit,
                n.male.eggs = Nests$n.male.eggs,
                H = Nests$H_inc,
                nest.hatch = Nests$Hatch,
                n.poults = Nests$Total_Offspring,
                n.clutch = Nests$Total.Eggs,
                nested = nested,
                H.brood = Nests$H_brood,
                brood.lived = Nests$Brood_lived,
                z.male = z.male1,
                harvested = harvested2,
                z.female = z.female1,
                nestingattempts = nestingattempts1,
                broodingattempts = broodingattempts1,
                broodstolive = broodstolive1,
                neststobrood2 = neststobroodsecond,
                neststobrood1 = neststobroodfirst
)


Nests$Attempt[is.na(Nests$Attempt)] <- 1
nest.site <-  ifelse(Nests$Site == 'SRS', 0, 1)
nim.consts <- list(n.nest = nrow(Nests),
                   doy = Nests$Date_Incubate_s,
                   age = Nests$Age,
                   n.tested= Nests$n.tested,
                   nest.site = nest.site,
                   attempt = Nests$Attempt,
                   n.birds = n.birds,
                   birds.age = birds.age,
                   male.age = male.age1,
                   first.m = first.male1,
                   last.m = last.male1,
                   female.age = female.age1,
                   first.f = first.female1,
                   last.f = last.female1,
                   n.females = length(last.female1),
                   n.males = length(last.male1),
                   nestingattempts_nodeath1 = nestingattempts_nodeath1b, #first time attempts
                   nestingattempts_nodeath2 = nestingattempts_nodeath1a, #second and later attempts
                   broodattempts_nodeath = broodattempts_nodeath1,
                   site = site,
                   fem.site = fem.site1,
                   malesite = malesite1,
                   precip = as.vector(scale(precip))
)

### Initial values ####
#Apologies if you're reading this part of the code, it's terrible
## need to make some initial objects for H, nest.hatch, n.poults, n.clutch, est.livepoults 
H.init <- nim.dat$H
H.init[is.na(H.init)] <- 1 #start all nests missing Hen info as alive (probably nests that never initiated)
H.init[!is.na(nim.dat$H)] <- NA #can't override real values
nest.hatch.init <- nim.dat$nest.hatch
nest.hatch.init[is.na(nest.hatch.init)] <- 0 #start all unknown fate nests as failures 
nest.hatch.init <- ifelse(is.na(nim.dat$H), nest.hatch.init, ifelse(nim.dat$H == 0, 0,nest.hatch.init))
nest.hatch.init[!is.na(nim.dat$nest.hatch)] <- NA #can't override real values
n.poults.init <- nim.dat$n.poults
n.poults.init[is.na(nim.dat$n.poults)] <- 0 #start with no hatch
n.poults.init[!is.na(nim.dat$n.poults)] <- NA #can't override real values
n.clutch.init <- nim.dat$n.clutch
n.clutch.init[is.na(nim.dat$n.clutch)] <- 16 #start with 16 eggs
n.clutch.init[!is.na(nim.dat$n.clutch)] <- NA #can't override real values
H.brood.init <- nim.dat$H.brood
H.brood.init[] <- 0 
H.brood.init[!is.na(nim.dat$H.brood)] <- NA
brood.lived.init <- nim.dat$brood.lived
brood.lived.init[] <- 0
brood.lived.init[!is.na(nim.dat$brood.lived)] <- NA
est.livepoults.init <- nim.dat$est.livepoults
est.livepoults.init[] <- 0
est.livepoults.init[!is.na(nim.dat$est.livepoults)] <- NA
nestingattempts.init <- nim.dat$nestingattempts
nestingattempts.init[] <- 1
nestingattempts.init[!is.na(nim.dat$nestingattempts)] <- NA
nim.inits <- list(mean.fail = rnorm(1, 1, 1),
                  b1 = rnorm(4,0, 1),
                  b2 = rnorm(1,0,1),
                  b3 = rnorm(1,0,1),
                  mean.attempt = rnorm(1, 1, 1),
                  a1 = rnorm(1,0, 1),
                  a2 = rnorm(1, 0, 1),
                  p2 = rnorm(1, 0, 1),
                  mean.phi.incubate = rnorm(1, 1, 1),
                  p1 = rnorm(1,0, 1),
                  p.hatch = rbeta(1,1,1),
                  mean.brood.fail = rnorm(1, 1, 1),
                  c1 = rnorm(1,0, 1),
                  c2 = rnorm(1,0, 1),
                  c3 = rnorm(1, 0, 1),
                  mean.phi.brood = rnorm(1, 1, 1),
                  d1 = rnorm(1,0, 1),
                  d2 = rnorm(1,0, 1),
                  AvgClutch = c(12, 12),
                  H = H.init,
                  H.brood = H.brood.init,
                  brood.lived = brood.lived.init,
                  nest.hatch = nest.hatch.init,
                  n.poults = n.poults.init,
                  n.clutch = n.clutch.init,
                  f1 = rnorm(3,1,1),
                  f2 = rnorm(3,1,1),
                  m1 = rnorm(4,1,1),
                  m2 = rnorm(4,1,1),
                  mean.f.phi = rep(1,3),
                  mean.m.phi = rep(1,4),
                  avg.nesting = array(2, c(2,2)),
                  avg.brooding = array(2, c(2,2)),
                  nestingattempts = nestingattempts.init,
                  notnestingphi = c(.7, .8),
                  p.avg.broodfail = array(.5, 2),
                  p.broodfail.site2 = rnorm(1),
                  p.age.broodfail = rnorm(1),
                  p.avg.fail = c(.5),
                  p.fail.site = rnorm(1),
                  p.harvest = array(.6, c(2,2, 5)),
                  p.age.fail = rnorm(1,0,1),
                  p.attempt2 = .1,
                  sex.ratio = c(.3, .3),
                  male.winter = c(.8, .8),
                  a3 = rnorm(1)
)


## Create RDS object ####
#saveRDS(list(nim.dat = nim.dat, nim.inits = nim.inits, nim.consts = nim.consts), 'ComboIPM.rds')


### Misc Graphs ####

### Hunting ####
hunting <- data.frame(Year = 2010:2018,
                      Jasper = c(0.7, 0.4, 0.9, 0.6, 0.4, 0.5, 0.7, 0.8 ,0.5),
                      Hampton = c(0.9, 1.0, 0.8, 0.5, 0.8, 0.9, 1.1, 1.0, 1.1))
hunting <- hunting %>% group_by(Year) %>% mutate(Avg = mean((Jasper + Hampton)/2))


## percent taken when? 
hunt2 <- data.frame(Year = 2018:2010, 
                    March20 = c(40, 40, 41, 0, 0, 0, 0 ,0 ,0),
                    April1 = c(18, 22, 20, 42, 44, 42, 50, 45, 48),
                    April8 = c(13, 14, 13, 22, 22, 23, 21, 22, 21),
                    April15 = c(11, 10, 11, 18, 18, 18, 18, 19, 18)
)
hunt2$May5 <- 100- rowSums(hunt2[,-1])
hunt3 <- hunt2 %>% group_by(Year) %>% summarize(All = sum(March20, April1, April8, April15),
                                                Pre_Lay = sum(March20, April1))


gghunt <- data.frame(Year = rep(hunt2$Year,5),
                     Harvested = c(hunt2$March20, hunt2$April1, hunt2$April8, hunt2$April15, hunt2$May5),
                     Group = rep(c('March 20', 'April 1', 'April 8', 'April 15', 'May 5'), each = 9))


gghunt <- gghunt %>%
  mutate(Group = as.Date(Group, format = "%B %d")) %>% 
  group_by(Year) %>%
  mutate(cum_harvest = cumsum(Harvested))

gghunt$Year <- as.character(gghunt$Year)
gghunt$Season <- ifelse(gghunt$Year < 2016, 'April 1st open', 'March 20 open')


aaaa <- ggplot(gghunt, aes(x = Group, y = cum_harvest, col  = Year, group = Year))+
  geom_line(lwd = 1)+
  geom_point(cex = 3)+
  theme_minimal()+
  #scale_color_manual(values = c('#461220', '#8C2F39', '#B23A48','#FCB9B2', '#FED0BB'))+
  scale_color_manual(values = c('#2E0B16', '#461220', '#6E1F2E', '#8C2F39', '#A73542', '#B23A48', '#D96C6F', '#FCB9B2', '#FED0BB', '#FFE5D9'))+
  geom_vline(xintercept = as.Date('04-22', format = '%m-%d'), lty = 3)+ #mean incubation
  geom_vline(xintercept = as.Date('04-10', format = '%m-%d'), lty = 2)+ #mean lay date
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        strip.text = element_text(size = 20),
        legend.position = 'bottom',
        plot.margin = margin(t = 40, r = 5, b = 5, l = 5))+
  ylab('Cummulative % harvest')+
  xlab('Date')

hunting2 <- data.frame(Year = rep(hunting$Year, 3),
                       Harvest = c(hunting$Jasper, hunting$Hampton, hunting$Avg),
                       County = rep(c('Jasper', 'Hampton', 'Combined'), each = 9))

hunting2$Harvest_sqkm <- hunting2$Harvest/2.58999

bbbb <- ggplot(hunting2, aes(x = Year, y = Harvest_sqkm, col  = County, group = County))+
  geom_line(lwd = 1)+
  geom_point(cex = 3)+
  theme_minimal()+
  ylim(0.1, 0.5)+
  scale_color_manual(values = c('#51A3A3','#CB904D', '#75485E'))+
  geom_vline(xintercept = 2014, lty = 1, col = 'grey30')+ 
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        strip.text = element_text(size = 20),
        plot.margin = margin(t = 40, r = 5, b = 5, l = 5),
        legend.position = 'bottom')+
  ylab(expression('Male harvest per km'^2))+
  xlab('Date')

ggarrange(bbbb, aaaa, nrow= 2, labels = 'AUTO', hjust = c(-2, -2), vjust= c(1.25, -1))

### Figure 1 ####
me1 <- Nests2 %>% filter(NestInit == 1 & Attempt == 1)
me2 <- Nests2 %>% filter(NestInit == 1 & Attempt == 1 & Age == 1)
mean(me1$Date_Incubate, na.rm = T)
mean(me2$Date_Incubate, na.rm = T)

median((Nests2 %>% filter(NestInit == 1 & Attempt == 1 & Site == "SRS"))$Date_Incubate, na.rm = T)
median((Nests2 %>% filter(NestInit == 1 & Attempt == 1 & Site != "SRS"))$Date_Incubate, na.rm = T)
range((Nests2 %>% filter(NestInit == 1 & Attempt == 1 & Site == "SRS"))$Date_Incubate, na.rm = T)
range((Nests2 %>% filter(NestInit == 1 & Attempt == 1 & Site != "SRS"))$Date_Incubate, na.rm = T)

# #when were opening dates public and private?
#2021 - no hunt (covid)
#April 30 2022- 15 people, 7 turkeys
#April 21-22 2023, 15 hunters, 19 turkeys
#April 19-20, 2024
#April 15, 2025 - 24 people, 16 turkeys

openingdays <- data.frame(OpenDay = c('04/01/2014','04/01/2015','03/20/2016',
                                      '03/20/2017','03/20/2018', '04/30/2022', '04/21/2023', '04/19/2024', '04/15/2025'),
                          Site = c(rep('Webb Complex', 5), rep('SRS', 4)))
openingdays$Day <- yday(as.Date(openingdays$OpenDay, format = "%m/%d/%Y"))
openingdays$Year <- year(as.Date(openingdays$OpenDay, format = "%m/%d/%Y"))

me1 %>% group_by(Site) %>% summarize(mean(Date_Incubate, na.rm = T), 
                                     max(Date_Incubate, na.rm = T), 
                                     min(Date_Incubate, na.rm = T))
me2 %>% group_by(Site) %>% summarize(mean(Date_Incubate, na.rm = T))
hh <- data.frame(Site = c('SRS', 'Webb Complex'), avg.date = c(108, 112))
h2 <- data.frame(Site = c('SRS', 'Webb Complex'), avg.date = c(107, 111))

hh$Site <- factor(hh$Site, levels = c('Webb Complex', 'SRS'))
h2$Site <- factor(h2$Site, levels = c('Webb Complex', 'SRS'))
me1$Site <- factor(me1$Site, levels = c('Webb Complex', 'SRS'))
me2$Site <- factor(me2$Site, levels = c('Webb Complex', 'SRS'))
openingdays$Site <- factor(openingdays$Site, levels = c('Webb Complex', 'SRS'))
aa <- ggplot(data = me1, aes(x = Year, group = Year))+
  geom_boxplot(aes(y = Date_Incubate), width = .5)+
  geom_point(data = openingdays, aes(y = Day, colour = 'Season opens'), cex = 3, pch = 19, alpha = .7)+
  theme_bw()+
  facet_wrap(~Site, scale = 'free_x')+
  ylab('Incubation date')+
  geom_hline(data = hh, aes(yintercept = avg.date), lty = 2, lwd = .2)+
  scale_y_continuous(breaks = c(50, 75, 100, 125, 150), limits = c(70, 165), labels = c('', 'March 16', 'April 10', 'May 5', 'May 30'))+
  scale_fill_brewer()+
  ggtitle('All first attempts')+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 0),
        strip.text = element_text(size = 20),
        legend.position = 'inside',
        legend.position.inside = c(.875, .9))

bb <- 
  ggplot(data = me2, aes(x = Year, group = Year))+
  geom_boxplot(aes(y = Date_Incubate), width = .5)+
  geom_point(data = openingdays, aes(y = Day, colour = 'Season opens'), cex = 3, pch = 19, alpha = .7)+
  theme_bw()+
  facet_wrap(~Site, scale = 'free_x')+
  ylab('Incubation date')+
  scale_y_continuous(breaks = c(50, 75, 100, 125, 150), limits = c(70, 165), labels = c('', 'March 16', 'April 10', 'May 5', 'May 30'))+
  geom_hline(data = h2, aes(yintercept = avg.date), lty = 2, lwd = .2)+
  scale_fill_brewer()+
  ggtitle('Adult first attempts')+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 0),
        strip.text = element_text(size = 20),
        legend.position = 'none')

ggarrange(aa, bb, nrow = 2)

## Figure S1 ####
firstonly <- subset(nesting_birds, nesting_birds$Attempt == 1)

nesting_info <- table(firstonly$Year, firstonly$Nested, firstonly$Age)
Juvis <- as.matrix(nesting_info[,,1])
Juvis <- cbind(Juvis, Juvis[,2]/rowSums(Juvis))
Adults <- as.matrix(nesting_info[,,2])
Adults <- cbind(Adults, Adults[,2]/rowSums(Adults))

Raw_success <- Nests %>% group_by(Year, Site) %>% filter(!is.na(Hatch)) %>% summarize(mean(Hatch))
raw_attempts <- firstonly %>% group_by(Year, Age, Site) %>% summarize(Initiated = sum(Nested, na.rm = T)/n())
raw_attempts <- raw_attempts[-10,]
raw_attempts2 <- firstonly %>% group_by(Year, Site) %>% summarize(Initiated = sum(Nested, na.rm = T)/n())
raw_attempts2$Age <- 'All'
raw_attempts_g <- rbind(as.data.frame(raw_attempts), as.data.frame(raw_attempts2))
raw_attempts_g$Site <- factor(raw_attempts_g$Site , levels = c('Webb Complex', 'SRS'))
raw_attempts_g$Age <- ifelse(raw_attempts_g$Age == 0, 'Juvenile', ifelse(raw_attempts_g$Age == '1', 'Adult', 'All'))
(a <- ggplot(raw_attempts_g[raw_attempts_g$Age != 'All',], aes(x = Year, y = Initiated, col = Age))+
    geom_point(cex = 4)+
    geom_line(lwd = 1)+
    theme_bw()+
    ylab('Nest initiation')+
    facet_wrap(~Site, scales = 'free_x')+
    scale_color_manual(values = c('#92140c', '#393e41'))+
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          strip.text = element_text(size = 20), 
          legend.position = 'none',
          panel.spacing.x = unit(1.5, "lines")))

Nests$Site <- factor(Nests$Site, c('Webb Complex', 'SRS'))
first_nest <- Nests %>% group_by(Year, Age, Site) %>% filter(!is.na(Hatch) & NestInit == 1 & Attempt == 1) %>% summarize(HatchRate = sum(Hatch, na.rm = T)/n())
first_nest2 <- Nests %>% group_by(Year, Site) %>% filter(!is.na(Hatch) & NestInit == 1 & Attempt == 1) %>% summarize(HatchRate = sum(Hatch, na.rm = T)/n())#, nn = n())
first_nest2$Age <- 'All'

first_nest$Age <- ifelse(first_nest$Age == 0, 'Juvenile', 'Adult')
(b <- ggplot(as.data.frame(first_nest), aes(x = Year, y = HatchRate, col = Age))+
    geom_point(cex = 4)+
    geom_line(lwd = 1)+
    theme_bw()+
    facet_wrap(~Site, scales = 'free_x')+
    ylab('Nest success')+
    scale_color_manual(values = c('#92140c', '#393e41'))+
    theme(axis.text = element_text(size = 18),
          axis.title = element_text(size = 22),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 20),
          strip.text = element_text(size = 20), 
          legend.position = 'none',  
          panel.spacing.x = unit(1.5, "lines")))



first_brood <- Nests %>% group_by(Year, Age, Site) %>% filter(!is.na(Brood_lived)) %>% summarize(BroodPhi = sum(Brood_lived, na.rm = T)/n())
first_brood2 <- Nests %>% group_by(Year, Site) %>% filter(!is.na(Brood_lived)) %>% summarize(BroodPhi = sum(Brood_lived, na.rm = T)/n())
first_brood2$Age <- 'All'

first_brood$Age <- ifelse(first_brood$Age == 0, 'Juvenile', 'Adult')

c <- ggplot(as.data.frame(first_brood), aes(x = Year, y = BroodPhi, col = Age))+
  geom_point(cex = 4)+
  geom_line(lwd = 1)+
  theme_bw()+
  ylim(0, 1)+
  facet_wrap(~Site, scales = 'free_x')+
  ylab('Brood success')+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        strip.text = element_text(size = 20),
        panel.spacing.x = unit(1.5, "lines"),
        legend.position = 'none')+
  #legend.position = c(.4, .65))+
  scale_color_manual(values = c('#92140c', '#393e41'))

ggarrange(a,b,c, nrow = 3, common.legend = T, legend = 'bottom')

### Precipitation vs initiation ####
precip <- data.frame(year = c(rep(2014:2018, each = 2), rep(2021:2025, each = 2)),
                     precip = c(5.06, 3.16, 5.34, 1.28, 1.87, 8.92, 3.16, 7.37,3.87,6.09,
                                1.66, 2.96, 4.58,  3.94, 6.74 , 4.78, 2.1, 3.11, 2.31, 8.61),
                     month = rep(4:5, 10)
)

precip2 <- precip %>% group_by(year) %>% summarize(totP = sum(precip),
                                                   meanP = mean(precip))
plot(precip2$year, precip2$totP, type = 'l')
plot(precip2$year, precip2$meanP, type = 'l')

AdultInt <- raw_attempts_g[raw_attempts_g$Age == 'Adult', ]
plot(AdultInt$Year, AdultInt$Initiated, type = 'l', ylim = c(0, 6))
lines(precip2$year, precip2$meanP)
cor(AdultInt$Initiated, precip2$totP)^2
cor(AdultInt$Initiated, precip2$meanP)^2

gg_precip <- data.frame(precip = precip2$totP*25.4,
                        year = precip2$year,
                        Initiation = AdultInt$Initiated,
                        Site = AdultInt$Site)

ggplot(gg_precip, aes(x = Initiation, y = precip, col = Site))+
  geom_point(pch = 19, cex= 3)+
  theme_classic()+
  ylab('Total precipitation (mm)')+
  xlab('% adult nest initation')

t.test(precip2$totP[1:5], precip2$totP[6:10])


### Prior for male winter survival ####

malewinter_literature <- 
  data.frame(phi = c(0.62, 0.71, 0.704,0.74, .917, .792, 0.97^5, 0.97^5, 0.90),
             se = c(0.12,0.11, .074,.05, .04, .03,0.01,0.01,0.015),
             age= c(0,1,1,1, 0, 1, 0, 1,1),
             reference = c('Norman 2004','Norman 2004', 
                           'Humberg 2009', 'Grisham 2008', 
                           'Holdstock 2010', 'Holdstock 2010',
                           'Collier 2010', 'Collier 2010',
                           'Norman 2022')
  )

(alpha3 <- malewinter_literature$phi^2 * ((1-malewinter_literature$phi)/malewinter_literature$se^2 - (1/malewinter_literature$phi)))
(beta3 <- alpha3 * (1/malewinter_literature$phi - 1))
nSims <- 5000
## Empty matrix to store the simulated survival probabilities
sim_phi_j <- matrix(NA, nrow = nSims, ncol = sum(malewinter_literature$age == 0))
sim_phi_a <- matrix(NA, nrow = nSims, ncol = sum(malewinter_literature$age == 1))
for(j in 1:ncol(sim_phi_j)){
  sim_phi_j[,j] <- rbeta(n = nSims, alpha3[malewinter_literature$age == 0][j], beta3[malewinter_literature$age == 0][j])
}
for(j in 1:ncol(sim_phi_a)){
  sim_phi_a[,j] <- rbeta(n = nSims, alpha3[malewinter_literature$age == 1][j], beta3[malewinter_literature$age == 1][j])
}

mu3_j <- mean(sim_phi_j)
var3_j <- sd(sim_phi_j) ^ 2
mu3_a <- mean(sim_phi_a)
var3_a <- sd(sim_phi_a) ^ 2

alpha_j <- mu3_j^2 * ((1-mu3_j)/var3_j - (1/mu3_j))
beta_j <- alpha_j * (1/mu3_j - 1)
alpha_a <- mu3_a^2 * ((1-mu3_a)/var3_a - (1/mu3_a))
beta_a <- alpha_a * (1/mu3_a - 1)

## Visualize briefly:
phi <- seq(0, 1, by = 0.001)
visual <- data.frame(phi = c(phi, phi),
                     density = c(dbeta(phi, alpha_j, beta_j),dbeta(phi, alpha_a, beta_a)),
                     Approach = rep(c('J', 'A'), each = length(phi)))

ggplot(visual, aes(x = phi, y= density, col = Approach))+
  geom_line()

#for adults: alpha = 13.89843, beta = 3.830493
#for juvis: alpha = 5.113408, beta = 1.285857
