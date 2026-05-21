#####################################################--
#Differences in wild turkey nest initiation and population growth rates between a hunted and nonhunted site in South Carolina, USA   
# This script loads cleaned data and runs the IPM
################################--
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to source file location

library(nimble)
library(coda)
library(parallel)
library(MCMCvis)

### Code ####
Turkey_IPM <- nimbleCode({
  ##* Nests ######
  ## Priors:
  mean.attempt ~ dnorm(1, 1) #first nest attempt probability
  a1 ~ dnorm(0, sd = 1) #impact of being an adult on initiation probability
  a2 ~ dnorm(0, sd = 1) #impact of being on site 2 on nesting probability 
  a3 ~ dnorm(0, sd = 1) #precipitation influence on nesting prob
  mean.phi.incubate ~ dnorm(1, 1) #surviving the incubation period pre-reg change
  p1 ~ dnorm(0, sd = 1) #effect of being an adult on surviving incubation 
  p2 ~ dnorm(0, sd = 1) #effect of site 2 on incubation survival 
  p.hatch ~ dbeta(1, 1) #P(hatching from surviving laid egg)
  sex.ratio[1] ~ dbeta(1, 1)  #sex ratio of tested eggs site 1; assuming no bias in testing
  sex.ratio[2] ~ dbeta(40, 62)  #sex ratio of Webb; 
  AvgClutch[1] ~ dgamma(shape = 6.2, scale = 2) #mean of 12 ish eggs
  AvgClutch[2] ~ dgamma(shape = 6.2, scale = 2) #mean of 12 ish eggs
  mean.fail ~ dnorm(0, 1) #average failure is around 80% 
  b1[1] ~ dnorm(0, 1) #impact of doy on nest failure, logit scale
  b1[2] ~ dnorm(0, 1) #impact of doy on nest failure #2, logit scale
  b1[3] <- b1[2] #impact of doy on 3rd nest failure, logit scale
  b1[4] <- b1[2] #impact of doy on 4th nest failure, logit scale
  b2 ~ dnorm(0, 1) #impact of site on nest failure
  b3 ~ dnorm(0, 1) #impact of age on nest
  mean.brood.fail ~ dnorm(1, 1)
  c1 ~ dnorm(0, 1) #impact of age on brood failure, logit scale
  c2 ~ dnorm(0, 1) #impact of site on brood failure, logit scale
  c3 ~ dnorm(0, 1) #impact of doy (for hatch) on brood failure, logit scale
  mean.phi.brood ~ dnorm(1, 1)  
  d1 ~ dnorm(0, 1) #impact of age on surviving brood stage, logit scale
  d2 ~ dnorm(0, 1) #impact of site on surviving brood stage, logit scale
  
  for(b in 1:n.birds){ #hens
    logit(p.attempt[b]) <- mean.attempt + a1*birds.age[b] + a2*site[b] + a3*precip[b] #site is binary
    nested[b] ~ dbern(p.attempt[b]) #did this bird make at least 1 nest attempt?
  }
  
  for(i in 1:n.nest){ #per nest
    logit(phi.incubate[i]) <- mean.phi.incubate + p1*age[i] + p2*(nest.site[i]) #site is binary
    H[i] ~ dbern(phi.incubate[i]) #hen survives from laying to hatching
    logit(p.fail[i]) <- mean.fail + b1[attempt[i]]*doy[i] + b2*(nest.site[i]) + b3*(age[i]) #day of year impacts likelihood of failure
    nest.hatch[i] ~ dbern(H[i]*(1-p.fail[i])*nest.init[i]) #survives to hatching? conditional on hen being alive [H], not being abandoned/predated, and nest being initiated
    n.clutch[i] ~ dpois(AvgClutch[age[i]+1])
    n.poults[i] ~ dbinom(size = n.clutch[i], prob = p.hatch*nest.hatch[i]) #can't hatch if your nest is dead
    n.male.eggs[i] ~ dbinom(size = n.tested[i], prob = sex.ratio[nest.site[i]+1])
    logit(brood.fail[i]) <- mean.brood.fail + c1*age[i] + c2*nest.site[i] + c3*doy[i] #based on age and site
    logit(phi.brood[i]) <- mean.phi.brood + d1*age[i] + d2*nest.site[i]
    H.brood[i] ~ dbern(phi.brood[i]) #hen survives from hatching to 28 days
    brood.lived[i] ~ dbern(H.brood[i]*(1-brood.fail[i])*nest.hatch[i]) #have to have mom live, brood can't fail, and nest had to have hatched for brood to succeed
  }
  
  ### *Hen part ####
  ## We can't re-use information from the nesting part so need to be cautious
  ## Priors for female part of model
  mean.f.phi[1] ~ dnorm(0, 1) #survival from beginning of year to March
  mean.f.phi[2] ~ dnorm(0, 1) #survival through nesting season (Aug) assuming you nest
  mean.f.phi[3] ~ dnorm(0, 1) #survival Aug to Jan
  notnestingphi[1] ~ dbeta(1, 1) #surviving nesting season given you don't nest, juvi
  notnestingphi[2] ~ dbeta(1, 1) #surviving nesting season given you don't nest, adult
  for(jj in 1:3){
    f1[jj] ~ dnorm(0, 1)
    f2[jj] ~ dnorm(0, 1)
  }
  avg.nesting[1,1] ~ dunif(0, 2) 
  avg.nesting[1,2] ~ dunif(0, 2) #adults site 1
  avg.nesting[2,1] ~ dunif(0, 2) 
  avg.nesting[2,2] ~ dunif(0, 2) #adults site 2
  avg.brooding[1,1] ~ dunif(0, 2)
  avg.brooding[2,1]  ~ dunif(0, 2) #adults site 1
  avg.brooding[1,2] ~ dunif(0, 2) 
  avg.brooding[2,2]  ~ dunif(0, 2) #adults site 2
  p.avg.broodfail[1] ~ dnorm(0, 1) #average probability of brood failing (not due to death); juvis
  p.avg.broodfail[2] ~ dnorm(0, 1) #average probability of brood failing (not due to death); adults
  p.broodfail.site2 ~ dnorm(0, 1)
  p.avg.fail ~ dnorm(0, 1) #intercept of nest failing (not due to death), logit scale
  p.fail.site ~ dnorm(0, 1) #intercept of nest failing (not due to death), logit scale, site 2
  p.attempt2 ~ dnorm(0, 1) #impact of this being a 2nd (or 3rd) attempt on nest failure prob
  p.age.fail ~ dnorm(0, 1) #effect of being an adult on average nest failure probability 
  p.age.broodfail ~ dnorm(0, 1) #effect of being an adult on brood failure 
  
  for(k in 1:n.females){
    logit(p.nest.fail1[k]) <- mean.fail + b3*female.age[k] + b2*fem.site[k] #average pre/post reg,
    logit(p.nest.fail2[k]) <- mean.fail + b3*female.age[k] + p.attempt2+  b2*fem.site[k] #average pre/post reg,
    logit(p.brood.fail[k]) <- mean.brood.fail + c1*female.age[k] + c2*fem.site[k]
    nestingattempts[k] ~ dpois(avg.nesting[fem.site[k]+1,female.age[k]+1])
    broodingattempts[k] ~ dpois(avg.brooding[fem.site[k]+1,female.age[k]+1]*(nestingattempts[k]>0)) #can't brood if you never nested
    neststobrood1[k] ~ dbinom(size = nestingattempts_nodeath1[k], p = 1- p.nest.fail1[k]) #how many nests succeeded, but only of the ones where the hen didn't die
    neststobrood2[k] ~ dbinom(size = nestingattempts_nodeath2[k], p = 1- p.nest.fail2[k]) #how many nests succeeded, but only of the ones where the hen didn't die
    broodstolive[k] ~ dbinom(size = broodattempts_nodeath[k], p = 1- p.brood.fail[k]) #how many brood attempts succeeded, but only of the ones where the hen didn't die
    zeronest[k] <- (nestingattempts[k] < 1)*1
    logit(femalephi[k, 1]) <- mean.f.phi[1] + f1[1]*female.age[k] + f2[1]*fem.site[k] #pre-nesting
    logit(femalephi[k, 3]) <- mean.f.phi[3] + f1[3]*female.age[k] 
    
    logit(femalenestingphi[k]) <- mean.phi.incubate + p1*female.age[k] + p2*(fem.site[k])
    logit(femalebroodingphi[k]) <-  mean.phi.brood + d1*female.age[k] + d2*fem.site[k]
    ## In order to survive nesting season, a bird needs to either:
    # 1. not nest and survive not nesting
    # 2. Nest some quantity of times and brood some quantity of times and survive each of those nesting and brooding attempts
    femalephi[k,2] <- (notnestingphi[female.age[k]+1]*zeronest[k]) + (1-zeronest[k])*((femalenestingphi[k])^nestingattempts[k])*((femalebroodingphi[k])^broodingattempts[k])
    for(t in (first.f[k]+1):last.f[k]){
      z.female[k,t] ~ dbern(femalephi[k, t-1]*z.female[k,t-1])
    }
  }
  
  #P surviving a year:
  for(r in 1:2){ #the two sites
    logit(juvi[r,1]) <- mean.f.phi[1] + f1[1]*0 + f2[1]*(r-1) #juvi
    logit(adult[r,1]) <- mean.f.phi[1] + f1[1]*1 + f2[1]*(r-1) #adult
    logit(juvi[r,3]) <- mean.f.phi[3] + f1[3]*0 
    logit(adult[r,3]) <- mean.f.phi[3] + f1[3]*1 
    logit(avgnestingphi.juvi[r]) <- mean.phi.incubate + p1*0 + p2*(r-1) #juvi surviving incubation
    logit(avgnestingphi.adult[r]) <- mean.phi.incubate + p1*1 + p2*(r-1) #adult
    logit(avgbroodingphi.juvi[r]) <- mean.phi.brood + d1*0 + d2*(r - 1)#juvi
    logit(avgbroodingphi.adult[r]) <- mean.phi.brood + d1*1 + d2*(r - 1)#adult
    
    logit(nestattempted[r, 1]) <- mean.attempt + a2*(r-1) #p(attempted a nest) juvis
    logit(nestattempted[r, 2]) <- mean.attempt + a1 + a2*(r-1)#p(attempted a nest) adults
    
    pTwoNest[r,1] <- ppois(1, avg.nesting[r,1], lower.tail = FALSE) #should give p(#nestsattempted >1)
    #should give p(#nestsattempted <=1) - p(0); aka p(1 attempt)
    pTwoNest[r,2] <- ppois(1, avg.nesting[r,2], lower.tail = FALSE) #should give p(#nestsattempted >1)
    
    logit(juvi[r,2]) <- notnestingphi[1]*(1-nestattempted[r,1]) + nestattempted[r,1]*((avgnestingphi.juvi[r])^avg.nesting[r,1])*((avgbroodingphi.juvi[r])^avg.brooding[1,r])
    logit(adult[r,2]) <- notnestingphi[2]*(1-nestattempted[r,2]) + nestattempted[r,2]*((avgnestingphi.adult[r])^avg.nesting[r,2])*((avgbroodingphi.adult[r])^avg.brooding[2,r])
    
    juvi.yr[r] <- juvi[r,1]*juvi[r,2]*juvi[r,3]
    adult.yr[r] <- adult[r,1]*adult[r,2]*adult[r,3]
    
    ## Some nest statistics for Leslie Matrix part
   
    logit(brood.fails.juvi[r]) <- mean.brood.fail + c2*(r-1)
    logit(brood.fails.adult[r]) <-mean.brood.fail + c2*(r-1) + c1 #adult brood fails but not b/c of death
    exp_f[r,1] <- AvgClutch[1]*p.hatch*(1-sex.ratio[r])
    exp_f[r,2] <- AvgClutch[2]*p.hatch*(1-sex.ratio[r])
    
    
    #p.avg.fail + p.age.fail*female.age[k] + p.fail.site*fem.site[k]
    logit(avg.fail1[r,1]) <- mean.fail + b3*(r-1)
    logit(avg.fail1[r,2]) <- mean.fail + b2 + b3*(r-1)
    logit(avg.fail2[r,1]) <- mean.fail+ p.attempt2 + b3*(r-1)#2nd attempt
    logit(avg.fail2[r,2]) <- mean.fail + b2 + p.attempt2 + b3*(r-1) #2nd attempt
    
    ### **Estimating Poult output for Leslie Matrix ####
    #hen has to survive through nesting season (time 1 and time 2)
    #then attempt a nest and succeed or fail first nest, attempt 2nd nest and succeed
    # then create some poults that survive 28 days
    nest.p.j[r] <- nestattempted[r,1]*(1-avg.fail1[r,1]) + nestattempted[r,1]*(avg.fail1[r,1])*pTwoNest[r,1]*(1-avg.fail2[r,1])
    nest.p.a[r] <- nestattempted[r,2]*(1-avg.fail1[r,2]) + nestattempted[r,2]*(avg.fail1[r,2])*pTwoNest[r,2]*(1-avg.fail2[r,2])
    juvi.poults[r] <- juvi[r,1]*juvi[r,2]*nest.p.j[r]*exp_f[r,1]*(1-brood.fails.juvi[r])
    adult.poults[r] <- adult[r,1]*adult[r,2]*nest.p.a[r]*exp_f[r,2]*(1-brood.fails.adult[r])
    
    ### for prob of a nesting juv producing a 28 day poult:
    npj[r] <- ((1-avg.fail1[r,1]) + (avg.fail1[r,1])*pTwoNest[r,1]*(1-avg.fail2[r,1]))*(1-(brood.fails.juvi[r]+1-avgbroodingphi.juvi[r]))
    npa[r] <- ((1-avg.fail1[r,2]) + (avg.fail1[r,2])*pTwoNest[r,2]*(1-avg.fail2[r,2]))*(1-(brood.fails.adult[r]+1-avgbroodingphi.adult[r]))
  }
  
  ### *Male bird part of model ####
  
  ## Priors for male part of model
  mean.m.phi[1] ~ dnorm(0, 1) #survival from beginning of year to hunting season
  mean.m.phi[2] ~ dnorm(0, 1) #survival through hunting season
  mean.m.phi[3] ~ dnorm(0, 1) #survival post hunting to Aug
  ## Since we lack this information (for now):
  mean.m.phi[4] ~ dnorm(0, 1) #survival Aug to Jan
  male.winter[1] ~ dbeta(18, 2)#prior for juveniles from winter survival literature
  male.winter[2] ~ dbeta(18, 2) #prior for adults from winter survival literature
  
  #gives us winter priors of ~0.86
  
  for(jj in 1:4){
    m1[jj] ~ dnorm(0, 1) #impact of age on survival 
    m2[jj] ~ dnorm(0, 1) #impact of site on survival
  }
  
  p.harvest[2,1, 2] ~ dbeta(1, 1) #proportion of young males that die to hunting vs. other causes during hunting season (for Webb)
  p.harvest[1,1, 2] <- 0 #proportion of young males that die to hunting vs. other causes during hunting season
  for(kk in 1:2){ #kk here is site
    p.harvest[kk,2, 2] ~ dbeta(1, 1) #proportion of old males that die to hunting vs. other causes during hunting season
    p.harvest[kk,1,1] <- 0
    p.harvest[kk,1,3] <- 0
    p.harvest[kk,1,4] <- 0
    p.harvest[kk,1,5] <- 0
    p.harvest[kk,2,1] <- 0
    p.harvest[kk,2,3] <- 0
    p.harvest[kk,2,4] <- 0
    p.harvest[kk,2,5] <- 0
  }
  
  for(j in 1:n.males){
    logit(malephi[j, 1]) <- mean.m.phi[1] + m1[1]*male.age[j] + m2[1]*malesite[j]
    logit(malephi[j, 2]) <- mean.m.phi[2] + m1[2]*male.age[j] + m2[2]*malesite[j] #Non-harvest survival
    logit(malephi[j, 3]) <- mean.m.phi[3] + m1[3]*male.age[j] + m2[3]*malesite[j]
    malephi[j,4] <- male.winter[male.age[j]+1]
    
    for(t in (first.m[j]+1):last.m[j]){ #stages of the year (1 = Jan to pre-hunt, 2 = hunting season, 3 = post-hunt to Aug, 4= Aug to Jan); we almost never observe #4
      z.male[j,t] ~ dbern(malephi[j, t-1]*z.male[j,t-1]*(1-p.harvest[malesite[j]+1, male.age[j]+1, t-1]))
      ### if males died during hunting season, was it hunting or natural?
      harvested[j,t] ~ dbern(malephi[j, 2]*z.male[j,t]*p.harvest[malesite[j]+1,male.age[j]+1,t]) #had to have been alive in time 2 first of all
    }
  }
  malsite[1] <- 1 #lazy coding sites
  malsite[2] <- 2
  for(r in 1:2){
    logit(juvi.m[r,1]) <- mean.m.phi[1] + m1[1]*0  + m2[1]*(malsite[r]-1) #juvi
    logit(adult.m[r,1]) <- mean.m.phi[1] + m1[1]*1 + m2[1]*(malsite[r]-1) #adult
    logit(juvi.m_p[r,2]) <- mean.m.phi[2] + m1[2]*0 + m2[2]*(malsite[r]-1) #juvi non-hunting survival
    logit(adult.m_p[r,2]) <- mean.m.phi[2] + m1[2]*1 + m2[2]*(malsite[r]-1) #adult non-hunting survival
    juvi.m[r,2] <- juvi.m_p[r,2]*(1-p.harvest[r,1,2])
    adult.m[r,2] <- adult.m_p[r,2]*(1-p.harvest[r,2,2])
    logit(juvi.m[r,3]) <- mean.m.phi[3] + m1[3]*0 + m2[3]*(malsite[r]-1) #juvi
    logit(adult.m[r,3]) <- mean.m.phi[3] + m1[3]*1 + m2[3]*(malsite[r]-1) #adult 
    juvi.m[r,4] <- male.winter[1]
    adult.m[r,4] <- male.winter[2]
    juvi.m.yr[r] <- juvi.m[r,1]*juvi.m[r,2]*juvi.m[r,3]*juvi.m[r,4]
    adult.m.yr[r] <- adult.m[r,1]*adult.m[r,2]*adult.m[r,3]*adult.m[r,4]
  }
  
  
  ### *Females Leslie Matrix ####
  ## Will calculate eigenvalues using posterior samples outside of model:
  #poult survival is unknown, at least by me, so we will put it as the same as juvenile survival 
  for(r in 1:2){ #pre and post regulation change
    poult_surv[r] <- juvi[r,3]
    LM[1,1,r] <- 0 # babies can't stay babies
    LM[2,1,r] <- poult_surv[r]  # female babies survive to become juvis
    LM[3,1,r] <- 0  # female babies can't become adults
    LM[1,2,r] <- juvi.poults[r]     # juvis making babies; survival of parents already baked in
    LM[2,2,r] <- 0 #juvis can't stay juvis
    LM[3,2,r] <- juvi.yr[r] #juvis survive to adults
    LM[1,3,r] <- adult.poults[r]    # adults making babies; survival of parents already baked in
    LM[2,3,r] <- 0 #can't get younger
    LM[3,3,r] <- adult.yr[r] #adults survive to adults next year
  }
  
})

### Grab objects ####
nim.stuff <- readRDS('./ModelObjects/ComboIPM.rds')
nim.dat <- nim.stuff$nim.dat
nim.consts <- nim.stuff$nim.consts
nim.inits <- nim.stuff$nim.inits

params <- c('p.hatch','AvgClutch', 'b1', 'mean.fail', 'mean.attempt', 'a1', 'a2', 'a3',
            'mean.phi.incubate', 'p1', 'p2', 'b2','b3', 'mean.brood.fail', 'c1', 'c2', 'c3', 'd2',
            'mean.phi.brood', 'd1', 'mean.f.phi', 'f1', 'f2', 'm1', 'm2',
            'mean.m.phi', 'juvi.poults', 'adult.poults', 'juvi', 
            'juvi.m', 'adult.m', 'adult', 'juvi.yr', 
            'adult.yr', 'juvi.m.yr', 'adult.m.yr', 'LM', 'avg.nesting', 
            'avg.brooding', 'exp_f', 'notnestingphi', 'avgnestingphi.juvi',
            'avgnestingphi.adult', 'avgbroodingphi.juvi', 'avgbroodingphi.adult', 
            'p.avg.broodfail', 'p.avg.fail', 'nestattempted', 'brood.fails.juvi', 
            'p.fail.site', 'p.broodfail.site2', 'brood.fails.adult', 'p.harvest', 
            'p.age.fail', 'p.age.broodfail', 'avg.fail1','avg.fail2', 
            'p.attempt2', 'pTwoNest', 'nest.p.a', 'sex.ratio',
            'nest.p.j', 'p.age.fail', 'npj', 'npa')


### Run Model ####
cl <- makeCluster(3)
clusterExport(cl = cl, varlist = c('nim.dat', "params", "Turkey_IPM", 'nim.consts', 'nim.inits'))
system.time(
  nim.out <- clusterEvalQ(cl = cl,{
    library(nimble)
    library(coda)
    
    prepnim <- nimbleModel(code = Turkey_IPM, constants = nim.consts,
                           data = nim.dat, inits = nim.inits, calculate = T)
    prepnim$initializeInfo()
    prepnim$calculate()
    mcmcnim <- configureMCMC(prepnim, monitors = params, print = T)
    nimMCMC <- buildMCMC(mcmcnim) #actually build the code for those samplers
    Cmodel <- compileNimble(prepnim) #compiling the model itself in C++;
    Compnim <- compileNimble(nimMCMC, project = prepnim) # compile the samplers next
    Compnim$run(niter = 10000, nburnin = 1000, thin = 1)
    return(as.mcmc(as.matrix(Compnim$mvSamples)))
  })
) #takes about a minute

ipm_SC <- mcmc.list(nim.out)
ipm_SC <- window(ipm_SC, start = 5000)
stopCluster(cl)

#saveRDS(ipm_SC, file = 'ipm_SC_out.rds')

### Grab Output ####
ipm_SC <- readRDS('./ModelObjects/ipm_SC_out.rds')
View(MCMCvis::MCMCsummary(ipm_SC, round = 2))


## Figure 2 ####
MCMCvis::MCMCplot(ipm_SC, 
                  ci = c(80, 95),
                  params = c('a2', 'b2', 'c2',
                             'm2[1]','m2[2]','m2[3]', 
                             'f2[1]', 'p2', 'd2'), 
                  ref_ovl = TRUE, guide_lines = T, exact = T, 
                  ISB = F, rank = F,
                  sz_labels = 2,
                  sz_ax_txt = 2,
                  labels = c(expression(alpha[2]^omega), #nest attempted 
                             expression(alpha[2]^rho), #nest failure
                             expression(alpha[3]^nu), #brood failure
                             expression(beta[list(1,2)]^"\u2642"), 
                             expression(beta[list(2,2)]^"\u2642"), 
                             expression(beta[list(3,2)]^"\u2642"),
                             expression(beta[list(1,2)]^phi), #pre-breeding survival
                             expression(beta[list(2,2)]^phi), #nesting survival
                             expression(beta[list(2,5)]^phi))) #brooding survival



### Leslie Matrices ####
LMS <- MCMCvis::MCMCsummary(ipm_SC, params = c('LM'), exact = F, Rhat = T, n.eff = T)
LeslieMat_SRS <- matrix(LMS$`50%`[1:9], nrow = 3, byrow = F)
LeslieMat_Webb <- matrix(LMS$`50%`[10:18], nrow = 3, byrow = F)

round(LeslieMat_SRS, digits = 3)
round(LeslieMat_Webb, digits = 3)

Re(eigen(LeslieMat_SRS)$values)[1]
Re(eigen(LeslieMat_Webb)$values)[1]
