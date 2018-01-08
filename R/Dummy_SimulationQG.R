# -----------------------------------------------------------------------------
# Dummy code for simulation of Quantitative genetics parameters
# Audrey Bourret
# 2018-01-05
# -----------------------------------------------------------------------------

#### Packages and datasets ####
if(!require(pedantics)){install.packages("pedantics")}
library(pedantics); data(gryphons)

#### Visualisation of the griffon pedigree ####
gryph.PED <- as.data.frame(cbind(gryphons$id, gryphons$dam, gryphons$sire))
   names(gryph.PED) <- c("animal", "dam", "sire")

gryph.PED <- fixPedigree(Ped=gryph.PED)

drawPedigree(gryph.PED) 


#### Simple simulation of body mass heritability ####

# Parameters
x.MASS     <- 300 # mean body mass
VP.MASS    <- 120 # phenotypic variance in body mass 
h2.MASS    <- 0.3 # body mass heritability

# Simulation
gryph.PHE <- phensim(gryph.PED, trait=1 , randomA=VP.MASS * h2.MASS  , randomE=VP.MASS * (1-h2.MASS), returnAllEffects = TRUE)$allEffects

names(gryph.PHE) <- c("animal", "dam", "sire", "MASS.G", "MASS.E", "MASS.P")
head(gryph.PHE)

# Add mean body mass to simulated phenotypic body mass (zero-centered)
gryph.PHE$MASS <- gryph.PHE$MASS.P + x.MASS 
head(gryph.PHE)

# Check if simulations went well
var(gryph.PHE$MASS.G) / var(gryph.PHE$MASS)

### Add a maternal (non-genetic) and a cohort effect ####

hMAT.MASS    <- 0.3 # maternal effect (proportion explained by)
hCOH.MASS    <- 0.1 # cohort effect (proportion explained by)

# Add infos on cohorts
gryph.PHE$cohort <- gryphons[match(gryphons$id, gryph.PHE$animal), "cohort"]

# remove unecessary column
gryph.PHE$MASS.E <- NULL
gryph.PHE$MASS.P <- NULL
gryph.PHE$MASS   <- NULL

# Simulation of cohort effect
cohort_effects <- as.data.frame(cbind(unique(gryph.PHE$cohort),
                  rnorm(length(unique(gryph.PHE$cohort)),0,sqrt(VP.MASS * hCOH.MASS ))))

gryph.PHE$MASS.C <- cohort_effects[match(gryph.PHE$cohort, cohort_effects[,1]),2]

gryph.PHE$cohort <- as.factor(gryph.PHE$cohort)

# Simulation of maternal effect
gryph.PHE$dam <- as.numeric(as.character(gryph.PHE$dam))

maternal_effects <- as.data.frame(cbind(unique(gryph.PHE$dam),
                  rnorm(length(unique(gryph.PHE$dam)),0,sqrt(VP.MASS * hMAT.MASS ))))

gryph.PHE$MASS.M <- maternal_effects[match(gryph.PHE$dam, maternal_effects[,1]),2]

gryph.PHE$dam <- as.factor(gryph.PHE$dam)

# Simulation of residual effect
gryph.PHE$MASS.R <- rnorm(nrow(gryph.PHE), 0, sqrt(VP.MASS * (1 - hCOH.MASS - hMAT.MASS - h2.MASS)))   

# Add all effect to obtain the body mass of each individual
gryph.PHE$MASS    <- gryph.PHE$MASS.G + gryph.PHE$MASS.C + gryph.PHE$MASS.M + gryph.PHE$MASS.R + x.MASS    
   
head(gryph.PHE)
hist(gryph.PHE$MASS)
var(gryph.PHE$MASS) # 120

# Check if simulation went well
var(gryph.PHE$MASS.G) / var(gryph.PHE$MASS) # 0.3
var(gryph.PHE$MASS.C) / var(gryph.PHE$MASS) # 0.1
var(subset(gryph.PHE, dam!="NA")$MASS.M) / var(subset(gryph.PHE, dam!="NA")$MASS) # 0.3
var(gryph.PHE$MASS.R) / var(gryph.PHE$MASS) # 0.3


#### Simulation of covariance between wing length and reproductive sucess ####

x.WING   <- 250 # mean wing length (cm)
VP.WING  <- 150 # phenotypic variance in wing length 
h2.WING  <- 0.4 # heritability of wing length 

x.RS     <- 5.5 # mean reproductive sucess
VP.RS    <- 2   # phenotypic variance in reproductive sucess 
h2.RS    <- 0.1 # heritability of reproductive sucess

ra       <- 0.1 # genetic correlation betwen wing length and reproductive sucess (should be between -1 and 1) 
re       <- 0.7 # environmental correlation betwen wing length and reproductive sucess (should be between -1 and 1) 

# 2 variance-covariance matrix for the simulation
# genetic matrix (G-matrix):
covGEN <- matrix(NA, 2, 2)
   covGEN[1,1] <- VP.WING * h2.WING                    # genetic variance of wing length
   covGEN[2,2] <- VP.RS * h2.RS                        # genetic variance of reproductive sucess 
   covGEN[1,2] <- ra * sqrt(covGEN[1,1] * covGEN[2,2]) # covariance
   covGEN[2,1] <- covGEN[1,2]                          # covariance (same as above)
covGEN

# environmental matrix
covENV <- matrix(NA, 2, 2)                             
   covENV[1,1] <- VP.WING * (1-h2.WING)                # environmental variance of wing length
   covENV[2,2] <- VP.RS * (1-h2.RS)                    # environmental variance of wing length
   covENV[1,2] <- re * sqrt(covENV[1,1] * covENV[2,2]) # covariance 
   covENV[2,1] <- covENV[1,2]                          # covariance (same as above)
covENV

# Simulation per se
gryph.PHE.bi <- phensim(gryph.PED, trait=2 , randomA=covGEN , randomE=covENV, returnAllEffects = TRUE)$allEffects
   names(gryph.PHE.bi) <- c("animal", "dam", "sire", "WING.G", "RS.G", "WING.E", "RS.E", "WING", "RS")
   gryph.PHE.bi$WING    <- gryph.PHE.bi$WING + x.WING        # add mean
   gryph.PHE.bi$RS    <- gryph.PHE.bi$RS + x.RS


# Add these new values to those already simulated
   
gryph.PHE <- merge(gryph.PHE,gryph.PHE.bi, by=c("animal","sire","dam"))
head(gryph.PHE)
   
var(gryph.PHE$WING.G) / var(gryph.PHE$WING) # 0.4
var(gryph.PHE$RS.G) / var(gryph.PHE$RS) # 0.1 

summary(gryph.PHE)

# Add the sex of each gryphon
gryph.PHE$sex <- "NA"
gryph.PHE$sex <- gryphons[match(gryphons$id, gryph.PHE$animal), "sex"]

gryph.PHE$sex <- as.factor(gryph.PHE$sex)

# Keep only the created dataset

gryph.data <- gryph.PHE[,c("animal","sire","dam","sex","cohort","MASS","WING","RS")]
names(gryph.data) <- c("animal","sire","dam","sex","cohort","mass","wing","rs")

rm(cohort_effects, covENV, covGEN, gryph.PHE, gryph.PHE.bi,gryphons, h2.MASS, h2.RS, h2.WING, hCOH.MASS,hMAT.MASS,maternal_effects,ra,re,VP.MASS, VP.RS,VP.WING, x.MASS, x.WING,x.RS)

summary(gryph.data)

# To do : add repeated values, add different mean depending on the sex of each gryphon.