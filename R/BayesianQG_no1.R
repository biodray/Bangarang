# -----------------------------------------------------------------------------
# Performing Quantitative genetics analysis with Baysian models
# No 1: Compute Gryffon body mass heritability 
# Audrey Bourret
# 2018-01-08
# -----------------------------------------------------------------------------

#### Packages and datasets ####
if(!require(MCMCglmm)){install.packages("MCMCglmm")}
library(MCMCglmm)


# Load dummy dataset
source("R/Dummy_SimulationQG.R", .GlobalEnv)

summary(gryph.data)
summary(gryph.PED)

#### Model 1.1 - heritability of body mass #1 ####

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002))

model1.1 <- MCMCglmm(mass ~ 1, random = ~animal, pedigree = gryph.PED, data = gryph.data, prior = prior1.1, verbose = TRUE) 
# Default parameters: burnin = 3000, nitt = 10000, thin = 10 (sampled every 10 iterations)

# Check the model
# Fixed effects 

plot(model1.1$Sol)

# Random effects
# residual = units
plot(model1.1$VCV)

# Test for autocorrelation: should be around 0, or < 0.1
autocorr(model1.1$Sol)
autocorr(model1.1$VCV)

# Effective size of the model: should be at least > 1000, ideally 10,000
effectiveSize(model1.1$Sol)
effectiveSize(model1.1$VCV)

# Convergence test: should be P > 0.05)
heidel.diag(model1.1$Sol)
heidel.diag(model1.1$VCV)

# Then check what is in this first model
summary(model1.1)

# Conclusion : small pattern of autocorrelation in random effects, to be corrected by increasing nitt / thin.

#### Model 1.2 - heritability of body mass #2 ####

model1.2 <- MCMCglmm(mass ~ 1, random = ~animal, pedigree = gryph.PED, data = gryph.data, prior = prior1.1, verbose = TRUE, burnin = 20000, nitt = 150000, thin = 150)

plot(model1.2$VCV)
autocorr(model1.2$VCV)

effectiveSize(model1.2$Sol)
effectiveSize(model1.2$VCV)

# Conclusion: this model seems better, so we can get the estimates of variance componants with 95%CI

posterior.mode(model1.2$VCV)
HPDinterval(model1.2$VCV)

# heritability

posterior.heritability1.1 <- model1.2$VCV[, "animal"] / (model1.2$VCV[,"animal"] + model1.2$VCV[,"units"]) 

posterior.mode(posterior.heritability1.1)
mean(posterior.heritability1.1)

HPDinterval(posterior.heritability1.1, 0.95)

plot(posterior.heritability1.1)

effectiveSize(posterior.heritability1.1)

#### Model 1.2 - Check if priors were adequate ####

# Does your first priors had an influence on results?
# Check by creating priors where genetics explained a lot

prior1.2.2 <- list(G = list(G1 = list(V = matrix(var(gryph.data$mass) * .95), nu = 1)), R = list(V = matrix(var(gryph.data$mass) * .05), nu = 1))

model1.2.2 <- MCMCglmm(mass ~ 1, random = ~animal, pedigree = gryph.PED, data = gryph.data, prior = prior1.2.2, verbose = TRUE, burnin = 20000, nitt = 100000, thin = 100)

# Then compare both models
posterior.mode(model1.2$VCV)   # non informatif priors
posterior.mode(model1.2.2$VCV) # informatif priors

# Conclusion: Similaur value so no impacts of your priors on your estimates

#### Model 1.3 - Add a fixed effect ####

model1.3 <- MCMCglmm(mass ~ sex, random = ~animal, pedigree = gryph.PED, data = gryph.data, prior = prior1.1, verbose = TRUE, burnin = 20000, nitt = 100000, thin = 100)

posterior.mode(model1.3$Sol)  
HPDinterval(model1.3$Sol, 0.95) 
plot(model1.3$Sol)  

#### Model 1.4 - Add other random effects ####

# Then we need to make new priors to include all random effects du nouvel effet aléatoire qu'on veut tester. Voici un exemple
# 1.4a : cohort + dam
prior1.4a <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, n = 0.002)), R = list(V = 1, nu = 0.002))

model1.4a <- MCMCglmm(mass ~ 1, random = ~ cohort + dam, pedigree = gryph.PED, data = gryph.data, prior = prior1.4a, verbose = TRUE, burnin = 20000, nitt = 100000, thin = 100)

# 1.4b: cohort + dam + animal
prior1.4b <- list(G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, n = 0.002), G3 = list(V = 1, n = 0.002)), R = list(V = 1, nu = 0.002))

model1.4b <- MCMCglmm(mass ~ 1, random = ~ cohort + dam + animal, pedigree = gryph.PED, data = gryph.data, prior = prior1.4b, verbose = TRUE, burnin = 20000, nitt = 100000, thin = 100)

# Now we can compare both model with DIC (lower = better; but take care: we do not really understand how DIC react)

model1.4a$DIC
model1.4b$DIC

posterior.mode(model1.4b$VCV)  
HPDinterval(model1.4b$VCV, 0.95) 
plot(model1.4b$VCV)  

posterior.heritability1.4b <- model1.4b$VCV[, "animal"] / (model1.4b$VCV[,"cohort"] + model1.4b$VCV[,"dam"] + model1.4b$VCV[,"animal"] + model1.4b$VCV[,"units"]) 

posterior.mode(posterior.heritability1.4b)
mean(posterior.heritability1.4b)

HPDinterval(posterior.heritability1.4b, 0.95)

# Conclusion: heritability of body mass = 0.38 [0.31 : 0.44]
