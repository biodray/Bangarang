# -----------------------------------------------------------------------------
# Dummy code for an explanation on Van de Pol et al. 2009
# Audrey Bourret
# 2018-01-04
# -----------------------------------------------------------------------------

#### 1. Create a dummy dataset ####

# Caracteristics (can be changed):
#  100 individuals, sampled between 2005 and 2010

N.id     <- 100  # Number of individuals (init value = 100)
year.min <- 2005 # First sampled year (init value = 2005)
year.max <- 2010 # Last sampled year (init value = 2010)

year.mast<- c(2005,2008,2009) # Years following a mast (init value = c(2005,2008,2009))

age.min <- 1 # Age at first reproduction (init value = 1)
age.max <- 5 # Max longevity (init value = 5)

id.effect   <- c(8, 0.2) # On mass: mean + sd in a normal distribution (init value = c(8, 0.2))
year.effect <- 0.2       # On mass: linear effect (init value = 0.2)
mast.effect <- 1         # On mass: add 1 g on Mast year (init value = 1)
age.effect  <- 0.4       # On mass: linear effect (init value = 0.4)
error.effect<- c(0, 0.2) # On mass: mean + sd in a normal distribution (init value = c(0, 0.2))

# Run the code below  (lines 27-74) to generate a dataset

data <- expand.grid(id = 1:N.id, year = year.min:year.max)
   data$mast <- 0
   data$mast[data$year%in%year.mast] <- 1
   data$age <- "NA"
   #data$repro <- "NA"
   data$mass <- "NA"

mass.id   <- data.frame(id= 1:N.id, age = round(runif(N.id, min=-(year.max-year.min), max=age.max)), x = rnorm(N.id, mean=id.effect[1], sd=id.effect[2]))
mass.year <- data.frame(year = year.min:year.max, x = c(0:(year.max-year.min))*year.effect)
mass.mast <- data.frame(mast = c(0, 1), x = c(0, mast.effect))
mass.age  <- data.frame(age = (-(year.max-year.min+age.max+1)):(year.max-year.min+age.max+1), x = c(rep(0,year.max-year.min+age.max+age.min+1), (1:age.max)*age.effect,rep(0, year.max-year.min+1)))

for(x in 1:nrow(data)){
   # add age
      age.obs <- mass.id$age[which(mass.id$id == data$id[x])]
      year.obs <- data$year[x]    
      data$age[x] <- year.obs - year.min + age.obs  
   # add a mass value
      id   <- mass.id$x[which(mass.id$id == data$id[x])]
      year <- mass.year$x[which(mass.year$year == data$year[x])]   
      mast <- mass.mast$x[which(mass.mast$mast == data$mast[x])]   
      age  <- mass.age$x[which(mass.age$age == data$age[x])]  
      data$mass[x] <- id + year + mast + age + rnorm(1, mean=error.effect[1], sd=error.effect[2])  
      #data$mass[x] <- round(id + year + mast + age) # for 0 or 1
   }

data$age    <- as.numeric(as.character(data$age))
data$mass   <- as.numeric(as.character(data$mass))
data$id     <- as.factor(data$id)
data$mast.f <- as.factor(data$mast)
data$mast   <- as.numeric(as.character(data$mast))

# Keep only possible values
data <- data[which(data$age%in%c(age.min:age.max) & data$year >= year.min),]

# Compute age.min
first.obs <- aggregate(data$year, list(data$id), min)
data$age <- NULL
for(x in 1:nrow(data)){ 
   data$age.min[x] <- data$year[x] - first.obs$x[which(first.obs$Group.1 == data$id[x])] + 1 
}

# Add year.r for year or sampling 1 to X
data$year.r <- data$year - year.min + 1

# End of the dataset creation!

#### 2. Dummy dataset ####

summary(data)

plot(y=data$mass, x=data$year)
plot(y=data$mass, x=data$age.min)

hist(data$mass)
hist(data$year)
hist(data$age.min)

#### 3. Normal LMM ####

library(lme4)
library(lmerTest) # to get p-values with lmer


m0 <- lm(mass ~ mast.f + year.r + age.min, data=data) 
summary(m0)

m1 <- lmer(mass ~ mast.f + year.r + age.min + (1|id), data=data)
summary(m1)

#### 4. Van de Pol LMM ####

# Step 0 - Keep only individual with 2 observations or more

n.obs <- aggregate(data$id, list(data$id), length) 
   names(n.obs) <- c("id", "n.obs")

data <- merge(data, n.obs, by = "id")   

data <- data[which(data$n.obs >=2),]
data$n.obs <- NULL

# Step 1 - Divide mast and year in mean + dev

# Mean of mast observed (mast.mean) and deviation (mast.dev)

mast.mean <- aggregate(data$mast, list(data$id), mean) 
   names(mast.mean) <- c("id", "mast.mean")

data <- merge(data, mast.mean, by = "id")

data$mast.dev <- data$mast - data$mast.mean

# Mean of year.r observed (year.r.mean) and deviation (year.r.dev)

year.r.mean <- aggregate(data$year.r, list(data$id), mean) 
   names(year.r.mean) <- c("id", "year.r.mean")

data <- merge(data, year.r.mean, by = "id")

data$year.r.dev <- data$year.r - data$year.r.mean

head(data)
summary(data)

# Step 2 - Model 

m2 <- lmer(mass ~ mast.mean + mast.dev + year.r.mean + year.r.dev + age.min + (1|id), data=data) 
summary(m2)

# Model with random slopes 

m2a <- lmer(mass ~ mast.mean + mast.dev + year.r.mean + year.r.dev + age.min + (mast.dev|id), data=data) 
summary(m2a)

m2b <- lmer(mass ~ mast.mean + mast.dev + year.r.mean + year.r.dev + age.min + (year.r.dev|id), data=data) 
summary(m2b)

anova(m2, m2a) #if significant you can keep the random slope
anova(m2, m2b)

# Step 3 - Graphic

# Compute individual slopes

pred.mast.b <- data.frame(id = names(tapply(data$mast.dev, list(data$id), min)),
                        mast.min = tapply(data$mast.dev, list(data$id), min),
                        mast.max = tapply(data$mast.dev, list(data$id), max))
pred.mast.b <- pred.mast.b[which(!is.na(pred.mast.b$mast.min)),]

pred.mast.a <- data.frame(coef(m2)$id) # or m2a/m2b depending if the random slope model is the best one

pred.mast <- data.frame(pred.mast.a, pred.mast.b) 

pred.mast$mass.mast.min <- pred.mast$X.Intercept. + (pred.mast$mast.min * pred.mast$mast.dev)
pred.mast$mass.mast.max <- pred.mast$X.Intercept. + (pred.mast$mast.max * pred.mast$mast.dev)

# Compute population slope

pop.INT <- fixef(m2)["(Intercept)"]
pop.SLP <- fixef(m2)["mast.mean"]
  
# graph (adjust manually xlim and ylim)

plot(x=0, y=0, xlab = "delta(Mast)", ylab="Mass (BLUP)", xlim = c(-0.7,0.7), ylim = c(8.5,11.5))
 # Add individual slopes
   for(x in unique(pred.mast$id)){
      lines(x=c(pred.mast[which(pred.mast$id == x), "mast.min"], pred.mast[which(pred.mast$id == x), "mast.max"]), y=c(pred.mast[which(pred.mast$id == x), "mass.mast.min"], pred.mast[which(pred.mast$id == x), "mass.mast.max"]),
         lwd=1, col="gray50")
   }
 # Add population slope
   abline(a=pop.INT, b=pop.SLP, lty="dashed", lwd=3, col="black")
