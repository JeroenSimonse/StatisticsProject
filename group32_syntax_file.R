### Title:          Statsistics & Methodology group project


###          ###
### Overview ###
###          ###


###--Data processing--###
##---Preliminaries---------------------------------------------------------------------------------------##

# Loading the packages
library(mice)
library(MASS)
library(dplyr) 
library(mitools)
library(wec)
library(miceadds)
library(MLmetrics)
library(car)

# Setting working directory
dataDir <- "../data/"
dat1 <- readRDS(paste0(dataDir, "wvs_data.rds"))

##---Data-cleaning--------------------------------------------------------------------------------------##

# All conservatism sub categories
tradition   <- c("V79","V100") 
hierarchy    <- c("V45","V47","V51","V53","V54")
authority   <- c("V69","V127")
nationalism <- c("V46","V107","V214")
religion    <- c("V145","V146","V152") 
pswb        <- c("V10","V23","V55","V216") 
pos_var_Satisf  <- c("V5","V6","V11","V59","V150","V153","V192","V239")
confounding_vars <- c("V240", "V2", "V57","V238","V242")


# Sub setting the data that will be used
realDat<-dat1[c(tradition,hierarchy,authority,nationalism,religion,pswb,pos_var_Satisf,confounding_vars)]


# All negative values to NA so that we can treat them later on
realDat[realDat<0]<- NA

# Recoding all variables so that high values and low values mean the same thing
realDat$V79<- recode(realDat$V79,"1 = 6; 2 = 5; 3 = 4; 4 = 3; 5 = 2; 6 = 1")
realDat$V100<- recode(realDat$V100, "1 = 10; 2 = 9; 3 = 8; 4 = 7; 5 = 6; 6 = 5; 7 = 4; 8 = 3; 9 = 2; 10 = 1")

realDat$V45<- recode(realDat$V45,"1 = 3; 2 = 2; 3 = 1")
realDat$V47<- recode(realDat$V47,"1 = 3; 2 = 2; 3 = 1")
realDat$V51<- recode(realDat$V51,"1 = 4; 2 = 3; 3 = 2; 4 = 1")
realDat$V53<- recode(realDat$V53,"1 = 4; 2 = 3; 3 = 2; 4 = 1")
realDat$V54<- recode(realDat$V54,"1 = 4; 2 = 3; 3 = 2; 4 = 1")

realDat$V69<- recode(realDat$V69,"1 = 3; 2 = 2; 3 = 1")
realDat$V127<- recode(realDat$V127,"1 = 4; 2 = 3; 3 = 2; 4 = 1")
realDat$V214<- recode(realDat$V214,"1 = 4; 2 = 3; 3 = 2; 4 = 1")

realDat$V107<- recode(realDat$V107,"1 = 4; 2 = 3; 3 = 2; 4 = 1")

realDat$V145<- recode(realDat$V145,"1 = 7; 2 = 6; 3 = 5; 4 = 4; 5 = 3; 6 = 2; 7 = 1")
realDat$V146<- recode(realDat$V146,"1 = 8; 2 = 7; 3 = 6; 4 = 5; 5 = 4; 6 = 3; 7 = 2; 8 = 1")

realDat$V10<- recode(realDat$V10,"1 = 4; 2 = 3; 3 = 2; 4 = 1")
realDat$V216<- recode(realDat$V216,"1 = 4; 2 = 3; 3 = 2; 4 = 1")


# realDat$tradition = rowMeans(subset(realDat, select = c(V79,V100)), na.rm = TRUE)
# realDat$hierarchy = rowMeans(subset(realDat, select = c(V45,V47,V51,V53,V54)), na.rm = TRUE)
# realDat$authority = rowMeans(subset(realDat, select = c(V69,V127)), na.rm = TRUE)
# realDat$nationalism = rowMeans(subset(realDat, select = c(V46,V107,V214)), na.rm = TRUE)
# realDat$religion = rowMeans(subset(realDat, select = c(V145,V146,V152)), na.rm = TRUE)
rowMeans(subset(realDat, select = c(V10,V23,V55,V216)), na.rm = TRUE)
# 
# realDat$tradition[realDat$tradition=="NaN"] <- NA
# realDat$hierarchy[realDat$hierarchy=="NaN"] <- NA
# realDat$authority[realDat$authority=="NaN"] <- NA
# realDat$nationalism[realDat$nationalism=="NaN"] <- NA
# realDat$religion[realDat$religion=="NaN"] <- NA
# realDat$pswb[realDat$pswb=="NaN"] <- NA
# 
# realDat[realDat$tradition>5 | realDat$tradition<7 ]<- "Average"

# Giving labels to each nominal/binary variable
realDat$V2 <- factor(realDat$V2,     levels = c(156, 276, 356, 643, 840), labels = c("China", "Germany", "India", "Russia", "US"))
realDat$V150 <- factor(realDat$V150, levels = c(1, 2, 3, 4), labels = c("follow religious norms and ceremonies", "do good to other people", " either of them,", " both"))
realDat$V240 <- factor(realDat$V240, levels = c(1, 2), labels= c("male", "female"))
realDat$V57 <- factor(realDat$V57,   levels = c(1, 2, 3, 4, 5, 6), labels = c("married", "living together as married", "divorced", "separated", "widowed", "single"))
realDat$V238 <- factor(realDat$V238, levels = c(1, 2, 3, 4, 5), labels = c("upper class", "upper middle class", "lower middle class", "working class", "lower class"))


# Proportion of missing data and covariance coverage rates for relevant variables
cm <- colSums(is.na(realDat)) #cases missing per variable
pm <- colMeans(is.na(realDat)) #percentage missing per variable
co <- colSums(!is.na(realDat)) #cases observed
cc <- md.pairs(realDat)$rr/nrow(realDat) # covariance coverage rate per variable
eps <- 0.80
all(cc > eps) # no variable has more then 20% missing cases

# Method for each ordinal variable = "pmm", method for each nominal variable = "polr" and method for each binary variable = "logreg"
meth        <- rep("pmm", ncol(realDat))
names(meth) <- colnames(realDat)
meth["V150"]<- "polr"
meth["V2"]  <- "polr"
meth["V240"]<- "logreg"
meth["V57"] <- "polr"
meth["V238"]<- "polr"

# Prediction matrix of variables
predMat <- quickpred(realDat, mincor = 0.2, include = "V240")


# Multiple imputation
miceOut<- mice(data        = realDat,
               m               = 10,
               maxit           = 20,
               method          = meth,
               predictorMatrix = predMat,
               seed            = 235711)

# We remove the nominal variables as we can't use them for our multivariate outlier detection.
miceOut2 <-
  subset_datlist(datlist = miceOut,
                 select  = setdiff(colnames(realDat), c("V150","V2","V240","V57","V238")),
                 toclass = "mids")
completedData <- complete(miceOut2, "all")

mdOutliers <-
  function(data, critProb, statType = "mcd", ratio = 0.75, seed = NULL)
  {
    ## Set a seed, if one is provided:
    if(!is.null(seed)) set.seed(seed)
    
    ## Compute (robust) estimates of the mean and covariance matrix:
    stats <- cov.rob(x             = data,
                     quantile.used = floor(ratio * nrow(data)),
                     method        = statType)
    
    ## Compute robust squared Mahalanobis distances
    md <- mahalanobis(x = data, center = stats$center, cov = stats$cov)
    
    ## Find the cutoff value:
    crit <- qchisq(critProb, df = ncol(data))
    
    ## Return row indices of flagged observations:
    which(md > crit)
  }

# For each imputed dataset we check for the outliers using the mahalnobis method
outliersPerSet <- lapply(completedData, mdOutliers, critProb = 0.99)

# How many times did each outlier get flagged as an outlier
occurancesAsOutlier <- table(unlist(outliersPerSet))

# At how many occurances are we deleting them? -> iteratons / 2
thresh <- 20/2

# Making a variable that consists of all row indices that are outliers
outs <- as.numeric(names(occurancesAsOutlier[occurancesAsOutlier >= thresh]))

# Deleting the multivariate outliers from our dataset
miceOut3 <- subset_datlist(datlist = miceOut, 
                           subset  = setdiff(1 : nrow(realDat), outs),
                           toclass = "mids")
completedData3 <- complete(miceOut3,"all")

miceOut4<-subset_datlist(datlist = miceOut3)
for (m in 1:length(miceOut4)){miceOut4[[m]][["tradition"]]<-rowMeans(subset(miceOut4[[m]], select = c(V79,V100)))
miceOut4[[m]][["hierarchy"]]<-rowMeans(subset(miceOut4[[m]], select = c(V45,V47,V51,V53,V54)))
miceOut4[[m]][["authority"]]<-rowMeans(subset(miceOut4[[m]], select = c(V69,V127)))
miceOut4[[m]][["nationalism"]]<-rowMeans(subset(miceOut4[[m]], select = c(V46,V107,V214)))
miceOut4[[m]][["religion"]]<-rowMeans(subset(miceOut4[[m]], select = c(V145,V146,V152)))
miceOut4[[m]][["pswb"]]<-rowMeans(subset(miceOut4[[m]], select = c(V10,V23,V55,V216)))
miceOut4[[m]][["conservatism"]]<-rowMeans(subset(miceOut4[[m]], select = c(hierarchy,authority,nationalism,religion)))
}

miceOut5<-subset_datlist(datlist = miceOut4,toclass = "mids")

###--Inference Analysis--###
##--------------------------------------------------------------------------------------------------------------##
# Question: Are conservative attitudes good or bad for your psychological well-being?

# Dependent (Satisfaction with Life)

#V10  ***   - Feeling of happiness                                            - 1: Not happy 4: Very Happy                 
#V23  **    - Satisfaction with your life                                     - 1: Dissatisfied 10: Satisfied
#V55  *     - Freedom and control of life                                     - 1: Not much 10: Much
#V216 *     - I see myself as an autonomous individual                        - 1: Strongly agree 4: Strongly disagree

# Independent/Predictors (Conservative attitudes)

# TRADITION
#V79  **    - Schwartz: Tradition is important to this person;                - 1: Not conservative     6 : Conservative
#V100 **    - Hard work brings success                                        - 1: Not conservative     10: Conservative

# HIERARCHY - 
#V45        - When jobs are scarce, men should have more right to a job than  - 1: Not conservative     3 : Conservative
#             a woman
#V47  *     - If a woman earns more money than her hushand, its almost        - 1: Not conservative     3 : Conservative
#             certain to cause problems
#V51        - On the whole, men make better political leaders than women do   - 1. Not conservative     4 : Conservative
#V54  *     - Being a housewife is just as fulfilling as working for pay      - 1: Not conservative     6 : Conservative

# AUTHORITY  -
#V69  **    - Future changes: Greater respect for authority                   - 1: Not conservative     3 : Conservative 
#V127 **    - Political system: Having a strong leader who does not have  to  - 1: Not conservative     4 : Conservative
#             bother with parliament and elections

# NATIONALISM
#V46        - When jobs are scarce, employers should give priority to people  - 1: Not conservative     3 : Conservative
#             of this country over immigrants           
#V107       - How much you trust: People of another nationality               - 1: Not conservative     4 : Conservative
#V214       - I see myself as part of the [country] nation                    - 1: Not conservative     4 : Conservative

# RELIGION
#V145       - How often do you attend religious services                      - 1: Not conservative     7 : Conservative
#V146       - How often do you play                                           - 1: Not conservative     8 : Conservative
#V152       - How important is God in your life                               - 1: Not conservative     10: Conservative

#EDA / Convergence checks to control the imputations by plotting the observed vs. imputed densities. 

densityplot(miceOut5, ~V10)
densityplot(miceOut5, ~V10|.imp)
densityplot(miceOut5, ~V23)
densityplot(miceOut5, ~V23|.imp)
densityplot(miceOut5, ~V69)
densityplot(miceOut5, ~V69|.imp)
densityplot(miceOut5, ~V79)
densityplot(miceOut5, ~V79|.imp)
densityplot(miceOut5, ~V100)
densityplot(miceOut5, ~V100|.imp)
densityplot(miceOut5, ~tradition)
densityplot(miceOut5, ~tradition|.imp)

# Now we've seen that our imputed data corresponds to our observed data from our original data set we can take a look at our DV's 
# We will do this by using a correlation matrix to see how strong our main DV's correlate with our chosen IV's.
# First we create a correlation vector, this will make our correlation matrix more flexible if our IV's and/or DV's will change in the future. 
CorVector1 <- c('pswb','tradition','hierarchy','authority','nationalism','religion')

# Now we've created our Correlation Vectors, we add these correlation vectors to our correlation matrixes, and create the correlation matrixes.

corMatrix1 <- micombine.cor(
  miceOut5,
  variables = CorVector1,
  method = c('pearson'))
corMatrix1.1<-  corMatrix1[order(corMatrix1$variable1,
                                 decreasing = FALSE),]
head(corMatrix1[c(30,29,27,24,20),],n=length(CorVector1)-1)


# Fitting our IV's one by one against our DV

fit1 <- lm.mids(pswb ~ tradition, data = miceOut5)
summary(pool(fit1))
pool.r.squared(fit1)[1] #0.02075772     #2

fit2 <- lm.mids(pswb ~ religion, data = miceOut5)
summary(pool(fit2))
pool.r.squared(fit2)[1] #0.007649452    #4

fit3 <- lm.mids(pswb ~ nationalism, data = miceOut5)
summary(pool(fit3))
pool.r.squared(fit3)[1] #0.02338606     #1

fit4 <- lm.mids(pswb ~ hierarchy, data = miceOut5)
summary(pool(fit4))
pool.r.squared(fit4)[1] #0.01555805     #3

fit5 <- lm.mids(pswb ~ authority, data = miceOut5)
summary(pool(fit5))
pool.r.squared(fit5)[1] #0.006055234    #5

fit6 <- lm.mids(pswb ~ conservatism, data = miceOut5)
summary(pool(fit6))
pool.r.squared(fit6)[1] #0.002142242    #6


# Two model to check if hierarchy add value to our model
fit1.4 <- lm.mids(pswb ~ nationalism + tradition + religion + hierarchy, data = miceOut5)
summary(pool(fit1.4))
pool.r.squared(fit1.4)[1] #0.06024238

fit1.5 <- lm.mids(pswb ~ nationalism + tradition + religion + hierarchy + authority, data = miceOut5)
summary(pool(fit1.5))
pool.r.squared(fit1.5)[1] #0.06340318

# Above we can see that every added subcategory explains more of our DV (Psychological well-being)
fTest <- pool.compare(fit1.5, fit1.4, method = ("wald"))
fTest$Dm      #36.59858 #F-value / Test statistic for table
fTest$pvalue  #2.065728e-09        #P-value  

# Fitting our Confounders one by one against our DV. 

confit1 <- lm.mids(pswb ~ V240, data = miceOut5)    #Controlling for gender
summary(pool(confit1))
pool.r.squared(confit1)[1]  #0.003235927

confit2 <- lm.mids(pswb ~ V2, data = miceOut5)      #Controlling for country
summary(pool(confit2))
pool.r.squared(confit2)[1]  #0.127336

confit3 <- lm.mids(pswb ~ V57, data = miceOut5)     #Controlling for maritial status
summary(pool(confit3))
pool.r.squared(confit3)[1]  #0.02557787

confit4 <- lm.mids(pswb ~ V238, data = miceOut5)    #Controlling for social class
summary(pool(confit4))
pool.r.squared(confit4)[1]  #0.06229852

confit5 <- lm.mids(pswb ~ V242, data = miceOut5)    #Controlling for age
summary(pool(confit5))
pool.r.squared(confit5)[1]  #0.003803789


# Controlling for all confounders
confitTot <- lm.mids(pswb ~ nationalism + tradition + hierarchy + religion + authority + V242 + V238 + V57 + V2 + V240, data = miceOut5)
summary(pool(confitTot))
pool.r.squared(confitTot)[1] #0.2125141 

confitTot2 <- lm.mids(pswb ~ conservatism + V242 + V238 + V57 + V2 + V240, data = miceOut5)
summary(pool(confitTot2))
pool.r.squared(confitTot2)[1] #0.1890759


# Identifying the difference with and without confounders
fTest2 <- pool.compare(confitTot, fit1.5, method = ("wald"))
fTest2$Dm          # 77.40966
fTest2$pvalue      # 0

# Hypothese H10: There is no effect of conservative attitudes on your psychological well-being
# Denied want:
summary(pool(fit1.5))

# Hypothese H1A: There is an effect of conservative attitudes on your psychological well-being
# RETAINED BECAUSE: 
summary(pool(fit1.5 ))

# Hypothese H20: There is no positive effect of our subgroup nationalism on your psychological well-being
# REJECTED BECAUSE:
summary(pool(fit3))

# Hypothese H2A: There is a positive effect of our subgroup nationalism on your psychological well-being
# RETAINED BECAUSE:
summary(pool(fit3))

# Hypothese H30: There is no positive effect of our subgroup tradition on your psychological well-being
#REJECTED BECAUSE: 
summary(pool(fit1))

# Hypothese H3A: There is a positive effect of our subgroup tradition on your psychological well-being
# RETAINED BECAUSE
summary(pool(fit1))

# Hypothese H40: There is no positive effect of our subgroup religion on your psychological well-being
# REJECTED BECAUSE:
summary(pool(fit2))

# Hypothese H4A: There is a positive effect of our subgroup religion on your psychological well-being
# RETAINED BECAUSE:
summary(pool(fit2))

# Hypothese H50: There is no positive effect of our subgroup authority on your psychological well-being
# RETAINED BECAUSE:
summary(pool(fit5))

# Hypothese H5A: There is a positive effect of our subgroup authority on your psychological well-being

# Hypothese H60: There is no negative effect of our subgroup hierarchy on your psychological well-being
# RETAINED BECAUSE: 
summary(pool(fit4))

# Hypothese H6A: There is a negative effect of our subgroup hierarchy on your psychological well-being

# Hypotheses H70: There is no effect of conservative attitudes on your psychological well-being, after controlling for your confounders:
# Denied want: 
summary(pool(confitTot))

# Hypothese H7A: There is an effect of conservative attitudes on your psychological well-being, after controlling for your confounders
# RETAINED BECAUSE: 
summary(pool(confitTot))


###--Predictive Analysis--###
##--------------------------------------------------------------------------------------------------------------##
# Task: Use multiple linear regression to build a model that predicts the outcome "Satisfaction with life"

source("studentFunctions.R")
source("miPredictionRoutines.R")

# Split the data into train, validation and test sets
n <- nrow(completedData3[[1]])
index <- sample(
  c(rep("train", 8000), rep("valid", 2000), rep("test", n - 10000))
)

impList <- splitImps(imps = completedData3, index = index)

# Define the models to test
mods <- c("V23 ~ V10",
          "V23 ~ V10 + V11",
          "V23 ~ V10 + V11 + V59",
          "V23 ~ V10 + V11 + V59 + V150 + V239",
          "V23 ~ V59 + V10 + V11 + V55 + V45 + V150 + V153 + V152 + V192 + V239  + V2 + V57")

# Conduct 20-fold cross-validation in each multiply imputed dataset:
tmp <- sapply(impList$train, cv.lm, K = 20, models = mods, seed = 235711)

# Aggregate the MI-based CVEs:
cve <- rowMeans(tmp)
cve

# Refit the winning model and compute test-set MSEs:
index2   <- gsub(pattern = "valid", replacement = "train", x = index)
impList2 <- splitImps(completedData3, index2)

fits <- lapply(X   = impList2$train,
               FUN = function(x, mod) lm(mod, data = x),
               mod = mods[which.min(cve)])
mse <- mseMi(fits = fits, newData = impList2$test)

mse
