rm(list = ls())

library('dplyr')
library('MASS')
library('plm')
library('xtable')

library('foreach')
library('doParallel')

#setwd('C:/Users/mayur/Dropbox/Industrial Organization/EC 536 Part 4 Production Function and Vertical Contracting/EC 536 PS7 Mayur')
setwd('C:/Users/mchoudhary/Dropbox/Industrial Organization/EC 536 Part 4 Production Function and Vertical Contracting/EC 536 PS7 Mayur')

## source all the required functions 
source('firstStage.R')
source('secondStage.R')
source('secondStageTab6.R')
source('locReg.R')
source('model.R')
source('modelTab8.R')
source('modelTab6.R')
source('gmmFun.R')
source('gmmFunPara.R')

## change command names for MASS and dplyr 

dselect <- dplyr::select
dfilter <- dplyr::filter

### read the data from finalSample.csv ### 

data.full <- read.csv('finalSample.csv', header = TRUE)

### Assign the residual quality measures (controlling for patient characteristics) 
q.reg <- lm(infect_rate ~ vat_cat1_ + avg_age_ + avg_comb_cond_ + avg_dur_esrd_ + avg_hemog_pct_ + female_all_pct_, data = data.full)
data.full$q.res <- -(q.reg$residuals)

## do the same for death rate ratio to be used as iv for measurement error 
q.iv.reg <- lm(death_ratio ~ vat_cat1_ + avg_age_ + avg_comb_cond_ + avg_dur_esrd_ + avg_hemog_pct_ + female_all_pct_, data = data.full)
data.full$q.iv.res <- -(q.iv.reg$residuals)

## Add the quality and factors product terms 
data.full <- data.full %>% mutate(q.k = q.res*lstations , q.l = q.res*lstaff, q.k.iv = q.iv.res*lstations , q.l.iv = q.iv.res*lstaff)


### subset data to get include 0 capital investment and non zero human capital investment 

data.full$hirepos <- ifelse(data.full$lhires > 0 , 1, 0) ## positive hires
data.subset <- data.full %>% filter(hires !=0 & invest == 0) ## subset data 


### We need a function that computes the first stage and the second stage of the original model Table 5
# Let that function be model() 

result.tab5 <- model(data.subset)
write.csv(result.tab5,'result_tab5.csv')

##### Bootstrapping to compute standard errors ##### 
detectCores()
registerDoParallel(cores=8)

n.boot <- 100

firm.unique <- unique(data.subset$cms_code)
firm.num <- length(firm.unique)

firm.id <- seq(1:firm.num)

## Bootstrap computations begin ## 
start.time <- Sys.time()
registerDoParallel(cores=8)

tab5.boot <- foreach(b.iter = 1:n.boot, .packages = c('MASS','dplyr', 'plm')) %dopar% {
  sample.firm <- sample(firm.unique, firm.num, replace = TRUE) ### Sample with replacement from the firms 
  sample.firm <- as.data.frame(cbind(sample.firm, firm.id))
  names(sample.firm) <- c('cms_code','firm.id')
  
  data.temp <- inner_join(sample.firm, data.subset, by = 'cms_code') ## get the sample using the draws 
  data.boot <- data.temp[order(data.temp$firm.id, data.temp$row_id),] ## set the data up in the form of the original file 
  model(data.boot)

  }

end.time <- Sys.time()

end.time - start.time  

tab5.boot.mat <- matrix(0, nrow(result.tab5), n.boot) 

for (i in 1:n.boot){
  tab5.boot.mat[,i] <- tab5.boot[[i]]
}

write.csv(tab5.boot.mat, 'tab5_boot.csv') 

### compute standard error ###
tab5.boot.temp <- (tab5.boot.mat - matrix(result.tab5, nrow(tab5.boot.mat),ncol(tab5.boot.mat)))^2
tab5.boot.se <- as.matrix(sqrt(rowMeans(tab5.boot.temp)))
write.csv(tab5.boot.se, 'tab5_boot_se.csv')


###################################################
#### Perform the same computation for Table 6 #####
###################################################

### First get the estimate ###
result.tab6 <- modelTab6(data.subset) 
write.csv(result.tab6, 'result_tab6.csv')

##### Bootstrapping to compute standard errors ##### 
registerDoParallel(cores=8)


## Bootstrap computations begin ## 
start.time <- Sys.time()
registerDoParallel(cores=8)

tab6.boot <- foreach(b.iter = 1:n.boot, .packages = c('MASS','dplyr', 'plm')) %dopar% {
  sample.firm <- sample(firm.unique, firm.num, replace = TRUE) ### Sample with replacement from the firms 
  sample.firm <- as.data.frame(cbind(sample.firm, firm.id))
  names(sample.firm) <- c('cms_code','firm.id')
  
  data.temp <- inner_join(sample.firm, data.subset, by = 'cms_code') ## get the sample using the draws 
  data.boot <- data.temp[order(data.temp$firm.id, data.temp$row_id),] ## set the data up in the form of the original file 
  modelTab6(data.boot)
  
}

end.time <-Sys.time()

end.time - start.time  

tab6.boot.mat <- matrix(0, nrow(result.tab6), n.boot) 

for (i in 1:n.boot){
  tab6.boot.mat[,i] <- tab6.boot[[i]]
}

write.csv(tab6.boot.mat, 'tab6_boot.csv') 

### compute standard error ###
tab6.boot.temp <- (tab6.boot.mat - matrix(result.tab6, nrow(tab6.boot.mat),ncol(tab6.boot.mat)))^2
tab6.boot.se <- as.matrix(sqrt(rowMeans(tab6.boot.temp)))
write.csv(tab6.boot.se, 'tab6_boot_se.csv')


###################################################
#### Perform the same computation for Table 8 #####
###################################################

### First get the estimate ###
result.tab8 <- modelTab8(data.subset) 
write.csv(result.tab8, 'result_tab8.csv')

##### Bootstrapping to compute standard errors ##### 
registerDoParallel(cores=8)

## Bootstrap computations begin ## 
start.time <- Sys.time()
registerDoParallel(cores=8)

tab8.boot <- foreach(b.iter = 1:n.boot, .packages = c('MASS','dplyr', 'plm')) %dopar% {
  sample.firm <- sample(firm.unique, firm.num, replace = TRUE) ### Sample with replacement from the firms 
  sample.firm <- as.data.frame(cbind(sample.firm, firm.id))
  names(sample.firm) <- c('cms_code','firm.id')
  
  data.temp <- inner_join(sample.firm, data.subset, by = 'cms_code') ## get the sample using the draws 
  data.boot <- data.temp[order(data.temp$firm.id, data.temp$row_id),] ## set the data up in the form of the original file 
  modelTab8(data.boot)
  
}

end.time <- Sys.time()

tab8.boot.time <- end.time - start.time  

tab8.boot.mat <- matrix(0, nrow(result.tab8), n.boot) 

for (i in 1:n.boot){
  tab8.boot.mat[,i] <- tab8.boot[[i]]
}

write.csv(tab8.boot.mat, 'tab8_boot.csv') 

### compute standard error ###
tab8.boot.temp <- (tab8.boot.mat - matrix(result.tab8, nrow(tab8.boot.mat),ncol(tab8.boot.mat)))^2
tab8.boot.se <- as.matrix(sqrt(rowMeans(tab8.boot.temp)))
write.csv(tab8.boot.se, 'tab8_boot_se.csv')


