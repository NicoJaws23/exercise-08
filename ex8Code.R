#Step 1, load and look at data
library(tidyverse)
library(skimr)
f <- "https://raw.githubusercontent.com/difiore/ada-datasets/main/Street_et_al_2017.csv"
d <- read_csv(f)
skim(d)

#Step 2: Plots ECV as a function of social group size,
#longevity, weaning, and reproductive lifespan
library(ggplot2)
library(cowplot)
ECV_GS <- ggplot(data = d, mapping = aes(Group_size, ECV)) +
  geom_point()
ECV_Long <- ggplot(data = d, mapping = aes(Longevity, ECV)) +
  geom_point()
ECV_Wean <- ggplot(data = d, mapping = aes(Weaning, ECV)) +
  geom_point()
ECV_Repro <- ggplot(data = d, mapping = aes(Repro_lifespan, ECV)) +
  geom_point()
plot_grid(ECV_GS, ECV_Long, ECV_Wean, ECV_Repro)
#Step 3: By hand get the beta1 and beta0
#for ECV (Y variable) as a function of Group_size (X variable)
#!!! Remove rows with missing data
d <- d |> 
  filter(!is.na(ECV) & !is.na(Group_size))

ECV_sd <- sd(d$ECV) #Y
GS_sd <- sd(d$Group_size) #X

b1 <- cor(d$Group_size, d$ECV) * (ECV_sd/GS_sd)
b0 <- mean(d$ECV) - (b1*mean(d$Group_size))
print(paste("Beta 1 =",b1, "; Beta 0 =",b0))

#Step 4: Use lm() to check results
ECVlm <- lm(ECV ~ Group_size, data = d) 
print(ECVlm) 

#Step 5, Repeat this process for each primate radiation
#catarrhines, platyrrhines, strepsirhines

#Calc b1 and b0 on filtered dataset by radiation
regCof <- function(df, rad, x, y){
  group <- df |>
    filter(Taxonomic_group == rad)
  sdX <- sd(group[[x]])
  sdY <- sd(group[[y]])
  b1 <- cor(group[[x]], group[[y]], use = "complete.obs") * (sdY/sdX)
  b0 <- mean(group[[y]]) - (b1*mean(group[[x]]))
  return(c(b1, b0))
}

#Use lm() on a filtered data set
filtLM <- function(df, rad, x, y){
  group <- df |>
    filter(Taxonomic_group == rad)
  ECVlm <- lm(group[[y]] ~ group[[x]])
  return(ECVlm)
}

cattys <- regCof(d, "Catarrhini", "Group_size", "ECV")
cattysLM <- filtLM(d, "Catarrhini", "Group_size", "ECV")
print(paste("Catarrhini Beta 1:", cattys[1], "and Beta 0:", cattys[2], "by hand"))
print(cattysLM)

platys <- regCof(d, "Platyrrhini", "Group_size", "ECV")
platysLM <- filtLM(d, "Platyrrhini", "Group_size", "ECV")
print(paste("Platyrrhini Beta 1:", platys[1], "and Beta 0:", platys[2], "by hand"))
print(platysLM)

strepys <- regCof(d, "Strepsirhini", "Group_size", "ECV")
strepysLM <- filtLM(d, "Strepsirhini", "Group_size", "ECV")
print(paste("Strepsirhini Beta 1:", strepys[1], "and Beta 0:", strepys[2]))
print(strepysLM)
#Step 6: From step 3, calc standard error for slope coefficient (b1)
#95% CI, and the p value associatd with this coefficient BY HAND
#Also get the same info from the results of the lm() function
#SEb1 <- sqrt(MSE/SSX), MSE is mean of remaining variance, SSX is how much
#variation there is
#SEb0 <- sqrt(MSE*sum(x^2))/(n * SSX) OR SEb1 * (sqrt(sum(x^2))/n)

#Step 6.1: Calculate MSE and SSX
SSE <- sum(ECVlm$residuals^2)
SSX <- sum((ECVlm$model$Group_size - mean(ECVlm$model$Group_size))^2)
df_error <- nrow(d) - 2
MSE <- SSE/df_error

#Step 6.2: Calc SE of b1
SEb1 <- sqrt(MSE/SSX)
print(paste("Slope coefficient standard error by hand:", SEb1))
#Step 6.3: 95% CI for b1
CIb1 <- b1 + qt(p=c(0.025, 0.975), ncp=0, df=150)*SEb1
print(paste("Upper and lower 95% confidence intervals for slope coefficient by hand:", CIb1[1], CIb1[2]))
#Step 6.4: p-value of b1, need t_stat (t = estimate/SE) to calc
#p-value (p = 2*pt(t, df = nrow - 2))
tB1 <- b1/SEb1
pB1 <- 2*pt(tB1, df = 149, lower.tail = FALSE)
print(paste("p-value for slope coefficient by hand:", pB1))
#Step 6.4: Compare results to those of the lm() function
summary(ECVlm)
lmCI <- confint(ECVlm)
lmCI

#Step 7: Run 1000 permutations to generate a null sampling distribution
#for the slope coefficient, need to permute Group_size since we want to break 
#any association it has with brain size/ECV

#Permutation test

library(mosaic)
library(broom)
reps <- 1000
slopes <- do(reps) * {
  d_new <- d
  d_new$Group_size <- sample(d_new$Group_size)
  m <- lm(data = d_new, ECV ~ Group_size)
  tidy(m) |>
    filter(term == "Group_size") |>
    pull(estimate)
}
hist(slopes$result)

permSE <- sd(slopes$result)

#Now we estimate the p value using the theory based method
#This method requires the standard error of our null sampling distribution (standard
#deviation), 
slopes_sd <- sd(slopes$result)
slopes_mean <- mean(slopes$result)
t <- (slopes_mean - b1)/slopes_sd
p_upper <- 1 - pt(abs(t), df = 150) #doing upper and lower 
p_lower <- pt(-1 * abs(t), df = 150)
p <- p_upper + p_lower
print(paste("p-value calulcated by hand useing theory-based method:", p)) #Slightly larger but still significant p-value
summary(ECVlm)

#Step 8: Use bootstrapping to get 95% CI for estimate of slope
#coefficient using both the quantile and theory based method
nb <- 1000
b <- do(nb) * {
  b_new <- d
  bs <- slice_sample(b_new, n = nb, replace = TRUE)
  bm <- lm(ECV ~ Group_size, data = bs)
  tidy(bm) |>
    filter(term == "Group_size") |>
    pull(estimate)
}

#CI by quantile
ci <- quantile(b$result, probs = c(0.025, 0.975))
print(paste("Upper and lower 95% confidence interval using quantile method:", ci[1], ci[2]))
#CI by theory based method
percent_ci <- 95
alpha <- 1 - percent_ci/100
bootMean <- mean(b$result)
bootSE <- sd(b$result)
lowerCI <- bootMean + qnorm(alpha/2)*bootSE
upperCI <- bootMean + qnorm(1 - alpha/2)*bootSE
ciTheory <- c(lowerCI, upperCI)
print(paste("Upper and lower 95% confidence interval using theory based method:", ciTheory[1], ciTheory[2]))
#Compare to original data
confint(ECVlm)
