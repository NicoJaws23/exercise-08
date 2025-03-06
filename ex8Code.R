#Step 1, load and look at data
library(tidyverse)
library(skimr)
f <- "https://raw.githubusercontent.com/difiore/ada-datasets/main/Street_et_al_2017.csv"
d <- read_csv(f)
skim(d)

#Step 2: Plots ECV as a function of social group size,
#longevity, weaning, and reproductive lifespan
ECV_GS <- plot(data = d, ECV ~ Group_size)
ECV_Long <- plot(data = d, ECV ~ Longevity)
ECV_Wean <- plot(data = d, ECV ~ Weaning)
ECV_Repro <- plot(data = d, ECV ~ Repro_lifespan)

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
  b1 <- cor(group[[x]], group[[y]]) * (sdY/sdX)
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

platys <- regCof(d, "Platyrrhini", "Group_size", "ECV")
platysLM <- filtLM(d, "Platyrrhini", "Group_size", "ECV")

strepys <- regCof(d, "Strepsirhini", "Group_size", "ECV")
strepysLM <- filtLM(d, "Strepsirhini", "Group_size", "ECV")

#Step 6: From step 3, calc standard error for slope coefficient (b1)
#95% CI, and the p value associatd with this coefficient BY HAND
#Also get the same info from the results of the lm() function


