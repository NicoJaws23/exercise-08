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
ECV_sd <- sd(d$ECV, na.rm = TRUE) #Y
GS_sd <- sd(d$Group_size, na.rm = TRUE) #X

b1 <- cor(d$Group_size, d$ECV, use = "complete.obs") * (ECV_sd/GS_sd)
b0 <- mean(d$ECV, na.rm = TRUE) - (b1*mean(d$Group_size, na.rm = TRUE))

#Step 4: Use lm() to check results
ECVlm <- lm(ECV ~ Group_size, data = d, na.action = na.omit) #Does not match, ask tony?

#Step 5, Repeat this process for each primate radiation
#catarrhines, platyrrhines, strepsirhines

regCof <- function(df, rad, x, y){
  group <- df |>
    filter(Taxonomic_group = rad)
  sdX <- group |>
    sd(x, na.rm = TRUE)
  sdY <- group |>
    sd(y, na.rm = TRUE)
  b1 <- cor(group[[x]], group[[y]]) * (sdY/sdX)
  b0 <- mean(group[[y]], na.rm = TRUE) - (b1*mean(group[[x]], na.rm = TRUE))
  return(b1, b0)
}



