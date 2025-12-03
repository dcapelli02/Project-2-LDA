## Libraries
library(haven)
library(sas7bdat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(GGally)
library(corrplot)
library(nlme)
library(lmtest)
library(DescTools)
library(Matrix)
library(MASS)
library(metafor)
library(geepack)
library(car)
library(lme4)
library(ORTH.Ord)
library(MuMIn)
library(patchwork)

## Import and fix the data
alz <- read_sas("C:/Users/Daniele/Desktop/2025 - 26 Primo Semestre/Longitudinal Data Analysis/Project 1 Alzheimer LDA/alzheimer25.sas7bdat")

head(alz)
summary(alz)

alz$trial <- as.factor(alz$trial)
alz$sex <- as.factor(alz$sex)
alz$edu <- as.factor(alz$edu)
alz$job <- as.factor(alz$job)
alz$wzc <- as.factor(alz$wzc)
alz$adl <- as.factor(alz$adl)
alz$adl_num <- as.numeric(alz$adl)
alz$n_obs_data <- rowSums(!is.na(alz[, c(18:24)]))

## Bins

alz$cdrsb_bin0 <- alz$cdrsb0
alz$cdrsb_bin0[alz$cdrsb0 > 10] <- 1
alz$cdrsb_bin0[alz$cdrsb0 <= 10] <- 0

alz$cdrsb_bin1 <- alz$cdrsb1
alz$cdrsb_bin1[alz$cdrsb1 > 10] <- 1
alz$cdrsb_bin1[alz$cdrsb1 <= 10] <- 0

alz$cdrsb_bin2 <- alz$cdrsb2
alz$cdrsb_bin2[alz$cdrsb2 > 10] <- 1
alz$cdrsb_bin2[alz$cdrsb2 <= 10] <- 0

alz$cdrsb_bin3 <- alz$cdrsb3
alz$cdrsb_bin3[alz$cdrsb3 > 10] <- 1
alz$cdrsb_bin3[alz$cdrsb3 <= 10] <- 0

alz$cdrsb_bin4 <- alz$cdrsb4
alz$cdrsb_bin4[alz$cdrsb4 > 10] <- 1
alz$cdrsb_bin4[alz$cdrsb4 <= 10] <- 0

alz$cdrsb_bin5 <- alz$cdrsb5
alz$cdrsb_bin5[alz$cdrsb5 > 10] <- 1
alz$cdrsb_bin5[alz$cdrsb5 <= 10] <- 0

alz$cdrsb_bin6 <- alz$cdrsb6
alz$cdrsb_bin6[alz$cdrsb6 > 10] <- 1
alz$cdrsb_bin6[alz$cdrsb6 <= 10] <- 0


## Create baseline values
alz$ab_base <- alz$abpet0
alz$tau_base <- alz$taupet0
alz$cdrsb_base <- alz$cdrsb0
alz$bprs_base <- alz$bprs0

summary(alz)

## Create longitudinal dataset
alz_df <- data.frame(alz)

alz_long <- alz_df %>%
  pivot_longer(
    
    cols = matches("^(bprs|cdrsb_bin|abpet|taupet)\\d+$"),
    
    
    names_to = c(".value", "year"),
    
    names_pattern = "(bprs|cdrsb_bin|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                          
    sample = factor(rep(1:nrow(alz_df), each = 7)) # ID per ogni paziente
  )

alz_long$cdrsb <- alz_long$cdrsb_bin
alz_long$cdrsb_bin <- as.factor(alz_long$cdrsb_bin)

## Discretize variables

## Maybe any 5 years?
alz_long$age_disc <- (alz_long$age %/% 5) * 5
alz_long$age_disc <- as.factor(alz_long$age_disc)

## bmi any 4
alz_long$bmi_disc <- (alz_long$bmi %/% 4) * 4
alz_long$bmi_disc <- as.factor(alz_long$bmi_disc)

## adl any 5
alz_long$adl_disc <- (as.numeric(alz_long$adl) %/% 5) * 5
alz_long$adl_disc <- as.factor(alz_long$adl_disc)

## abpet any 0.2
alz_long$abpet_disc <- (alz_long$abpet %/% 1) * 1
alz_long$abpet_disc <- as.factor(alz_long$abpet_disc)

## taupet any 0.2
alz_long$taupet_disc <- (alz_long$taupet %/% 1) * 1
alz_long$taupet_disc <- as.factor(alz_long$taupet_disc)

## adl any 5
alz$adl_disc <- (as.numeric(alz$adl) %/% 5) * 5
alz$adl_disc <- as.factor(alz$adl_disc)


# 4 gruppi (quartili)
alz_long <- alz_long %>%
  mutate(
    age_q4 = cut(age, breaks = quantile(age, probs = seq(0,1,0.25), na.rm = TRUE),
                 include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4")),
    bmi_q4  = cut(bmi,  breaks = quantile(bmi,  probs = seq(0,1,0.25), na.rm = TRUE),
                  include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4")),
    adl_q4  = cut(adl_num,  breaks = quantile(adl_num,  probs = seq(0,1,0.25), na.rm = TRUE),
                  include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))
  )

# con ntile (interi 1..4)
alz_long <- alz_long %>%
  mutate(
    age_ntile4 = ntile(age, 4),
    adl_ntile4 = ntile(adl, 4)
  )

## year discrete
alz_long$year_seq <- ave(alz_long$year, alz_long$sample,
                         FUN = function(x) as.integer(factor(x)))

#####################################################################
#####################################################################

#### INITIAL EXPLORATORY DATA ANALYSIS ####

ggplot(alz_long, aes(x = year, y = cdrsb)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun.data = mean_sdl, geom = "ribbon", fun.args = list(mult = 1),
               alpha = 0.2, fill = "lightblue") +
  theme_minimal() 
# + labs(title = "Fraction of CDRSB 1 over year")

## Sembra una funzione sigmoide --> la logistic potrebbe fittarla bene

p_sex <- ggplot(alz_long, aes(x = year, y = cdrsb, group = sex, color = sex)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  theme_minimal() +
  labs(title = "Mean trajectory by SEX")

p_age <- ggplot(alz_long, aes(x = year, y = cdrsb, group = age_q4, color = age_q4)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  theme_minimal() +
  labs(title = "Mean trajectory by AGE")

p_bmi <- ggplot(alz_long, aes(x = year, y = cdrsb, group = bmi_q4, color = bmi_q4)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  theme_minimal() +
  labs(title = "Mean trajectory by BMI")

p_adl <- ggplot(alz_long, aes(x = year, y = cdrsb, group = adl_q4, color = adl_q4)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  theme_minimal() +
  labs(title = "Mean trajectory by ADL")

p_ab <- ggplot(alz_long, aes(x = year, y = cdrsb, group = abpet_disc, color = abpet_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  theme_minimal() +
  labs(title = "Mean trajectory by Aβ")

p_tau <- ggplot(alz_long, aes(x = year, y = cdrsb, group = taupet_disc, color = taupet_disc)) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  theme_minimal() +
  labs(title = "Mean trajectory by TAU")

# Composition 2 × 3
final_plot <- (p_sex | p_age | p_bmi) /
  (p_adl | p_ab | p_tau)

final_plot

# Saving
ggsave("composizione.png", final_plot, width = 12, height = 8, dpi = 300)

### EMPIRICAL VARIANCE
var_by_year <- alz_long %>%
  group_by(year) %>%
  summarise(var_cdrsb = var(cdrsb, na.rm = TRUE))

ggplot(var_by_year, aes(x = year, y = var_cdrsb)) +
  geom_line(size = 1) + geom_point(size=2.5) +
  theme_minimal() 
#labs(title = "Variance of CDRSB over time")


### CORRELATION
cor_matrix_cdrsb <- cor(alz[, c(41:47)], use = "pairwise.complete.obs")
round(cor_matrix_cdrsb, 2)

heatmap(cor_matrix_cdrsb, Rowv = NA, Colv = NA)  ## --> come interpretare??
print(cor_matrix_cdrsb)


### SPAGHETTI PLOT
## For the spaghetti plot we need a random sample

set.seed(1)
casual1 <- sample(1:length(alz$patid), 20)

## Now we can start by looking at random values for the mean and see
## if we can work on the mean and so on

alz_rist1 <- alz[casual1, ]
alz_rist1_df <- data.frame(alz_rist1)

alz_long_rist <- alz_rist1_df %>%
  pivot_longer(
    
    cols = matches("^(bprs|cdrsb_bin|abpet|taupet)\\d+$"),
    
    
    names_to = c(".value", "year"),
    
    names_pattern = "(bprs|cdrsb_bin|abpet|taupet)(\\d+)"
  ) %>%
  mutate(
    year = as.numeric(year),                          
    sample = factor(rep(1:nrow(alz_rist1_df), each = 7)) # ID per ogni paziente
  )

alz_long_rist$cdrsb <- alz_long_rist$cdrsb_bin
alz_long_rist$cdrsb_bin <- as.factor(alz_long_rist$cdrsb_bin)

ggplot(alz_long_rist, aes(x = year, y = cdrsb, group = patid, 
                          color = sample, show.legend = FALSE)) + 
  geom_line(alpha = 1, show.legend = FALSE, size = 0.5) +
  theme_bw() +
  labs(title = "Time Evolution of CDRSB")



