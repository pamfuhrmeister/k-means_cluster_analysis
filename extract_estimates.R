#!/usr/bin/env Rscript

# Data wrangling and fitting models to extract individual estimates for k-means-cluster analysis

# setwd() for directory where file is----
this_dir <- function(directory)
  setwd( file.path(getwd(), directory) )

# load packages----
library(tidyverse)
library(nlme)
library(ggthemes)
library(data.table)
library(lme4)

# fit 3-parameter logistic model for ba-da continuum, extract participant estimates for slope----
data <- read.csv("native_speech_data_raw.csv", skip = 1, header = TRUE)

bada <- data %>%
  filter(Procedure == "bada") %>%
  select("Subject", "Sound_file", "Slide1.Slider1.Value")

colnames(bada) <- c("Subject", "ContinuumPoint", "SliderResponse")
bada$ContinuumPoint <- gsub("ba-da-ga_00","", bada$ContinuumPoint)
bada$ContinuumPoint <- as.numeric(bada$ContinuumPoint)
bada <- bada %>%
  mutate(ContinuumPoint = (ContinuumPoint-1))

# logistic
bada$ContinuumPoint <- as.numeric(as.character(bada$ContinuumPoint))
bada$SliderResponse <- as.numeric(as.character(bada$SliderResponse))
bada$Subject <- as.factor(as.character(bada$Subject))

m1 <- nlme(SliderResponse ~ max*((exp(slope*(ContinuumPoint-boundary)))/(1+exp(slope*(ContinuumPoint-boundary)))),
           data = bada,
           fixed = max + boundary + slope ~ 1,
           random = list(Subject = pdDiag(max + boundary + slope ~ 1)),
           start = c(7,3.4,.7))

bada$ba_da_resid <- resid(m1)
bada$ba_da_resid2 <- bada$ba_da_resid^2

speech_data <- as.data.frame(coef(m1))
speech_data <- mutate(speech_data, Subject = rownames(speech_data))
speech_data <- speech_data[,c(4,3)]
colnames(speech_data)[2] <- "slope_ba_da"

temp <- bada %>%
  group_by(Subject) %>%
  summarize(consistency_ba_da = mean(ba_da_resid2)) %>%
  mutate(consistency_ba_da = -1*consistency_ba_da)

colnames(temp)[1] <- "Subject"

speech_data <- merge(speech_data, temp, by = "Subject")


# fit 3-parameter logistic model for s-sh continuum, extract participant estimates for slope----

signshine <- data %>%
  filter(Procedure == "signshine") %>%
  select("Subject", "Sound_file", "Slide2.Slider1.Value") %>%
  mutate(Slide2.Slider1.Value=8-Slide2.Slider1.Value) %>%
  mutate(ContinuumPoint = dplyr::recode(Sound_file, 
                                        "SIGNblend_20s80h" = 1,
                                        "SIGNblend_30s70h" = 2,
                                        "SIGNblend_40s60h" = 3,
                                        "SIGNblend_50s50h" = 4,
                                        "SIGNblend_60s40h" = 5,
                                        "SIGNblend_70s30h" = 6,
                                        "SIGNblend_80s20h" = 7)) %>%
  select(-Sound_file)

signshine <- signshine[,c(1,3:2)]
colnames(signshine)[3] <- "SliderResponse"

# 3-parameter logistic
signshine$ContinuumPoint <- as.numeric(as.character(signshine$ContinuumPoint))
signshine$SliderResponse <- as.numeric(as.character(signshine$SliderResponse))
signshine$Subject <- as.factor(as.character(signshine$Subject))

m2 <- nlme(SliderResponse ~ max*((exp(slope*(ContinuumPoint-boundary)))/(1+exp(slope*(ContinuumPoint-boundary)))),
            data = signshine,
            control = nlmeControl(maxIter = 200000),
            fixed = max + boundary + slope ~ 1,
            random = list(Subject = pdDiag(max + boundary + slope ~ 1)), 
            start = c(7,4,.5))

signshine$s_sh_resid <- resid(m2)
signshine$s_sh_resid2 <- signshine$s_sh_resid^2

s_sh <- as.data.frame(coef(m2))
s_sh <- mutate(s_sh, Subject = rownames(s_sh))
s_sh <- s_sh[,c(4,3)]
colnames(s_sh)[2] <- "slope_s_sh"

temp2 <- signshine %>%
  group_by(Subject) %>%
  summarize(consistency_s_sh = mean(s_sh_resid2)) %>%
  mutate(consistency_s_sh = -1*consistency_s_sh)

s_sh <- merge(s_sh, temp2, by = "Subject")
speech_data <- merge(speech_data, s_sh)


# non-native speech sound discrimination task, fit model and extract estimates----
setwd("./non-native_speech_task/")
file_names <- list.files(pattern = "*.csv") #where you have your files

read <- function(x){
  fread(x, header = TRUE, select = c("subject_nr", "block", "condition", "correct", "ID_training_1_soundfile", "posttest_discrim_soundfile", "pretest_discrim_soundfile"))
}

df <- lapply(file_names, read) %>% bind_rows()

df_pre <- df %>%
  filter(block == "discrim_pretest") %>%
  select(-c(ID_training_1_soundfile, posttest_discrim_soundfile, block))
colnames(df_pre) <- c("Subject", "condition", "correct", "item")

# change 1s and 0s to numeric
df_pre$correct <- as.numeric(as.character(df_pre$correct))

df_pre <- df_pre %>%
  mutate(ans_dif = ifelse(correct == 1 & condition == "Different", 1,
                          ifelse(correct == 0 & condition == "Same", 1, 0)))

df_pre$Subject <- as.factor(as.character(df_pre$Subject))

df_pre$condition <- as.factor(df_pre$condition)
contrasts(df_pre$condition) <- c(.5, -.5)

m3 <- glmer(ans_dif ~ condition + (condition|Subject) + (1|item),
            data = df_pre,
            family = binomial(link = "probit"))

discrim <- as.data.frame(coef(m3)$Subject)
discrim <- discrim %>%
  rownames_to_column(var="Subject") %>%
  select(c(Subject, condition1))
colnames(discrim)[2] <- "dprime"

speech_data <- merge(speech_data, discrim, by = "Subject")
write.csv(speech_data, "../speech_data.csv")


