#saving the directroty
setwd("C:/Users/debar/Desktop/R prog/Take Home exam")

library(tidyverse)
library(magrittr)
library(stargazer)
library(readxl)
library(foreign)
library(memisc)
library(dplyr)
library(labelled)
library(expss)
library(Hmisc)
library(haven)
library(rio)
library(ggplot2)
library(lessR)
library(janitor)
library(gtools)
library(psych)

rm(list=ls())

load("C:/Users/debar/Desktop/R prog/Take Home exam/Take_home_Exam.RData")
mydata=take_home_exam
attach(mydata)
names(mydata)

#ANSWERS:
########################################################################################################################

#i) Checking for unique identifier in the dataset
var_label(mydata$ID_P)
var_label(mydata$SURVEY)

get_dupes(mydata, ID_P)
get_dupes(mydata, c(ID_P,SURVEY))

########################################################################################################################

#ii) Summary statistics for key variables
sub=na.omit(dplyr::select(mydata, c("widow","widow_support","hh_assets","std_edu","hh_children","ln_cons","girl","cons_capita")))

summary_stats <- summary(sub)
write.csv(summary_stats, file = "summary_statistics.csv")

#stargazer(sub, type="text", title="Descriptive statistics", digits=3, out="summary_stat.rtf")
#This command is not working for the given dataset

########################################################################################################################

#iii) Is it true the households headed by widowed mothers are poorer (in the sense that they earn less)?
reg5=lm(cons_capita~widow, data=mydata)
reg1=lm(log(cons_capita)~widow, data=mydata)
stargazer(reg5, reg1, type = "text",title="Regression Table (Table 2)", digits=3, out="Regression Table q3.rtf")

########################################################################################################################
#(iv) In studying the association between widowhood and educational attainments, why does
#it make sense to use standardized years of education rather than completed years of
#education? 
########################################################################################################################

#v) Comment on the association between widowhood and children’s educational attainment.
#What background characteristics of the household would you like to control for?
reg2=lm(std_edu~widow, data=mydata)
reg3=lm(std_edu~widow + hh_children + hh_assets, data=mydata)

stargazer(reg2,reg3, type = "text",title="Regression Table (Table 3)", digits=3,out="Regression Table q5.rtf")


########################################################################################################################

#vi)Does widowhood appear to affect girls differently than it affects boys? 
reg4=lm(std_edu~ widow*girl  + widow + girl, data = mydata)
stargazer(reg4, type = "text", title="Regression Table (Table 4)", digits=3,out="Regression Table q6.rtf")

########################################################################################################################