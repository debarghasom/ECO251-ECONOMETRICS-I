#saving the directroty
direct = "C:/Users/debar/Desktop/R prog"
setwd(direct)
getwd()

install.packages("tidyverse")
install.packages("magrittr")
install.packages("stargazer")
install.packages("readxl")
install.packages("foreign")
install.packages("memisc")
install.packages("dplyr")
install.packages("labelled")
install.packages("expss")
install.packages("Hmisc")
install.packages("haven")
install.packages("rio")
install.packages("ggplot2")
install.packages("lessR")
install.packages("janitor")
install.packages("gtools")
install.packages("psych")
install.packages("wooldridge")


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
library(wooldridge)

###########################################################

x=5

abs(x)
sqrt(x)
factorial(34)
exp(x)
logb(exp(x), base=exp(1))

#generate objects immediately display using ():
(y=-3)
(z = x + y)

ls() #all the variables in the environment
rm(x) #clears the variable x
rm(list = ls()) #clears the environment


#define a vector
(a=c(1,2,3,4,5,6,7,8,9))
(b = a+1)
(c = a+b)
(d = b*c)
(e = sqrt(a))

(a=c(1,3,9,8,6,2,5,4,7))
sort(a)
length(a)
min(a)
max(a)
sum(a)
prod(a)

#character wala vector
countries = c("India","USA","Canada")
countries

#logical 
0 == 1

###########################################################

# Matrix
a = c(1,5,6,8,5,-2)
(A = matrix(a, nrow=2)) #data length must be a sub-multiple or multiple of the number of rows/cols
(Ac = matrix(a, ncol=2))
#binding elements
r1 = c(1,5,7); r2 = c(6,-2,9); r3=c(2,-4,8) #same number of elements
(B = rbind(r1, r2, r3))
(C = cbind(r1, r2, r3))
#naming the rows and cols
colnames(A) = c("x","y","z")
rownames(A) = c("a","b")
A
#Diagonal and identity matrices:
D = diag(c(1,-5,9)) #Diagoanl matrix
E = diag(4) #Identity matrix
#Indexing for extracting elements from A:
A[2,3]
A["a","y"]
#Matrix operators
A = matrix( c(4,-3,-2,1,6,0), nrow=2)
B = matrix( c(5,2,-3,5,-4,2), nrow=2)
A
B
A*B #this is not matrix multiplication(simple element multiplication)
#Transpose:
(C = t(B))
#Matrix multiplication by %*%:
(D = A %*% C)
#Inverse:
(E=solve(D))

###########################################################

#Data Frames: collection of variables with observations
A = matrix( c(4,-3,-2,1,6,0,5,6,8,9,-5,-7), nrow=4)
A
#Create a data frame and display it:
data_1 = as.data.frame(A)
data_1
view(data_1)
#attach database to the R search path
attach(data_1)#later
names(data_1) #what are the names of the variables in data_1 or dataframe

#different ways to create new variables
#Generating a new variable in the data frame:
data_1$V4 = data_1$V1 + data_1$V2 #$-sign is used to refer to the var name
#The same but using "with":
data_1$V5 = with(data_1, V1+V2) #shortens the code
#The same but using "attach":
attach(data_1)
data_1$V6 = V1+V2
detach(data_1)
#Result:
data_1
#subset of data_1 with V3>0 (will delete all other rows)
data_2 = subset(data_1, V3>0) #subset is the command

#Saving data in .RData
setwd("C:/Users/debar/Desktop/R prog/Aux_data")
save(data_1, file = "my_data_1.RData") #will save in in the current working directory

#remove data frame data_1 from memory
rm(data_1)
# Does data_1 exist?
exists("data_1")

#Load data set(from current working directory):
load("my_data_1.RData")
# Does data_1 exist?
exists("data_1")

#Export and import this small dataframe in .txt, csv,.dta format
install.packages("rio")
library(rio)
#export
rio::export(data_1, "my_data_2.txt") #.txt format
rio::export(data_1, "my_data_3.csv") #.csv format
rio::export(data_1, "my_data_4.dta") #.dta format for stata

#import
data_3<-rio::import("my_data_2.txt") #.txt format
data_4<-rio::import("my_data_3.csv") #.csv format
data_5<-rio::import("my_data_4.dta") #.dta format

###########################################################

#(2)Basic of Graphs
#(a) by using curve command: curve(function(x), xmin, xmax )
curve(x^2, -3,3)
curve(x^3, -2,2)

#(b) by using plot command: plot(x,y, type=" ") type= p l b o s h
x = c(0,2,3,5,8,9)
y = c(0,4,7,3,9,8)
#different ways of same plot
plot(x,y)        #default is p
plot(x,y, type="l") #this is solid line
plot(x,y, type="b") #line with dots but spaces
plot(x,y, type="o") #line with dots connected
plot(x,y, type="s") #horizontal and vertical connected
plot(x,y, type="h") #vertical lines


#7/03/2024
#customizing graphs
#change markers by pch
plot(x,y, pch="d") #marker by letter d (you can choose any letter)
plot(x,y, pch=0) #hollow squares
plot(x,y, pch=1) #hollow circles
plot(x,y, pch=15) #solid squares
plot(x,y, type="o", pch=16) #solid circle

#change the line type by "lty". Helps us to distinguish between the lines.
plot(x,y, type="o", pch=16, lty=1) #default lty=1 solid line
plot(x,y, type="o", pch=16, lty=2)
plot(x,y, type="o", pch=16, lty=3)

#size of points by cex=1,2,...
plot(x,y, type="o", pch=16, cex=0.5) #default is 1
#width of lines by lwd=1,2,...
plot(x,y, type="o", pch=16, lty=1, cex=1, lwd=3) 

#colors
colors() #list of colors
plot(x,y, type="o", pch=16, lty=1, cex=1, lwd=1, col="red") #col=colorname
#title and subtitle
plot(x,y, type="o", pch=16, lty=1, cex=1, lwd=1, col="black", 
     main="My Graph", sub="My Sub Graph")
#x label and y label
plot(x,y, type="o", pch=16, lty=1, cex=1, lwd=1, col="magenta", main="My Graph",
     sub="My Sub Graph", xlab="x var", ylab="y var")

#Overlaying Several Plots, use the "add" command = TRUE
#(b) plot the normal density with dnorm(x,0,n) with mean 0 and sd=n
curve(dnorm(x,0,1), -10, 10, lty=1, cex=1, lwd=1, col="red", main="My Graph",
      sub="My Sub Graph", xlab="x var", ylab="y var") #mean 0 and sd=1
curve(dnorm(x,0,2), -10, 10, add=TRUE, lty=2, cex=1, lwd=2, col="black") #mean 0 and sd=2 
curve(dnorm(x,0,3), -10, 10, add=TRUE, lty=3, cex=1, lwd=3, col="green") #mean 0 and sd=3
#legend()
legend("topright",c("SD=1","SD=2","SD=3"), lty=1:3, cex=1, lwd=1, col= c("red", "black", "green"))

###########################################################

#(2) tidyverse: set of packages for manipulation and visualization
rm(list = ls())
#set the wd
setwd("C:/Users/debar/Desktop/R prog/Raw_Data")
#import the raw data "mpg.csv"
mydata=rio::import("mpg.csv") #.csv format

install.packages("ggplot2")
library(ggplot2)

help(mpg) #this mpg is part of the ggplot2, so help will give you details
view(mydata)
attach(mydata)#to attach the dataset online for this session
names(mydata)


#ggplot geoms:can be points, lines, or other objects
#geom_point(): points
#geom_line(): lines
#geom_smooth(): nonparametric regression
#geom_area(): ribbon
#geom_boxplot(): boxplot

# Generate ggplot2 graph:
#scatter plot
ggplot() + geom_point(data=mydata, mapping=aes(x=displ, y=cty)) #aes means aesthetics we have to mention the x and y variable
#line plot
ggplot() + geom_line(data=mydata, mapping=aes(x=displ, y=cty))
#smooth(non parametric)
ggplot() + geom_smooth(data=mydata, mapping=aes(x=displ, y=cty)) #95% CI by default
ggplot() + geom_smooth(data=mydata, se=FALSE, mapping=aes(x=displ, y=cty)) #remove CI, se=FALSE and se=F are the same


#00/3/24
#Overlay the two graphs
ggplot() + 
  geom_point(data=mydata, mapping=aes(x=displ, y=cty)) +
  geom_smooth(data=mydata, mapping=aes(x=displ, y=cty))

#Describe aesthetic: aes and data at start in ggplot() only once
ggplot(mydata, aes(displ,cty)) + 
  geom_point() +
  geom_smooth(se=F)   #a more concise code 

#colors and shapes
#different colors
ggplot(mydata, aes(displ,cty) ) + 
  geom_point(color="red") +
  geom_smooth(color="yellow")
#remove CI
ggplot(mydata, aes(displ,cty) ) + 
  geom_point(color="green") +
  geom_smooth(se=FALSE, color="red")

#overlay line graphs
ggplot() + 
  geom_smooth(data=mydata,se=FALSE, color="yellow3", mapping=aes(x=displ, y=cty)) +
  geom_smooth(data=mydata,se=FALSE, color="red", mapping=aes(x=displ, y=hwy))


#Fine graphs by ggplot #to remove background colour we use theme_light
ggplot(mydata, aes(displ,cty) ) + 
  geom_point(color="green") +
  geom_smooth(se=FALSE, color="red") +
  theme_light() +
  labs(title="My ggplot graph",
       subtitle="sub ggplot",
       caption="Source: ECO251",
       x = "x var",
       y = "y var") +
  coord_cartesian(xlim=c(1,7), ylim=c(5,35))

#add legend details
ggplot() + 
  geom_smooth(data=mydata,se=FALSE, mapping=aes(x=displ, y=cty, color="cty")) +
  geom_smooth(data=mydata,se=FALSE, mapping=aes(x=displ, y=hwy, color="hwy")) +
  theme_light() +
  labs(title="My ggplot graph",
       subtitle="sub ggplot",
       caption="Source: ECO251",
       x = "x var",
       y = "y var", color = "my legend") +
  coord_cartesian(xlim=c(0,7), ylim=c(0,35) + 
                    scale_color_manual(values = c("red","green")))


#save the graph  
ggsave("C:/Users/debar/Desktop/R prog/Results/my_ggplot.png", width = 7, height = 5)

###########################################################
#12/3/2024

mydat =  rio::import("C:/Users/debar/Desktop/R prog/Raw_data/mpg.csv")
help(mpg)

#(2)Data Manipulation: dplyr package
rm(list = ls())
#set the wd
setwd("C:/Users/debar/Desktop/R prog")
getwd()
#import the raw data "ihds_women_2.dta"
mydata=rio::import("Raw_Data/ihds_women_2.dta") #.csv format
org_data=mydata
attach(mydata)
names(mydata)

#label of variables 
var_label(mydata$EW8)
var_label(mydata$EW6)

#label the variables
var_label(mydata) = list(
  EW8 = "educ of woman (yrs)",
  EW6 = "age of woman (yrs)")

#rename the variables
mydata = rename(mydata,  "EW8" = "educ_women", "EW6" = "age_women")
#label of variables 
var_label(mydata$educ_women)
var_label(mydata$age_women)

#drop variables from the data
mydata= dplyr::select(mydata, c(-SURVEY,-HHID)) #select and - sign removes the variables
#abive one is not working so use the next one!!!
#or you can use this command
mydata$PSUID<- NULL  #drop one var
mydata[,c("GE11", "GE12")]= NULL #drop selected vars
mydata[,c("SURVEY","HHID")] = NULL

#tabulate a variable
table(mydata$STATEID) #here it gives state i (1,2,etc.) has number of observations (eg.state 34 should have 98, when ran in class)
#or more precise is to count the obs(a frequency table)
count(mydata, STATEID)

#filter data only for some values or characters
mydata= filter(mydata, STATEID==1|STATEID==2)

#sort the data by STATEID codes or any other category
mydata = arrange(mydata, STATEID)

#Generate new variables: mutate function
mydata=mutate(mydata, educ_womensq=educ_women^2, 
                ln_age_women=log(age_women))

#keep only selected variables. Either one works
mydata = dplyr::select(mydata, c("educ_women", "age_women"))
mydata = mydata[, c("educ_women", "age_women")]


#plot for this clean data set
ggplot(mydata, aes(age_women, educ_women)) +
  geom_smooth(se=FALSE) + 
  theme_light() +
  labs(title="Educ vs Age",
       subtitle="IHDS-12",
       x = "age",
       y = "Education")


#keep only selected variables
mydata_1<- select(mydata, c(educ_women, age_women, ID13))

#(3) pipe operator: %>%
a=9
sqrt(9)
a %>% sqrt()

rm(list = ls())
#set the wd
setwd("C:/Users/debar/Desktop/R prog")
getwd()
#import the raw data "ihds_women_2.dta"
mydata=rio::import("Raw_Data/ihds_women_2.dta") #.csv format
org_data=mydata
attach(mydata)
names(mydata)

#label of variables 
var_label(mydata$EW8)
var_label(mydata$EW6)

#label the variables
var_label(mydata) = list(
  EW8 = "educ of woman (yrs)",
  EW6 = "age of woman (yrs)")

#pipe operator can do all the following four functions in one line
#rename the variables
mydata_1<-dplyr::rename(mydata, "educ_women"="EW8", "age_women"="EW6")
#drop variables from the data
mydata_1<- dplyr::select(mydata_1, c(-SURVEY,-HHID))
#filter data only for some values or characters
mydata_1<- dplyr::filter(mydata_1, STATEID==1|STATEID==2)
#sort the data by STATEID codes or any other category
mydata_1<- dplyr::arrange(mydata_1, STATEID)

#all the above 4 in here by %>%
mydata_2<- mydata %>% rename("EW8"="educ_women", "EW6"="age_women") %>%
  dplyr::select(c(-SURVEY,-HHID)) %>% 
  dplyr::filter(STATEID==1|STATEID==2) %>% 
  arrange(STATEID)


#(3) Missing values: NA (not available); NaN(not a number)

# Missing for all variables in data (count)
colSums(is.na(mydata_2))

# missing for particular variable
mis_educ = is.na(mydata_2$SPED6)
table(mis_educ)    #miss for 121 obs

#mean of SPED6(spouse education) variable
mean(mydata_2$SPED6) #this contains missing values so no mean is calculated
mean(mydata_2$SPED6, na.rm=TRUE) #this will give us mean of SPED6 dropping missing values


#(4)outliers and cleaning
#count without pipe
count(mydata_2,INCOME)
#count with pipe
mydata_2%>%count(INCOME)
#to understand where the outliers are, we calculate the min and max
#tab the min of Income
min(mydata_2$INCOME)
#tab the max of Income
max(mydata_2$INCOME)

#now we can cal the mean(with negative income)
mean(mydata_2$INCOME)

#replace the negative values of INCOME by missing (NA)
mydata_2$INCOME[mydata_2$INCOME<0]=NA
count(mydata_2,INCOME)

#now we can cal the new mean
mean(mydata_2$INCOME) #no mean is calculated as we have NA
mean(mydata_2$INCOME, na.rm=TRUE)


#(5) Descriptive Statistics
#count of obs for educ_women
count(mydata_2, educ_women)
#alternative way of count
table(mydata_2$educ_women)
#proportion of each category of educ_women
prop.table(table(mydata_2$educ_women))

#Fundamental Statistics
#mean
mean(mydata_2$educ_women)
#median
median(mydata_2$educ_women)
#SD
sd(mydata_2$educ_women)

#summary stat of selected variables
summary(mydata_2[, c("educ_women","age_women")])

#export a table of descriptive statistics for selected variables 
mydata_2 %>% select(educ_women,age_women,SPED6) %>% stargazer(type="text",
                                                              title="Descriptive statistics", 
                                                              digits=2, out="Results/summary_stat.rtf")

#########################################################################################

#14/3/2024

#import the raw data "wage.csv"
mydata=rio::import("Raw_Data/wage.csv") #.csv format
attach(mydata)
names(mydata)
head(mydata,10)


#(1) Descriptive Statistics
#count of obs for educ
count(mydata,educ)
#alternative way of count
table(mydata$educ)
#proportion of each category of educ
prop.table(table(mydata$educ))

#Fundamental Statistics
#mean
mean(mydata$educ)
#median
median(mydata$educ)
#SD
sd(mydata$educ)

#summary stat of selected variables
summary(mydata[, c("wage","educ")])
summary(mydata)

#export a table of descriptive statistics for selected variables 
mydata %>% dplyr::select(wage,educ,exper) %>% stargazer(type="text",
                                                 title="Descriptive statistics", 
                                                 digits=2, out="Results/summary_stat.rtf")
#data visualization
#histogram: graph to show frequency distributions
hist((mydata$educ), main = "Histogram: Education",
     xlab="education", ylab="frequency", col = "blue")

#kernel density: more sophisticated version of a histogram
plot(density(mydata$educ), main = "Density: Education",
     xlab="education", ylab="density", col = "red")
##############################################

#(2) Simple OLS
#run the simple ols: wages on education
reg1=lm(wage~educ, data=mydata)

#view the results
summary(reg1)
stargazer(reg1, type = "text")

#scatter plot
with(mydata, plot(educ,wage, main = "Scatter Wage",
                  ylab= "Wage", xlab="Educ", col="brown"))

#plot the regression line
abline(reg1, col="blue")


#Export the results
stargazer(reg1,
          omit.stat = c("f", "ser"),
          type = "text",
          out = "Results/reg_wages_educ.rtf")


#gen new variable
mydata=mydata %>% mutate(educsq=educ^2)

#run the simple ols: wages on educ square
reg2=lm(wage~educsq, data=mydata)

#run the simple ols: wages on experience
reg3<-lm(wage~exper, data=mydata)


#Export the results: col1 for reg1, col2 for reg2 and col3 for reg3
stargazer(reg1, reg2, reg3, 
          omit.stat = c("f", "ser"),
          type = "text",
          out = "Results/reg_wages_educ.rtf")



#(3) Coefficients, Fitted Values, and Residuals
names(reg1)
reg1$coefficients
reg1$fitted.values
reg1$residuals

#to store these
bhat=coef(reg1)
yhat=fitted(reg1)
uhat=resid(reg1)

#check the list in the object
list(bhat)

#(4) Non-linearity: change y to log(y)
reg4<-lm(log(wage)~educ, data=mydata)

stargazer(reg4, type = "text")


#Export the results: reg4 with reg1-reg3
stargazer(reg1, reg2, reg3, reg4, 
          omit.stat = c("f", "ser"),
          type = "text",
          out = "Results/reg_wages_educ.rtf")

#(5) Regression through origin and constant
#through origin (no intercept) we add zero
reg5=lm(wage~0+educ, data=mydata)
reg5

#on constant(no slope)
reg6=lm(wage~1, data=mydata)
reg6
#the estimate above is simple mean of y
mean(mydata$wage)

#scatter plot
with(mydata, plot(educ,wage, main = "Wage and Educ",
                  ylab= "Wage", xlab="Educ", col="black"))

#plot the regression line
abline(reg1, col="blue", lty=1)
abline(reg5, col="red", lty=2)
abline(reg6, col="green", lty=3)

#add a legend
legend("topleft", c("full","origin","constant"), lty=1:3, col= c("blue", "red", "green"))

############################################################################################
#19/3/24

mydata=rio::import("Raw_Data/wage.csv") #.csv format
attach(mydata)
names(mydata)
head(mydata,10)

mydata = apply_labels(mydata, wage = "Wage Rate ($/hr)",
                       educ = "Years of Education",
                       exper = "Years of Experience",
                       tenure = "Years of Tenure",
                       nonwhite = "Person is a Non-White",
                       female = "Person is a Female",
                       lwage = "Log of Wage Rate",
                       expersq = "Square of Years of Experience",
                       tenursq = "Square of Years of Tenure")

var_label(mydata$educ)


#(1) Multiple regression
#run the simple ols: wages on education
reg1=lm(wage~educ, data=mydata)
summary(reg1)

#Multiple Regression analysis: Let's add more covariates to this simple regression
reg2<-lm(wage~educ + exper + tenure, data=mydata)

reg3<-lm(wage~educ + exper + tenure +  expersq  + tenursq, data=mydata)


#Export the results
stargazer(reg1, reg2, reg3,
          omit.stat = c("f", "ser"),
          type = "text",
          out = "Results/mult_reg_wages.rtf")

#details about stargazer
help(stargazer)

#(2) Dummies
count(mydata, female)
count(mydata, nonwhite)
table(mydata$female)

#Dummy variable as a regressor: lets see how returns to education vary by gender and race
reg4=lm(lwage ~ educ + female + nonwhite, data=mydata)
## See the results below
stargazer(reg4, type = "text")
#interpret coeff on female dummy: wage per hour of a female is 36.1% lower than wage per hour of a male with same educ and race
#or a woman makes 36.1% less than a man with same educ and race

#Dummy interactions: married*female (four outcomes: single female, single male, married female, married male)
reg5=lm(lwage ~ educ + exper + tenure + married*female, data=mydata)
#the ref category (omitted category) for married:female is unmarried-males, i.e., when 0:0

## See the results below
stargazer(reg5, type = "text")
summary(reg5)

#reference category= single male
#married=> female=0 married=1 => married males earning 29.2%
#female=> female=1 married=0 => single females earning -9.7%

#interpret coeff on married*female:
#married-female makes=(coeff on married + coeff on female + coeff on married*female)% less/more than single-males 


#(3) Categorical variables
#Add categorical variable to the data frame
mydata$educ_cat = as.factor(ifelse(mydata$educ<=5, 'Prim', 
                                    ifelse(mydata$educ<=8, 'Mid', 
                                           ifelse(mydata$educ<=12, 'High', 
                                                  ifelse(mydata$educ<=18, 'Higher')))))

#Pie and bar plots : Distributions of categorical variable
#pie chart from 'lessR'
PieChart(educ_cat, hole = 0, values = "%", data=mydata, main="my pie chart")
#a bar plot
barplot(table(mydata$educ_cat), las=1, col=c("red", "blue" , "green", "orange")) #las=1 means lables for x-axis
#bar plot of one cat var by other cat var
barplot(table(mydata[,c("female", "educ_cat")]), las=1, col=c("red", "blue"),
        legend=TRUE, args.legend=c(x="topright"), main="my bar-plot", xlab="Education")
#side by side bars
barplot(table(mydata[,c("female", "educ_cat")]), las=1, las=1, col=c("red", "blue"),
        legend=TRUE, args.legend=c(x="topright"), main="my bar-plot", xlab="Education", beside=TRUE)


#label the female dummy options
mydata$female <- factor(mydata$female, levels=c(1,0), labels=c("Female", "Male")) 

#bar plot with label names for gender dummy
barplot(table(mydata[,c("female", "educ_cat")]), las=1, las=1, col=c("red", "blue"),
        legend=TRUE, args.legend=c(x="topright"), main="my bar-plot", xlab="Education")
############################################################################################################################################

#(A)MERGE
#left right addition of datasets
# There are three type of merging -> 1:1; 1:m; m:1
#(1)import first dataset(ihds-2 household data)
hhd_2_raw<-rio::import("Raw_Data/ihds_2_hhd.dta")
attach(hhd_2_raw)
names(hhd_2_raw)

#the unique identifier in data (check for duplicates).Since there are 42152 hh, so we use 
# state, district, village, hh and seperated IDs 
# to make them a unique ID
get_dupes(hhd_2_raw, c(STATEID,DISTID,PSUID,HHID,HHSPLITID))

#get the var label
var_label(hhd_2_raw$INCOME)   #income is at hhd level

#keep the unique identifier vars and one additional var (for analysis)
hhd_2<- dplyr::select(hhd_2_raw, c(STATEID,DISTID,PSUID,HHID,HHSPLITID,INCOME))


#(2)import first dataset(ihds-2 individual data)
ind_2_raw<-rio::import("Raw_Data/ihds_2_ind.dta")
attach(ind_2_raw)
names(ind_2_raw)

#the unique identifier in data (check for duplicates)
get_dupes(ind_2_raw, c(STATEID,DISTID,PSUID,HHID,HHSPLITID)) #not included PersonID

get_dupes(ind_2_raw, c(STATEID,DISTID,PSUID,HHID,HHSPLITID,PERSONID))  #included PersonID. Now this is the unique identifier

#get the var label
var_label(ind_2_raw$ED6)   #educ is at individual level

#keep the unique identifier vars and one additional var (for analysis)
ind_2<- dplyr::select(ind_2_raw, c(STATEID,DISTID,PSUID,HHID,HHSPLITID,PERSONID,ED6))


####Merge hhd_2 with ind_2 dataset####

hhd_ind_2<- merge(hhd_2, ind_2, by=c("STATEID","DISTID","PSUID","HHID","HHSPLITID"), all=TRUE) 
#all=TRUE will merge all the variables. Or else make a list of variables you want to merge

#ols(individual educ vs hhd income)
reg1<-lm(ED6~INCOME, data=hhd_ind_2)


#(B)APPEND
# joining from below
#(3)import first dataset(ihds-1 household data)
hhd_1_raw<-rio::import("Raw_Data/ihds_1_hhd.dta")
attach(hhd_1_raw)
names(hhd_1_raw)

#the unique identifier in data (check for duplicates)
get_dupes(hhd_1_raw, c(STATEID,DISTID,PSUID,HHID,HHSPLITID))

#get the var label
var_label(hhd_1_raw$INCOME)   #income is at hhd level

#keep the unique identifier vars and one additional var (for analysis)
hhd_1<- dplyr::select(hhd_1_raw, c(STATEID,DISTID,PSUID,HHID,HHSPLITID,INCOME))


#Append the hhd_1 and hhd_2
hhd_apend=smartbind(hhd_1,hhd_2) #no. of obs. will increase, since they are being added


#(C)RESHAPE
rm(list = ls())
getwd()


#(1)import first dataset(ihds-2 individual data)
ind_2_raw<-rio::import("Raw_Data/ihds_2_ind.dta")
attach(ind_2_raw)
names(ind_2_raw)

#get the var label
var_label(ind_2_raw$ED6)   #educ is at individual level

#keep required variables for reshape
ind_2<- dplyr::select(ind_2_raw, c(IDHH,PERSONID,ED6)) 
#IDHH is all state, village, etc all combined to make a unique id. Assignment will be to create this ID from the given data.

view(ind_2)

#tabulate the time variable(this will give us total number of new columns generated for new var)
count(ind_2,PERSONID)

#Reshape from long to wide
wide=reshape(data = ind_2,
             idvar = "IDHH",
             timevar= "PERSONID",
             v.names = "ED6",
             direction = "wide")
#idvar: is the var that will be in row
#use help to get clear idea


#Reshape back (from wide to long)
long=reshape(data = wide,
             v.names = "ED6",
             idvar = c("IDHH"),
             direction = "long")

