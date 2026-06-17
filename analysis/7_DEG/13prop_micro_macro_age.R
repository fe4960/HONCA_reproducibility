data=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Micro_Macro_prop_age_anno",header=T,row.names=1,sep=",")

data1=data[(data$majorclass=="Microglia")&(data$tissue=="ON"),]
summary(lm(age_year~age+race+gender,data=data1))

Call:
lm(formula = age_year ~ age + race + gender, data = data1)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.73595 -0.08035  0.03068  0.12968  0.22697 

Coefficients:
                   Estimate Std. Error t value Pr(>|t|)    
(Intercept)       0.9167550  0.1451837   6.314 4.37e-08 ***
age               0.0002753  0.0012889   0.214    0.832    
raceAsian-Indian -0.2652688  0.2158847  -1.229    0.224    
raceBlack        -0.1293011  0.1131998  -1.142    0.258    
raceHispanic     -0.0582719  0.1113743  -0.523    0.603    
raceWhite        -0.0965729  0.1114905  -0.866    0.390    
raceunknown       0.0745301  0.1230052   0.606    0.547    
genderMale       -0.0373550  0.0541186  -0.690    0.493    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.1928 on 57 degrees of freedom
Multiple R-squared:  0.126,	Adjusted R-squared:  0.01863 
F-statistic: 1.174 on 7 and 57 DF,  p-value: 0.3324

summary(lm(age_year~age,data=data1))

Call:
lm(formula = age_year ~ age, data = data1)

Residuals:
     Min       1Q   Median       3Q      Max
-0.78324 -0.07690  0.09199  0.12968  0.16353

Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept) 8.334e-01  7.443e-02  11.197   <2e-16 ***
age         7.212e-05  1.264e-03   0.057    0.955
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Residual standard error: 0.1962 on 63 degrees of freedom
Multiple R-squared:  5.167e-05,	Adjusted R-squared:  -0.01582
F-statistic: 0.003255 on 1 and 63 DF,  p-value: 0.9547
