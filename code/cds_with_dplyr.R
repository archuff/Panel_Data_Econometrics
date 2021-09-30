# Récupération des données
library(RCurl)
urlfile<- getURL('https://raw.githubusercontent.com/archuff/Panel_Data_Econometrics/main/code/CDS.csv')
df <- read.csv(text = urlfile)

# Chargement de tidyverse

library(tidyverse)
head(df)

# Filter, select, mutate, slice, arrange
filter(df, ID=='BRAZIL')
select(df, CDS)
select(df, ID,SM)
mutate(df, mean(CDS))
slice(df, 1:157)
arrange(df,CDS)

# La puissance du pipe

slice(select((filter(df,ID=='BRAZIL')),ID,CDS), 5:10)
df %>% filter(ID=='BRAZIL') %>% select(ID,CDS) %>% slice(5:10)

# Estimation within à la main 
df <- df %>% group_by(ID) %>% mutate(DCDS = c(NA,diff(log(CDS))))
df <- df%>% drop_na(DCDS)
df <- df %>% group_by(ID) %>% mutate(meanDCDS = mean(DCDS))
df <- df %>% group_by(ID) %>% mutate(meanSM = mean(SM))
df <- df %>% mutate(WCDS = DCDS-meanDCDS)
df <- df %>% mutate(WSM = SM - meanSM)
library(plm)
model <- lm(WCDS ~ 0 + WSM, data = df)
summary(model)
model <- plm(DCDS ~ SM, data = df)
summary(model)
