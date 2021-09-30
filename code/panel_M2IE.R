########################################################
##### Code R pour les exemples du cours de  Panel  #####
#####   Universtié de Franche-Comté - Octobre 2018 #####
########################################################

## Chapitre 1: example estimation Fixed Effects, Random Effects et Pooled après avoir 
## crée un Panel à l'aide du package PSIDR et le Panel Study of Income Dynamics

rm(list = ls())
load('panel_psid.RData')

# Build PSID panel data
# Chargement des packages utiles

library(plm)
library(psidR)
library(knitr)
library(data.table)
# Part 1 - Non obligatoire si RData load
r = system.file(package="psidR")
# Les fichiers suivants doivent être placés dans le dossiers psid-lists du package psidR
f = fread(file.path(r,"psid-lists","famvars.txt")) # Les variables fam que l'on souhaite dans le PSID. 
i = fread(file.path(r,"psid-lists","indvars.txt")) # Les variables indiv que l'on souhaite dans le PSID

f[1:4,vgroup := "own"]
f[5:8,vgroup := "hvalue"]
f[9:12,vgroup := "mortgage"]
f[13:16,vgroup := "wage"]
f[17:20,vgroup := "educ"]
f[21:24,vgroup := "experience"]
f[25:28,vgroup := "experience_month"]
setkey(f,vgroup)

i[1:4,vgroup := "weight"]
i[5:8,vgroup := "age"]
i[9:12,vgroup := "empstat"]
i[13:16,vgroup := "educ"]
setkey(i,vgroup)

ind = cbind(i[J("age"),list(year,age=variable)],
            i[J("educ"),list(educ=variable)],
            i[J("weight"),list(weight=variable)])

fam = cbind(f[J("wage"),list(year,wage=variable)],
            f[J("educ"),list(educ=variable)],
            f[J("experience"),list(experience=variable)],
            f[J("experience_month"),list(experience_month=variable)])

# ATTENTION, pour lancer cette commande il vous faut un login et un mot de passe sur le site du PSID
# https://psidonline.isr.umich.edu/

d = build.panel(fam.vars=fam,
                ind.vars = ind,
                SAScii = TRUE, 
                heads.only = TRUE,
                sample="SRC",
                design=2)

# Part 2 - Clean and regress
# Pour le cours, j'ai supprimé quelques observations. 
# ATTENTION, pour une analyse correcte, il faudrait nettoyer la base plus rigoureusememt

df <- d$data[!(d$data$wage == 0),]
df <- df[!(df$educ == 99),]
df <- df[!(df$experience == 99),]
df <- df[!(df$experience == 98),]
df <- df[!(df$wage == 9999999),]
df <- df[!(df$wage == 9999998),]
df <- pdata.frame(df, index = c("pid","year"))
df = make.pbalanced(df,balance.type = c("shared.individuals"))
df$exptot = df$experience*12 + df$experience_month +1
df$exptot2 = df[, "exptot"]^2

summary(plm(log(wage) ~ educ  + exptot, data = df, index = c("pid","year"), 
            model = "pooling"))
summary(plm(log(wage) ~ educ  + exptot , data = df, index = c("pid","year"), 
            model = "within"))
summary(plm(log(wage) ~ 0+ educ  + exptot , data = df, index = c("pid","year"), 
            model = "random"))
x <- plm(log(wage) ~ educ  + exptot , data = df, index = c("pid","year"), 
       model = "within") 
x2 <- plm(log(wage) ~ 0+ educ  + exptot , data = df, index = c("pid","year"), 
          model = "random")
phtest(x,x2)

## Chapitre 2 - Simulation du biais dynamique
## Programme principal

g = (0:7)/10
t = c(10,30,50,100)
R = 500 ## Nombre de réplications
g_hat <- rep(NA,R)
biais_g <- matrix(0, nrow = length(g), ncol = length(t))
biais_gam_hat <- rep(NA,R)
mg_hat <- matrix(0, nrow = length(g), ncol = length(t))
for (l in 1:length(t)){
  for (k in 1:length(g)){
    G <- g[k] ## Paramètre autorégressif
    T <- t[l]  ## Nombre de périodes
    N <- 10 ## Nombre d'individus
    year <- get_year(T,N) ## Construction de l'index temporel
    id <- get_id(T,N)
    for (r in 1:R){ ## Boucle sur les réplications
      g_hat[r] <- coeff_lsdv_arsim(T,N,G)
      biais_gam_hat[r] <- g_hat[r] - G
    }
    mg_hat[k,l] <- mean(g_hat)
    biais_g[k,l] <- mean(biais_gam_hat)
  }
}
col_set <- rainbow(4)
matplot(g,biais_g, type = 'l', col = col_set)
legend("bottomleft", c("T=10", "T=30","T=50", "T=100"), col = col_set, lty = 1:4)

## Fonction permettant de créer l'index temporel
get_year <- function(t,n){
  return(rep(1:t,n))
}

## Fonction permettant de créer l'index individuel
get_id <- function(t,n){
  id <- rep(0,(t*n))
  for (i in 1:n){
    id[(1+(t*(i-1))):(t*i)] <- rep(i,t)
  }
  return(id)
}

## Fonction simulant les données
coeff_lsdv_arsim <- function(t,n,g){
  alpha <- runif(n,-1,1)  ## Simulation des paramètres non-observés
  y <- array(rep(0, (t+1)*n), dim=c(t+1, n)) ## Initialisation de la variable dépendante
  e <- array(rnorm((t+1)*n), dim=c(t+1, n)) ## Simulation des erreurs
  for (t in 2:(t+1)){ ## On simule la variable expliqué
    y[t,] <- alpha + g*y[t-1,] + e[t,]
  }
  y0 <- y[2:t,] ## y0 est la variable dépendante
  y1 <- y[1:(t-1),] ## y1 est le lag de la variable dépendante
  y0 <- c(y0)
  y1 <- c(y1)
  df <- data.frame(id,year,y0,y1) ## Construction du dataframe
  ## Estimateur LSDV
  lsdv <- plm(y0 ~ y1, index = c("id","year"), data = df, model = "within") 
  gam_hat <- lsdv$coefficients
  return(gam_hat)
}
