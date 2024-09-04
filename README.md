# survival-analysis-Getug-12
Analyse de survie avancée en utilisant la modélisation multi-états sur les données de l'essai clinique randomisé GETUG 12. L'objectif est d'évaluer l’effet du  traitement expérimental sur les transitions entre les différents états de la maladie, afin d'améliorer les approches thérapeutiques en oncologie.
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
# rm(list=ls(all=T))
# options(digits = 3)
library(haven)
library(dplyr)
library(tableone)
library(knitr)
library(mstate)
library(survival)
library(RColorBrewer)
library(tidyverse)
library(DT)
library(ggplot2)
library(survminer)
library(gridExtra)
library(etm)
library(RColorBrewer)
library(reshape2)
library(cowplot)
library("colorspace")
# install.packages("caret")
library(caret)

```

# Datamanagement

## Import du jeu de données

```{r}
# getug_i <- read_sas("C:/Users/I_BELHECHMI/Desktop/nouvelle base de donnees/multi_stage.sas7bdat", NULL)
getug_i <- read_sas("C:/Users/Khaled/Documents/multi_stage_nvl.sas7bdat", NULL)
getug <- getug_i
names(getug) <-  tolower(names(getug)) # Mettre les noms des variables en minuscule
getug <- getug[, c("numtas","glearec", "trec", "psarec", "pnrec", "trait", "age", "age_cl",
                   "nb_crit", "nb_crit_cl", "evt", "dattas", "phou", "datpsah", "pc", "datpc", "pm",
                   "datpm", "etat", "datdc", "fin_multi_etat", "ddn")]
#head(getug)
#names(getug)
View(getug)

```

## Correction d'une erreur de code

On note que : un patient (numtas = 5) a une date de décès différente de (<) la date fin_multi_etat => Remplacer date fin_multi_etat = date de décès


```{r}

# Supprimer les patients qui ont des trajets diffrents des transitions de notre modèle
# Les patients 126 et 335 ont fait Rx->PC->PB->DC => Ils ont des trajets diffrents des transitions de notre modèle
getug <- getug %>%
            filter(!numtas %in% c(4, 71, 140, 223, 360, 126,335 ))

# Les patients qui ont fait Rx -> PM->DC => On a éliminé la transition Rx->PM donc ces patients sont considérés avec les patients qui sont décédés directement
        
getug <- getug %>%
            filter(!numtas %in% c(36, 37, 72, 175, 199, 209,412)) 

# Les patients qui ont fait Rx->PC/PM->censuré => Sont considérés avec les patients censurés après Rando sans faire aucun evenement puisqu'on a ignoré la transition Rx -> PC/PM => On va supprimer ces patients
getug <- getug %>%
            filter(!numtas %in% c(109,122, 174, 201, 206, 222, 227, 245,  250, 251, 258,276, 298, 318, 329, 403, 409))

# Les patients qui ont fait une progression biologique puis décédés => Supprimer ces patients car on a supprimer la transition pb->dc
getug <- getug %>% 
            filter(!numtas %in% c(5, 12, 15, 32, 44, 53, 54, 79, 83, 143, 159, 166, 184, 243, 247, 264, 302, 353, 362, 386))

# => On a supprimé 51 patients + 4 patients(eliminé à la préparation de la base de données): reste 358 patients

```

## Caractéristiques initiales des patients et traitement

```{r}
class(getug$age_cl)
getug$age_c <- ifelse(getug$age_cl==1, 0, 1)
#table(getug$age_c, getug$age_cl, useNA = "always")

getug$age_cf <- factor(getug$age_c, 0:1, c("<65",">=65"))
# table(getug$age_cl, getug$age_c, useNA = "always")

getug$glearec_f <- factor(getug$glearec, 0:1, c("<8",">=8"))
#table(getug$glearec_f, getug$glearec, useNA = "always")

getug$pnrec_f <- factor(getug$pnrec, 0:1, c("N0","N1"))
#table(getug$pnrec_f, getug$pnrec, useNA = "always")

getug$trec_f <- factor(getug$trec, 0:1, c("T1-T2","T3-T4"))
#table(getug$trec_f, getug$trec, useNA = "always")

getug$psarec_f <- factor(getug$psarec, 0:1, c("=<20",">20"))
#table(getug$psarec_f, getug$psarec, useNA = "always")

getug$trt <- ifelse(getug$trait==2, 0, 1)  # ADT=2 et ADT+DE=1  ==> ADT=0 et ADT+DE=1
#table(getug$trt, getug$trait, useNA = "always")
getug$trt_f <- factor(getug$trt, 0:1, c("TAA","TAA + DE"))
#table(getug$trt_f, getug$trt, useNA = "always")

getug$nb_crit_f <- as.factor(getug$nb_crit)
# table(getug$nb_crit_f, getug$nb_crit, useNA = "always")
```

Crée une table descriptives des caractéristiques des patients lors de la randomisation, stratifiée par la variable traitement.

```{r}
var_base <- c("age_cf","glearec_f","trec_f", "psarec_f", "pnrec_f","nb_crit_f")
baseline_char <- CreateTableOne(vars = var_base, strata = c("trt_f"), addOverall = T, data = getug)
print(baseline_char, showAllLevels = TRUE, test=FALSE, formatOptions = list(big.mark = ","), nonnormal ="age")

```

## Evénements

- Progression biochimique, Progression local/clinique + Progression à distance/métastase, Décès

```{r}
# Valeurs manquantes
# sum(is.na(getug$ddn))

getug$pb <- ifelse(is.na(getug$phou), 0, getug$phou)
getug$pc2 <- getug$pc
getug$pc <- ifelse(is.na(getug$pc2), 0,  getug$pc2)
getug$pm2 <- getug$pm
getug$pm <- ifelse(is.na(getug$pm2), 0,  getug$pm2)
getug$dc <- ifelse(is.na(getug$etat), 0, getug$etat)


# la variable "pc_pm" est un indicateur binaire (0 ou 1) qui indique si le patient a connu un événement de progression clinique ou métastatique 
getug$pc_pm <- ifelse(getug$pc == 1 | getug$pm == 1, 1, 0)

getug$pb_f <- as.factor(getug$pb)
getug$pc_f <- as.factor(getug$pc)
getug$pm_f <- as.factor(getug$pm)
getug$pc_pm_f <- as.factor(getug$pc_pm)
getug$dc_f <- as.factor(getug$dc)


# Evénéments (sans le traitement de rattrapage)
getug$event <- ifelse(getug$pb==1 | getug$pc_pm==1 | getug$dc==1, 1, 0)
# table(getug$event, useNA="always")
getug$event_f <- as.factor(getug$event)

View(getug)
```

## Délai d'apparition des événements

```{r}

# En jours
getug$futimed <- as.numeric(difftime(getug$fin_multi_etat, getug$dattas, units = "days"))
getug$pbtimed <- as.numeric(difftime(getug$datpsah, getug$dattas, units = "days"))
getug$pctimed <- as.numeric(difftime(getug$datpc, getug$dattas, units = "days"))
getug$pmtimed <- as.numeric(difftime(getug$datpm, getug$dattas, units = "days"))
getug$dctimed <- as.numeric(difftime(getug$datdc, getug$dattas, units = "days"))

# En années
getug$futime <- as.numeric(getug$futimed / 365.25) # follow-up time
getug$pbtime <- as.numeric(getug$pbtimed / 365.25)
getug$pctime <- as.numeric(getug$pctimed / 365.25)
getug$pmtime <- as.numeric(getug$pmtimed / 365.25)
getug$dctime <- as.numeric(getug$dctimed / 365.25)

# Créer la variable pc_pmtimed  pc_pmtimed
getug$pc_pmtimed <- pmin(getug$pctimed, getug$pmtimed, na.rm = TRUE)
# summary(getug[,c("futimed","pbtimed","pctimed","pmtimed","dctimed", "pc_pmtimed")])
getug$pc_pmtime <- as.numeric(getug$pc_pmtimed / 365.25 )
# summary(getug[,c("futimed","pbtimed","pctimed","pmtimed","dctimed", "pc_pmtimed")])


# Supprimer des colonnes
# getug <- getug %>% select(-pc2, -phou, -pm2, -etat, -age, -nb_crit_cl, -age_cl, -pc, -pm, -pctime, -pmtime, -pctimed, -pmtimed,
#                           -datpsah, -datdc, -datpm, -datpc, -dattas, -ddn, -fin_multi_etat,
#                           -nb_crit_cl, -nb_crit_f, -nb_crit, -trait)
# View(getug)

```

## Kaplan Meier

```{r} 
# fit_tvim_os <- survfit(Surv(futime, status) ~ 1, data = getug)
#  
# fit_tvim_os
#  
# ggsurvplot(fit_tvim_os,
#            pval = TRUE, conf.int = TRUE,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk http://127.0.0.1:47319/graphics/plot_zoom_png?width=1920&height=1017table color by groups
#            linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c("#E7B800", "#2E9FDF"))
```

## Délai d'apparition du premier événement 

```{r}
getug$time_to_prog <- ifelse(getug$pb==1 | getug$pc==1 | getug$pm==1 | getug$dc==1, 
                             pmin(getug$pbtime, getug$pctime, getug$pmtime, getug$dctime, na.rm = T), 
                             getug$futime)
#getug[,c("numtas","pbtime","pctime","pmtime","dctime","ratime","futime","event_n","time_to_prog")]
```

## Modèle multi-état Simplifié

Structure du modèle: définir les transitions possibles

```{r, fig.height=9}

figmodel <- function(what, ...) {
    sname <- c("1: RX", "2: PB", "3: PC + PM", "4: DC")
    connect <- matrix(0, 4, 4, dimnames=list(sname, sname))
    connect[1, c(2, 4)] <- c(1, 1.1)
    connect[2, 3] <- 1
    connect[3, 4] <- 1
    statefig(matrix(c(1,2,1)), connect, cex=.8, ...)
}

figmodel()

mtext("RX: Randomisation", side=2, adj=0.5, cex=.75, font=1, line=4)
mtext("PB: Progression biochimique", side=2, adj=0.5, cex=.75, font=1, line=3)
mtext("PC + PM: Progression locale/Progression à distance", side=2, adj=0.5, cex=.75, font=1, line=2)
mtext("DC: Décès", side=2, adj=0.5, cex=.75, font=1, line=1)

mtext("Rx -> Progression biochimique : transition 1", side=1, adj=0.5, cex=0.5, font=2, line=1)
mtext("Rx -> Décès :  transition 2",side=1, adj=0.5, cex=0.5, font=2, line=2)
mtext("PB -> PC + PM : transition 3",side=1, adj=0.5, cex=0.5, font=2, line=3)
mtext("PC + PM -> Décès : transition 4",side= 1, adj=0.5, cex=0.5, font=2, line=4)

```

6 transitions possibles : 

1.  Randomisation -> Progression biochimique (1 -> 2) 
2.  Randomisation -> Décès (1 -> 4) 
3.  Progression biochimique -> Progression locale ou à distance (2 -> 3)   
4.  Progression locale ou à distance -> Décès (3 -> 4)  

## Praparation des données

<span style="color:#e40000; font-weight:bold; font-size: larger;"> Remarque:</span> 

<span style="color:#e40000; font-weight:bold; font-size: larger;"> Pour les patients qui n'ont pas eu d'événement, soit pb = 0 , pc_pm = 0 ou dc = 0, leur temps d'événement est indiqué comme NA. </span>


```{r}
getug %>% filter(pb == 0) %>% select(numtas, pb, pbtime, pbtimed)
getug %>% filter(pc_pm == 0) %>% select(numtas, pc_pm, pc_pmtime, pc_pmtimed)
getug %>% filter(dc == 0) %>% select(numtas, dc, dctime, dctimed)

# Create interactive HTML tables
DT::datatable(getug)

```

### Traitement des valeurs manquantes

```{r}
getug$pbtimed <- ifelse(is.na(getug$pbtimed), getug$futimed, getug$pbtimed)
getug$pbtime  <- ifelse(is.na(getug$pbtime), getug$futime, getug$pbtime)

getug$pc_pmtimed <- ifelse(is.na(getug$pc_pmtimed), getug$futimed, getug$pc_pmtimed)
getug$pc_pmtime  <- ifelse(is.na(getug$pc_pmtime), getug$futime, getug$pc_pmtime)

getug$dctimed <- ifelse(is.na(getug$dctimed), getug$futimed, getug$dctimed)
getug$dctime <- ifelse(is.na(getug$dctime), getug$futime, getug$dctime)

# sum(is.na(getug[,c("pbtimed", "pbtime", "pc_pmtimed", "pc_pmtime", "dctimed", "dctime" )]))
# View(getug)

```

### Verification de la nouvelle base

```{r}
# Figure du modèle:

# Les patients qui ont eu une progression biologique pb en premier lieu
unique_numtas <- getug$numtas[getug$pb==1 & (getug$pc_pm==0 | getug$pbtimed < getug$pc_pmtimed)]
print(length(unique_numtas))

# Les patients qui sont décédés sans faire aucun evènement
unique_numtas_3 <- getug$numtas[getug$dc==1 & getug$pb==0 & getug$pc_pm==0]
print(length(unique_numtas_3))

# Les patients qui n'ont pas fait aucun evènement
unique_numtas_4 <- getug$numtas[getug$dc==0 & getug$pb==0 & getug$pc_pm==0]
print(length(unique_numtas_4))


# Les patients qui ont fait une progression biologique puis une progression à distance pb -> pc_pm
# View(getug[which(getug$pb==1 & getug$pc_pm==1 & getug$pbtimed < getug$pc_pmtimed), c("numtas", "pb", "pm", "dc","pbtimed", "pmtimed", "dctimed")])


# Les patients qui ont fait une progression local puis décédés pc_pm -> dc
getug[which(getug$pc_pm==1 & getug$dc==1 & (getug$pb == 0 | (getug$pc_pmtimed > getug$pbtimed))),  
              c("numtas", "pb", "pm", "pc", "dc", "pbtimed","pctimed","pmtimed", "dctimed")]


# Nombre total des patients qui son décédés
length(getug$numtas[which(getug$dc==1)])

# Nombre de patients censurés après une progression biologique pb -> censuré
getug[which(getug$pb==1 & getug$dc==0 & (getug$pc_pm==0 | getug$pc_pmtime < getug$pbtime)),]
getug$numtas[which(getug$pb==1 & getug$dc==0 & (getug$pc_pm==0 | getug$pc_pmtime < getug$pbtime))]

# Nombre de patients censurés après une progression local pc_pm-> censuré
# View(getug[which(getug$pc_pm==1 & getug$dc==0 & (getug$pb==0 | getug$pbtime < getug$pc_pmtime)),])

```

#### Kaplan Meier 
#####Survie globale

```{r}
fit_tvim_os <- survfit(Surv(futime, dc) ~ 1, data = getug)
 
fit_tvim_os
summary_fit <- summary(fit_tvim_os)


ggsurvplot(fit_tvim_os,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk http://127.0.0.1:47319/graphics/plot_zoom_png?width=1920&height=1017table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

# Calculer l'estimateur de Kaplan-Meier pour chaque groupe
surv_obj <- Surv(getug$futime, getug$dc)
km_fit_group <- survfit(surv_obj ~ trt_f, data = getug)

# Créer le graphique log-negative-log
ggsurvplot(km_fit_group, data = getug, 
           fun = "cloglog",
           title = "Vérification graphique de l'hypothèse des risques proportionnels",
           xlab = "log(Temps)", 
           ylab = "log(-log(Survie))")

```

##### Survie global par bras de trt

```{r}
fit_tvim_os <- survfit(Surv(futime, dc) ~ trt_f, data = getug)
 
fit_tvim_os
 
ggsurvplot(fit_tvim_os,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk http://127.0.0.1:47319/graphics/plot_zoom_png?width=1920&height=1017table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"),
           legend.title = "Treatments",
           legend.labs = c("TAA ", "TAA + DE "))


surv_diff <- survdiff(Surv(futime, dc) ~ trt_f, data = getug)

surv_diff

```

#####Survie sans progression

```{r}
fit_getug_pss <- survfit(Surv(time_to_prog, evt) ~ 1, data = getug)
 
fit_getug_pss
 
ggsurvplot(fit_getug_pss,
           pval = FALSE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk http://127.0.0.1:47319/graphics/plot_zoom_png?width=1920&height=1017table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           ylab = "Relapse-free survival %",
           palette = c("darkorange2"))

```

##### Survie sans progression par bras de traitement

```{r}
fit_getug_pss <- survfit(Surv(time_to_prog, evt) ~ trt_f, data = getug)
 
fit_getug_pss
 
ggsurvplot(fit_getug_pss,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk http://127.0.0.1:47319/graphics/plot_zoom_png?width=1920&height=1017table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("darkorange2", "darkorchid4"),
           legend.title = "Treatments",
           legend.labs = c("TAA ", "TAA + DE "))

surv_diff <- survdiff(Surv(time_to_prog, evt) ~ trt_f, data = getug)
surv_diff

```

#### Modèles de Cox
#####Survie globale

```{r}
mfit0 <- coxph(Surv(futime, dc) ~ trt + age + glearec + trec + pnrec + psarec, data = getug)
mfit0

```


#####Survie sans progression

```{r}
mfit0_prog <- coxph(Surv(time_to_prog, evt) ~ trt + age + glearec + trec + pnrec + psarec, data = getug)
mfit0_prog
```

### Matrice de transitions

```{r}

tmat <- transMat(x = list(c(2, 4), c(3), c(4), c()), names = c("RX", "PB", "PC_PM", "DC"))
tmat

paths(tmat)

```

### Format long

```{r}
covs <- c("numtas","trt","age_c","glearec","trec","psarec","pnrec", "nb_crit_cl")

msgetug <-  msprep(time = c(NA,"pbtime","pc_pmtime","dctime"), status = c(NA, "pb", "pc_pm", "dc"), data = getug, trans = tmat, keep=covs)

head(msgetug)
# View(msgetug)
# events(msgetug)
 
```

```{r}
# censored = subset(msgetug, !id %in% unique(msgetug$id[msgetug$status != 0]))
# 
# unique_patients <- length(unique(censored$numtas))
# unique_patients
# 
# unique_numtas <- unique(censored$numtas)
# print(unique_numtas)

```

## Modelisation
### Modèle non paramétrique
#### Modèle 0 : Modèle stratifié sans covariables 

```{r}
c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msgetug, method = "breslow")
c0
summary(c0)

# strata(trans) : permet au risque de base de varier entre les différents niveaux de trans
```


```{r}
# transition-specific covariates
# msgetug$trt <- as.character(msgetug$trt)
msgetug <- expand.covs(msgetug, covs = covs[!covs %in% c("numtas")])
head(msgetug)  
View(msgetug)

```

### Modèle semi paramétrique
#### Modèle 1 : Modèle stratifié "Markov stratified hazards model"

```{r}
c1 <- coxph(Surv(Tstart, Tstop, status) ~  
               trt.1 + age_c.1  + trec.1  + glearec.1  + pnrec.1  + psarec.1 +
               trt.2 + age_c.2  + trec.2  + glearec.2  + pnrec.2  + psarec.2 +
               trt.3 + age_c.3  + trec.3  + glearec.3  + pnrec.3  + psarec.3 +
               trt.4 + age_c.4  + trec.4  + glearec.4  + pnrec.4  + psarec.4 +
               strata(trans),
               data = msgetug, method = "breslow")
c1
summary(c1)
cox.zph(c1)

```
 ==> La valeur p globale (GLOBAL) est de 0.465 : le modèle ne montre pas de violation significative de l'hypothèse de proportionnalité des risques. Nous pouvons supposer les risques proportionnel.

#### Modèle 2 : Modèle à intensités proportionnelles "Markov PH model" 

```{r}
# Le modèle suivant considéré est le modèle de Markov dans lequel les risques de transition vers le décès (qui correspondent aux transitions 3, 5 et 6) sont supposés être proportionnels; ces transitions ont le même état de réception, donc la même valeur de to.
# Afin de distinguer les transitions 2 et 4, nous introduisons deux covariables temporelles ilns.1 qui indiquent si la PB ou PC+PM ont déjà eu lieu ou non.

msgetug$ilns.1 <- 0

msgetug$ilns.1[msgetug$trans == 4] <- 1


# On suppose que les risques de bases des transitions 3, 5 et 6 sont proportionnelles

c2 <- coxph(Surv(Tstart, Tstop, status) ~  
               trt.1 + age_c.1  + trec.1  + glearec.1  + pnrec.1  + psarec.1 +
               trt.2 + age_c.2  + trec.2  + glearec.2  + pnrec.2  + psarec.2 +
               trt.3 + age_c.3  + trec.3  + glearec.3  + pnrec.3  + psarec.3 + 
               trt.4 + age_c.4  + trec.4  + glearec.4  + pnrec.4  + psarec.4 + ilns.1 + 
               strata(to),
               data = msgetug, method = "breslow")
c2

```

- Les coefficients de régression estimés dans le modèle 2 pour la **transition 1** sont identiques à ceux estimés par le modèle 1


## Test de l'hypothèse de la proportionnalité des risques

```{r}
# Tester l'hypothèse des risques proportionnels 
test.ph <- cox.zph(c2)
test.ph

pvalg <- test.ph$table[rownames(test.ph$table)=="GLOBAL",3]
#pvalg
pvali.1 <- test.ph$table[rownames(test.ph$table)=="ilns.1",3]
#pvali.1

```

==> La valeur p globale (GLOBAL) est de 0.00085, inférieure à 0.05, ce qui indique une violation globale de l'hypothèse des risques proportionnels pour le modèle.

- ilns.1 (p = 2.2e-07) : La p-value indique que l'hypothèse de proportionnalité des risques n'est pas violée de manière significative pour ilns.1. Cela signifie que l'effet de ilns.1 sur le risque de décès est relativement constant au fil du temps. Ainsi, même si ilns.1 est une covariable temporelle, son influence sur la transition vers le décès ne varie pas de manière significative dans le temps.


```{r}
#ggcoxzph(test.ph)
# plot(test.ph)
# str(test.ph)

# y : Matrice des résidus de Schoenfeld pour chaque covariable à chaque temps d'événement.
# plot(test.ph$time, test.ph$y[,2])

```

### Vérification graphique de la validité de l'hypothèse de proportionnalité des risques de base pour les transitions 2 et 4 

```{r}
newd <- data.frame(trt = rep(0, 4), age_c = rep(0, 4), trec = rep(0, 4), glearec = rep(0, 4), pnrec = rep(0, 4), psarec = rep(0, 4),trans = 1:4)

# newd$trt <- factor(newd$trt, levels = 0:1, labels = levels(getug$trt_f))
# newd$age_c <- factor(newd$age_c, levels = 0:1, labels = levels(getug$age_cf))
# newd$trec <- factor(newd$trec, levels = 0:1, labels = levels(getug$trec_f))
# newd$glearec <- factor(newd$glearec, levels = 0:1, labels = levels(getug$glearec_f))
# newd$pnrec <- factor(newd$pnrec, levels = 0:1, labels = levels(getug$pnrec_f))
# newd$psarecpnrec <- factor(newd$psarec, levels = 0:1, labels = levels(getug$psarec_f))

attr(newd, "trans") <- tmat
class(newd) <- c("msdata", "data.frame")
newd <- expand.covs(newd, covs[2:7], longnames = FALSE)
# strata : specifying to which stratum in the coxph object each transition belongs. Here each transition corresponds to a separate stratum, so we specify 1:5
newd$strata = 1:4
newd

msf1 <- msfit(c1, newdata = newd, trans = tmat)
summary(msf1)

# Haz est l'estimation du risque de base calculée pour chaque moment spécifique (time) et pour chaque transition (trans)

# Pour le modèle 2 : 
newd2 <- newd

# strata spécifie à quelle stratum de l'objet coxph appartient chaque transition
newd2$strata = c(1, 3, 2, 3)

# ilns.1 needs to be included, taking the value 0 for transitions  2 and 5, and 1 for transition 4
# ilns.2 needs to be included, taking the value 0 for transitions  2 and 4, and 1 for transition 5

newd2$ilns.1 <- c(0, 0, 0, 1)

# Calcul des risques de tranbsition ou  risques cumulés
msf2 <- msfit(c2, newdata = newd2, trans = tmat)  
# summary(msf2)
# Estimations des risques (hazards) de transition. (time : le temps auquel le risque est estimé. Haz : l'estimation du risque (hazard) à ce temps)

```

```{r, fig.width= 8}
# Baseline cumulative hazard curves 
par(mfrow = c(1, 2))
plot(msf1, cols = 1:4, lwd = 2, lty = 1:4 , xlab = "Temps en année", ylab = "Risques de base stratifiés",
     legend.pos = "topleft", main = "Modèle 1")

plot(msf2, cols = 1:4, lwd = 2, lty = 1:4, xlab = "Temps en année", ylab = "Risques de base cumulatifs",
     legend.pos = "topleft", main = "Modèle 2")


```

## Estimation non paramétrique 
## Modèle non paramétrique ****
#### Intensités cumulées de transition 

```{r, fig.width= 8}

msf0 <- msfit(object=c0, vartype = "aalen", trans = tmat)
summary(msf0)
plot0 <- plot(msf0, las=1, lwd=1, xlab="Time", ylab="Cumulative Hazard of prostate cancer diagnosis", use.ggplot = T)
plot0 +  theme_bw()

# Graphes séparés pour chaque transition
plot(
  msf0, 
  type = "separate", 
  use.ggplot = TRUE,
  conf.int = 0.95,
  scale_type = "fixed"
)

```

#### Intensités cumulées de transition en fonction du traitement

```{r, fig.width= 8}

msgetug_0 <- msgetug[msgetug$trt == 0,]
msgetug_1 <- msgetug[msgetug$trt == 1,]

c0_0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msgetug_0, method = "breslow")
c0_1 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msgetug_1, method = "breslow")

msf0_0 <- msfit(object = c0_0, vartype = "aalen", trans = tmat)
# summary(msf0_0)
msf0_1 <- msfit(object = c0_1, vartype = "aalen", trans = tmat)
# summary(msf0_1)

# Intensités cumulées de transition trt 0
plot0_0 <- plot(msf0_0, las=1, lwd=0.5, xlab="Follow up time (years)", ylab="Cumulative Hazard TAA", use.ggplot = T)
plot0_0 +  theme_bw() 

# Intensités cumulées de transition trt 1
plot0_1 <- plot(msf0_1, las=1, lwd=0.5, xlab="Follow up time (years)", ylab="Cumulative Hazard TAA+DE", use.ggplot = T)
plot0_1 +  theme_bw() 


# Plots des risques cumulés en fonction du traitement sur le meme graphe (pour chaque transition et chaque traitement) 

df <- data.frame(temps = msf0_0[[1]]$time, 
                 intensites = msf0_0[[1]]$Haz,
                 trans = msf0_0[[1]]$trans,
                 Group = "TAA") 

df$trans<- factor( df$trans, 1:4, c("RX -> PB","RX -> DC", "PB -> PC_PM", "PC_PM -> DC"))


df2 <- data.frame(temps = msf0_1[[1]]$time, 
                 intensites = msf0_1[[1]]$Haz,
                 trans = msf0_1[[1]]$trans,
                 Group = "TAA + DE")

df2$trans<- factor( df2$trans, 1:4, c("RX -> PB","RX -> DC", "PB -> PC_PM", "PC_PM -> DC"))

D <- rbind(df, df2)


# Définir les couleurs personnalisées 
couleurs_personnalisees <- c( "RX -> PB" = "blueviolet", 
                              "RX -> DC" = "#ff7f0e", 
                              "PB -> PC_PM" = "#2ca02c", 
                              "PC_PM -> DC" = "darkred") 



# Plot des intensités cumulées en fonction du traitement 
ggplot(D, aes(x = temps, y = intensites, color = trans, linetype = Group)) +
  geom_line(linewidth = 0.5) +
  scale_linetype_manual(values = c("TAA" = "dashed", "TAA + DE" = "solid")) +
  scale_x_continuous(breaks = seq(1, max(D$temps), by = 1)) +
  scale_color_manual(values = couleurs_personnalisees) +
  # scale_x_continuous(breaks = seq(0, 17, by = 1), limits = x_limits) +
  labs(x = "Follow-up time (years)", y = "Cumulative Hazard ") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1))

# Intensite de transition en fonction du traitement plot pour chaque transition
ggplot(D, aes(x = temps, y = intensites, color = Group)) +
  geom_line() +
  facet_wrap(~trans, scales = "free_y") + # Un graphe pour chaque transition
  labs(title = "Intensités de transition par traitement",
       x = "Temps",
       y = "Intensité de transition",
       color = "Groupe de traitement") +
  theme_minimal()

```


#### Intensités cumulées de transition en fonction de l'age

```{r, fig.width= 8}

# Séparation des données selon l'âge
msgetug_c0 <- msgetug[msgetug$age_c == 0,]
msgetug_c1 <- msgetug[msgetug$age_c== 1,]

# Modèles de Cox pour chaque groupe d'âge
c0_c0 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msgetug_c0, method = "breslow")
c0_c1 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msgetug_c1, method = "breslow")

# calcule les hasards cumulatifs pour chacune des transitions 
msf0_c0 <- msfit(object = c0_c0, vartype = "aalen", trans = tmat)
# summary(msf0_c0)
msf0_c1 <- msfit(object = c0_c1, vartype = "aalen", trans = tmat)
# summary(msf0_c1)

# Création des graphiques des hasards cumulatifs 
plot0_c0 <- plot(msf0_c0, las=1, lwd=0.5, xlab="Time", ylab="Cumulative Hazard < 65", use.ggplot = T)
plot0_c0 +  theme_bw()

plot0_c1 <- plot(msf0_c1, las=1, lwd=0.5, xlab="Time", ylab="Cumulative Hazard > 65", use.ggplot = T)
plot0_c1 +  theme_bw() 


# Plot des risques cumulés pour chaque transition et chaque traitement 

df_c1 <- data.frame(temps = msf0_c0[[1]]$time, 
                 intensites = msf0_c0[[1]]$Haz,
                 trans = msf0_c0[[1]]$trans,
                 groupe = "<65") 

df_c1$trans<- factor( df_c1$trans, 1:4, c("RX -> PB (<65)","RX -> DC (<65)", "PB -> PC_PM (<65)", "PC_PM -> DC (<65)"))


df_c2 <- data.frame(temps = msf0_c1[[1]]$time, 
                 intensites = msf0_c1[[1]]$Haz,
                 trans = msf0_c1[[1]]$trans,
                 groupe = ">=65")

df_c2$trans<- factor( df_c2$trans, 1:4, c("RX -> PB (>=65)","RX -> DC (>=65)", 
                                      "PB -> PC_PM (>=65)", "PC_PM -> DC (>=65)"))

D_age <- rbind(df_c1, df_c2)


# Définir les couleurs personnalisées 
couleurs_personnalisees <- c( "RX -> PB (<65)" = "blueviolet", 
                              "RX -> DC (<65)" = "#ff7f0e", 
                              "PB -> PC_PM (<65)" = "#2ca02c", 
                              "PC_PM -> DC (<65)" = "darkred", 
                              "RX -> PB (>=65)" = "blueviolet", 
                              "RX -> DC (>=65)" = "#ff7f0e",
                              "PB -> PC_PM (>=65)" = "#2ca02c", 
                              "PC_PM -> DC (>=65)" = "darkred" ) 


# Plot des intensités cumulées en fonction du traitement 
ggplot(D_age, aes(x = temps, y = intensites, color = trans, linetype = groupe)) + 
  geom_line(linewidth = 0.5) + 
  scale_linetype_manual(values = c("<65" = "dashed", ">=65" = "solid")) +
  scale_color_manual(values = couleurs_personnalisees) + 
  labs(x = "Time", y = "Cumulative Hazard ") +
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1))

```


#### Intensités cumulées de transition en fonction du "score de risque"


```{r}
msgetug_r1 <- msgetug[msgetug$nb_crit_cl == 1,]
msgetug_r2 <- msgetug[msgetug$nb_crit_cl == 2,]
msgetug_r3 <- msgetug[msgetug$nb_crit_cl == 3,]

c0_1 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msgetug_r1, method = "breslow")
c0_2 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msgetug_r2, method = "breslow")
c0_3 <- coxph(Surv(Tstart, Tstop, status) ~ strata(trans), data = msgetug_r3, method = "breslow")

msf0_1 <- msfit(object = c0_1, vartype = "aalen", trans = tmat)
# summary(msf0_1)
msf0_2 <- msfit(object = c0_2, vartype = "aalen", trans = tmat)
# summary(msf0_2)
msf0_3 <- msfit(object = c0_3, vartype = "aalen", trans = tmat)
# summary(msf0_2)

plot0_1 <- plot(msf0_1, las=1, lwd=0.5, xlab="Time", ylab="Cumulative Hazard ", use.ggplot = T)
plot0_1 +  theme_bw() 
plot0_2 <- plot(msf0_2, las=1, lwd=0.5, xlab="Time", ylab="Cumulative Hazard ", use.ggplot = T)
plot0_2 +  theme_bw() 
plot0_3 <- plot(msf0_3, las=1, lwd=0.5, xlab="Time", ylab="Cumulative Hazard ", use.ggplot = T)
plot0_3 +  theme_bw() 
# Ajouter des titres aux plots
plot0_1 <- plot0_1 + ggtitle("risk score = 1")
plot0_2 <- plot0_2 + ggtitle("risk score = 2")
plot0_3 <- plot0_3 + ggtitle("risk score >= 3")

plot_grid(plot0_1, plot0_2, plot0_3, ncol = 2, nrow = 2)
'''
