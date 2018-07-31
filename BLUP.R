library(lattice)
library(lme4)
library(pedigree)
library(MCMCglmm)
library(ibdreg)
library(pastecs)
require(ggplot2)
require(reshape2)
require(gridExtra)


###wczytanie i przekszta³cenie zbioru danych
dane <- read.table(file = "C:/Users/Admin/Desktop/danex.csv", 
                   sep = ";", header = TRUE)
dane2 = dane[!(dane$NRMATKI %in% dane$NROJCA),]
dane_final = dane2[!(dane2$NROJCA %in% dane2$NRMATKI),]
dane_final = dane_final[!duplicated(dane_final), ]
attach(dane_final)
#podstawowe statystyki opisowe
stat.desc(dane_final$KGMLE) 


###wstepna analiza graficzna 
#podzial RAHF i ROHF na kategorie
kat_RAHF <- cut(RAHF,breaks=4)
levels(kat_RAHF) <- c('0-23.5%', '23.5-47%', '47-70.5%', '70.5%+')
kat_ROHF <- cut(dane_final$ROHF,breaks=4)
levels(kat_ROHF) <- c('0-25%', '25-50%', '50-75%', '75%+')
licz = table(kat_RAHF)
licz2 = table(kat_ROHF)
licz_df = data.frame(RAHF = rownames(licz),
                     liczebnoœæ = as.vector(licz))
licz_df2 = data.frame(ROHF = rownames(licz2),
                     liczebnoœæ = as.vector(licz2))
#wykres - liczebnosci podgrup RAHF i ROHF
p_RAHF <- ggplot(data=licz_df, aes(x=RAHF, y=liczebnoœæ, ymax=500)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=liczebnoœæ), vjust=-1, color="black", size=3.5) +
  ggtitle("A") +
  theme_minimal()
p_ROHF <- ggplot(data=licz_df2, aes(x=ROHF, y=liczebnoœæ, ymax=800)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=liczebnoœæ), vjust=-1, color="black", size=3.5) +
  ggtitle("B") +
  theme_minimal()
#wykres - srednie cech prodykcyjnych w latach
ustaw_tytul <- list(par.main.text = list(font = 2,
                                         just = "left", 
                                         x = grid::unit(20, "mm")))
œrednie = aggregate(dane_final[,11:15], list(dane_final$ROURO), mean)
œr_RAHF = aggregate(dane_final[,4], list(dane_final$ROURO), mean)
œrednie <- cbind(œrednie, œr_RAHF[,2])
names(œrednie)[7] <- "RAHF"
names(œrednie)[1] <- "ROURO"
p1 <- xyplot(RAHF ~ ROURO, data = œrednie, pch = 19, type = "o",
             scales=list(y=list(at=seq(0.4,0.8,0.1), limits=c(0.4,0.8))),
             par.settings = ustaw_tytul,
             main = "A")
p2 <- xyplot(KGMLE ~ ROURO, data = œrednie, pch = 19, type = "o",
        xlab = "ROURO",
        ylab = "KGMLE",
        scales=list(y=list(at=seq(2500,6000,500), limits=c(2500,6000))),
        par.settings = ustaw_tytul,
        main = "B")
p3 <- xyplot(KGT£U + KGBIA£ ~ ROURO, data = œrednie, pch = 19, type = "o",
        xlab = "ROURO",
        ylab = "WARTOŒÆ",
        scales=list(y=list(at=seq(50,250,50), limits=c(50,250))),
        par.settings = ustaw_tytul,
        main = "C",
        auto.key = TRUE)
p4 <- xyplot(PRT£U + PRBIA£ ~ ROURO, data = œrednie, pch = 19, type = "o",
       xlab = "ROURO",
       ylab = "WARTOŒÆ",
       scales=list(y=list(at=seq(2.5,5,0.3), limits=c(2.5,4.6))),
       par.settings = ustaw_tytul,
       main = "D",
       auto.key = TRUE)
grid.arrange(p1,p2,p3,p4, ncol = 2)
grid.arrange(p_RAHF,p_ROHF, ncol = 2)
#wykres - srednie cech w kategoriach RAHF
danex <- data.frame(dane_final[,11:15])
danex <- cbind(danex, kat_RAHF)
layout(matrix(c(1,2,3,4,4,4), ncol = 3, byrow=TRUE), heights=c(4, 1))
par(mai=rep(0.5, 4))
boxplot(KGMLE~kat_RAHF, data = danex, ylab = "KGMLE", main = "A")
boxplot(KGT£U~kat_RAHF, data = danex, ylab = "KGT£U", main = "B")
boxplot(KGBIA£~kat_RAHF, data = danex, ylab = "KGBIA£", main = "C")
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=4,legend=c("0-23.5%","23.5-47%","47-70.5%","70.5%+"),
       fill=c("black","black","black", "black"), title="Kategorie RAHF")
layout(matrix(c(1,2,3,3), ncol = 2, byrow=TRUE), heights=c(6, 1))
par(mai=rep(0.8, 4))
boxplot(PRT£U~kat_RAHF, data = danex, ylab = "PRT£U", main = "D")
boxplot(PRBIA£~kat_RAHF, data = danex, ylab = "PRBIA£", main = "E")


###analiza modeli
#wybor uogolnionego liniowego modelu
dane_final$ROURO <- as.factor(dane_final$ROURO)
dane_final$ROKWY <- as.factor(dane_final$ROKWY)
dane_final$MIEWY <- as.factor(dane_final$MIEWY)
dane_final$MURO <- as.factor(dane_final$MURO)
m1 <- glm(KGMLE ~ ROURO + MURO + RAHF + RMHF + ROHF + ROKWY + MIEWY, data = dane_final)
step(m1,direction="backward")
m2 <- glm(KGT£U ~ ROURO + MURO + RAHF + RMHF + ROHF + ROKWY + MIEWY, data = dane_final)
step(m2,direction="backward")
m3 <- glm(PRT£U ~ ROURO + MURO + RAHF + RMHF + ROHF + ROKWY + MIEWY, data = dane_final)
step(m3,direction="backward")
m4 <- glm(KGBIA£ ~ ROURO + MURO + RAHF + RMHF + ROHF + ROKWY + MIEWY, data = dane_final)
step(m4,direction="backward")
m5 <- glm(PRBIA£ ~ ROURO + MURO + RAHF + RMHF + ROHF + ROKWY + MIEWY, data = dane_final)
step(m5,direction="backward")
#estymacja komponentow wariancyjnych modeli mieszanych metodami REML i ML
lmer_1 <- lmer(KGMLE ~ ROURO + RAHF + MIEWY + (1|NrKR), data = dane_final,
               control=lmerControl(check.nobs.vs.nlev="ignore",
                                   check.nobs.vs.rankZ="ignore",
                                   check.nobs.vs.nRE="ignore"))
summary(lmer_1)
lme_1 <- lmer(KGMLE ~ ROURO + RAHF + MIEWY + (1|NrKR), REML = FALSE, data = dane_final,
              control=lmerControl(check.nobs.vs.nlev="ignore",
                                  check.nobs.vs.rankZ="ignore",
                                  check.nobs.vs.nRE="ignore"))
summary(lme_1)
lmer_2 <- lmer(KGT£U ~ ROURO + RAHF + (1|NrKR), data = dane_final,
               control=lmerControl(check.nobs.vs.nlev="ignore",
                                   check.nobs.vs.rankZ="ignore",
                                   check.nobs.vs.nRE="ignore"))
summary(lmer_2)
lme_2 <- lmer(KGT£U ~ ROURO + RAHF + (1|NrKR), REML = FALSE, data = dane_final,
              control=lmerControl(check.nobs.vs.nlev="ignore",
                                  check.nobs.vs.rankZ="ignore",
                                  check.nobs.vs.nRE="ignore"))
summary(lme_2)
lmer_3 <- lmer(PRT£U ~ ROURO + RAHF + ROHF + ROKWY + MIEWY + (1|NrKR), data = dane_final,
               control=lmerControl(check.nobs.vs.nlev="ignore",
                                   check.nobs.vs.rankZ="ignore",
                                   check.nobs.vs.nRE="ignore"))
summary(lmer_3)
lme_3 <- lmer(PRT£U ~ ROURO + RAHF + ROHF + ROKWY + MIEWY + (1|NrKR), REML = FALSE, data = dane_final,
              control=lmerControl(check.nobs.vs.nlev="ignore",
                                  check.nobs.vs.rankZ="ignore",
                                  check.nobs.vs.nRE="ignore"))
summary(lme_3)
lmer_4 <- lmer(KGBIA£ ~ ROURO + RAHF + MIEWY + (1|NrKR), data = dane_final,
               control=lmerControl(check.nobs.vs.nlev="ignore",
                                   check.nobs.vs.rankZ="ignore",
                                   check.nobs.vs.nRE="ignore"))
summary(lmer_4)
lme_4 <- lmer(KGBIA£ ~ ROURO + RAHF + MIEWY + (1|NrKR), REML = FALSE, data = dane_final,
              control=lmerControl(check.nobs.vs.nlev="ignore",
                                  check.nobs.vs.rankZ="ignore",
                                  check.nobs.vs.nRE="ignore"))
summary(lme_4)
lmer_5 <- lmer(PRBIA£ ~ ROURO + ROHF + (1|NrKR), data = dane_final,
               control=lmerControl(check.nobs.vs.nlev="ignore",
                                   check.nobs.vs.rankZ="ignore",
                                   check.nobs.vs.nRE="ignore"))
summary(lmer_5)
lme_5 <- lmer(PRBIA£ ~ ROURO + ROHF + (1|NrKR), REML = FALSE, data = dane_final,
              control=lmerControl(check.nobs.vs.nlev="ignore",
                                  check.nobs.vs.rankZ="ignore",
                                  check.nobs.vs.nRE="ignore"))
summary(lme_5)


###test ilorazu wiarygodnosci dla efektu losowego 
#model1
IEL <- logLik(lmer_1 <- lmer(KGMLE ~ ROURO + RAHF + MIEWY + (1|NrKR), data = dane_final,
                             control=lmerControl(check.nobs.vs.nlev="ignore",
                                                 check.nobs.vs.rankZ="ignore",
                                                 check.nobs.vs.nRE="ignore")))
IELb <- logLik(lm_1 <- lm(KGMLE ~ ROURO + RAHF + MIEWY, data = dane_final))
roznica = as.numeric(IELb - IEL)
pchisq(-2*roznica, 1, lower.tail = FALSE) 
#model2
IEL <- logLik(lmer_2 <- lmer(KGT£U ~ ROURO + RAHF + (1|NrKR), data = dane_final,
                             control=lmerControl(check.nobs.vs.nlev="ignore",
                                                 check.nobs.vs.rankZ="ignore",
                                                 check.nobs.vs.nRE="ignore")))
IELb <- logLik(lm_2 <- lm(KGT£U ~ ROURO + RAHF, data = dane_final))
roznica = as.numeric(IELb - IEL)
pchisq(-2*roznica, 1, lower.tail = FALSE) 
#model3
IEL <- logLik(lmer_3 <- lmer(PRT£U ~ ROURO + RAHF + ROHF + ROKWY + MIEWY + (1|NrKR), data = dane_final,
                             control=lmerControl(check.nobs.vs.nlev="ignore",
                                                 check.nobs.vs.rankZ="ignore",
                                                 check.nobs.vs.nRE="ignore")))
IELb <- logLik(lm_3 <- lm(PRT£U ~ ROURO + RAHF + ROHF + ROKWY + MIEWY, data = dane_final))
roznica = as.numeric(IELb - IEL)
pchisq(2*roznica, 1, lower.tail = FALSE)
#model4
IEL <- logLik(lmer_4 <- lmer(KGBIA£ ~ ROURO + RAHF + MIEWY + (1|NrKR), data = dane_final,
                             control=lmerControl(check.nobs.vs.nlev="ignore",
                                                 check.nobs.vs.rankZ="ignore",
                                                 check.nobs.vs.nRE="ignore")))
IELb <- logLik(lm_4 <- lm(KGBIA£ ~ ROURO + RAHF + MIEWY, data = dane_final))
roznica = as.numeric(IELb - IEL)
pchisq(-2*roznica, 1, lower.tail = FALSE) 
#model5
IEL <- logLik(lmer_5 <- lmer(PRBIA£ ~ ROURO + ROHF + (1|NrKR), data = dane_final,
                             control=lmerControl(check.nobs.vs.nlev="ignore",
                                                 check.nobs.vs.rankZ="ignore",
                                                 check.nobs.vs.nRE="ignore")))
IELb <- logLik(lm_5 <- lm(PRBIA£ ~ ROURO + ROHF, data = dane_final))
roznica = as.numeric(IELb - IEL)
pchisq(2*roznica, 1, lower.tail = FALSE)


###diagnostyka graficzna modeli
#wykres kwantylowy dla reszt
par(mfrow=c(2,3))
e1 = residuals(lmer_1)
qqnorm(e1, pch = 16, main = "MODEL_1")
qqline(e1)
e2 = residuals(lmer_2)
qqnorm(e2, pch = 16, main = "MODEL_2")
qqline(e2)
e3 = residuals(lmer_3)
qqnorm(e3, pch = 16, main = "MODEL_3")
qqline(e3)
e4 = residuals(lmer_4)
qqnorm(e4, pch = 16, main = "MODEL_4")
qqline(e4)
e5 = residuals(lmer_5)
qqnorm(e5, pch = 16, main = "MODEL_5")
qqline(e5)
#zaleznosc miedzy predykcjami a ocenami reszt 
pe1 <- ggplot(lmer_1, aes(.fitted, .resid)) + geom_point() +
  stat_smooth(method="glm")+
  geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values")+ylab("Residuals") +
  ggtitle("MODEL_1") +
  theme_bw()
pe2 <- ggplot(lmer_2, aes(.fitted, .resid)) + geom_point() +
  stat_smooth(method="glm")+
  geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values")+ylab("Residuals") +
  ggtitle("MODEL_2") +
  theme_bw()
pe3 <- ggplot(lmer_3, aes(.fitted, .resid)) + geom_point() +
  stat_smooth(method="glm")+
  geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values")+ylab("Residuals") +
  ggtitle("MODEL_3") +
  theme_bw()
pe4 <- ggplot(lmer_4, aes(.fitted, .resid)) + geom_point() +
  stat_smooth(method="glm")+
  geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values")+ylab("Residuals") +
  ggtitle("MODEL_4") +
  theme_bw()
pe5 <- ggplot(lmer_5, aes(.fitted, .resid)) + geom_point() +
  stat_smooth(method="glm")+
  geom_hline(yintercept=0, col="red", linetype="dashed") +
  xlab("Fitted values")+ylab("Residuals") +
  ggtitle("MODEL_5") +
  theme_bw()
grid.arrange(pe1,pe2,pe3,pe4,pe5, ncol = 3)


###analiza BLUP
#tworzenie pedigree i odwrotnej macierzy spokrewnien
ID <- dane_final$NrKR
DAM <- dane_final$NRMATKI
SIRE <- dane_final$NROJCA
rodowod <- data.frame(ID,DAM,SIRE)
rodowod2 <- add.Inds(rodowod)
ord <- orderPed(rodowod2)
rodowod2 <- rodowod2[order(ord),]
Ainv = inverseA(rodowod2)$Ainv 
#przygotowanie ramki danych do tworzenia macierzy ukladu rownan modelu mieszanego 
kat_ROHF <- cut(dane_final$ROHF,breaks=4)
levels(kat_ROHF) <- c('0-25%', '25-50%', '50-75%', '75%+')
kat_RAHF <- cut(dane_final$RAHF,breaks=4)
levels(kat_RAHF) <- c('0-23.5%', '23.5-47%', '47-70.5%', '70.5%+')
dane_BLUP1 <- dane_final[1:3]
dane_BLUP1 <- cbind(dane_BLUP1, RAHF = kat_RAHF, ROHF = kat_ROHF, ROKWY = dane_final$ROKWY, MIEWY = dane_final$MIEWY)
dane_BLUP1$ROURO <- as.numeric(as.factor(dane_BLUP1$ROURO))
dane_BLUP1$ROKWY <- as.numeric(as.factor(dane_BLUP1$ROKWY))
dane_BLUP1$RAHF <- as.numeric(as.factor(dane_BLUP1$RAHF))
dane_BLUP1$ROHF <- as.numeric(as.factor(dane_BLUP1$ROHF))
dane_BLUP1$NrKR <- as.character(as.factor(dane_BLUP1$NrKR))
#deklarowanie zmiennych MME
Z1 = matrix(0, nrow = 895, ncol = 917)
Z2 = diag(x = 1, nrow = 895, ncol = 895)
Z = cbind(Z1,Z2)
k1 = 1.005562
k2 = 1
k3 = 1.563943
k4 = 0.9922083
k5 = 0.992674
y1 <- as.vector(dane_final$KGMLE)
y2 <- as.vector(dane_final$KGT£U)
y3 <- as.vector(dane_final$PRT£U)
y4 <- as.vector(dane_final$KGBIA£)
y5 <- as.vector(dane_final$PRBIA£)
#BLUP
#model 1
n1 = length(unique(dane_BLUP1$ROURO))
n2 = length(unique(dane_BLUP1$RAHF))
n3 = length(unique(dane_BLUP1$MIEWY))
nX = n1 + n2 + n3
X <- matrix(0, nrow = nrow(dane_final), ncol = nX+1)
X[,1] <- c(1)
for (i in 1:nrow(dane_final)){
  for (j in 1:(n1+1)){
    if (dane_BLUP1$ROURO[i] == j){X[i,j+1] <- 1
    }
  }
  for (j in 1:n2){
    if (dane_BLUP1$RAHF[i] == j){X[i,j+(n1+1)] <- 1
    }
  }
  for (j in 1:n3){
    if (dane_BLUP1$MIEWY[i] == j){X[i,j+(n1+n2+1)] <- 1
    }
  }
}
Xpy <- t(X)%*%y1
Zpy <- t(Z)%*%y1
XpX <- t(X)%*%X
XpZ <- t(X)%*%Z
ZpX <- t(Z)%*%X
ZpZ <- t(Z)%*%Z
LHS <- rbind(cbind(XpX, XpZ), cbind(ZpX, ZpZ+Ainv*k1))
RHS <- rbind(Xpy, Zpy)
LHSinv = Ginv(LHS)
b = LHSinv$Ginv%*%as.matrix(RHS)
ebv1 <- data.frame(NR = dane_final$NrKR,
                  ROURO = dane_final$ROURO,
                  RAHF = kat_RAHF,
                  MIEWY = dane_final$MIEWY,
                  ODCH = b[945:nrow(b),],
                  EBV = b[1]+b[945:nrow(b),])

L2 = LHSinv$Ginv[945:nrow(b), 945:nrow(b)]
dokl = sqrt(1-(diag(L2)*k1))
ebv1 <- cbind(ebv1, r = dokl)
ebv_sort1 <- ebv1[order(-ebv1$EBV),] 
head(ebv_sort1, n = 10)
#model 2
n1 = length(unique(dane_BLUP1$ROURO))
n2 = length(unique(dane_BLUP1$RAHF))
nX = n1 + n2
X <- matrix(0, nrow = nrow(dane_final), ncol = nX+1)
X[,1] <- c(1)
for (i in 1:nrow(dane_final)){
  for (j in 1:(n1+1)){
    if (dane_BLUP1$ROURO[i] == j){X[i,j+1] <- 1
    }
  }
  for (j in 1:n2){
    if (dane_BLUP1$RAHF[i] == j){X[i,j+(n1+1)] <- 1
    }
  }
}
Xpy <- t(X)%*%y2
Zpy <- t(Z)%*%y2
XpX <- t(X)%*%X
XpZ <- t(X)%*%Z
ZpX <- t(Z)%*%X
ZpZ <- t(Z)%*%Z
LHS <- rbind(cbind(XpX, XpZ), cbind(ZpX, ZpZ+Ainv*k2))
RHS <- rbind(Xpy, Zpy)
LHSinv = Ginv(LHS)
b = LHSinv$Ginv%*%as.matrix(RHS)
ebv2 <- data.frame(NR = dane_final$NrKR,
                  ROURO = dane_final$ROURO,
                  RAHF = kat_RAHF,
                  ODCH = b[933:nrow(b),],
                  EBV = b[1]+b[933:nrow(b),])
L2 = LHSinv$Ginv[933:nrow(b), 933:nrow(b)]
dokl2 = sqrt(1-(diag(L2)*k2))
ebv2 <- cbind(ebv2, r = dokl2)
ebv_sort2 <- ebv2[order(-ebv2$EBV),] 
head(ebv_sort2, n = 10)
#model 3
n1 = length(unique(dane_BLUP1$ROURO))
n2 = length(unique(dane_BLUP1$RAHF))
n3 = length(unique(dane_BLUP1$ROHF))
n4 = length(unique(dane_BLUP1$ROKWY))
n5 = length(unique(dane_BLUP1$MIEWY))
nX = n1 + n2 + n3 + n4 + n5
X <- matrix(0, nrow = nrow(dane_final), ncol = nX+1)
X[,1] <- c(1)
for (i in 1:nrow(dane_final)){
  for (j in 1:(n1+1)){
    if (dane_BLUP1$ROURO[i] == j){X[i,j+1] <- 1
    }
  }
  for (j in 1:n2){
    if (dane_BLUP1$RAHF[i] == j){X[i,j+(n1+1)] <- 1
    }
  }
  for (j in 1:n3){
    if (dane_BLUP1$ROHF[i] == j){X[i,j+(n1+n2+1)] <- 1
    }
  }
  for (j in 1:n4){
    if (dane_BLUP1$ROKWY[i] == j){X[i,j+(n1+n2+n3+1)] <- 1
    }
  }
  for (j in 1:n5){
    if (dane_BLUP1$MIEWY[i] == j){X[i,j+(n1+n2+n3+n4+1)] <- 1
    }
  }
}
Xpy <- t(X)%*%y3
Zpy <- t(Z)%*%y3
XpX <- t(X)%*%X
XpZ <- t(X)%*%Z
ZpX <- t(Z)%*%X
ZpZ <- t(Z)%*%Z
LHS <- rbind(cbind(XpX, XpZ), cbind(ZpX, ZpZ+Ainv*k3))
RHS <- rbind(Xpy, Zpy)
LHSinv = Ginv(LHS)
b = LHSinv$Ginv%*%as.matrix(RHS)
ebv3 <- data.frame(NR = dane_final$NrKR,
                  ROURO = dane_final$ROURO,
                  RAHF = kat_RAHF,
                  ODCH = b[958:nrow(b),],
                  EBV = b[1]+b[958:nrow(b),])
L2 = LHSinv$Ginv[958:nrow(b), 958:nrow(b)]
dokl3 = sqrt(1-(diag(L2)*k3))
ebv3 <- cbind(ebv3, r = dokl3)
ebv_sort3 <- ebv3[order(-ebv3$EBV),] 
head(ebv_sort3, n = 10)
#model 4
n1 = length(unique(dane_BLUP1$ROURO))
n2 = length(unique(dane_BLUP1$RAHF))
n3 = length(unique(dane_BLUP1$MIEWY))
nX = n1 + n2 + n3
X <- matrix(0, nrow = nrow(dane_final), ncol = nX+1)
X[,1] <- c(1)
for (i in 1:nrow(dane_final)){
  for (j in 1:(n1+1)){
    if (dane_BLUP1$ROURO[i] == j){X[i,j+1] <- 1
    }
  }
  for (j in 1:n2){
    if (dane_BLUP1$RAHF[i] == j){X[i,j+(n1+1)] <- 1
    }
  }
  for (j in 1:n3){
    if (dane_BLUP1$MIEWY[i] == j){X[i,j+(n1+n2+1)] <- 1
    }
  }
}
Xpy <- t(X)%*%y4
Zpy <- t(Z)%*%y4
XpX <- t(X)%*%X
XpZ <- t(X)%*%Z
ZpX <- t(Z)%*%X
ZpZ <- t(Z)%*%Z
LHS <- rbind(cbind(XpX, XpZ), cbind(ZpX, ZpZ+Ainv*k4))
RHS <- rbind(Xpy, Zpy)
LHSinv = Ginv(LHS)
b = LHSinv$Ginv%*%as.matrix(RHS)
ebv4 <- data.frame(NR = dane_final$NrKR,
                  ROURO = dane_final$ROURO,
                  RAHF = kat_RAHF,
                  ODCH = b[945:nrow(b),],
                  EBV = b[1]+b[945:nrow(b),])
L2 = LHSinv$Ginv[945:nrow(b), 945:nrow(b)]
dokl4 = sqrt(1-(diag(L2)*k4))
ebv4 <- cbind(ebv4, r = dokl4)
ebv_sort4 <- ebv4[order(-ebv4$EBV),] 
head(ebv_sort4, n = 10)
#model 5
n1 = length(unique(dane_BLUP1$ROURO))
n2 = length(unique(dane_BLUP1$ROHF))
nX = n1 + n2
X <- matrix(0, nrow = nrow(dane_final), ncol = nX+1)
X[,1] <- c(1)
for (i in 1:nrow(dane_final)){
  for (j in 1:(n1+1)){
    if (dane_BLUP1$ROURO[i] == j){X[i,j+1] <- 1
    }
  }
  for (j in 1:n2){
    if (dane_BLUP1$ROHF[i] == j){X[i,j+(n1+1)] <- 1
    }
  }
}
Xpy <- t(X)%*%y5
Zpy <- t(Z)%*%y5
XpX <- t(X)%*%X
XpZ <- t(X)%*%Z
ZpX <- t(Z)%*%X
ZpZ <- t(Z)%*%Z
LHS <- rbind(cbind(XpX, XpZ), cbind(ZpX, ZpZ+Ainv*k5))
RHS <- rbind(Xpy, Zpy)
LHSinv = Ginv(LHS)
b = LHSinv$Ginv%*%as.matrix(RHS)
ebv5 <- data.frame(NR = dane_final$NrKR,
                  ROURO = dane_final$ROURO,
                  RAHF = kat_RAHF,
                  ODCH = b[933:nrow(b),],
                  EBV = b[1]+b[933:nrow(b),])
L2 = LHSinv$Ginv[933:nrow(b), 933:nrow(b)]
dokl5 = sqrt(1-(diag(L2)*k5))
ebv5 <- cbind(ebv5, r = dokl5)
ebv_sort5 <- ebv5[order(-ebv5$EBV),] 
head(ebv_sort5, n = 10)

#wykres - wartosc hodowlana w latach + trend genetyczny
œr_BLUP = aggregate(ebv1$EBV, list(ebv1$ROURO), mean)
b1 <- xyplot(x ~ Group.1, data = œr_BLUP, pch = 19,
             panel=function(x, y, col, ...) {
               panel.xyplot(x, y, col="black",...)
               panel.abline(lm(y~x), col='blue')},
             type='o',
             scales=list(y=list(at=seq(2500,3500,200), limits=c(2500,3500))), 
             xlab = "ROURO",
             ylab = "EBV",
             main = "KGMLE")
œr_BLUP2 = aggregate(ebv2$EBV, list(ebv2$ROURO), mean)
b2 <- xyplot(x ~ Group.1, data = œr_BLUP2, pch = 19,
             panel=function(x, y, col, ...) {
               panel.xyplot(x, y, col="black",...)
               panel.abline(lm(y~x), col='blue')},
             type='o',
             scales=list(y=list(at=seq(100,150,10), limits=c(100,150))),
             xlab = "ROURO",
             ylab = "EBV",
             main = "KGT£U")
œr_BLUP3 = aggregate(ebv3$EBV, list(ebv3$ROURO), mean)
b3 <- xyplot(x ~ Group.1, data = œr_BLUP3, pch = 19,
             panel=function(x, y, col, ...) {
               panel.xyplot(x, y, col="black",...)
               panel.abline(lm(y~x), col='blue')},
             type='o',
             scales=list(y=list(at=seq(2.2,2.35,0.03), limits=c(2.2,2.35))),
             xlab = "ROURO",
             ylab = "EBV",
             main = "PRT£U")
œr_BLUP4 = aggregate(ebv4$EBV, list(ebv4$ROURO), mean)
b4 <- xyplot(x ~ Group.1, data = œr_BLUP4, pch = 19,
             panel=function(x, y, col, ...) {
               panel.xyplot(x, y, col="black",...)
               panel.abline(lm(y~x), col='blue')},
             type='o',
             scales=list(y=list(at=seq(75,105,5), limits=c(75,105))),
             xlab = "ROURO",
             ylab = "EBV",
             main = "KGBIA£")
œr_BLUP5 = aggregate(ebv5$EBV, list(ebv5$ROURO), mean)
b5 <- xyplot(x ~ Group.1, data = œr_BLUP5, pch = 19,
             panel=function(x, y, col, ...) {
               panel.xyplot(x, y, col="black",...)
               panel.abline(lm(y~x), col='blue')},
             type='o',
             scales=list(y=list(at=seq(2.20,2.4,0.04), limits=c(2.20,2.4))),
             xlab = "ROURO",
             ylab = "EBV",
             main = "PRBIA£")
grid.arrange(b1,b2,b3,b4,b5, ncol = 3)

