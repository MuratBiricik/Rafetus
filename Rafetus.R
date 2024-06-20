
##   Rafetus euphraticus
##   Sand habitat comparison 2018/2022, Tigris TR

#-------------------------------------------------------------------------------

# Sandbank data 2018

setwd("C:/.../Rafetus")

sand_alldata <- read.csv("SandAll.csv")
table(sand_alldata["sand"])

#-------------------------------------------------------------------------------

# Explanatory variables 2018

library(raster)

# Sentinel-2 tiles (sensed 2018) covering entire area clipped by Tigris 2018

{
  Bd02 <- raster('./rs18/Bd02.tif')
  Bd07 <- raster('./rs18/Bd07.tif')
  Bd08 <- raster("./rs18/Bd08.tif")
  Bd09 <- raster("./rs18/Bd09.tif")
  Bd12 <- raster("./rs18/Bd12.tif")
}

# create a raster stack
rs18 <- stack(Bd02,	Bd07,	Bd08,	Bd09,	Bd12)

# set coordinate reference system
crs(rs18) <- "+init=epsg:4326 +WGS84"

rm(Bd02, Bd07,	Bd08,	Bd09,	Bd12)

#-------------------------------------------------------------------------------

# Modeling

library(sdm)
library(dplyr)
library(usdm)

# create "sdmData" required by the sdm package
d <- sdmData(sand ~ Bd02 + Bd07 +	Bd08 +	Bd09 + Bd12 +	
                    coords(x + y), train = sand_alldata)
d

# model building

modelSDM <- sdm(sand ~ Bd02 + Bd07 +	Bd08 +	Bd09 + Bd12 +	coords(x + y), d, 
                methods         = c("GAM", "GLM", "GBM", "RF"), 
                replication     = c("boot", "subs"),
                test.p          = 30,
                n               = 20,
                parallelSetting = list(ncore = 4, method = "parallel"))
modelSDM

rm(d)

#-------------------------------------------------------------------------------

# Evaluation

library(shiny)
gui(modelSDM)

# to see details (e.g.)
# modelSDM@models$sand$GAM$`2`

# average variable importance for all models
vi <- getVarImp(modelSDM)
vi

# plot variable importance (using R standard plot instead of ggplot)
par(mfrow = c(1, 1))
plot(vi, gg = F)

rm(vi)

#-------------------------------

# plot ROC curves
roc(modelSDM, 
    method      = c("glm", "gam"),
    replication = NULL,
    run         = NULL,
    wtest       = NULL,
    smooth      = T,
    legend      = T)

# response curve
rcurve(modelSDM, gg = T, main = "...")

getEvaluation(modelSDM, stat=c('TSS','Kappa','AUC'),opt=2)
getEvaluation(modelSDM,w=1,wtest='training')
getEvaluation(modelSDM,w=200,wtest='training',stat=1) # stat=1 (threshold_independent)
ev <- getEvaluation(modelSDM,w=5,wtest=NULL,stat=2) # stat=2 (threshold_based)
getEvaluation(modelSDM,w=1:3,wtest='test.dep',stat=c('AUC','TSS'),opt=2)
getEvaluation(modelSDM,opt=1) # all models
getEvaluation(modelSDM,stat=c('TSS','Kappa','AUC'),opt=1) # all models

#-------------------------------------------------------------------------------

# Ensemble model 2018

ens18 <- ensemble(modelSDM,
                  newdata  = rs18,
                  filename = "ens18_01.tif",
                  setting  = list(method = "weighted", stat = "TSS"))

# plot(ens18)

#-------------------------------------------------------------------------------
 
# Ensemble model 2018 evaluation

ens18 <- raster("ens18_01.tif")
sand_coord <- sand_alldata[, c("x", "y")]

p <- extract(ens18, sand_coord)

(ev <- evaluates(sand_alldata[, "sand"], p))
plot(ev, main = "...")

rm(sand_coord, ev)

#-------------------------------------------------------------------------------

# Explanatory variables 2022

# Sentinel-2 tiles (sensed 2022) covering entire area clipped by Tigris 2022

{
  Bd02 <- raster("./rs22/Bd02.tif")
  Bd07 <- raster("./rs22/Bd07.tif")
  Bd08 <- raster("./rs22/Bd08.tif")
  Bd09 <- raster("./rs22/Bd09.tif")
  Bd12 <- raster("./rs22/Bd12.tif")
}

rs22 <- stack(Bd02, Bd07, Bd08,	Bd09,	Bd12)

# set coordinate reference system
crs(rs22) <- "+init=epsg:4326 +WGS84"

rm(Bd02, Bd07, Bd08,	Bd09,	Bd12)

#-------------------------------------------------------------------------------

# Ensemble model 2022

ens22 <- ensemble(modelSDM,
                  newdata  = rs22,
                  filename = "ens22_01.tif",
                  setting  = list(method = "weighted", stat = "TSS"))
plot(ens22)

#-------------------------------------------------------------------------------

# Comparison & Visualization

# load potential habitat maps prepared in GIS environment
pothab18 <- raster("ens18_01pothab.tif")
pothab22 <- raster("ens22_01pothab.tif")

# get cell values & remove NAs
val18 <- values(pothab18)
val18[val18 <= 0] <- NA
val18 <- na.omit(val18)
summary(val18)
sd(val18)

val22 <- values(pothab22)
val22[val22 <= 0] <- NA
val22 <- na.omit(val22)
summary(val22)
sd(val22)

# standard error
std.error <- function(x) sd(x) / sqrt(length(x))
std.error(val18)
std.error(val22)

rm(pothab18, pothab22, std.error)

#-------------------------------

# boxplot
par(mfrow = c(1, 1))

boxplot(val18, val22,
        names = c("2018", "2022"),
        outline = T,
        ylab  = "Cell values",
        main  = "Cell no. in potential habitat")

means <- c(mean(val18), mean(val22))
#add means as circles to each boxplot
points(means, pch = 15, cex = 1.5)

# inspect outliers
a <- summary(val18)
# a[5]
a <- unname(a[5])
val18_outliers <- val18[val18 > a]
length(val18_outliers)
sum(val18_outliers)
mean(val18_outliers)
sd(val18_outliers)
std.error(val18_outliers)
# summary(val18_outliers)

a <- summary(val22)
# a[5]
a <- unname(a[5])
val22_outliers <- val22[val22 > a]
length(val22_outliers)
sum(val22_outliers)
mean(val22_outliers)
sd(val22_outliers)
std.error(val22_outliers)

# summary(val22_outliers)

#-------------------------------

# histogram
par(mfrow = c(1, 2))

(sand18 <- hist(val18,
                xlim   = c(0, 0.6), 
                ylim   = c(0, 80000), 
                breaks = 10,
                xlab   = "Cell values",
                ylab   = "Cell numbers",
                main   = "2018"))

(sand22 <- hist(val22,
                xlim   = c(0, 0.6),
                ylim   = c(0, 80000),
                breaks = 10,
                xlab   = "Cell values",
                ylab   = "Cell numbers",
                main   = "2022"))

# number & percentage of cells in groups
sand18 <- data.frame(sand18$mids, sand18$counts)
names(sand18) <- c("mids", "counts")
sand18
perc18 <- sand18$counts / sum(sand18$counts)
sand18 <- cbind(sand18, perc18)
names(sand18) <- c("mids", "counts", "percent")
year <- "2018"
sand18 <- cbind(sand18,year)
sand18

sand22 <- data.frame(sand22$mids, sand22$counts)
names(sand22) <- c("mids", "counts")
sand22
perc22 <- sand22$counts / sum(sand22$counts)
sand22 <- cbind(sand22, perc22)
names(sand22) <- c("mids", "counts", "percent")
year <- "2022"
sand22 <- cbind(sand22,year)
sand22

sandval <- rbind(sand18, sand22)
sandval

#-------------------------------

# categories of values (1: "very poor" ... 6: "best")

catg18 <- cut(val18, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 
                     labels = c(1:6))
table(catg18)
df18 <- data.frame(catg18, val18)
# plot(df18, ylim = c(0, 1))

# set variables for the table (cn: cell numbers, cv: cell values)
cn18 <- 1:6
cv18 <- 1:6

for (i in 1:6) {
  df <- subset(df18, catg18 == i, select = c(catg18, val18))
  cv18[i] <- sum(df$val18)
  cn18[i] <- nrow(df)
}
cv18
sum(cv18)
cn18
sum(cn18)

catg22 <- cut(val22, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6), 
                     labels = c(1:6))
table(catg22)
df22 <- data.frame(catg22, val22)
# plot(df22, ylim = c(0, 1))

cn22 <- 1:6
cv22 <- 1:6

for (i in 1:6) {
  df <- subset(df22, catg22 == i, select = c(catg22, val22))
  cv22[i]  <- sum(df$val22)
  cn22[i] <- nrow(df) 
}
cv22
sum(cv22)
cn22
sum(cn22)

# evaluation table building
# catg <- c("very poor", "poor", "moderate", "good", "best")
catg <- 1:6
df <- data.frame(catg, cn18, round((cn18 / sum(cn18) * 100),1),
                       round(cv18), round((cv18 / sum(cv18) * 100),1),
                       cn22, round((cn22 / sum(cn22) * 100),1),
                       round(cv22), round((cv22 / sum(cv22) * 100),1),
                       (cn22-cn18), round((cv22-cv18))
)
# replace column names (ch_cn: changes in cell numbers, ch_cv: changes in cell values)
names(df) <- c("category", "cn18", "cn18%", "cv18", "cv18%", 
                           "cn22", "cn22%", "cv22", "cv22%",
                            "ch_cn", "ch_cv")
df
# write.csv(df, "eval.csv")

#-------------------------------

# barplot
par(mfrow = c(1, 1))

library(ggplot2)
library(tidyverse)
library(dplyr)

# check data
# (histogram setting above should be: "breaks = 10" for both years)
sandval
# (a small correction to gain equal size of bins)
sandval[22,] <- c(0.525, 0, 0, 2022)
sandval

years <- as.factor(sandval$year)
bins  <- as.numeric(sandval$mids)
perc  <- as.numeric(sandval$percent)

# pivot data

df <- data.frame(years, bins, perc) %>%
  pivot_longer(cols      = bins,
               names_to  = "percent", 
               values_to = "bins")
df

ggplot(df, aes(x = as.factor(bins), y = perc, fill = years)) +
  geom_bar(stat = "identity", position = "dodge2", show.legend = T) +
#  ggtitle("Percentual distribution of cell values by years") +
  xlab("Groups of cell values") +
  ylab("Percent") +
  scale_fill_manual(values = c("black", "darkgray")) +
  theme(panel.background      = element_rect(fill = "transparent"),
        plot.background       = element_rect(fill = "transparent", color = NA),
        panel.grid.major      = element_blank(),
        panel.grid.minor      = element_blank(),
        legend.background     = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        axis.line             = element_line(color = "black"), 
        axis.ticks            = element_line(color = "black"))

rm(df, sandval, years, bins, perc)

#-------------------------------------------------------------------------------

# Mann-Whitney-Wilcoxon Test

# data arrangement
yr                 <- val18
yr[]               <- 2018
allvalues18        <- data.frame(yr, val18)
names(allvalues18) <- c("yr", "val")
rm(yr)

yr                 <- val22
yr[]               <- 2022
allvalues22        <- data.frame(yr, val22)
names(allvalues22) <- c("yr", "val")

allvalues <- rbind(allvalues18, allvalues22)
head(allvalues)
tail(allvalues)

rm(yr, allvalues18, allvalues22)

# MWU test
(MWU <- wilcox.test(val ~ yr, data = allvalues))
MWU[["statistic"]][["W"]]
MWU[["p.value"]]

rm(MWU)

################################################################################

# Selecting threshold of occurrence

library(SDMTools)

# get observations
obs <- sand_alldata$sand

# predictions (already made!!)
# ens18_1 <- raster("ens18_01.tif")
# sand_coord <- sand_alldata[, c("x", "y")]
# p <- extract(ens18_1, sand_coord)

# get predictions (2018), remove NAs
df <- na.omit(data.frame(obs, p))

# sort data
df <- df[order(df$obs, df$p),]

# calculate optimum thresholds
(optth <- optim.thresh(df$obs, df$p))

# select mean occurrence prediction as threshold metric
optth <- optth$mean.occurence.prediction
# optth = 0.2751128

# get rid of unnecessaries
rm(obs, df, p)

#-------------------------------

# Comparison & Visualization

# read data
sandvalall <- read.csv("SandValAll.csv")

#Create a function to generate a continuous color palette
ColorPal <- colorRampPalette(c("blue", "red"))

# add a column of color values based on the cell values
sandvalall$ColorPal <- ColorPal(length(sandvalall))[as.numeric(cut(sandvalall$cellvalue,
                                                                   breaks = length(sandvalall)))]

layout(matrix(1:2, ncol = 1))
par(mar = c(0.3, 5, 4, 2), mgp = c(3, 1, 0), cex.axis = 1.5)
plot(sandvalall$cellvalue, sandvalall$X2018,
     xlim = c(0, 0.5), ylim = c(0, 30),
     cex.lab = 1.5,
     xlab = "Cell value", ylab = "Proportion cells (%)",
     title("2018", cex.main = 1.5, line = -1.5),
     pch = 20, cex = 2, col = sandvalall$ColorPal,
     xaxt = "n")
abline(v = optth, lty = "dotted", lwd = 2, col = "darkred")

par(mar = c(4, 5, 0.3, 2))
plot(sandvalall$cellvalue, sandvalall$X2022, 
     xlim = c(0, 0.5), ylim = c(0, 30),
     cex.lab = 1.5,
     xlab = "Cell value", ylab = "Proportion cells (%)",
     title("2022", cex.main = 1.5, line = -1.5),
     pch = 20, cex = 2, col = sandvalall$ColorPal)
abline(v = optth, lty = "dotted", lwd = 2, col = "darkred")

################################################################################


