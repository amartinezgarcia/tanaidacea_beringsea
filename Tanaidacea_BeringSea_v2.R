##########################################################################################
### 
### Title: Depth-structured diversity: high number of undescribed species 
###           across depth gradient in subarctic regions of the NE Pacific 
###
### Authors: Głuchowska K., Błażewicz M., Fontaneto D., Gellert M., Pourebrahimi S., Martínez A. 
###
### Script written by Alejandro Martínez
### v.1.0 Pallanza 10/10/2024
### Last version: Pallanza, 22/07/2025 [during peer review process]
###
### Requirements: R.4.x
### packages: dplyr, nlme, vegan, BAT, mvabund, psych
##########################################################################################


#### Set working environment  ---------

type_your_wd <- "/Users/amartinez/Dropbox/_Papers/_READY/submitted - Kamilla_Tanaidacea_BeringSea/01 analyses"

setwd(type_your_wd)



#### Load packages ---------

library(dplyr)
library(BAT)
library(nlme)
library(lme4)
library(vegan)


#### Load functions --------

source_path <- file.path(getwd(), "functions", "Tanaidacea_Bering_functions.R")
if (file.exists(source_path)) {
  source(source_path)
} else {
  stop("Function file not found at: ", source_path)
}



############ Part 0: Load and clean data -------------------------------

data <- read.csv2("Tanaids North Pacific.csv")

data$species.code <- paste0(data$Genus,"_",data$MorphoSpecies)
data$species.code <- gsub(". ", "", data$species.code, fixed = T)


# Environmental data table

ecol <- data %>%
  select(Station, Area, depth_start, depth_end, Current.velocity,
         Current.direction, Salinity, Silicate, Phosphate, pH,
         Dissolved.oxygen, Nitrate, Iron, Temperature,
         Lat_start_dec, Lat_end_dec, Long_start_dec, Long_end_dec) %>%
  distinct()


# Calculate richness by station

richness <- data %>% count(Station, name = "Richness")
ecol <- left_join(ecol, richness, by = "Station")


## Preparing community matrix

comm <- abspres(data, sites.col = "Station", sp.col = "species.code")
row.names(comm) <- comm[,1]
comm <- comm[,-1]


## Checking for correlation amongst ecological variables
psych::pairs.panels(ecol %>% select(depth_start, Current.velocity, Current.direction,
                                    Salinity, Silicate, Phosphate, pH,
                                    Dissolved.oxygen, Nitrate, Iron, Temperature))

# Silicate and Iron are correlated (I keep silicate, cause Iron also correlate with pH)
# Dissolved.oxygen correlates with pH (I keep O2, cause it is easier to discuss)
# "Phosphate, "Nitrate" and "Temperature" are correlated, I keep temperature
# "Salinity" and "Temperature" are correlated, I keep temperature
# Depth and temperature are correlated... I keep depth, cause it is one of hte variables you want to test.

ecol00 <- ecol[c("Station",
                 "Area",
                  "depth_start",
                  "Current.velocity",
                  "Current.direction",
                  "Silicate",
                  "Dissolved.oxygen",
                 "Iron",
                  "Lat_start_dec",
                   "Lat_end_dec",
                   "Long_start_dec",
                  "Long_end_dec",
                  "Richness")]


psych::pairs.panels(ecol00 %>% select(depth_start, Current.velocity,
                                    Current.direction, Silicate,
                                    Dissolved.oxygen))



## I check the structure of the table
str(ecol00)

# Final cleaned ecological dataframe
ecol00 <- ecol %>%
  select(Station, Area, depth_start, Current.velocity, Current.direction,
         Silicate, Dissolved.oxygen, Iron, Lat_start_dec, Lat_end_dec,
         Long_start_dec, Long_end_dec, Richness) %>%
  mutate(
    across(c(depth_start, Current.velocity, Current.direction,
             Silicate, Dissolved.oxygen, Lat_start_dec, Long_start_dec), as.numeric), # numeric variables
    Lat = Lat_start_dec + runif(n(), -0.0001, 0.0001),
    Long = Long_start_dec + runif(n(), -0.0001, 0.0001) # add a fuzzy decimal to coordinates
  )


############ Part 1: RICHNESS ANALYSES  -------------------------------

##### 1.1. Modelling species richness version environmental parameters ----------

f.tax <- formula(Richness ~ depth_start + Current.velocity +
                   Current.direction + Silicate + Dissolved.oxygen)


# Spatial GLS models
Tax.B0 <- nlme::gls(f.tax, data = ecol00)
Tax.B0A <- nlme::gls(f.tax, correlation = corSpher(form = ~ Long + Lat, nugget = TRUE), data = ecol00)
Tax.B0B <- nlme::gls(f.tax, correlation = corLin(form = ~ Long + Lat, nugget = TRUE), data = ecol00)
Tax.B0C <- nlme::gls(f.tax, correlation = corRatio(form = ~ Long + Lat, nugget = TRUE), data = ecol00)
Tax.B0D <- nlme::gls(f.tax, correlation = corGaus(form = ~ Long + Lat, nugget = TRUE), data = ecol00)
Tax.B0E <- nlme::gls(f.tax, correlation = corExp(form = ~ Long + Lat, nugget = TRUE), data = ecol00)

AIC(Tax.B0, Tax.B0A, Tax.B0B, Tax.B0C, Tax.B0D, Tax.B0E) ## We keep model B0A

summary(Tax.B0A) # Richness is uniform.



### Let's avoid spatial structure and instead, include the areas in the model.

Tax.B1 <- glm(Richness ~ depth_start + Current.velocity +
                Current.direction + Silicate + Dissolved.oxygen + Area,
              data = ecol00,
              family = "poisson")

performance::check_model(Tax.B1)
summary(Tax.B1)
car::Anova(Tax.B1)
## Silicate, Dissolved.oxygen, and depth are correlated with Area, but
# more seriously, the variance have a lot of structure (probably correspodning with the areas). We need to use the 
# spatially explicit model

rm(Tax.B0,Tax.B0A,Tax.B0B,Tax.B0C,Tax.B0D,Tax.B0E,Tax.B1,f.tax) #Clean the room





############ Part 2: SPECIES COMPOSITION ANALYSES //-------------------------------

###### 2.1. Distance-based method using adonis ------------------------------------


beta.tax <-  BAT::beta(comm, abund=F) 

model.Btotal <- beta.tax$Btotal ~ depth_start + Current.velocity + Current.direction + 
  Silicate + Dissolved.oxygen + Area
adonis2(model.Btotal, data = ecol00, permutations = 9999)

model.Bnest <- beta.tax$Brich ~ depth_start + Current.velocity + Current.direction + 
  Silicate + Dissolved.oxygen + Area
adonis2(model.Bnest, data = ecol00, permutations = 9999)

model.Bturn <- beta.tax$Brepl ~ depth_start + Current.velocity + Current.direction + 
  Silicate + Dissolved.oxygen + Area
adonis2(model.Bturn, data = ecol00, permutations = 9999)


rm(model.Btotal,model.Bnest,model.Bturn)


###### 2.2. Mvabund composition analysis ------------------------------------

## See discussion here: https://environmentalcomputing.net/statistics/mvabund/


comm.mv <- mvabund::mvabund(comm)
mvabund::meanvar.plot(comm.mv)

mod.multi <- mvabund::manyglm(comm.mv ~ depth_start + Current.velocity + Current.direction + 
                                Silicate + Dissolved.oxygen + Area,
                     family = "negative_binomial", data = ecol00)
plot(mod.multi)


anova_mod.multi_adj <-  anova(mod.multi, p.uni = "adjusted")
saveRDS(anova_mod.multi_adj,"anova_mod.multi_adj.rds")
#anova <- readRDS("anova_mod.multi_adj.rds")


uni.p <- as.data.frame(t(anova$uni.p))
uni.p$species <- row.names(uni.p)
uni.test <- as.data.frame(t(anova$uni.test))
uni.test$species <- row.names(uni.test)

uni.analyses <- merge(uni.test,uni.p, by = "species")

colnames(uni.analyses) <- c("species",
                            "Dev.intercept",
                            "Dev.depth.start",
                            "Dev.current.velocity",
                            "Dev.current.direction",
                            "Dev.silicate",
                            "Dev.dissolved.oxygen",
                            "Dev.Area",
                            "p.intercept",
                            "p.depth.start",
                            "p.current.velocity",
                            "p.current.direction",
                            "p.silicate",
                            "p.dissolved.oxygen",
                            "p.Area")

write.csv2(uni.analyses, "mvabund.species.csv")


##### Cluster showing the differences in stations per area, but using beta diversity

dist_matrix <- as.dist(beta.tax$Btotal)
clustering <- hclust(dist_matrix, method = "average") # Using average linkage
plot(clustering, main = "Clustering of Sites Based on Beta Diversity")

     