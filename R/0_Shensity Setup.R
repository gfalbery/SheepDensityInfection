
# Shensity ####

rm(list = ls())

library(tidyverse); library(readxl); library(cowplot); library(ggregplot); library(magrittr)
library(colorspace); library(INLA); library(gsheet); library(patchwork); library(adehabitatHR)

Root <- paste0("Data/")

Sheep <- read.csv(paste0(Root, "/SheepIndividuals.csv"), header = T)
SheepCensuses <- read.csv(paste0(Root, "/SheepCensuses.csv"), header = T)
SheepParasites <- read.csv(paste0(Root, "/SheepParasites.csv"), header = T)
SheepBirths <- read.csv(paste0(Root, "/SheepBirths.csv"), header = T)
HirtaVillPop <- read.csv(paste0(Root, "/HirtaVillPop.csv"), header = T)

Dates <- read.csv(paste0(Root, "/Dates.csv"), header = T)

# Sorting out censuses ####

SheepCensuses <- merge(SheepCensuses, Dates[,c("DateUK","Ndate")], all.x = T, by.x = "Date", by.y = "DateUK")
SheepCensuses$ID <- as.factor(SheepCensuses$ID)

YearStarts <- Dates[substr(Dates$DateUK,1,5) == "01/01","Ndate"]
SpringStarts <- Dates[substr(Dates$DateUK,1,5) == "01/05","Ndate"][2:52]

Timing <- data.frame(Year = 1969:2019, YearStarts, SpringStarts)

SheepCensuses$Year <- as.numeric(as.character(cut(SheepCensuses$Ndate,breaks=YearStarts,labels=1969:2018)))
SheepCensuses$SheepYear <- as.numeric(as.character(cut(SheepCensuses$Ndate,breaks=SpringStarts,labels=1969:2018)))

YearLocations <- with(SheepCensuses, cbind(reshape2::melt(tapply(Easting, list(ID, SheepYear), function(m) mean(m, na.rm=T))),
                                           reshape2::melt(tapply(Northing, list(ID, SheepYear), function(m) mean(m, na.rm=T)))$value))

names(YearLocations) <- c("ID", "SheepYear", "Easting", "Northing")

YearLocations <- na.omit(YearLocations[YearLocations$Northing>975,])

# Sorting out parasites ####

SheepParasites <- merge(SheepParasites, Dates[,c("DateUK","Ndate")], all.x = T, by.x = "DateCollect", by.y = "DateUK")

SheepParasites$Year <- as.numeric(as.character(cut(SheepParasites$Ndate,breaks=YearStarts,labels=1969:2018)))
SheepParasites$SheepYear <- as.numeric(as.character(cut(SheepParasites$Ndate,breaks=SpringStarts,labels=1969:2018)))
SheepParasites$August <- ifelse(substr(SheepParasites$DateCollect, 4, 5) == "08", "Y", "N")

SheepParasites %<>% separate(DateCollect, sep = "/", into = paste0(c("Day", "Month", "Year"), "Collect"), remove = F)

Shpace <- merge(SheepParasites[SheepParasites$August == "Y", ], YearLocations, by = c("ID", "SheepYear"))

Shpace$ID <- as.factor(Shpace$ID)
Shpace$fYear <- as.factor(Shpace$SheepYear)

Shpace %>% rename(Coccidia = Coccidea) -> Shpace

Shpace %>% mutate_at("Coccidia", ~ifelse(is.na(.x), 0, .x))
Shpace[Shpace$Year%in%1988:1992, "Coccidia"] <- NA

# Sorting out individuals #####

Sheepvar <- c("ID", "BirthRef", "Sex", "Coat", "Horn")

Shpace %>% left_join(Sheep[,Sheepvar] %>% mutate(ID = as.factor(ID)), by = "ID") %>%
  left_join(SheepBirths[,c("BirthRef", "BirthYear")], by ="BirthRef") %>%
  mutate(Age = SheepYear - BirthYear,
         AgeCat = cut(Age, breaks = c(-1, 0.5, 1.5, Inf), 
                      labels = c("Lamb", "Yearling", "Adult"), include.lowest = T),
         Sex = factor(c("F","M")[Sex])) %>%
  left_join(HirtaVillPop[,c("Year", "VillTotal")], by = c("SheepYear" = "Year")) %>%
  rename(NPop = VillTotal) %>% 
  left_join(HirtaVillPop[,c("Year", "VillTotal")] %>% mutate_at("Year", ~.x + 1) %>% 
              rename(NPop.t0 = VillTotal), by = c("SheepYear" = "Year")) %>% 
  left_join(HirtaVillPop[,c("Year", "VillTotal")] %>% mutate_at("Year", ~.x + 2) %>% 
              rename(NPop.t_1 = VillTotal), by = c("SheepYear" = "Year")) -> 
  Shpace

Shpace %<>% rename(X = Easting, Y = Northing)

# Adding Keds ####

KedCounts <- read.csv(paste0(Root, "/KedCounts.csv"), header = T) %>%
  rename(SheepYear = Day.Month.Year,
         Month = txtCapmonth,
         Day = txtCapDay) %>% 
  filter(!SheepYear %in% 1985:1987, 
         Month == 8, Keds<15)

KedCounts %>%
  left_join(YearLocations, by = c("ID", "SheepYear")) %>%
  left_join(Sheep[,Sheepvar], by = "ID") %>%
  left_join(SheepBirths[,c("BirthRef", "BirthYear")], by ="BirthRef") %>%
  mutate(Age = SheepYear - BirthYear,
         AgeCat = cut(Age, breaks = c(-1, 1, 2, Inf),
                      labels = c("Lamb", "Yearling", "Adult"), include.lowest = T),
         Sex = factor(c("F","M")[Sex])) %>%
  left_join(HirtaVillPop[,c("Year", "VillTotal")] %>% mutate_at("Year", ~.x + 1) %>% 
              rename(NPop.t0 = VillTotal), by = c("SheepYear" = "Year")) -> 
  
  KedSheep

KedSheep %<>% rename(X = Easting, Y = Northing)

# Adding density ####

FocalYears <- Shpace$Year %>% unique %>% sort

spdf <- SpatialPointsDataFrame(data = Shpace[,c("X", "Y", "Year")], coords = Shpace[,c("X", "Y")])

spdf <- spdf[,"Year"]

KUDLYear <- kernelUD(spdf, same4all = TRUE, grid = 500)

Shpace2 <- 
  FocalYears[-1] %>% map(function(x){
    
    SubDF <- Shpace %>% filter(Year == x)
    
    KUDLYear[[as.character(x)]] %>% 
      raster::raster() %>% raster::extract(SubDF[,c("X", "Y")]) ->
      SubDF$AnnualDensity
    
    SubDF %>% return
    
  }) %>% bind_rows()

SaveShpace <- Shpace

Shpace <- Shpace2

# Adding Keds ####

KedSheep %<>% rename(Year = SheepYear)

KedSheep %<>% filter(!is.na(X))

FocalYears <- KedSheep$Year %>% unique %>% sort

spdf <- SpatialPointsDataFrame(data = KedSheep[,c("X", "Y", "Year")], coords = KedSheep[,c("X", "Y")])

spdf <- spdf[,"Year"]

KUDLYear <- kernelUD(spdf, same4all = TRUE, grid = 500)

KedSheep2 <- 
  FocalYears[-1] %>% map(function(x){
    
    SubDF <- KedSheep %>% filter(Year == x)
    
    KUDLYear[[as.character(x)]] %>% 
      raster::raster() %>% raster::extract(SubDF[,c("X", "Y")]) ->
      SubDF$AnnualDensity
    
    SubDF %>% return
    
  }) %>% bind_rows()

SaveKedSheep <- KedSheep

KedSheep <- KedSheep2

KedSheep$SheepYear <- KedSheep$Year

# Reproduction ####

ReproDF <- left_join(SheepBirths, Sheep, by = "BirthRef") 

Shpace <- 
  ReproDF %>% 
  filter(!is.na(MumID)) %>% 
  filter(!is.na(BirthYear)) %>% 
  group_by(MumID, BirthYear) %>% 
  count() %>% filter(n>0) %>% 
  mutate(Reproductive = 1) %>% mutate_at("MumID", as.factor) %>% 
  left_join(Shpace, ., by = c("ID" = "MumID", "SheepYear" = "BirthYear")) %>% 
  mutate(Status=case_when(AgeCat=="Lamb" ~ "Lamb",
                          AgeCat=="Yearling" ~ "Yearling", 
                          AgeCat=="Adult" & Reproductive==1 ~ "ReproAdultF", 
                          AgeCat=="Adult" & Sex=="F" & is.na(Reproductive) ~ "NonReproAdultF", 
                          AgeCat=="Adult" & Sex=="M" ~ "AdultM")) 


Shpace %<>% # ADDING CHECK IT HASN'T FUCKED ANYTHING
  mutate_at("ID", ~as.numeric(as.character(.x))) %>% 
  left_join(
    Sheep %>%
      mutate(DeathSheepYear = ifelse(DeathMonth < 5, DeathYear - 1, DeathYear)) %>%
      mutate_at("ID", ~as.numeric(as.character(.x))) %>%
      dplyr::select(ID, setdiff(colnames(.), c(colnames(Shpace), "Keds"))), by = "ID") %>%
  mutate(Survived = as.numeric(DeathSheepYear > SheepYear))

Shpace %>% count(AgeCat, #SexRepro, 
                 Reproductive)

Shpace %>% 
  mutate_at("ID", ~.x %>% as.factor %>% as.numeric) %>% 
  dplyr::select(-c(Key, PERSON, n, BirthRef, Coat, Horn)) %>% 
  saveRDS("Data/Shpace.rds")

KedSheep <- 
  ReproDF %>% 
  filter(!is.na(MumID)) %>% 
  filter(!is.na(BirthYear)) %>% 
  group_by(MumID, BirthYear) %>% 
  count() %>% filter(n>0) %>% 
  mutate(Reproductive = 1) %>% mutate_at("MumID", as.factor) %>% 
  left_join(KedSheep %>% mutate_at("ID", as.factor), ., 
            by = c("ID" = "MumID", "SheepYear" = "BirthYear")) %>% 
  mutate(Status=case_when(AgeCat=="Lamb" ~ "Lamb",
                          AgeCat=="Yearling" ~ "Yearling", 
                          AgeCat=="Adult" & Reproductive==1 ~ "ReproAdultF", 
                          AgeCat=="Adult" & Sex=="F" & is.na(Reproductive) ~ "NonReproAdultF", 
                          AgeCat=="Adult" & Sex=="M" ~ "AdultM")) 

KedSheep %>% 
  mutate_at("ID", ~.x %>% as.factor %>% as.numeric) %>% 
  dplyr::select(-c(Tag, BirthRef, Coat, Horn, New.Tag)) %>% 
  saveRDS("Data/KedSheep.rds")

# Adding boundary? ####

Boundary <- 
  SheepCensuses %>%
  filter(Northing > 900) %>%
  # mutate_at("Northing", ~(.x - 800000)/100) %>%
  dplyr::select(X = Easting,
                Y = Northing) %>%
  count(X, Y) %>% 
  arrange(X, Y) %>% mutate(N = 1:n()) %>% 
  slice(c(
    
    1:4, 14, 26, 217, 230, 227, 218, 
    183, 173, 132, 121, 119, 118, 101, 61, 27, 5
    
  )) %>% dplyr::select(X, Y)

BoundaryMatrix <- 
  Boundary %>% 
  slice(n():1) %>% 
  as.matrix

BoundaryMatrix %>% saveRDS("Data/BoundaryMatrix.rds")

