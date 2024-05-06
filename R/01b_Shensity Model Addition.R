
# 01b_Shensity Model Addition ####

rm(list = ls())

library(tidyverse); library(readxl); library(cowplot); library(ggregplot); library(magrittr)
library(colorspace); library(INLA); library(gsheet); library(patchwork); library(fs)

dir_create("Output")

if(!file.exists("Data/Shpace.rds")) source("0_Shensity Setup.R") else{
  
  Shpace <- readRDS("Data/Shpace.rds")
  
  KedSheep <- readRDS("Data/KedSheep.rds")
  
  BoundaryMatrix <- readRDS("Data/BoundaryMatrix.rds")
  
}

{
  
  theme_set(theme_cowplot())
  
  Resps <- names(Shpace)[which(names(Shpace) == "Strongyles"):which(names(Shpace) == "Monezia")]
  
  Resps <- Resps %>% setdiff(c("Monezia", "Trichuris", 
                               "Strongyloides", "Capillaria")[1:2])
  
  ParasiteOutliers <- c(Strongyles = 2500, 
                        Strongyloides = 400, 
                        Coccidia = 20000, 
                        Nematodirus = 1,
                        Capillaria = 1)
  
  Shpace %>% mutate_at(vars(Resps[ParasiteOutliers == 1]), AsBinary) ->
    Shpace
  
  YearlingResps <- Resps %>% setdiff("Nematodirus")
  AdultResps <- Resps %>% setdiff(c("Nematodirus", "Capillaria"))
  
  Covar <- c("Status", "Sex")
  NonCovar <- c("ID", "X", "Y", "SheepYear")
  AddCovar <- c(#"NPop", 
    "NPop.t0", #"NPop.t_1", 
    "AnnualDensity") #, "AnnualDensityt0")
  
  FamilyList <- c(rep("nbinomial", 3), rep("binomial", 2), "nbinomial")
  names(FamilyList) <- Resps[1:5] %>% c("Keds")
  
}

# Overall Models ####

IMList <- list()

r <- 1

for(r in r:length(Resps)){
  
  print(Resps[r])
  
  TestDF <- 
    Shpace %>% 
    mutate(Status = fct_relevel(Status, "Lamb")) %>% 
    dplyr::select(all_of(c(Covar, 
                           NonCovar, 
                           AddCovar,
                           Resps[r])), AnnualDensity, AgeCat) %>% 
    na.omit %>% droplevels
  
  TestDF <- 
    TestDF[TestDF[,Resps[r]] <= ParasiteOutliers[[Resps[r]]], ]
  
  TestDF %>% nrow %>% print
  
  IM1 <- 
    INLAModelAdd(Data = TestDF, 
                 Response = Resps[r],
                 Family = FamilyList[[Resps[r]]],
                 Explanatory = Covar, # %>% c("AnnualDensity"), # %>% setdiff("AgeCat"), 
                 Random = c("ID", "SheepYear"), RandomModel = 'iid',
                 AllModels = T, Base = T, Delta = -Inf,
                 Add = AddCovar,
                 Boundary = BoundaryMatrix,
                 AddSpatial = T
    )
  
  IMList[[Resps[r]]] <- IM1
  
}

# Keds ####

print("Keds")

TestDF <- 
  KedSheep %>% 
  # filter(AgeCat == "Lamb") %>% 
  mutate(Status = fct_relevel(Status, "Lamb")) %>% 
  dplyr::select(all_of(c(Covar, 
                         AddCovar,
                         NonCovar)), Keds, AnnualDensity, AgeCat) %>% 
  na.omit %>% droplevels

TestDF %>% nrow %>% print

IM1 <- 
  INLAModelAdd(Data = TestDF, 
               Response = "Keds",
               Family = "nbinomial",
               Explanatory = Covar, # %>% c("AnnualDensity"), # %>% setdiff("AgeCat"), 
               Random = c("ID", "SheepYear"), RandomModel = 'iid',
               AllModels = T, Base = T, Delta = -Inf,
               Add = AddCovar,
               Boundary = BoundaryMatrix,
               AddSpatial = T
  )

IMList[["Keds"]] <- IM1

saveRDS(IMList, "Output/AdditionFullModels.rds")

# lapply(1:length(IMList), function(x){
#   Efxplot(IMList[[x]]$FinalModel) + 
#     ggtitle(names(IMList)[x])
# }) -> AllAgePlots 

# Lamb Only Models ####

CovarLamb <- "Sex" 

IMList <- list()

r <- 1

for(r in r:length(Resps)){
  
  print(Resps[r])
  
  TestDF <- 
    Shpace %>% 
    filter(AgeCat == "Lamb") %>%
    dplyr::select(all_of(c(CovarLamb, 
                           AddCovar,
                           NonCovar, Resps[r])), AnnualDensity) %>% 
    na.omit %>% droplevels
  
  TestDF %>% nrow %>% print
  
  IM1 <- 
    INLAModelAdd(Data = TestDF, 
                 Response = Resps[r],
                 Family = FamilyList[[Resps[r]]],
                 Explanatory =  CovarLamb, # %>% c("AnnualDensity") ,#%>% setdiff("AgeCat"), 
                 Random = c("ID", "SheepYear"), RandomModel = 'iid',
                 Boundary = BoundaryMatrix,
                 Add = AddCovar,
                 AllModels = T, Base = T, Delta = -Inf,
                 AddSpatial = T
    )
  
  IMList[[Resps[r]]] <- IM1
  
}

# Keds ####

print("Keds")

TestDF <- 
  KedSheep %>% 
  filter(AgeCat == "Lamb") %>%
  dplyr::select(all_of(c(CovarLamb, 
                         AddCovar,
                         NonCovar)), Keds, AnnualDensity) %>% 
  na.omit %>% droplevels

TestDF %>% nrow %>% print

IM1 <- 
  INLAModelAdd(Data = TestDF, 
               Response = "Keds",
               Family = "nbinomial",
               Explanatory = CovarLamb, # %>% c("AnnualDensity"), #%>% setdiff("AgeCat"), 
               Random = c("ID", "SheepYear"), RandomModel = 'iid',
               Boundary = BoundaryMatrix,
               Add = AddCovar,
               AllModels = T, Base = T, Delta = -Inf,
               AddSpatial = T
  )

IMList[["Keds"]] <- IM1

saveRDS(IMList, "Output/AdditionLambModels.rds")

# Yearling Only Models ####

CovarLamb <- "Sex" 

IMList <- list()

r <- 1

for(r in (r:length(YearlingResps))){
  
  print(YearlingResps[r])
  
  TestDF <- 
    Shpace %>% 
    filter(AgeCat == "Yearling") %>%
    dplyr::select(all_of(c(CovarLamb, 
                           AddCovar,
                           NonCovar, YearlingResps[r])), AnnualDensity) %>% 
    na.omit %>% droplevels
  
  TestDF %>% nrow %>% print
  
  IM1 <- 
    INLAModelAdd(Data = TestDF, 
                 Response = YearlingResps[r],
                 Family = FamilyList[[YearlingResps[r]]],
                 Explanatory = CovarLamb, # %>% c("AnnualDensity"), #%>% setdiff("AgeCat"), 
                 Random = c("ID", "SheepYear"), RandomModel = 'iid',
                 Add = AddCovar,
                 Boundary = BoundaryMatrix,
                 AllModels = T, Base = T, Delta = -Inf,
                 AddSpatial = T
    )
  
  IMList[[YearlingResps[r]]] <- IM1
  
}

# Keds ####

print("Keds")

TestDF <- 
  KedSheep %>% 
  filter(AgeCat == "Yearling") %>%
  dplyr::select(all_of(c(CovarLamb, 
                         AddCovar,
                         NonCovar)), Keds, AnnualDensity) %>% 
  na.omit %>% droplevels

TestDF %>% nrow %>% print

IM1 <- 
  INLAModelAdd(Data = TestDF, 
               Response = "Keds",
               Family = "nbinomial",
               Explanatory = CovarLamb, # %>% c("AnnualDensity"), #%>% setdiff("AgeCat"), 
               Random = c("ID", "SheepYear"), RandomModel = 'iid',
               Boundary = BoundaryMatrix,
               Add = AddCovar,
               AllModels = T, Base = T, Delta = -Inf,
               AddSpatial = T
  )

IMList[["Keds"]] <- IM1

saveRDS(IMList, "Output/AdditionYearlingModels.rds")

# Adult Only Models ####

IMList <- list()

r <- 1

for(r in r:length(AdultResps)){
  
  print(AdultResps[r])
  
  TestDF <- 
    Shpace %>% 
    filter(AgeCat == "Adult") %>%
    dplyr::select(all_of(c(Covar, "Age", 
                           NonCovar,
                           AddCovar,
                           AdultResps[r]))) %>% 
    na.omit %>% droplevels
  
  TestDF %>% nrow %>% print
  
  IM1 <- 
    INLAModelAdd(Data = TestDF, 
                 Response = AdultResps[r],
                 Family = FamilyList[[AdultResps[r]]],
                 Explanatory = Covar %>% c("Age") %>% # c("AnnualDensity", "Age") %>% 
                   setdiff("AgeCat") %>% 
                   setdiff("Sex"), 
                 Random = c("ID", "SheepYear"), RandomModel = 'iid',
                 AllModels = T, Base = T, Delta = -Inf,
                 Add = AddCovar,
                 Boundary = BoundaryMatrix,
                 Beep = F,
                 AddSpatial = T
    )
  
  IMList[[AdultResps[r]]] <- IM1
  
}

# Keds ####

print("Keds")

TestDF <- 
  KedSheep %>% 
  filter(AgeCat == "Adult") %>%
  dplyr::select(all_of(c(Covar, "Age", 
                         AddCovar,
                         NonCovar)), Keds, AnnualDensity) %>% 
  na.omit %>% droplevels

TestDF %>% nrow %>% print

IM1 <- 
  INLAModelAdd(Data = TestDF, 
               Response = "Keds",
               Family = "nbinomial",
               Explanatory = Covar %>% c("Age") %>% #c("AnnualDensity", "Age") %>% 
                 setdiff("AgeCat") %>% 
                 setdiff("Sex"), 
               Random = c("ID", "SheepYear"), RandomModel = 'iid',
               Boundary = BoundaryMatrix,
               AllModels = T, Base = T, Delta = -Inf,
               Add = AddCovar,
               Beep = F,
               AddSpatial = T
  )

IMList[["Keds"]] <- IM1

saveRDS(IMList, "Output/AdditionAdultModels.rds")

# Overall Spatial Models ####

IMList <- list()

r <- 1

for(r in r:length(Resps)){
  
  print(Resps[r])
  
  TestDF <- 
    Shpace %>% 
    # filter(AgeCat == "Lamb") %>% 
    mutate(Status = fct_relevel(Status, "Lamb")) %>% 
    dplyr::select(all_of(c(Covar, NonCovar, Resps[r])), AnnualDensity, AgeCat) %>% 
    na.omit %>% droplevels
  
  TestDF <- 
    TestDF[TestDF[,Resps[r]] <= ParasiteOutliers[[Resps[r]]], ]
  
  TestDF %>% nrow %>% print
  
  IM1 <- 
    INLAModelAdd(Data = TestDF, 
                 Response = Resps[r],
                 Family = FamilyList[[Resps[r]]],
                 Explanatory = Covar, # %>% c("AnnualDensity"), # %>% setdiff("AgeCat"), 
                 Random = c("ID", "SheepYear"), RandomModel = 'iid',
                 Boundary = BoundaryMatrix,
                 AllModels = T, Base = T, Delta = -Inf,
                 Groups = T, GroupVar = "Status",
                 AddSpatial = T
    )
  
  IMList[[Resps[r]]] <- IM1
  
}

# Keds ####

print("Keds")

TestDF <- 
  KedSheep %>% 
  mutate(Status = fct_relevel(Status, "Lamb")) %>% 
  dplyr::select(all_of(c(Covar, NonCovar)), Keds, AnnualDensity, AgeCat) %>% 
  na.omit %>% droplevels

TestDF %>% nrow %>% print

IM1 <- 
  INLAModelAdd(Data = TestDF, 
               Response = "Keds",
               Family = "nbinomial",
               Explanatory = Covar, # %>% c("AnnualDensity"), # %>% setdiff("AgeCat"), 
               Random = c("ID", "SheepYear"), RandomModel = 'iid',
               Boundary = BoundaryMatrix,
               AllModels = T, Base = T, Delta = -Inf,
               AddSpatial = T
  )

IMList[["Keds"]] <- IM1

saveRDS(IMList, "Output/AdditionFullSpatialModels.rds")
