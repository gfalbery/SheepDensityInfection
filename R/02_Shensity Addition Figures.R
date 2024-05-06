
# 02_Shensity Addition Figures ####

{
  
  library(tidyverse); library(readxl); library(cowplot); library(ggregplot); library(magrittr)
  library(colorspace); library(INLA); library(gsheet); library(patchwork); library(fs)
  
  theme_set(theme_cowplot())
  
  OutputFiles <- "Output" %>% dir_ls(regex = ".rds") %>% setdiff("Output/FullSpatialModels.rds")
  
  OutputFiles <- OutputFiles[str_detect(OutputFiles, "Addition")]
  
  AdultIMList <- OutputFiles[[1]] %>% readRDS
  FullIMList <- OutputFiles[[2]] %>% readRDS
  LambIMList <- OutputFiles[[3]] %>% readRDS
  YearlingIMList <- OutputFiles[[4]] %>% readRDS
  
  IMList <- OutputFiles %>% map(readRDS)
  
  names(IMList) <- OutputFiles
  
  ParasiteColours2 <- ParasiteColours %>% c(AlberColours[["Pink"]])
  
  names(ParasiteColours2) <- LambIMList %>% names
  
  names(ParasitePalettes) <- c("Strongyles", "Nematodirus", "Coccidia", "Keds", "Capillaria")
  
}

# Figure 1_Maps ####

ParasiteColours2 <- ParasiteColours %>% c(AlberColours[["Pink"]])

names(ParasiteColours2) <- LambIMList %>% names

SpatialParasites <- LambIMList %>% names %>% extract(c(1, 3, 4, 6))

SpatialList <- readRDS("Output/FullSpatialModels.rds")

SpatialList %>% 
  map(~list(.x$FinalModel, .x$Spatial$Model) %>% MDIC %>% unlist %>% diff) %>% 
  unlist %>% round(2)

FourMaps <- 
  SpatialParasites %>% 
  map(~ggField(SpatialList[[.x]]$Spatial$Model,
               Boundary = Boundary,
               SpatialList[[.x]]$Spatial$Mesh) +
        scale_fill_discrete_sequential(palette = ParasitePalettes[[.x]]) +
        geom_point(data = Shpace, inherit.aes = F, aes(X, Y), alpha = 0.01) +
        labs(fill = .x))

FourMaps %>% 
  ArrangeCowplot() + 
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2)

ggsave("Figures/FourMaps.jpeg", units = "mm", height = 250, width = 250)

# Figure 2_Tiles ####

Prevalences <- 
  Shpace %>% summarise_at(Resps, ~Prev(na.omit(.x))) %>% 
  cbind(KedSheep %>% summarise_at("Keds", ~Prev(.x))) %>% 
  sort %>% round(2)*100

OrderParasites <- 
  Shpace %>% summarise_at(Resps, ~Prev(na.omit(.x))) %>% 
  cbind(KedSheep %>% summarise_at("Keds", ~Prev(.x))) %>% 
  sort %>% #rev %>% 
  names

OrderParasiteLabels <- 
  OrderParasites %>% 
  paste0("\n(", Prevalences, "%)")

TileDF <- 
  IMList %>% 
  map(~map(.x, "FinalModel") %>% 
        map(~extract2(.x, "summary.fixed") %>% 
              rename(Estimate = 1, Lower = 3, Upper = 5) %>% 
              rownames_to_column("Var")) %>% 
        bind_rows(.id = "Parasite")) %>% 
  bind_rows(.id = "Model") %>% 
  mutate_at("Model", ~str_remove(.x, "Output/") %>% str_remove("Models.rds$") %>% str_remove("Addition")) %>% 
  filter(Var != "(Intercept)") %>% filter(Var == "AnnualDensity") %>% 
  # filter(Model != "Full") %>% 
  mutate_at("Model", ~factor(.x, levels = c("Lamb", "Yearling", "Adult", "Full"))) %>% 
  mutate_at("Parasite", ~factor(.x, levels = (OrderParasites)))

TileDF <- 
  IMList %>% 
  map(~map(.x, "FinalModel") %>% 
        map(~INLAPValue(.x, c("AnnualDensity"))) %>% 
        bind_rows(.id = "Parasite")) %>% 
  bind_rows(.id = "Model") %>% rename(P = AnnualDensity) %>% 
  mutate_at("Model", ~str_remove(.x, "Output/") %>% str_remove("Models.rds$") %>% str_remove("Addition")) %>% 
  # filter(Model != "Full") %>% 
  mutate_at("Model", ~factor(.x, levels = c("Lamb", "Yearling", "Adult", "Full"))) %>% 
  mutate_at("Parasite", ~factor(.x, levels = (OrderParasites))) %>%
  left_join(TileDF, ., by = c("Model", "Parasite"))

TileDF %<>%
  mutate_at(c("Estimate", "Lower", "Upper"), ~round(.x, 2)) %>% 
  mutate_at("P", ~round(.x, 3)) %>% 
  mutate_at("P", ~ifelse(.x == 0, "<0.001", paste0("=", .x))) %>% 
  # mutate(Label = paste0(Estimate, "\n(", Lower, ", ", Upper, ")")) %>% 
  mutate(Label = paste0(Estimate, "\n(", Lower, ", ", Upper, ")\nP", P)) %>% 
  mutate(Sig = as.numeric(Lower*Upper > 0))

(TilePlot <- 
    TileDF %>% filter(Model != "Full") %>% 
    # mutate_at("Label", ~ifelse(Sig, expression(bold(.x)), .x)) %>% 
    ggplot(aes(Model, Parasite)) +
    # ggplot(aes(Parasite, Model)) +
    geom_tile(aes(fill = Estimate)) +
    # geom_text(aes(label = Label, alpha = Sig, size = Sig)) +
    geom_text(aes(label = Label, alpha = Sig, size = Sig)) +
    # coord_fixed() +
    scale_fill_continuous_diverging(palette = "Tropic") +
    scale_alpha_continuous(range = c(0.4, 1)) +
    scale_size_continuous(range = c(1.5, 2)*2) +
    guides(alpha = F, size = F) + labs(x = NULL, y = NULL) +
    scale_y_discrete(labels = OrderParasiteLabels) +
    # scale_x_discrete(labels = rev(OrderParasiteLabels)) +
    # theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    # coord_flip() +
    NULL)

# Density map

(DensityMap <-
    KUDL %>% GetKUDL() %>% 
    filter(X > 95, X < 105, Y > 985, Y < 998) %>% 
    # filter(Density>0.005) %>%
    filter(Density>0.001) %>%
    ggplot(aes(X, Y, fill = Density)) +
    geom_tile() +
    # theme(panel.background = element_rect(fill = "light blue")) +
    geom_polygon(fill = "light grey", alpha = 0,
                 colour = "black",
                 data = Boundary, inherit.aes = F, aes(X, Y)) +
    coord_fixed() +
    # geom_tile() +
    scale_fill_continuous_sequential(palette = AlberPalettes[[2]], limits = c(0, NA)) +
    geom_contour(aes(z = Density), colour = "white") +
    labs(x = "Easting", y = "Northing") +
    # geom_point(data = Shpace, aes(x = X, y = Y), inherit.aes = F, alpha = 0.01) +
    NULL)

((DensityMap + theme(legend.position = c(0.7, 0.05), legend.justification = c(1, 0)))|
    (TilePlot + theme(legend.position = "none"))) + 
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure2.jpeg", units = "mm", height = 200, width = 320)

EffectDF <- 
  IMList %>% 
  map(~map(.x, "FinalModel") %>% 
        map(~extract2(.x, "summary.fixed") %>% 
              rename(Estimate = 1, Lower = 3, Upper = 5) %>% 
              rownames_to_column("Var")) %>% 
        bind_rows(.id = "Parasite")) %>% 
  bind_rows(.id = "Model") %>% 
  mutate_at("Model", ~str_remove(.x, "Output/") %>% str_remove("Models.rds$") %>% str_remove("Addition")) %>% 
  filter(Var != "(Intercept)") %>% filter(Var %in% c("AnnualDensity", "NPop.t0")) %>% 
  mutate_at("Var", ~str_replace_all(.x, c("AnnualDensity" = "Local", "NPop.t0" = "Global"))) %>% 
  # filter(Model != "Full") %>% 
  mutate_at("Model", ~factor(.x, levels = c("Lamb", "Yearling", "Adult", "Full"))) %>% 
  mutate_at("Parasite", ~factor(.x, levels = (OrderParasites))) %>% 
  mutate(Significant = as.numeric(Lower*Upper > 0))

(EffectPlot <- 
  EffectDF %>% filter(Model != "Full") %>% 
  ggplot(aes(Var, Estimate)) +
  geom_hline(yintercept = 0, lty = 2, alpha = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, 
                    alpha = Significant,
                    colour = Parasite), position = position_dodge(w = 0.8)) +
  geom_point(aes(ymin = Lower, ymax = Upper, 
                 shape = Var,
                 group = Parasite), 
             colour = "black", size = 2,
             position = position_dodge(w = 0.8)) +
  geom_point(aes(ymin = Lower, ymax = Upper, 
                 shape = Var,
                 colour = Parasite), 
             # colour = "black", 
             # size = 2,
             position = position_dodge(w = 0.8)) +
  coord_flip() +
  facet_wrap(Model ~ ., #scales = "free", 
             nrow = 3) +
    scale_shape_manual(values = (c(1:2)+15), limits = c("Local", "Global")) +
  scale_colour_manual(values = rev(ParasiteColours2[OrderParasites]), limits = rev(OrderParasites)) +
  labs(x = NULL, shape = "Density type") +
    theme(strip.background = element_rect(colour = "dark grey", fill = NA)) +
  guides(alpha = "none"))

((DensityMap + theme(legend.position = c(0.7, 0.05), legend.justification = c(1, 0)))|
    (TilePlot + theme(legend.position = "none"))) /
    EffectPlot + 
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1.5, 1))

# ggsave("Figures/Figure2.jpeg", units = "mm", height = 200, width = 350)

ggsave("Figures/Figure2.jpeg", units = "mm", height = 220, width = 300)


# Figure 3_Juvenile Figures ####

ParasiteColours2 <- ParasiteColours %>% c(AlberColours[["Pink"]])

names(ParasiteColours2) <- c("Strongyles", "Nematodirus", "Coccidia", "Keds", "Capillaria", "Strongyloides")

ParasiteOutliers$Keds <- Inf

ParasiteOutliers %<>% unlist

LambPlotList <- 
  
  LambIMList %>% 
  names %>% #setdiff("Strongyloides") %>% 
  map(function(a){
    SmoothOutput(Data = LambIMList[[a]]$Data[LambIMList[[a]]$Data[,a] <= ParasiteOutliers[a],],
                 Model = LambIMList[[a]]$FinalModel,
                 Covariates = c(Covar %>% setdiff("AgeCat") %>% setdiff("Status"), "AnnualDensity"), Response = a,
                 OutputCovariates = c("AnnualDensity"),
                 Output = "Data", 
                 Family = ifelse(a %in% Resps[c(4:5)], "Binomial", "NBinomial"),
                 AddPoints = T, TestDF = LambIMList[[a]]$Data[LambIMList[[a]]$Data[,a] <= ParasiteOutliers[a],], PointAlpha = 0.1, 
                 TextColour = ParasiteColours2[[a]], PointColour = ParasiteColours2[[a]],
                 AddP = T, AddEstimate = T, LimitClip = F, #(a %in% c("Strongyles", "Coccidia", "Keds")),
                 ReturnPlot = F)[[1]]
    
  })

YearlingPlotList <- 
  
  YearlingIMList %>% 
  names %>% 
  map(function(a){
    SmoothOutput(Data = YearlingIMList[[a]]$Data[YearlingIMList[[a]]$Data[,a] <= ParasiteOutliers[a],],
                 Model = YearlingIMList[[a]]$FinalModel,
                 Covariates = c(Covar %>% setdiff("AgeCat") %>% setdiff("Status"), "AnnualDensity"), Response = a,
                 OutputCovariates = c("AnnualDensity"),
                 Output = "Data", 
                 Family = ifelse(a %in% Resps[c(4:5)], "Binomial", "NBinomial"),
                 LineAlpha = 0.1,
                 AddPoints = T, TestDF = YearlingIMList[[a]]$Data[YearlingIMList[[a]]$Data[,a] <= ParasiteOutliers[a],], 
                 PointAlpha = 0.1, 
                 TextColour = ParasiteColours2[[a]], PointColour = ParasiteColours2[[a]],
                 AddP = T, AddEstimate = T, LimitClip = F, #(a %in% c("Strongyles", "Coccidia", "Keds")),
                 ReturnPlot = F)[[1]]
    
  })

AdultPlotList <- 
  
  AdultIMList %>% 
  names %>% 
  map(function(a){
    SmoothOutput(Data = AdultIMList[[a]]$Data[AdultIMList[[a]]$Data[,a] <= ParasiteOutliers[a],],
                 Model = AdultIMList[[a]]$FinalModel,
                 Covariates = c(Covar %>% setdiff("AgeCat") %>% setdiff("Status"), "AnnualDensity"), Response = a,
                 OutputCovariates = c("AnnualDensity"),
                 Output = "Data", 
                 Family = ifelse(a %in% Resps[c(4:5)], "Binomial", "NBinomial"),
                 LineAlpha = 0.1,
                 AddPoints = T, TestDF = AdultIMList[[a]]$Data[AdultIMList[[a]]$Data[,a] <= ParasiteOutliers[a],], PointAlpha = 0.1, 
                 TextColour = ParasiteColours2[[a]], PointColour = ParasiteColours2[[a]],
                 AddP = T, AddEstimate = T, LimitClip = F, #(a %in% c("Strongyles", "Coccidia", "Keds")),
                 ReturnPlot = F)[[1]]
    
  })

PlotList <- 
  LambPlotList[c(1, 4)] %>% 
  append(YearlingPlotList[4]) %>% 
  append(LambPlotList[3]) %>% 
  append(YearlingPlotList[3]) %>% 
  append(AdultPlotList[3])

PlotList[c(1, 4:6)] %<>% map(function(a) a + scale_y_log10())

Titles <- c("Lamb", "Lamb", "Yearling", "Lamb", "Yearling", "Adult")

PlotList <- 
  Titles %>% seq_along %>% 
  map(function(x){
    
    PlotList[[x]] + ggtitle(Titles[x]) + labs(x = "Density")
    
  })

PlotList %>% #map(1) %>% map(~.x + labs(x = "Density")) %>% 
  ArrangeCowplot() + 
  plot_layout(nrow = 2) +
  plot_annotation(tag_levels = "A")

# ggsave("Figures/Figure3.jpeg", units = "mm", height = 220, width = 220)
ggsave("Figures/Figure3.jpeg", units = "mm", height = 200, width = 300)

# Figure 4_Keds ####

KedPlotList <- LambPlotList %>% last %>% list %>% 
  append(YearlingPlotList %>% last %>% list) %>% 
  append(AdultPlotList %>% last %>% list)

KedPlotList[[1]] <- KedPlotList[[1]] + ggtitle("Lambs") + scale_y_log10()

KedPlotList[[2]] <- KedPlotList[[2]] + ggtitle("Yearlings") + scale_y_log10()

KedPlotList[[3]] <- KedPlotList[[3]] + ggtitle("Adults") + scale_y_log10()

KedPlotList %>% 
  map(~.x + labs(x = "Density")) %>% 
  ArrangeCowplot() + 
  plot_layout(nrow = 1) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure4.jpeg", units = "mm", height = 100, width = 250)


# Supplement ####

# Effects Tables ####

TileDF <- 
  IMList %>% 
  map(~map(.x, "FinalModel") %>% 
        map(~extract2(.x, "summary.fixed") %>% 
              rename(Estimate = 1, Lower = 3, Upper = 5) %>% 
              rownames_to_column("Var")) %>% 
        bind_rows(.id = "Parasite")) %>% 
  bind_rows(.id = "Model") %>% 
  mutate_at("Model", ~str_remove(.x, "Output/") %>% str_remove("Models.rds$") %>% str_remove("Addition")) %>% 
  filter(Var != "(Intercept)") %>% #filter(Var == "AnnualDensity") %>% 
  # filter(Model != "Full") %>% 
  mutate_at("Model", ~factor(.x, levels = c("Lamb", "Yearling", "Adult", "Full"))) %>% 
  mutate_at("Parasite", ~factor(.x, levels = (OrderParasites)))

TileDF <- 
  IMList %>% 
  map(~map(.x, "FinalModel") %>% 
        map(~INLAPValue(.x, c("AnnualDensity"))) %>% 
        bind_rows(.id = "Parasite")) %>% 
  bind_rows(.id = "Model") %>% rename(P = AnnualDensity) %>% 
  mutate_at("Model", ~str_remove(.x, "Output/") %>% str_remove("Models.rds$") %>% str_remove("Addition")) %>% 
  # filter(Model != "Full") %>% 
  mutate_at("Model", ~factor(.x, levels = c("Lamb", "Yearling", "Adult", "Full"))) %>% 
  mutate_at("Parasite", ~factor(.x, levels = (OrderParasites))) %>%
  left_join(TileDF, ., by = c("Model", "Parasite"))

TileDF %<>%
  mutate_at(c("Estimate", "Lower", "Upper"), ~round(.x, 2)) %>% 
  # mutate_at("P", ~round(.x, 3)) %>% 
  # mutate_at("P", ~ifelse(.x == 0, "<0.001", paste0("=", .x))) %>% 
  # # mutate(Label = paste0(Estimate, "\n(", Lower, ", ", Upper, ")")) %>% 
  # mutate(Label = paste0(Estimate, "\n(", Lower, ", ", Upper, ")\nP", P)) %>% 
  mutate(Sig = as.numeric(Lower*Upper > 0))

TileDF %<>% mutate_at("Var", ~str_replace_all(.x, c("NPop.t0" = "GlobalDensity", "AnnualDensity" = "LocalDensity")))

TileDF %>% dplyr::select(1:2, Variable = Var, Estimate, Lower, Upper, Sig) %>% 
  write.csv("SupplementaryEstimates.csv")

# Supplementary Tables ####

FullSpatialModels <- readRDS("Output/FullSpatialModels.rds")

FullModels <- 
  OutputFiles %>% 
  map(readRDS) %>% map(~map(.x, "FinalModel")) %>% 
  append(list(FullSpatialModels %>% map(c("Spatial", "Model"))), .)

OutputDF <- 
  FullModels %>% 
  map(function(a) 
    
    SubDF <- 
      map(a, ~extract2(.x, "summary.fixed") %>% 
            rename(Estimate = 1, Lower = 3, Upper = 5) %>% 
            rownames_to_column("Var") %>% 
            mutate_at("Var", ~str_replace_all(.x, c("NPop.t0" = "GlobalDensity", "AnnualDensity" = "LocalDensity"))) %>% 
            mutate(Significant = as.numeric(Lower*Upper > 0))
          
      ) %>% 
      bind_rows(.id = "Parasite"))

1:length(OutputDF) %>% 
  map(function(a) OutputDF[[a]] %>% 
        dplyr::select(1:2, Variable = Var, Estimate, Lower, Upper, Significant) %>% 
        mutate(Model = c("Spatial", OutputFiles %>% str_remove("Output/") %>% str_remove(".rds"))[a]) %>% 
        write.csv(paste0("Output/SupplementaryTable", a, ".csv"), row.names = F))

bind_rows(.id = "Model") %>% 
  mutate_at("Model", ~str_remove(.x, "Output/") %>% str_remove("Models.rds$")) %>% 
  filter(Var != "(Intercept)") %>% filter(Var == "AnnualDensity") %>% 
  # filter(Model != "Full") %>% 
  mutate_at("Model", ~factor(.x, levels = c("Lamb", "Yearling", "Adult", "Full"))) %>% 
  mutate_at("Parasite", ~factor(.x, levels = rev(OrderParasites)))

# New spatial ###'

FullSpatialModels <- readRDS("Output/FullSpatialModels.rds")

names(ParasitePalettes) <- 
  FullSpatialModels[c(1, 4, 3, 6)] %>% names 

FullSpatialModels[c(1, 4, 3, 6)] %>% names %>% 
  map(~ggField(FullSpatialModels[[.x]]$Spatial$Model, 
               FullSpatialModels[[.x]]$Spatial$Mesh) +
        scale_fill_discrete_sequential(palette = ParasitePalettes[[.x]]) +
        geom_density_2d(data = Shpace, aes(X, Y))) %>% 
  ArrangeCowplot() + 
  plot_annotation(tag_levels = "A") +
  plot_layout(ncol = 2)

ggsave("Figures/FourMaps.jpeg", units = "mm", height = 250, width = 180)

# Tiles from spatial models ####

Prevalences <- 
  Shpace %>% summarise_at(Resps, ~Prev(na.omit(.x))) %>% 
  cbind(KedSheep %>% summarise_at("Keds", ~Prev(.x))) %>% 
  sort %>% round(2)*100

OrderParasites <- 
  Shpace %>% summarise_at(Resps, ~Prev(na.omit(.x))) %>% 
  cbind(KedSheep %>% summarise_at("Keds", ~Prev(.x))) %>% 
  sort %>% #rev %>% 
  names

OrderParasiteLabels <- 
  OrderParasites %>% 
  paste0("\n(", Prevalences, "%)")

TileDF <- 
  IMList %>% 
  map(~map(.x, c("Spatial", "Model")) %>% 
        map(~extract2(.x, "summary.fixed") %>% 
              rename(Estimate = 1, Lower = 3, Upper = 5) %>% 
              rownames_to_column("Var")) %>% 
        bind_rows(.id = "Parasite")) %>% 
  bind_rows(.id = "Model") %>% 
  mutate_at("Model", ~str_remove(.x, "Output/") %>% str_remove("Models.rds$")) %>% 
  filter(Var != "(Intercept)") %>% filter(Var == "AnnualDensity") %>% 
  # filter(Model != "Full") %>% 
  mutate_at("Model", ~factor(.x, levels = c("Lamb", "Yearling", "Adult", "Full"))) %>% 
  mutate_at("Parasite", ~factor(.x, levels = rev(OrderParasites)))

TileDF <- 
  IMList %>% 
  map(~map(.x, c("Spatial", "Model")) %>% 
        map(~INLAPValue(.x, "AnnualDensity")) %>% 
        bind_rows(.id = "Parasite")) %>% 
  bind_rows(.id = "Model") %>% rename(P = AnnualDensity) %>% 
  mutate_at("Model", ~str_remove(.x, "Output/") %>% str_remove("Models.rds$")) %>% 
  # filter(Model != "Full") %>% 
  mutate_at("Model", ~factor(.x, levels = c("Lamb", "Yearling", "Adult", "Full"))) %>% 
  mutate_at("Parasite", ~factor(.x, levels = rev(OrderParasites))) %>%
  left_join(TileDF, ., by = c("Model", "Parasite"))

TileDF %<>%
  mutate_at(c("Estimate", "Lower", "Upper"), ~round(.x, 2)) %>% 
  mutate_at("P", ~round(.x, 3)) %>% 
  mutate_at("P", ~ifelse(.x == 0, "<0.001", paste0("=", .x))) %>% 
  # mutate(Label = paste0(Estimate, "\n(", Lower, ", ", Upper, ")")) %>% 
  mutate(Label = paste0(Estimate, "\n(", Lower, ", ", Upper, ")\nP", P)) %>% 
  mutate(Sig = as.numeric(Lower*Upper > 0))

(TilePlot <- 
    TileDF %>% 
    # ggplot(aes(Model, Parasite)) +
    ggplot(aes(Parasite, Model)) +
    geom_tile(aes(fill = Estimate)) +
    geom_text(aes(label = Label, alpha = Sig, size = Sig)) +
    coord_fixed() +
    scale_fill_continuous_diverging(palette = "Tropic") +
    scale_alpha_continuous(range = c(0.4, 1)) +
    scale_size_continuous(range = c(1.5, 2)*2) +
    guides(alpha = F, size = F) + labs(x = NULL, y = NULL) +
    # scale_y_discrete(labels = OrderParasiteLabels) +
    scale_x_discrete(labels = rev(OrderParasiteLabels)) +
    # theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    # coord_flip() +
    NULL)

ggsave("Figures/Figure2b.jpeg", units = "mm", height = 200, width = 220)

# Full addition process ####

"Output" %>% 
  dir_ls(regex = "Addition") %>% 
  map(readRDS) ->
  FullAdditionList

"Output" %>% 
  list.files(pattern = "Addition") %>% str_remove("Addition") %>% str_remove("Models.rds") ->
  names(FullAdditionList)

FullAdditionList[[3]] %>% 
  map("AllModels") %>% map(2) %>% 
  map(~Efxplot(.x, VarOrder = c("NPop", "NPop.t0", "NPop.t_1", "AnnualDensity", "AnnualDensityt0"))) %>% 
  ArrangeCowplot() + plot_layout(guides = "collect")

FullAdditionList %>% names

FullAdditionList %>% 
  map(function(Model){
    
    Model %>% 
      map("FinalModel") %>%
      # map(c("Spatial", "Model")) %>% 
      Efxplot(VarOrder = c("NPop", "NPop.t0", "NPop.t_1", "AnnualDensity", "AnnualDensityt0")) +
      scale_colour_manual(#breaks = c("Strongyles", "Coccidia", "Nematodirus", "Keds"),
        limits = c("Strongyles", "Coccidia", "Nematodirus", "Keds"),
        values = c(AlberColours[[1]], AlberColours[[2]], AlberColours[[3]], AlberColours[[4]]))
    
  }) %>% ArrangeCowplot + plot_layout(guides = "collect")

FullAdditionList[[1]] %>% 
  map("dDIC")

FullAdditionList[[2]] %>% 
  map("dDIC")

FullAdditionList[[3]] %>% 
  map("dDIC")

FullAdditionList[[4]] %>% 
  map("dDIC")

}

# Age plots for conference ####

AgePlotList <- 
  
  AdultIMList[c(1, 3, 4)] %>% 
  names %>% #setdiff("Strongyloides") %>% 
  map(function(a){
    
    print(a)
    
    SmoothOutput(Data = AdultIMList[[a]]$Data[AdultIMList[[a]]$Data[,a] <= ParasiteOutliers[a],],
                 Model = AdultIMList[[a]]$FinalModel,
                 Covariates = c(Covar %>% setdiff(c("Sex", "AgeCat")) %>% c("Age"), "AnnualDensity"), Response = a,
                 OutputCovariates = c("Age"),
                 Output = "Data", 
                 Family = ifelse(a %in% Resps[c(4:5)], "Binomial", "NBinomial"),
                 AddPoints = T, TestDF = AdultIMList[[a]]$Data[AdultIMList[[a]]$Data[,a] <= ParasiteOutliers[a],], PointAlpha = 0.1, 
                 TextColour = ParasiteColours2[[a]], PointColour = ParasiteColours2[[a]],
                 AddP = T, AddEstimate = T, LimitClip = F, #(a %in% c("Strongyles", "Coccidia", "Keds")),
                 ReturnPlot = F)[[1]]
    
  })



AgePlotList %>% 
  ArrangeCowplot()
