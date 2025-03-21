#### Import training data and script ####
source("Import/Script/BRT_script.R")
Meco <- readRDS("Import/Training/Meco.Rds")
Mclim <- readRDS("Import/Training/Mclim.Rds")
Msurf.mean <- readRDS("Import/Training/MbrGDGT.Rds")
M.br.GDGT <- Msurf.mean[,grep("^f.I", colnames(Msurf.mean))]

#### Fuzzy C-mean clustering ####
Cmean = T
if(Cmean == T){
  #### Subsets des FA ####
  To.remove <- grep("[7,8]", colnames(M.br.GDGT))
  To.keep <- seq(1, ncol(M.br.GDGT))
  To.keep <- setdiff(To.keep, To.remove)
  M.br.GDGT.km <- M.br.GDGT[,To.keep]
  M.br.GDGT.km <- data.frame(t(M.br.GDGT.km), check.names = F)
  M.br.GDGT.km <- apply(M.br.GDGT.km, 2, MESS::round_percent)
  M.br.GDGT.km <- data.frame(t(M.br.GDGT.km/100), check.names = F)
  
  Keep.samptype <- c("Soil", "Moss", "Lacustrine")
  M.br.GDGT.km <- M.br.GDGT.km[row.names(M.br.GDGT.km) %in% row.names(Meco[Meco$Sample.type %in% Keep.samptype,]),]
  Mclim.km <- Mclim[row.names(Mclim) %in% row.names(Meco[Meco$Sample.type %in% Keep.samptype,]),]
  Meco.km <- Meco[row.names(Meco) %in% row.names(Meco[Meco$Sample.type %in% Keep.samptype,]),]
  
  #### Silhouet test (PAM) ####
  PAM.test = T
  if(PAM.test == T){
    silhouette_score <- function(k){
      km <- kmeans(M.br.GDGT.km, centers = k, nstart=25)
      ss <- silhouette(km$cluster, dist(M.br.GDGT.km))
      mean(ss[, 3])
    }
    k <- 2:10
    avg_sil <- sapply(k, silhouette_score)
    plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
    pHAC <- fviz_nbclust(M.br.GDGT.km, kmeans, method='silhouette')
    
    W = 1000; H = 500; Save.path = "Figures/Silhouette_cluster.pdf"
    ggsave(filename = Save.path, plot = pHAC, width = W*0.026458333, height = H*0.026458333, units = "cm")
  }
  
  #### Clustering (Fuzzy C-means) ####
  cm <- cmeans(M.br.GDGT.km, 2)
  H <- get_clust_tendency(M.br.GDGT.km, n = 50)
  hopkins <- round(H$hopkins_stat, digits = 2)
  
  cm$cluster[cm$cluster == 1] <- "K-warm/arid"
  cm$cluster[cm$cluster == 2] <- "K-cold/wet"
  M.br.GDGT.km$Cluster.km <- cm$cluster
  
  #### K-means evaluation and stats ####
  Training.set <- M.br.GDGT.km
  
  SMOTE = T
  if(SMOTE == T){
    smote_result <- SMOTE(Training.set[-ncol(Training.set)], Training.set[ncol(Training.set)], K = 5, dup_size = 1)
    C1 <- smote_result$data[smote_result$data$class == "K-warm/arid",]
    C2 <- smote_result$data[smote_result$data$class == "K-cold/wet",]
    set.seed(42)
    if(nrow(C2) > nrow(C1)){C2 <- C2[sample(nrow(C2), nrow(C1)),]}
    else{C1 <- C1[sample(nrow(C1), nrow(C2)),]}
    Training.set <- rbind(C1, C2)
    names(Training.set)[names(Training.set) == "class"] <- "Cluster.km"
  }
  Training.set$Cluster.km <-  factor(Training.set$Cluster.km)
  
  Cluster.prediction.ACADB.brGDGT <- randomForest(Cluster.km ~ ., data = Training.set, ntree = 500)
  saveRDS(Cluster.prediction.ACADB.brGDGT, "Results/Cluster.prediction.ACADB.brGDGT.Rds")
  
  #### Export clusters ####  
  Meco <- full_join(rownames_to_column(Meco), rownames_to_column(M.br.GDGT.km["Cluster.km"]), by = "rowname")
  row.names(Meco) <- Meco$rowname
  Meco <- Meco[setdiff(names(Meco), "rowname")]
  
  Meco$Cluster.km[is.na(Meco$Cluster.km)] <- "K-lacustrine"
  Meco$Cluster.km <- factor(Meco$Cluster.km, c("K-warm/arid", "K-cold/wet", "K-lacustrine"), ordered = T)
  Mfull <- cbind(Msurf.mean, Mclim, subset(Meco, select = -c(Latitude, Longitude)))
  saveRDS(Meco, "Results/Meco_kmean.Rds")
  # saveRDS(Meco3, "Resultats/ACA/GDGT/Meco3_kmean.Rds")
  Mbar <- full_join(rownames_to_column(Meco), rownames_to_column(Mclim), by = join_by(rowname))
  
  M.br.GDGT.km <- full_join(rownames_to_column(M.br.GDGT), rownames_to_column(Meco["Cluster.km"]), by = join_by(rowname))
  M.br.GDGT.km <- subset(M.br.GDGT.km, select = -c(rowname))
  M.br.GDGT.km <- M.br.GDGT.km[,c(To.keep, ncol(M.br.GDGT.km))]
  
  Table.sample.by.km <- Meco %>%
    mutate(Sample.type = if_else(Sample.type == "Moss", "Soil", Sample.type)) %>%
    group_by(Cluster.km) %>%
    count(Sample.type)
  pLake.C1 <- round(Table.sample.by.km[1,3]/(Table.sample.by.km[1,3] + Table.sample.by.km[2,3]), digits = 2)*100
  pLake.C2 <- round(Table.sample.by.km[3,3]/(Table.sample.by.km[3,3] + Table.sample.by.km[4,3]), digits = 2)*100
  
  #### Plots clustering (multivariate analysis) ####
  MyColors <- c("K-warm/arid" = "darkorange", "K-cold/wet" = "royalblue", "K-lacustrine" = "#0073C2")
  p1 <- fviz_cluster(list(data = subset(M.br.GDGT.km, select = -c(Cluster.km)), cluster = M.br.GDGT.km$Cluster.km), 
                     ellipse.type = "norm", repel = F,
                     ellipse.level = 0.80,
                     # palette = "jco", 
                     geom = "point", main = "(A) Fuzzy C-Means Clustering",
                     ggtheme = theme_minimal())
  p1 <- p1 + scale_color_manual(values = MyColors) + scale_fill_manual(values = MyColors)
  
  # if(exists("Eurasia_map") == F){source("Scripts/Pollen_fun_trans.R")}
  
  #### Plot k-mean map ####
  p2 <- ggplot(Meco, aes(y = Latitude, x = Longitude, colour = Cluster.km))+
    ggtitle("(C) Clusters location")+
    geom_polygon(data = Eurasia_map, aes(x=long, y=lat, group = group), alpha = 1, fill = "grey85", color = "grey30", size = 0.3)+
    geom_polygon(data = ACA.bo.proj, aes(x=long, y=lat),colour="black", fill = NA, size = .5) +
    geom_point(size = 1.5, alpha = .4)+
    scale_color_manual(values = MyColors) + 
    coord_quickmap(xlim = c(45, 125), ylim = c(25, 63)) +
    labs(x = c("Longitude (°)"), y = c("Latitude (°)"))
  
  Mano <- summary.aov(manova(cbind(AI, MAAT, MPWAQ, MAF) ~ Cluster.km, Mbar))
  Mano.lab <- c(paste("F = ", round(Mano$` Response AI`$`F value`[1], digits = 1), "***", sep = ""),
                paste("F = ", round(Mano$` Response MAAT`$`F value`[1], digits = 1), "***", sep = ""),
                paste("F = ", round(Mano$` Response MPWAQ`$`F value`[1], digits = 1), "***", sep = ""),
                paste("F = ", round(Mano$` Response MAF`$`F value`[1], digits = 1), "***", sep = "")
  )
  Mano.lab <- data.frame(variable = c("AI", "MAAT", "MPWAQ", "MAF"), 
                         Fval = Mano.lab, x = 1, y = c(7640, 14.05, 283, 16.2), Cluster.km = "white")
  
  Mbar <- subset(Mbar, select = c(AI, MAAT, MPWAQ, MAF, Cluster.km))
  Mbar <- melt(Mbar, "Cluster.km")
  Mbar2 <- melt(M.br.GDGT.km, "Cluster.km")
  Mbar$variable <- factor(Mbar$variable, levels = c("AI", "MAAT", "MPWAQ", "MAF"), ordered = T)
  Mbar <- Mbar[!is.na(Mbar$value),]
  
  #### Plot k-mean Climate param ####
  PlotBar <- ggplot(data = Mbar, aes(x = as.factor(variable), y = value, fill = Cluster.km))+
    geom_boxplot(outliers = F)+facet_wrap(vars(as.factor(variable)), scales = "free")+
    ggtitle("(B) Climate parameters")+
    geom_label(data = Mano.lab, aes(x = 1, y = Inf, label = Fval), size = 3.5,vjust = 1.3, hjust = 0.5, show.legend = F, fill = "white")+
    scale_fill_manual(values = MyColors)+
    theme(axis.title = element_blank(), axis.text.x = element_text(size = 12), strip.placement = "none", strip.background = element_blank(), strip.text = element_blank())
  
  Plot.kmean <- ((p1 / p2 + plot_layout(heights = c(0.66,0.33))) | PlotBar) & 
    theme(legend.position = "none",  panel.border = element_rect(fill = NA), panel.grid = element_blank(), panel.background = element_blank())
  
  #### Plot k-mean FA histogram ####
  Boxplot.ACADB <- GDGT.histo.plot.surf.core(Msurf = Msurf.mean, 
                                             Mtype = Meco, Select.type = "Cluster.km", Leg.pos = c(0.28,0.7), Leg.iso = F,
                                             Iso.GDGT = F, Remove.8Me = T, Remove.7Me = T, Overlap.OK = F, Show.Plotly = F,
                                             # Color.choice = c("#0073C2", "#80cec1", "#8c510a"),
                                             Color.choice = c(MyColors[[2]], MyColors[[1]]),
                                             # Color.choice = c("#8c510a", "#80cec1"),
                                             Return.plot = T, Global.box = T,
                                             Boxplot.title = "(D) brGDGT fractional adundances", Annot.size = 4, Leg.size = 7, Leg.box = T, Ymax = 45)
  
  #### Linear plots ####
  PIR <- ggplot(Mfull, aes(y = IR, x = MBTp5Me, color = Cluster.km))+geom_point(alpha = .5)+ylim(0,1)+ xlim(0,1)+ geom_smooth(method = "lm", se = F, formula = 'y ~ x')+ stat_poly_eq(label.y = "bottom", label.x = "right", size = 2)+ xlab(expression(paste(MBT,"'"[5~Me])))+ggtitle("(E)")+ scale_color_manual(values = MyColors) + theme(legend.position = "none",  panel.border = element_rect(fill = NA), panel.grid = element_blank(), panel.background = element_blank())
  PIR2 <- ggplot(Mfull, aes(y = MBTp5Me, x = MAAT, color = Cluster.km))+geom_point(alpha = .5)+ylim(0,1)+ xlim(-10, 20)+ geom_smooth(method = "lm", se = F, formula = 'y ~ x')+ stat_poly_eq(label.y = "top", size = 2) + ylab(expression(paste(MBT,"'"[5~Me]))) +xlab("MAAT (°C)")+ggtitle("(F)")+ scale_color_manual(values = MyColors) + theme(legend.position = "none",  panel.border = element_rect(fill = NA), panel.grid = element_blank(), panel.background = element_blank())
  
  #### Export ####
  PIR <- (PIR|PIR2) 
  
  Plot.kmean <- Plot.kmean/ ((Boxplot.ACADB)/(PIR))  + plot_layout(heights = c(0.5, 0.5))
  
  W = 700*.82; H = 1400*.82; Save.path = "Figures/ACA_kmean_presentation.pdf"
  ggsave(filename = Save.path, plot = Plot.kmean, width = W*0.026458333, height = H*0.026458333, units = "cm")
} else{
  Meco <- readRDS("Results/Meco_kmean.Rds")
  Mfull <- cbind(Msurf.mean, Mclim, subset(Meco, select = -c(Latitude, Longitude)))
}

#### Machine learning ####
Machine.learning = T
if(Machine.learning == T){
  #### Import and clean data ####
  # source("Scripts/Pollen_fun_trans.R")
  M.MAAT.GDGT <- Msurf.mean[,grep("^MAAT", colnames(Msurf.mean))]
  M.br.GDGT <- Msurf.mean[,grep("^f.I", colnames(Msurf.mean))]
  To.remove <- grep("[7,8]", colnames(M.br.GDGT))
  To.keep <- seq(1, ncol(M.br.GDGT))
  To.keep <- setdiff(To.keep, To.remove)
  M.br.GDGT <- M.br.GDGT[,To.keep]
  M.br.GDGT <- M.br.GDGT[,grep("^f.I", colnames(M.br.GDGT))]
  M.br.GDGT <- M.br.GDGT[row.names(M.br.GDGT) %in% row.names(Meco[Meco$Sample.type %in% c("Soil", "Moss"),]),]
  M.MAAT.GDGT <- M.MAAT.GDGT[row.names(M.MAAT.GDGT) %in% row.names(Meco[Meco$Sample.type %in% c("Soil", "Moss"),]),]
  Mclim <- Mclim[row.names(Mclim) %in% row.names(Meco[Meco$Sample.type %in% c("Soil", "Moss"),]),]
  Meco <- Meco[row.names(Meco) %in% row.names(Meco[Meco$Sample.type %in% c("Soil", "Moss"),]),]
  
  Mcoord <- subset(Meco, select = c(Latitude, Longitude))
  Mcoord$ALTI <- NA
  names(Mcoord) <- c("LAT", "LONG", "ALTI")
  MC <- subset(Mclim, select = c(MPCOQ, MAAT, MAF, AI))
  
  Condition1 <- row.names(Meco)[which(Meco$Cluster.km %in% c("K-warm/arid"))] #n = 44, r2 = 0.47 (temp)
  Condition2 <- row.names(Meco)[which(Meco$Cluster.km %in% c("K-cold/wet"))] #n = 44, r2 = 0.47 (temp)
  MC.karid <- MC[Condition1,]
  Mcoord.karid <- Mcoord[Condition1,]
  MC.kwet <- MC[Condition2,]
  Mcoord.kwet <- Mcoord[Condition2,]
  M.br.GDGT.karid <- M.br.GDGT[Condition1,]
  M.br.GDGT.kwet <- M.br.GDGT[Condition2,]
  saveRDS(M.br.GDGT, "Results/M_brGDGT_ACADB.Rds")
  saveRDS(M.br.GDGT.karid, "Results/M_brGDGT_ACADB_karid.Rds")
  saveRDS(M.br.GDGT.kwet, "Results/M_brGDGT_ACADB_kwet.Rds")
  saveRDS(MC, "Results/M_clim_ACADB.Rds")
  saveRDS(MC.karid, "Results/M_clim_ACADB_karid.Rds")
  saveRDS(MC.kwet, "Results/M_clim_ACADB_kwet.Rds")
  
  #### Calibration surface (WAPLS, MAT, RF et BRT) ####
  Calculate.FT = T
  if(Calculate.FT == T){
    RF.brACA.bool = T
    if(RF.brACA.bool == T){
      RF.brACA <- FT.quantif(M.br.GDGT, MC, Mcoord = Mcoord, Model = "RF", Save.RDS = T,
                             Save.path = "Results/RF_brACA.csv")
      RF.brACA.karid <- FT.quantif(M.br.GDGT.karid, MC.karid, Mcoord = Mcoord.karid, Model = "RF", Save.RDS = T,
                                   Save.path = "Results/RF_brACA_karid.csv")
      RF.brACA.kwet <- FT.quantif(M.br.GDGT.kwet, MC.kwet, Mcoord = Mcoord.kwet, Model = "RF", Save.RDS = T,
                                  Save.path = "Results/RF_brACAkwet.csv")}
    
    BRT.brACA.bool = T
    if(BRT.brACA.bool == T){
      BRT.brACA <- FT.quantif(M.br.GDGT, MC, Mcoord = Mcoord, Model = "BRT", Save.RDS = T, 
                              Save.path = "Results/BRT_brACA.csv")
      BRT.brACA.karid <- FT.quantif(M.br.GDGT.karid, MC.karid, Mcoord = Mcoord.karid, Model = "BRT", Save.RDS = T,
                                    Save.path = "Results/BRT_brACA_karid.csv")
      BRT.brACA.kwet <- FT.quantif(M.br.GDGT.kwet, MC.kwet, Mcoord = Mcoord.kwet, Model = "BRT", Save.RDS = T,
                                   Save.path = "Results/BRT_brACAkwet.csv")}
  }
  else{
    BRT.brACA       <- readRDS("Results/BRT_brACA.Rds")
    BRT.brACA.karid <- readRDS("Results/BRT_brACA_karid.Rds")
    BRT.brACA.kwet  <- readRDS("Results/BRT_brACAkwet.Rds")
    RF.brACA        <- readRDS("Results/RF_brACA.Rds")
    RF.brACA.karid  <- readRDS("Results/RF_brACA_karid.Rds")
    RF.brACA.kwet   <- readRDS("Results/RF_brACAkwet.Rds")
  }
  
  #### Tableaux des résultats ####
  Table.ML = T
  if(Table.ML == T){
    Table.res.ML <- function(M1 = NULL){
      A <- M1[[1]]
      B <- M1[[2]]
      C <- M1[[3]]
      A$Calib <- "ACADB"
      B$Calib <- "k-cold/wet"
      C$Calib <- "k-warm/arid"
      A <- rbind(A,B,C)
      A$Param.clim <- gsub("[0-9]", "", row.names(A))
      A <- A[order(A$Param.clim), c(10,9,4:7)]
      A$I <- ifelse(as.numeric(A$p.val)<0.001, paste(A$I, "*", sep = ""), A$I)
      A <- A[-c(6)]
      
      
      Position.clim <- match(unique(A$Param.clim), A$Param.clim)
      A$Param.clim[setdiff(seq(1:nrow(A)),Position.clim)] <- ""
      Position.clim <- c(Position.clim, nrow(A))
      
      colnames(A) <- c("Climate parameters$^{(a)}$", "Database", "$\\mathrm{R^2}$", "RMSE", "Moran's I$^{(b)}$")
      row.names(A) <- NULL
      return(list(A,Position.clim))
      
    }
    
    T.BRT <- Table.res.ML(M1 = list(BRT.brACA$Best.Param,
                                    BRT.brACA.karid$Best.Param,
                                    BRT.brACA.kwet$Best.Param))
    T.RF <- Table.res.ML(M1 = list(RF.brACA$Best.Param,
                                   RF.brACA.karid$Best.Param,
                                   RF.brACA.kwet$Best.Param))
    Position.clim <- T.BRT[[2]] 
    T.BRT <- T.BRT[[1]]
    T.RF <- T.RF[[1]]
    
    T.all <- cbind(T.BRT, T.RF) 
    T.all <- T.all[-c(6,7)]
    T.all <- rbind(names(T.all),T.all)
    names(T.all) <- c("", "", "\\textbf{BRT}", "", "", "\\textbf{RF}", "", "")
    T.all[c(1),] <- gsub("\\.1", "", T.all[c(1),]) 
    
    Save.path <- "Results/Tableau_BRT_brGDGT.csv"
    LateX.caption <- "Statistical results of the BRT calibration for four different climate parameters (AI, MAAT, MAF and MAP) based on the XXXX brGDGT compounds."
    Path.to.create <- gsub("(.*/).*\\.csv.*","\\1", Save.path)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    write.table(T.all, file=Save.path, row.names=T, col.names=NA, sep=",", dec = ".")
    
    library(xtable)
    Save.path.tex <- gsub("\\.csv", "\\.tex", Save.path)
    Tlatex <- xtable(T.all, caption = LateX.caption, type = "latex", label = "FT_table")
    print(Tlatex, file = Save.path.tex, booktabs = T, include.rownames = F, comment = F,
          caption.placement = "top", sanitize.text.function = function(x){x},
          hline.after = c(-1,1,Position.clim[-c(1, length(Position.clim))],nrow(T.all)))
  }
  
  #### Plot residuals ####
  Plot.res.ok = T
  if(Plot.res.ok == T){
    #### Function ####
    Residual.plot <- function(MML = NULL, FT = "BRT", Msurf = NULL, 
                              Mclim = NULL, Annot = NULL, model = NULL){
      #### Residual calculation ####
      if(is.null(MML) == F){
        if(FT == "BRT"){
          Mconfu1 <- FT.core(Model.BRT = MML, MCore = Msurf, Only.fit = F, Save.RDS = F, Displot = F)
        }
        if(FT == "RF"){
          Mconfu1 <- FT.core(Model.RF = MML, MCore = Msurf, Only.fit = F, Save.RDS = F, Displot = F)
        }
        if(FT == "MAT"){
          Mconfu1 <- FT.core(Model.MAT = MML, MCore = Msurf, Only.fit = F, Save.RDS = F, Displot = F)
        }
        if(FT == "WAPLS"){
          Mconfu1 <- FT.core(Model.WAPLS = MML, MCore = Msurf, Only.fit = F, Save.RDS = F, Displot = F)
        }
        
        Mconfu1 <- do.call(rbind, lapply(Mconfu1, as.data.frame))
        Mconfu1$Sites <- row.names(Mconfu1)
        Mconfu1 <- setNames(data.frame(reshape2::melt(Mconfu1, id = "Sites")), c("Sites", "Pred.clim.lab", "Pred.clim.val"))
        
        Mconfu1.clim <- Mclim[match(levels(Mconfu1$Pred.clim.lab), names(Mclim))]
      }
      else{
        Mconfu1 <- Msurf
        Mconfu1$Sites <- row.names(Msurf)
        Mconfu1 <- setNames(data.frame(reshape2::melt(Mconfu1, id = "Sites")), c("Sites", "Pred.clim.lab", "Pred.clim.val"))
        Mconfu1$Calib <- Mconfu1$Pred.clim.lab
        Mconfu1$Pred.clim.lab <- gsub("_.*", "", Mconfu1$Pred.clim.lab)
        
        Mconfu1.clim <- Mclim[match(unique(Mconfu1$Pred.clim.lab), names(Mclim))]
      }
      
      Mconfu1.clim$Sites <- unique(Mconfu1$Sites)
      Mconfu1.clim <- setNames(data.frame(reshape2::melt(Mconfu1.clim, id = "Sites")), c("Sites", "Obs.clim.lab", "Obs.clim.val"))
      Mconfu1 <- cbind(Mconfu1, Mconfu1.clim[3])
      Mconfu1$Residuals <- (Mconfu1$Obs.clim.val - Mconfu1$Pred.clim.val)
      
      if(is.null(MML) == T){
        Mconfu1 <- subset(Mconfu1, select = -c(Pred.clim.lab))
        names(Mconfu1)[names(Mconfu1) == "Calib"] <- "Pred.clim.lab"
        Mconfu1 <- Mconfu1[c(1,3,2,4,5)]
      }
      
      #### Residual scaling ####
      Mconfu1_scaled <- Mconfu1 %>%
        dplyr::group_by(Pred.clim.lab) %>%
        mutate(Residuals.sc = (Residuals - mean(Residuals)) / sd(Residuals))
      
      #### Lim settlement ####
      Mlims <- Mconfu1[names(Mconfu1) %in% c("Pred.clim.lab", "Pred.clim.val", "Obs.clim.val")]
      Mlims <- rbind(setNames(Mlims[c(1,2)], c("Pred.clim.lab", "B")), setNames(Mlims[c(1,3)], c("Pred.clim.lab", "B")))
      Mlims <- Mlims %>%
        dplyr::group_by(Pred.clim.lab) %>%
        dplyr::summarize(
          xlim = min(B),
          ylim = max(B)
        )
      Mlims <- reshape2::melt(Mlims, id = "Pred.clim.lab")
      Mlims$Sites <- "Dummy.dot"
      Mlims <- subset(Mlims, select = -c(variable))
      Mlims$Obs.clim.val <- Mlims$value
      Mlims$Obs.clim.lab <- Mlims$Pred.clim.lab
      Mlims$Residuals <- 0
      names(Mlims)[names(Mlims) == "value"] <- "Pred.clim.val"
      
      Mlims <- Mlims[match(names(Mconfu1), names(Mlims))]
      Mconfu1 <- rbind(Mconfu1, Mlims)
      Mlims$Residuals.sc <- 0
      Mconfu1_scaled <- rbind(Mconfu1_scaled, Mlims)
      
      #### Global graphical param ####
      New.lab1 <- New.lab[match(levels(Mconfu1$Pred.clim.lab), names(New.lab))]
      
      if(all(is.na(unique(New.lab1)))){
        New.lab1 <- levels(Mconfu1$Pred.clim.lab)
        New.lab1 <- paste(sub("_", "~(", New.lab1), ")", sep = "")
        New.lab1 <- sub("_", "~", New.lab1)
        New.lab1 <- gsub("5Me", "[5*Me]", New.lab1)
        New.lab1 <- gsub("~\\[", "[", New.lab1)
      }
      
      if(is.null(model) == F){
        Ytitle <- paste(model, "predictions")
      }
      else{Ytitle <- "Predicted \nclimate parameters"}
      
      #### Plotting loops ####
      plots <- list()
      for(i in 1:length(unique(Mconfu1$Pred.clim.lab))){
        #### Param i ####
        Param <- unique(Mconfu1$Pred.clim.lab)[i]
        New.lab1.i <- New.lab1[i]
        if(is.null(Annot) == F){New.lab1.i <- paste("(", Annot[i], ")~", New.lab1.i, sep = "")}
        Mconfu.i <- Mconfu1[Mconfu1$Pred.clim.lab == Param,]
        Mconfu1_scaled.i <- Mconfu1_scaled[Mconfu1_scaled$Pred.clim.lab == Param,]
        
        if(i == 1){Ylab <- element_text()}
        if(i > 1){Ylab <- element_blank()}
        
        #### Plot (main) ####
        p <- ggplot(Mconfu.i, aes(x = Obs.clim.val, y = Pred.clim.val)) +
          geom_hex(bins = 30, aes(fill = ..count..)) +
          geom_smooth(method = "lm", se = F, linewidth = 0.7, formula = 'y ~ x') +
          scale_fill_gradientn(name = "Nb. obs.",
                               colors = c("grey90", "royalblue", "darkorange", "darkred"),  # Discrete color steps
                               values = scales::rescale(c(0, 1, 10, 100)),  # Scale the breaks
                               breaks = c(0, 1, 10, 100)) +
          geom_abline(slope = 1, intercept = 0, color = "grey20", linetype = "dashed") + 
          labs(x = NULL, 
               y = Ytitle) + 
          ggtitle(parse(text = New.lab1.i[1]))+
          stat_poly_eq(size = 3, vstep = 0.07, label.y = "bottom", label.x = "right")+
          theme(panel.background = element_rect(fill = NA, colour = 'black', linewidth = .5), panel.border = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 13),
                legend.position = "none", panel.grid = element_line(colour = "grey80", linetype = 2, linewidth = .1),
                axis.text.x = element_blank(),
                axis.title.y = Ylab,
                axis.line.x = element_blank())
        
        #### Plot (residuals) ####
        presidu <- ggplot(Mconfu1_scaled.i, aes(x = Obs.clim.val, y = Residuals)) +
          geom_point(aes(color = Residuals)) +
          geom_hline(yintercept = 0, lty = "dashed")+
          labs(x = "Observed climate parameters",
               y = "Residuals") +
          scale_color_gradient2(name = "Residual (z-scores)",
                                low = "#963327",
                                mid = "white",
                                high = "#963327",
                                midpoint = 0,
                                limit = c(-max(abs(Mconfu1_scaled.i$Residuals)), max(abs(Mconfu1_scaled.i$Residuals))))+
          theme(
            plot.background = element_blank(),  plot.margin = unit(c(0,0,0,0), 'pt'),
            axis.title.y = Ylab,
            legend.position = "none", panel.grid = element_line(colour = "grey80", linetype = 2, linewidth = .1),
            panel.background = element_rect(fill = NA, colour = 'black'), panel.border = element_blank(),
            strip.text = element_blank(), strip.background = element_blank())
        
        plots[[i]] <- (p / presidu)+ plot_layout(nrow = 2, heights = c(.75,.25))
      }
      
      plots <- wrap_plots(plots, nrow = 1)
      return(plots)
    }
    
    #### Plot all values BRT ####
    p.all.BRT.full <- Residual.plot(MML = BRT.brACA, Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("A", seq(1:4), sep = ''), model = c("BRT (ACADB)"))
    p.karid.BRT.full <- Residual.plot(MML = BRT.brACA.karid, Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("B", seq(1:4), sep = ''), model = c("BRT (K-arid)"))
    p.kwet.BRT.full <- Residual.plot(MML = BRT.brACA.kwet, Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("C", seq(1:4), sep = ''), model = c("BRT (K-wet)"))
    p.all.RF.full <- Residual.plot(MML = RF.brACA, Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("D", seq(1:4), sep = ''), model = c("RF (ACADB)"), FT = "RF")
    p.karid.RF.full <- Residual.plot(MML = RF.brACA.karid, Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("E", seq(1:4), sep = ''), model = c("RF (K-arid)"), FT = "RF")
    p.kwet.RF.full <- Residual.plot(MML = RF.brACA.kwet, Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("F", seq(1:4), sep = ''), model = c("RF (K-wet)"), FT = "RF")
    
    p.full <- wrap_plots(p.all.BRT.full, p.karid.BRT.full, p.kwet.BRT.full, p.all.RF.full, p.karid.RF.full, p.kwet.RF.full, nrow = 6)
    
    H = 2000; W = 1000; Save.plot ="Figures/Obs_pred_scatterplot_ML.pdf"
    ggsave(filename = Save.plot, p.full, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Plot MAAT for all calibs ####
    p.all.BRT <- Residual.plot(MML = BRT.brACA[c(2,5)], Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("A", 1, sep = ''), model = c("BRT (ACADB)"))
    p.karid.BRT <- Residual.plot(MML = BRT.brACA.karid[c(2,5)], Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("A", 2, sep = ''), model = c("BRT (K-arid)"))
    p.kwet.BRT <- Residual.plot(MML = BRT.brACA.kwet[c(2,5)], Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("A", 3, sep = ''), model = c("BRT (K-wet)"))
    
    p.all.RF <- Residual.plot(MML = RF.brACA[c(2,5)], Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("B", 1, sep = ''), model = c("RF (ACADB)"), FT = "RF")
    p.karid.RF <- Residual.plot(MML = RF.brACA.karid[c(2,5)], Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("B", 2, sep = ''), model = c("RF (K-arid)"), FT = "RF")
    p.kwet.RF <- Residual.plot(MML = RF.brACA.kwet[c(2,5)], Msurf = M.br.GDGT, Mclim = Mclim, Annot = paste("B", 3, sep = ''), model = c("RF (K-wet)"), FT = "RF")
    
    p.calib <- Residual.plot(MML = NULL, Msurf = M.MAAT.GDGT[c(2,8,9)], Mclim = Mclim, Annot = paste("C", seq(1:4), sep = ''))
    
    k <- list(p.all.BRT, p.karid.BRT, p.kwet.BRT, 
              p.all.RF, p.karid.RF, p.kwet.RF, 
              p.calib[[1]], p.calib[[2]], p.calib[[3]])
    p <- wrap_plots(k, nrow = 3)
    H = 1100; W = 900; Save.plot ="Figures/Obs_pred_scatterplot_ML_MAAT.pdf"
    ggsave(filename = Save.plot, p, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
  }
  
  #### brGDGT contributions ####
  GDGT.contrib = T
  if(GDGT.contrib == T){
    #### Import contributions ####
    BRT.brACA.c <- data.frame(read.csv(file="Results/BRT_brACA_Taxa_importance.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
    BRT.brACA.karid.c <- data.frame(read.csv(file="Results/BRT_brACA_karid_Taxa_importance.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
    BRT.brACA.kwet.c <- data.frame(read.csv(file="Results/BRT_brACAkwet_Taxa_importance.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
    RF.brACA.c <- data.frame(read.csv(file="Results/RF_brACA_Taxa_importance.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
    RF.brACA.karid.c <- data.frame(read.csv(file="Results/RF_brACA_karid_Taxa_importance.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
    RF.brACA.kwet.c <- data.frame(read.csv(file="Results/RF_brACAkwet_Taxa_importance.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
    
    #### Clean data ####
    BRT.brACA.c$Models <-"BRT ACADB"
    BRT.brACA.karid.c$Models <-"BRT K-warm/arid"
    BRT.brACA.kwet.c$Models <-"BRT K-cold/wet"
    RF.brACA.c$Models <-"RF ACADB"
    RF.brACA.karid.c$Models <-"RF K-warm/arid"
    RF.brACA.kwet.c$Models <-"RF K-cold/wet"
    
    BRT.brACA.c$Compounds <- row.names(BRT.brACA.c)
    BRT.brACA.karid.c$Compounds <- row.names(BRT.brACA.karid.c)
    BRT.brACA.kwet.c$Compounds <- row.names(BRT.brACA.kwet.c)
    RF.brACA.c$Compounds <- row.names(RF.brACA.c)
    RF.brACA.karid.c$Compounds <- row.names(RF.brACA.karid.c)
    RF.brACA.kwet.c$Compounds <- row.names(RF.brACA.kwet.c)
    
    BRT.brACA.c <- melt(BRT.brACA.c, id = c("Models", "Compounds"))
    BRT.brACA.karid.c <- melt(BRT.brACA.karid.c, id = c("Models", "Compounds"))
    BRT.brACA.kwet.c <- melt(BRT.brACA.kwet.c, id = c("Models", "Compounds"))
    RF.brACA.c <- melt(RF.brACA.c, id = c("Models", "Compounds"))
    RF.brACA.karid.c <- melt(RF.brACA.karid.c, id = c("Models", "Compounds"))
    RF.brACA.kwet.c <- melt(RF.brACA.kwet.c, id = c("Models", "Compounds"))
    
    
    df_list <- list(BRT.brACA.c, BRT.brACA.karid.c, BRT.brACA.kwet.c, RF.brACA.c, RF.brACA.karid.c, RF.brACA.kwet.c)
    MTB <- Reduce(function(x, y) merge(x, y, all = T, sort = F), df_list)
    
    MTB$Calibration <- gsub(".* ","", MTB$Models)
    MTB$Models <- gsub(" .*","", MTB$Models)
    MTB$Models[MTB$Models == "BRT"] <- "(A) BRT"
    MTB$Models[MTB$Models == "RF"] <- "(B) RF"
    MTB$Compounds <- gsub("_5Me", "", MTB$Compounds)
    MTB$Compounds <- gsub("_6Me", "'", MTB$Compounds)
    MTB$Compounds <- gsub("f.", "", MTB$Compounds)
    
    TAB <- MTB %>%
      group_by(Compounds) %>%
      summarise_at(vars(value), mean)
    TAB <- TAB[order(TAB$value, decreasing = T),]
    
    #### Plots + export ####
    MTB$Compounds <- factor(MTB$Compounds, levels = TAB$Compounds, ordered = T)
    Pcontrib <- ggplot(data = MTB, aes(x = Compounds, y = value, fill = Calibration, group = Calibration))+
      geom_vline(xintercept = seq(0, nlevels(MTB$Compounds), by = 1), colour = "grey80", size = 0.1) +
      geom_hline(yintercept = c(0,5,10,20,30), colour = "grey60", size = 0.2, linetype = "dashed") +
      geom_bar(stat = "identity", position = "dodge")+
      coord_polar()+
      
      ylim(-10,max(MTB$value)) +
      scale_fill_manual(values = c("ACADB" = "bisque3", "K-cold/wet" = "royalblue", "K-warm/arid" = "darkorange"))+
      facet_grid(variable ~ Models, scale = "free")+
      theme_minimal()+
      theme(legend.position = "bottom", axis.title = element_blank(), legend.background = element_blank(),
            panel.border = element_blank(), legend.key = element_blank(), axis.ticks = element_blank(),
            axis.text.y = element_blank(), strip.text = element_text(size = 12), panel.grid  = element_blank()
      )
    H = 900; W = H*.5; Save.plot ="Figures/BRT_RF_contribs.pdf"
    ggsave(filename = Save.plot, Pcontrib, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
  }
  
  #### Clim param independance ####
  Mat.cor.clim = T
  if(Mat.cor.clim == T){
    Keep.clim <- c("AI", "MAAT", "MAF", "MPCOQ")
    Mmatcor.ACADB <- cbind(Meco, Mclim)
    Mmatcor.ACADB.karid <- Mmatcor.ACADB[Mmatcor.ACADB$Cluster.km == "K-warm/arid", Keep.clim]
    Mmatcor.ACADB.kwet <- Mmatcor.ACADB[Mmatcor.ACADB$Cluster.km == "K-cold/wet", Keep.clim]
    Mmatcor.ACADB <- Mmatcor.ACADB[Keep.clim]
    
    H = 300*.8; W = 800*.8; Save.matcor.trait = "Figures/Clim_param_indep.pdf"
    pdf(file = Save.matcor.trait, width = W*0.01041666666667, height = H*0.01041666666667)
    par(mfrow = c(1,3))
    A <- Mat.corel(Mmatcor.ACADB, Mmatcor.ACADB, Display.pval = "pch", Label = F, I.confiance = 0.99, Label.simple = T, Disp.R = "number", Display = "lower", Title = "(A) ACADB")
    B <- Mat.corel(Mmatcor.ACADB.kwet, Mmatcor.ACADB.kwet, Display.pval = "pch", Label = F, I.confiance = 0.99, Label.simple = T, Disp.R = "number", Display = "lower", Title = "(B) K-cold/wet")
    C <- Mat.corel(Mmatcor.ACADB.karid, Mmatcor.ACADB.karid, Display.pval = "pch", Label = F, I.confiance = 0.99, Label.simple = T, Disp.R = "number", Display = "lower", Title = "(C) K-warm/arid")
    dev.off()
    
  }
}

#### Paleo applications ####
Paleo = T
if(Paleo == T){
  #### Import data ####
  GDGT.Van       <- readRDS("Import/Paleo/Vanevan.Rds")
  GDGT.XRD       <- readRDS("Import/Paleo/XRD.Rds")
  GDGT.NRX       <- readRDS("Import/Paleo/NRX.Rds")
  Actual.val     <- readRDS("Import/Paleo/Core_metadata.Rds")
  BRT.brACA      <- readRDS("Results/BRT_brACA.Rds")
  BRT.brACA.karid<- readRDS("Results/BRT_brACA_karid.Rds")
  BRT.brACA.kwet <- readRDS("Results/BRT_brACAkwet.Rds")
  M.br.GDGT      <- readRDS("Results/M_brGDGT_ACADB.Rds")
  GDGT.Van.conv  <- GDGT.Van[which(names(GDGT.Van) %in%  names(M.br.GDGT))]
  GDGT.XRD.conv  <- GDGT.XRD[which(names(GDGT.XRD) %in% names(M.br.GDGT))]
  GDGT.NRX.conv  <- GDGT.NRX[which(names(GDGT.NRX) %in% names(M.br.GDGT))]
  
  #### Histogramme f[br-GDGT] ####
  br.GDGT.ATM = T
  if(br.GDGT.ATM == T){
    GDGT.histo.plot.surf.core(Mcore = GDGT.Van, Iso.GDGT = T, Remove.8Me = F, Remove.7Me = F, W = 1600, H = 700, Save.path = "Figures/Van_Hist_br-GDGT.pdf")
    GDGT.histo.plot.surf.core(Mcore = GDGT.XRD, Iso.GDGT = T, Remove.8Me = F, Remove.7Me = F, W = 1600, H = 700, Save.path = "Figures/NRX_Hist_br-GDGT.pdf")
    GDGT.histo.plot.surf.core(Mcore = GDGT.NRX, Iso.GDGT = T, Remove.8Me = F, Remove.7Me = F, W = 1600, H = 700, Save.path = "Figures/XRD_Hist_br-GDGT.pdf")
  }
  
  #### BRT predictions ####
  plot.Vanevan.ft.COST = T
  if(plot.Vanevan.ft.COST == T){
    Vanevan.brACA <- FT.core(Model.BRT = BRT.brACA,
                           MCore = GDGT.Van.conv, MAge = GDGT.Van$Age, GDGT = T,
                           LakeName = "Vanevan", Only.fit = T, Save.RDS = T, Displot = F, 
                           Save.path = "Results/Vanevan.csv")
    
    Vanevan.brACA.karid <- FT.core(Model.BRT = BRT.brACA.karid,
                                 MCore = GDGT.Van.conv, MAge = GDGT.Van$Age, GDGT = T,
                                 LakeName = "Vanevan", Only.fit = T, Save.RDS = T, Displot = F,
                                 Save.path = "Results/Vanevan.csv")
    
    Vanevan.brACA.kwet <- FT.core(Model.BRT = BRT.brACA.kwet,
                                MCore = GDGT.Van.conv, MAge = GDGT.Van$Age, GDGT = T,
                                LakeName = "Vanevan", Only.fit = T, Save.RDS = T, Displot = F,
                                Save.path = "Results/Vanevan.csv")
    
    NRX.brACA <- FT.core(Model.BRT = BRT.brACA,
                         MCore = GDGT.NRX.conv, MAge = GDGT.NRX$Age,
                         LakeName = "NRX", Only.fit = T, Save.RDS = T, Displot = F, GDGT = T,
                         Save.path = "Results/NRX.csv")
    
    NRX.brACA.karid <- FT.core(Model.BRT = BRT.brACA.karid,
                               MCore = GDGT.NRX.conv, MAge = GDGT.NRX$Age,
                               LakeName = "NRX", Only.fit = T, Save.RDS = T, Displot = F, GDGT = T,
                               Save.path = "Results/NRX.csv")
    
    NRX.brACA.kwet <- FT.core(Model.BRT = BRT.brACA.kwet,
                              MCore = GDGT.NRX.conv, MAge = GDGT.NRX$Age,
                              LakeName = "NRX", Only.fit = T, Save.RDS = T, Displot = F, GDGT = T,
                              Save.path = "Results/NRX.csv")
    
    XRD.brACA <- FT.core(Model.BRT = BRT.brACA,
                         MCore = GDGT.XRD.conv, MAge = GDGT.XRD$Age,
                         LakeName = "XRD", Only.fit = T, Save.RDS = T, Displot = F, GDGT = T,
                         Save.path = "Results/XRD.csv")
    
    XRD.brACA.karid <- FT.core(Model.BRT = BRT.brACA.karid,
                               MCore = GDGT.XRD.conv, MAge = GDGT.XRD$Age,
                               LakeName = "XRD", Only.fit = T, Save.RDS = T, Displot = F, GDGT = T,
                               Save.path = "Results/XRD.csv")
    
    XRD.brACA.kwet <- FT.core(Model.BRT = BRT.brACA.kwet,
                              MCore = GDGT.XRD.conv, MAge = GDGT.XRD$Age,
                              LakeName = "XRD", Only.fit = T, Save.RDS = T, Displot = F, GDGT = T,
                              Save.path = "Results/XRD.csv")
    
    
  }
  
  #### Cluster prediction ####
  Cm.predict.Van = T
  if(Cm.predict.Van == T){
    #### Import + settings ####
    Cluster.prediction.ACADB.brGDGT <- readRDS("Results/Cluster.prediction.ACADB.brGDGT.Rds")
    Vanevan.brACA       <- readRDS("Results/Vanevan_brACA.Rds")[[2]]
    Vanevan.brACA.karid <- readRDS("Results/Vanevan_karid.Rds")[[2]]
    Vanevan.brACA.kwet  <- readRDS("Results/Vanevan_kwet.Rds")[[2]]
    NRX.brACA         <- readRDS("Results/NRX_brACA.Rds")[[2]]
    NRX.brACA.karid   <- readRDS("Results/NRX_karid.Rds")[[2]]
    NRX.brACA.kwet    <- readRDS("Results/NRX_kwet.Rds")[[2]]
    XRD.brACA         <- readRDS("Results/XRD_brACA.Rds")[[2]]
    XRD.brACA.karid   <- readRDS("Results/XRD_karid.Rds")[[2]]
    XRD.brACA.kwet    <- readRDS("Results/XRD_kwet.Rds")[[2]]
    
    #### MAAT plot ####
    pBRT.Vanevan <- Combine.ML.cluster(
      List.models = list(M1 = Vanevan.brACA, M2 = Vanevan.brACA.karid, M3 = Vanevan.brACA.kwet),
      Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
      Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
      GDGT.paleo = GDGT.Van,
      Surf.val = Actual.val$MAAT[row.names(Actual.val) == "Vanevan"],
      Compare.curve = c("MAAT_mr_DJ", "MAAT_soil_Naaf", "MAAT_LSun"), Core.name = "(B) Vanevan (Armenia)",
      Plot.y = "Age", Plot.y.lab = NULL, Param.clim = "MAAT", Cluster.prob = "K-warm/arid",
      Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T,
      Save.path = "Results/Vanevan_brACA_combined.Rds")
    
    pBRT.NRX <- Combine.ML.cluster(
      List.models = list(M1 = NRX.brACA, M2 = NRX.brACA.karid, M3 = NRX.brACA.kwet),
      Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
      Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
      GDGT.paleo = GDGT.NRX, Time.lim = c(0, 7200), Cluster.prob = "K-warm/arid",
      Surf.val = Actual.val$MAAT[row.names(Actual.val) == "NRX"],
      Compare.curve = c("MAAT_mr_DJ", "MAAT_DJ_5Me", "MAAT_NMSDB_mr5"), Core.name = "(A) NRX (Altai)",
      Plot.y = "Age",  Plot.y.lab = NULL, Param.clim = "MAAT", 
      Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T,
      Save.path = "Results/NRX_brACA_combined.Rds"
    )
    
    pBRT.XRD <- Combine.ML.cluster(
      List.models = list(M1 = XRD.brACA, M2 = XRD.brACA.karid, M3 = XRD.brACA.kwet),
      Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
      Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
      GDGT.paleo = GDGT.XRD, Time.lim = c(0, 7200), 
      Surf.val = Actual.val$MAAT[row.names(Actual.val) == "XRD"],
      Compare.curve = c("MAAT_mr_DJ", "MAAT_DJ_5Me", "MAAT_NMSDB_mr5"), Core.name = "(C) XRD (Qaidam)",
      Plot.y = "Age", Plot.y.lab = "Age (yr cal BP)", Param.clim = "MAAT", Cluster.prob = "K-cold/wet",
      Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T,
      Save.path = "Results/XRD_brACA_combined.Rds"
    )
    
    pCal <- pBRT.NRX[[1]]   / pBRT.NRX[[2]]   / pBRT.NRX[[3]] /
      pBRT.Vanevan[[1]] / pBRT.Vanevan[[2]] / pBRT.Vanevan[[3]] /
      pBRT.XRD[[1]]   / pBRT.XRD[[2]]   / pBRT.XRD[[3]]+
      plot_layout(guides = 'collect', heights = c(1,3,3,1,3,3,1,3,3))&
      theme(legend.position = "bottom")
    H = 1300; W = 600; Save.plot = "Figures/BRT_valid_China.pdf"
    ggsave(filename = Save.plot, pCal, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### AI plot ####
    pBRT.Vanevan.AI <- Combine.ML.cluster(
      List.models = list(M1 = Vanevan.brACA, M2 = Vanevan.brACA.karid, M3 = Vanevan.brACA.kwet),
      Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
      Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
      GDGT.paleo = GDGT.Van,
      Surf.val = Actual.val$AI[row.names(Actual.val) == "Vanevan"],
      # Compare.curve = c("MAAT_mr_DJ", "MAAT_soil_Naaf", "MAAT_LSun"), Core.name = "(B) Vanevan (Armenia)",
      Plot.y = "Age", Plot.y.lab = NULL, Param.clim = "AI", Cluster.prob = "K-warm/arid",
      Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T)
    
    pBRT.NRX.AI <- Combine.ML.cluster(
      List.models = list(M1 = NRX.brACA, M2 = NRX.brACA.karid, M3 = NRX.brACA.kwet),
      Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
      Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
      GDGT.paleo = GDGT.NRX, Time.lim = c(0, 7200), Cluster.prob = "K-warm/arid",
      Surf.val = Actual.val$AI[row.names(Actual.val) == "NRX"],
      # Compare.curve = c("MAAT_mr_DJ", "MAAT_DJ_5Me", "MAAT_NMSDB_mr5"),
      Core.name = "(A) NRX (Altai)",
      Plot.y = "Age",  Plot.y.lab = NULL, Param.clim = "AI", 
      Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T)
    
    pBRT.XRD.AI <- Combine.ML.cluster(
      List.models = list(M1 = XRD.brACA, M2 = XRD.brACA.karid, M3 = XRD.brACA.kwet),
      Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
      Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
      GDGT.paleo = GDGT.XRD, Time.lim = c(0, 7200), 
      Surf.val = Actual.val$AI[row.names(Actual.val) == "XRD"],
      # Compare.curve = c("MAAT_mr_DJ", "MAAT_DJ_5Me", "MAAT_NMSDB_mr5"), 
      Core.name = "(C) XRD (Qaidam)",
      Plot.y = "Age", Plot.y.lab = "Age (yr cal BP)", Param.clim = "AI", Cluster.prob = "K-cold/wet",
      Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T)
    
    pCal.AI <- pBRT.NRX.AI[[1]]   / pBRT.NRX.AI[[2]]   /
      pBRT.Vanevan.AI[[1]] / pBRT.Vanevan.AI[[2]] /
      pBRT.XRD.AI[[1]]   / pBRT.XRD.AI[[2]] +
      plot_layout(guides = 'collect', heights = c(1,3,1,3,1,3))&
      theme(legend.position = "bottom")
    H = 1300; W = 600; Save.plot = "Figures/BRT_valid_China_AI.pdf"
    ggsave(filename = Save.plot, pCal.AI, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### MAF plot ####
    pBRT.Vanevan.MAF <- Combine.ML.cluster(
      List.models = list(M1 = Vanevan.brACA, M2 = Vanevan.brACA.karid, M3 = Vanevan.brACA.kwet),
      Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
      Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
      GDGT.paleo = GDGT.Van,
      Surf.val = Actual.val$MAF[row.names(Actual.val) == "Vanevan"],
      Compare.curve = c("MAF_MSosa"),
      Core.name = "(B) Vanevan (Armenia)",
      Plot.y = "Age", Plot.y.lab = NULL, Param.clim = "MAF", Cluster.prob = "K-warm/arid",
      Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T)
    
    pBRT.NRX.MAF <- Combine.ML.cluster(
      List.models = list(M1 = NRX.brACA, M2 = NRX.brACA.karid, M3 = NRX.brACA.kwet),
      Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
      Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
      GDGT.paleo = GDGT.NRX, Time.lim = c(0, 7200), Cluster.prob = "K-warm/arid",
      Surf.val = Actual.val$MAF[row.names(Actual.val) == "NRX"],
      Compare.curve = c("MAF_MSosa"),
      Core.name = "(A) NRX (Altai)", Plot.y = "Age",  Plot.y.lab = NULL, Param.clim = "MAF", 
      Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T)
    
    pBRT.XRD.MAF <- Combine.ML.cluster(
      List.models = list(M1 = XRD.brACA, M2 = XRD.brACA.karid, M3 = XRD.brACA.kwet),
      Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
      Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
      GDGT.paleo = GDGT.XRD, Time.lim = c(0, 7200), 
      Surf.val = Actual.val$MAF[row.names(Actual.val) == "XRD"],
      Compare.curve = c("MAF_MSosa"), Core.name = "(C) XRD (Qaidam)",
      Plot.y = "Age", Plot.y.lab = "Age (yr cal BP)", Param.clim = "MAF", Cluster.prob = "K-cold/wet",
      Facet = F, Only.best = F, Show.proba = T, Highlight.combined = T)
    
    pCal.MAF <- pBRT.NRX.MAF[[1]]   / pBRT.NRX.MAF[[2]]   / pBRT.NRX.MAF[[3]] /
      pBRT.Vanevan.MAF[[1]] / pBRT.Vanevan.MAF[[2]] / pBRT.Vanevan.MAF[[3]] /
      pBRT.XRD.MAF[[1]]   / pBRT.XRD.MAF[[2]]   / pBRT.XRD.MAF[[3]]+
      plot_layout(guides = 'collect', heights = c(1,3,3,1,3,3,1,3,3))&
      theme(legend.position = "bottom")
    H = 1300; W = 600; Save.plot = "Figures/BRT_valid_China_MAF.pdf"
    ggsave(filename = Save.plot, pCal.MAF, width = W*0.026458333, height = H*0.026458333, units = "cm")    
  }
  
  #### Test (randomTF) #### 
  test.randomTF = F
  if(test.randomTF == T){
    Vanevan.rTF.ACADB <- Plot.randomTF(MPsurf = M.br.GDGT, MPpaleo = GDGT.Van.conv,  Mclim = MC, Database = "ACADB", Lake = "Vanevan", 
                                     Save.path = "Results/Vanevan_randomFT_ACADB.Rds", 
                                     H = 500, W = 2000, Save.plot = "Figures/Vanevan_randomTF_ACADB.pdf")
    
    Vanevan.rTF.ACADB.comb <- Plot.randomTF(MPsurf = M.br.GDGT.kwet, MPpaleo = GDGT.Van.conv,  Mclim = MC.kwet, Database = "ACADB-combined", Lake = "Vanevan",
                                          Save.path = "Results/Vanevan_randomFT_ACADB_comb.Rds", 
                                          H = 500, W = 2000, Save.plot = "Figures/Vanevan_randomTF_ACADB_comb.pdf")
  
    XRD.rTF.ACADB <- Plot.randomTF(MPsurf = M.br.GDGT, MPpaleo = GDGT.XRD.conv,  Mclim = MC, Database = "ACADB", Lake = "XRD", 
                                   Save.path = "Results/XRD_randomFT_ACADB.Rds", 
                                   H = 500, W = 2000, Save.plot = "Figures/XRD_randomTF_ACADB.pdf")
    
    XRD.rTF.ACADB.comb <- Plot.randomTF(MPsurf = M.br.GDGT.karid, MPpaleo = GDGT.XRD.conv,  Mclim = MC.karid, Database = "ACADB-combined", Lake = "XRD",
                                        Save.path = "Results/XRD_randomFT_ACADB_comb.Rds",
                                        H = 500, W = 2000, Save.plot = "Figures/XRD_randomTF_ACADB_comb.pdf")
  
    NRX.rTF.ACADB <- Plot.randomTF(MPsurf = M.br.GDGT, MPpaleo = GDGT.NRX.conv,  Mclim = MC, Database = "ACADB", Lake = "NRX",
                                   Save.path = "Results/NRX_randomFT_ACADB.Rds",
                                   H = 500, W = 2000, Save.plot = "Figures/NRX_randomTF_ACADB.pdf")
    
    NRX.rTF.ACADB.comb <- Plot.randomTF(MPsurf = M.br.GDGT.kwet, MPpaleo = GDGT.NRX.conv,  Mclim = MC.kwet, Database = "ACADB-combined",
                                        Save.path = "Resultats/NRX_randomFT_ACADB_comb.Rds",
                                        H = 500, W = 2000, Save.plot = "Figures/NRX_randomTF_ACADB_comb.pdf")}
  
  #### Table randomTF() ####
  Full.table = F
  if(Full.table == T){
    XRD.rTF.ACADB          <- readRDS("Results/XRD_randomFT_ACADB.Rds")
    XRD.rTF.ACADB.comb     <- readRDS("Results/XRD_randomFT_ACADB_comb.Rds")
    NRX.rTF.ACADB          <- readRDS("Results/NRX_randomFT_ACADB.Rds")
    NRX.rTF.ACADB.comb     <- readRDS("Results/NRX_randomFT_ACADB_comb.Rds")
    Vanevan.rTF.ACADB      <- readRDS("Results/Vanevan_randomFT_ACADB.Rds")
    Vanevan.rTF.ACADB.comb <- readRDS("Results/Vanevan_randomFT_ACADB_comb.Rds")
    
    XRD.rTF.ACADB$lake <- "XRD" 
    XRD.rTF.ACADB.comb$lake <- "XRD" 
    NRX.rTF.ACADB$lake <- "NRX" 
    NRX.rTF.ACADB.comb$lake <- "NRX" 
    Vanevan.rTF.ACADB$lake <- "Vanevan" 
    Vanevan.rTF.ACADB.comb$lake <- "Vanevan" 
    
    Table.rTF.full <- rbind(XRD.rTF.ACADB,
                            XRD.rTF.ACADB.comb,
                            NRX.rTF.ACADB,
                            NRX.rTF.ACADB.comb,
                            Vanevan.rTF.ACADB,
                            Vanevan.rTF.ACADB.comb) 
    
    Table.rTF.full <- Table.rTF.full[c(ncol(Table.rTF.full),1:(ncol(Table.rTF.full)-1))]
    Table.rTF.full <- Table.rTF.full[Table.rTF.full$Model == "BRT",]
    Table.rTF.full <- Table.rTF.full[order(Table.rTF.full$lake, Table.rTF.full$Database),]
    
    Position.clim <- match(unique(Table.rTF.full$lake), Table.rTF.full$lake)
    Table.rTF.full$lake[setdiff(seq(1:nrow(Table.rTF.full)),Position.clim)] <- ""
    names(Table.rTF.full)[ncol(Table.rTF.full)] <- "$\\mathrm{PCA_{1} (var. \\%)^{(b)}}$"
    names(Table.rTF.full)[ncol(Table.rTF.full)-1] <- "Threshold$^{(a)}$"
    names(Table.rTF.full)[1] <- "Archive name"
    
    Save.path.tex <- "Results/Results_randomTF.tex"
    LateX.caption <- "Statistical results for the significance test conducted with \texttt{randomTF()}, designed to assess the reliability of the BRT climate reconstruction based on various calibration sets (ACADB and the dominant cluster for the combined models), were compared to a transfer function trained on randomly selected calibration sets, using 999 permutations. The three archives (NRX, XRD and Vanevan) are tested."
    Tlatex <- xtable(Table.rTF.full, caption = LateX.caption, type = "latex", label = "Table_randomTF")
    print(Tlatex, file = Save.path.tex, booktabs = T, include.rownames = F, comment = F,
          caption.placement = "top", sanitize.text.function = function(x){x},
          hline.after = c(-1,0,Position.clim[-c(1)]-1,nrow(Table.rTF.full)))
  }
}

