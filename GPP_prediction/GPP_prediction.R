library(tidyverse)
library(randomForest)
library(VSURF)
library(tictoc)
library(caret)
library(spm)
library(viridis)
library("gplots")
library(UBL)
library(ggpubr)
library(doParallel)
library(lme4)
library(MuMIn)
library(latex2exp)
library(lubridate)
library(CAST)

# input dataset ----

predictors.v1 <- c("arvi_mean",
                   "dvi_mean",
                   "gemi_mean",
                   "gndvi_mean",
                   "ipvi_mean",
                   "ireci_mean",
                   "mcari_mean",
                   "msavi2_mean",
                   "msavi_mean",
                   "mtci_mean",
                   "ndi45_mean",
                   "pssra_mean",
                   "pvi_mean",
                   "rvi_mean",
                   "s2rep_mean",
                   "savi_mean",
                   "tndvi_mean",
                   "tsavi_mean",
                   "wdvi_mean",
                   "cir_mean",
                   "cig_mean",
                   "ndvi_mean",
                   "nirv_mean",
                   "B1_mean", 
                   "B2_mean", 
                   "B3_mean", 
                   "B4_mean", 
                   "B5_mean", 
                   "B6_mean", 
                   "B7_mean", 
                   "B8_mean", 
                   "B8A_mean", 
                   "B9_mean", 
                   "B11_mean", 
                   "B12_mean")

sentinel2.efps <- read_csv("~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_iqr.csv")
sentinel2.efps %>% nrow()

sentinel2.efps <- sentinel2.efps %>% 
  mutate(igbp = as_factor(igbp))

# data structure ---

# number of observations 
sentinel2.efps %>% 
  group_by(site_name) %>% 
  summarise(n = length(site_name)) %>% 
  print(n = 58)

# number of sites per vegetation type
sentinel2.efps %>% 
  group_by(igbp) %>% 
  summarise(n = length(unique(site_name)))


# number of observations per vegetation type

sentinel2.efps %>% 
  group_by(igbp) %>% 
  summarise(n = length(igbp))

# sentinel2.efps %>% 
#   ggplot(aes(x = igbp)) +
#   geom_bar() +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 550)) +
#   ylab("Number of observations") +
#   xlab("Vegetation type") +
#   theme_bw(base_size = 13)

# ggsave("~/dpabon/results_EFPs_ML/plots/gpp_sentinel2_imbalanced_vt.png", width = 7, height = 5)

# ggsave("~/dpabon/results_EFPs_ML/plots/gpp_sentinel2_imbalanced_vt.pdf", width = 7, height = 5)

# dataset balancing ----


sentinel2.efps %>% 
  group_by(igbp) %>% 
  summarise(n = length(igbp))
# per vegatation type (Kind of artesanal balancing!!)

# undersampling

# OSH is the category with less observations

# sentinel2.efps.undersampling <- sentinel2.efps %>%
#   group_by(igbp) %>%
#   sample_n(158, replace = F) %>%
#   ungroup()
# 
# write_csv(sentinel2.efps.undersampling, "~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_under.csv")

sentinel2.efps.undersampling <- read_csv("~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_under.csv")


#oversampling

# each category is filled with random samples of the same category until complete the number of observations of the higher category

# categor <- sentinel2.efps %>%
#   group_by(igbp) %>%
#   summarise(n = length(igbp))
# 
# 
# temp.list <- list()
# for (i in 1:nrow(categor)){
#   temp <- sentinel2.efps %>%
#     filter(igbp == categor$igbp[i])
#   if(max(categor$n) != nrow(temp)){
#     n.to.replace <- max(categor$n)-nrow(temp)
#     if(n.to.replace < nrow(temp)){
#       replace.tmp <- F
#     }else{
#       replace.tmp <- T
#     }
#     to.fill <- temp %>%
#       sample_n(n.to.replace, replace = replace.tmp)
#     temp.list[[i]] <- bind_rows(temp, to.fill)
#   }else{
#     temp.list[[i]] <- sentinel2.efps %>%
#       filter(igbp == categor$igbp[i])
#   }
# }
# 
# temp.list
# 
# sentinel2.efps.oversampling <- bind_rows(temp.list)
# write_csv(sentinel2.efps.oversampling, "~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_over.csv")

sentinel2.efps.oversampling <- read_csv("~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_over.csv")

# Smote Regression balancing (A fancy one, doesn't base on categories instead in the occurence of observations) (Torgo et al., 2013)

# sentinel2.efps.temp <- as.data.frame(sentinel2.efps)
# sentinel2.efps.temp$site_name <- as.factor(sentinel2.efps.temp$site_name)
# sentinel2.efps.temp$anomaly <- NULL
# sentinel2.efps.temp$manual_outlier <- NULL
# sentinel2.efps.temp$date_index <- as.factor(sentinel2.efps$date_index)
# set.seed(10)
# sentinel2.efps.smoter <- as_tibble(SmoteRegress(form = GPP.day.smooth ~ ., dat = sentinel2.efps.temp, dist = "HEOM", k = 5))
# sentinel2.efps.smoter$date_index <- as.Date(sentinel2.efps.smoter$date_index)
# 
# View(sentinel2.efps.smoter)
# write_csv(sentinel2.efps.smoter, "~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_smoter.csv")

sentinel2.efps.smoter <- read_csv("~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_smoter.csv")

# SMOTER GPPsat
# 
# sentinel2.efps.temp <- sentinel2.efps %>%
# select(gppsat, igbp, predictors.v1) %>%
# as.data.frame()
# 
# 
# sentinel2.efps.smoter.gppsat <- as_tibble(SmoteRegress(form = gppsat ~ ., dat = sentinel2.efps.temp, dist = "HEOM", k = 5))
# 
# write_csv(sentinel2.efps.smoter.gppsat, "~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_smoter_gppsat.csv")
# 

# sentinel2.efps.smoter.gppsat <- read_csv("~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_smoter_gppsat.csv")

# sentinel2.efps.smoter.gppsat <- inner_join(sentinel2.efps.smoter.gppsat, sentinel2.efps)

# kNDVI estimation ----

sentinel2.efps$kndvi_mean <- kNDVI(NIR = sentinel2.efps$B8_mean, RED = sentinel2.efps$B4_mean, sigma = "average", site = sentinel2.efps$site_name)

sentinel2.efps.oversampling$kndvi_mean <- kNDVI(NIR = sentinel2.efps.oversampling$B8_mean, RED = sentinel2.efps.oversampling$B4_mean, sigma = "average", site = sentinel2.efps.oversampling$site_name)

sentinel2.efps.undersampling$kndvi_mean <- kNDVI(NIR = sentinel2.efps.undersampling$B8_mean, RED = sentinel2.efps.undersampling$B4_mean, sigma = "average", site = sentinel2.efps.undersampling$site_name)

sentinel2.efps.smoter$kndvi_mean <- kNDVI(NIR = sentinel2.efps.smoter$B8_mean, RED = sentinel2.efps.smoter$B4_mean, sigma = "average", site = sentinel2.efps.smoter$site_name)


# Red-edge VIs vs Non Red-edge VIs -----

# Cross validation 
summary_cv <- data.frame()
datasets <- c("Imbalanced", "Undersampling", "Oversampling", "SMOTER")


# creating folds ----


# imbalanced dataset
sentinel2.efps$doy <- yday(sentinel2.efps$date_index)
indices.imbalanced <- list(list(list()), list(list()))
names(indices.imbalanced) <- c("index", "indexOut")

resampling <- 50
k <- 10
count <- 1
for(i in 1:resampling){
  indices.temp <- CreateSpacetimeFolds(sentinel2.efps, spacevar = "site_name", timevar = "doy", seed = i, k = k)
  for(j in 1:k){
    indices.imbalanced$index[[count]] <- indices.temp$index[[j]]
    indices.imbalanced$indexOut[[count]] <- indices.temp$indexOut[[j]]
    count <- count + 1
  }
}
str(indices.imbalanced)

# undersampling

sentinel2.efps.undersampling$doy <- yday(sentinel2.efps.undersampling$date_index)

indices.undersampling <- list(list(list()), list(list()))
names(indices.undersampling) <- c("index", "indexOut")
resampling <- 50
k <- 10
count <- 1
for(i in 1:resampling){
  indices.temp <- CreateSpacetimeFolds(sentinel2.efps.undersampling, spacevar = "site_name", timevar = "doy", seed = i, k = k)
  for(j in 1:k){
    indices.undersampling$index[[count]] <- indices.temp$index[[j]]
    indices.undersampling$indexOut[[count]] <- indices.temp$indexOut[[j]]
    count <- count + 1
  }
}
str(indices.undersampling)

# oversampling

sentinel2.efps.oversampling$doy <- yday(sentinel2.efps.oversampling$date_index)

indices.oversampling <- list(list(list()), list(list()))
names(indices.oversampling) <- c("index", "indexOut")
resampling <- 50
k <- 10
count <- 1
for(i in 1:resampling){
  indices.temp <- CreateSpacetimeFolds(sentinel2.efps.oversampling, spacevar = "site_name", timevar = "doy", seed = i, k = k)
  for(j in 1:k){
    indices.oversampling$index[[count]] <- indices.temp$index[[j]]
    indices.oversampling$indexOut[[count]] <- indices.temp$indexOut[[j]]
    count <- count + 1
  }
}
str(indices.oversampling)

# SMOTER

sentinel2.efps.smoter$doy <- yday(sentinel2.efps.smoter$date_index)

indices.smoter <- list(list(list()), list(list()))
names(indices.smoter) <- c("index", "indexOut")
resampling <- 50
k <- 10
count <- 1
for(i in 1:resampling){
  indices.temp <- CreateSpacetimeFolds(sentinel2.efps.smoter, spacevar = "site_name", timevar = "doy", seed = i, k = k)
  for(j in 1:k){
    indices.smoter$index[[count]] <- indices.temp$index[[j]]
    indices.smoter$indexOut[[count]] <- indices.temp$indexOut[[j]]
    count <- count + 1
  }
}
str(indices.smoter)


## running folds imbalanced dataset -----
# tic()
# for (i in 1:length(datasets)){
#   if(datasets[i] == "Imbalanced") {
# 
#     indices <- indices.imbalanced
#     input <- sentinel2.efps
# 
#   }else if(datasets[i] == "Undersampling") {
# 
#     indices <- indices.undersampling
#     input <- sentinel2.efps.undersampling
# 
#   }else if(datasets[i] == "Oversampling") {
# 
#     indices <- indices.oversampling
#     input <- sentinel2.efps.oversampling
# 
#   }else if(datasets[i] == "SMOTER") {
# 
#     indices <- indices.smoter
#     input <- sentinel2.efps.smoter
# 
#   }
#   input$doy <- yday(input$date_index)
# 
#   set.seed(10)
# 
#   model.cir <- train(input[,"cir_mean"],input$GPP.day.smooth,
#                      method="lm",
#                      trControl=trainControl(method="cv",
#                                             index = indices$index, indexOut = indices$indexOut))
# 
#   set.seed(10)
#   model.ireci <- train(input[,"ireci_mean"],input$GPP.day.smooth,
#                        method="lm",
#                        trControl=trainControl(method="cv",
#                                               index = indices$index, indexOut =  indices$indexOut))
# 
#   set.seed(10)
#   model.nirv <- train(input[,"nirv_mean"],input$GPP.day.smooth,
#                       method="lm",
#                       trControl=trainControl(method="cv",
#                                              index = indices$index, indexOut = indices$indexOut))
#   set.seed(10)
#   model.ndvi <- train(input[,"ndvi_mean"],input$GPP.day.smooth,
#                       method="lm",
#                       trControl=trainControl(method="cv",
#                                              index = indices$index, indexOut = indices$indexOut))
#   set.seed(10)
#   model.kndvi <- train(input[,"kndvi_mean"],input$GPP.day.smooth,
#                        method="lm",
#                        trControl=trainControl(method="cv",
#                                               index = indices$index, indexOut = indices$indexOut))
# 
#   ##
# 
#   output <- data.frame(VI = rep(c("CIR", "IRECI", "NIRV", "NDVI", "kNDVI"), each = 500*3), metric_name = rep(c("RMSE", "R2", "MAE"), each = 500), metric_value = c(as.vector(as.matrix(model.cir$resample[,-4])), as.vector(as.matrix(model.ireci$resample[,-4])), as.vector(as.matrix(model.nirv$resample[,-4])), as.vector(as.matrix(model.ndvi$resample[,-4])), as.vector(as.matrix(model.kndvi$resample[,-4]))), dataset = datasets[i], technique = "lm")
# 
#   summary_cv <- rbind(summary_cv, output)
#   print(i)
# 
# }
# toc()
# 
# 
# summary_cv <- as_tibble(summary_cv)
# 
# saveRDS(summary_cv, file = "~/dpabon/results_EFPs_ML/linear_comparison_50.rds")
summary_cv <- read_rds("~/dpabon/results_EFPs_ML/linear_comparison_50.rds")

my_comparisons <- list( c("CIR", "IRECI"), c("CIR", "kNDVI"), c("CIR", "NDVI"), c("CIR", "NIRV"), c("IRECI", "kNDVI"), c("IRECI", "NDVI"), c("IRECI", "NIRV"), c("kNDVI", "NDVI"), c("kNDVI", "NIRV"), c("NDVI", "NIRV"))

# Two plots:

# Only Imbalanced data set 

summary_cv$metric_name <- as.factor(summary_cv$metric_name)

# levels(summary_cv$metric_name)[2] <- bquote('R^2')

summary_cv %>% 
  filter(dataset == "Imbalanced", VI != "kCIR") %>%
ggplot(aes(x = VI, y = metric_value, col = VI)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  facet_wrap(~metric_name, scales = "free_y", ncol = 3)+
  xlab("Vegetation Index") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  ylab("Value") +
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label = "p.signif", paired = T) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave("~/dpabon/results_EFPs_ML/plots/gpp_sentinel2_red-edge_vs_non_red_edge.png", width = 7, height = 5)

ggsave("~/dpabon/results_EFPs_ML/plots/gpp_sentinel2_red-edge_vs_non_red_edge.pdf", width = 7, height = 5)

# short-test relationship between bands

# All except Imbalanced dataset 
summary_cv %>% 
  filter(dataset != "Imbalanced",  VI != "kCIR") %>% 
  ggplot(aes(x = VI, y = metric_value, col = VI)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  facet_wrap(vars(dataset, metric_name), scales = "free_y", ncol = 3)+
  xlab("Vegetation Index") +
  ylab("Value") +
  theme(legend.position = "none") +
  theme_bw() +
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label = "p.signif")
  

ggsave("~/dpabon/results_EFPs_ML/plots/gpp_sentinel2_red-edge_vs_non_red_edge_balanced.png", width = 7.89, height = 10.89)

ggsave("~/dpabon/results_EFPs_ML/plots/gpp_sentinel2_red-edge_vs_non_red_edge_balanced.pdf", width = 7.89, height = 10.89)

# comparison imbalanced vs balanced datasets

summary_cv %>% 
  ggplot(aes(x = dataset, y = metric_value, col = dataset)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  facet_wrap(vars(metric_name, VI), scales = "free_y", ncol = 6)+
  xlab("Vegetation Index") +
  ylab("Value") +
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  theme(legend.position = "none") +
  stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label = "p.signif", paired = T)



# Some important numbers

summary_cv %>% 
  group_by(dataset, VI, metric_name) %>% 
  summarise(mean_value = mean(metric_value, na.rm = T)) %>%
  mutate(mean_value = round(mean_value, 2)) %>% 
  print(n = 72)

AIC(model.cir$finalModel)
AIC(model.nirv$finalModel)
summary_cv %>% 
  group_by(VI, metric_name) %>% 
  summarize(mean.v = mean(metric_value))
  

pft <- unique(as.character(sentinel2.efps$igbp))

predictors <- c("ireci_mean", "cir_mean", "nirv_mean", "ndvi_mean")

group_variables <- c(pft, "global_imba", "global_under", "global_over", "global_smoter")

metrics <- (c("r2", "RMSE", "MEF"))

gpp.linear.table <- matrix(data = NA, nrow = length(group_variables), ncol = length(predictors))
colnames(gpp.linear.table) <- predictors
rownames(gpp.linear.table) <- group_variables

metrics.index <- list()

response <- c("GPP.day.smooth", "gppsat")

for(r in 1:length(response)){
  for (m in 1:length(metrics)){
    for (i in 1:length(predictors)){
      for (j in 1:length(group_variables)) {
        if(is.element(group_variables[j], pft)){
          tempora <- sentinel2.efps %>% 
            filter(igbp == group_variables[j]) %>% 
            as.data.frame()
        }
        if (group_variables[j] == "global_imba"){
          tempora <- sentinel2.efps %>% 
            as.data.frame()
        }
        if (group_variables[j] == "global_under"){
          tempora <- sentinel2.efps.undersampling %>% 
            as.data.frame()
        }
        if(group_variables[j] == "global_over"){
          tempora <- sentinel2.efps.oversampling %>% 
            as.data.frame()
        }
        if(group_variables[j] == "global_smoter"){
          if(response[r] == "GPP.day.smooth") {
            tempora <- sentinel2.efps.smoter %>% 
              as.data.frame()
          }else{
            tempora <- sentinel2.efps.smoter.gppsat %>% 
              as.data.frame()
          }
        }
        lm.model <- lm(tempora[,response[r]] ~ tempora[,predictors[i]])
        tempora.predicted <- predict(lm.model)
        model.eval <- modelEval(obs = tempora[,response[r]], mod = tempora.predicted)
        
        gpp.linear.table[j,i] <- model.eval[which(names(model.eval) == metrics[m])]
      }
    }
    write.csv(as.data.frame(gpp.linear.table), file = paste0("~/dpabon/results_EFPs_ML/prediction_linear_VT/", response[r], "_", metrics[m], ".csv"))
    gpp.linear.table <- matrix(data = NA, nrow = length(group_variables), ncol = length(predictors))
    colnames(gpp.linear.table) <- predictors
    rownames(gpp.linear.table) <- group_variables
  }
}


# GLMM ----
# sentinel2.efps <- sentinel2.efps %>% mutate(across(where(is_character),as_factor))

gpp.glm <- lmer(GPP.day.smooth ~ cir_mean + (cir_mean|igbp), data = sentinel2.efps)
gpp.glm.2 <- lmer(GPP.day.smooth ~ cir_mean + (1|igbp), data = sentinel2.efps)

AIC(gpp.glm.2, gpp.glm)

# predict(gpp.glm)
# summary(gpp.glm)
r.squaredGLMM(gpp.glm)
r.squaredGLMM(gpp.glm.2)

gpp.glm <- lmer(GPP.day.smooth ~ ireci_mean + (ireci_mean|igbp), data = sentinel2.efps)

# summary(gpp.glm)
r.squaredGLMM(gpp.glm)


gpp.glm <- lmer(GPP.day.smooth ~ nirv_mean + (nirv_mean|igbp), data = sentinel2.efps)
# summary(gpp.glm)
r.squaredGLMM(gpp.glm)

summary(gpp.glm)
r.squaredGLMM(gpp.glm)

gpp.glm.oversampling <- lmer(GPP.day.smooth ~ cir_mean + (1|igbp/site_name), data = sentinel2.efps.oversampling)

gpp.glm.undersampling <- lmer(GPP.day.smooth ~ cir_mean + (1|igbp/site_name), data = sentinel2.efps.undersampling)


plot(gpp.glm)
qqnorm(residuals(gpp.glm))
qqline(residuals(gpp.glm))
qqplot(gpp.glm)

r.squaredGLMM(gpp.glm)[2] - r.squaredGLMM(gpp.glm)[1]

r.squaredGLMM(gpp.glm.oversampling)
r.squaredGLMM(gpp.glm.oversampling)[2] - r.squaredGLMM(gpp.glm.oversampling)[1]


# variables selection analysis per site and time ----

predictors.v1 <- c("arvi_mean",
                   "dvi_mean",
                   "gemi_mean",
                   "gndvi_mean",
                   "ipvi_mean",
                   "ireci_mean",
                   "mcari_mean",
                   "msavi2_mean",
                   "msavi_mean",
                   "mtci_mean",
                   "ndi45_mean",
                   "pssra_mean",
                   "pvi_mean",
                   "rvi_mean",
                   "s2rep_mean",
                   "savi_mean",
                   "tndvi_mean",
                   "tsavi_mean",
                   "wdvi_mean",
                   "cir_mean",
                   "cig_mean",
                   "ndvi_mean",
                   "nirv_mean",
                   "B1_mean",
                   "B2_mean",
                   "B3_mean",
                   "B4_mean",
                   "B5_mean",
                   "B6_mean",
                   "B7_mean",
                   "B8_mean",
                   "B8A_mean",
                   "B9_mean",
                   "B11_mean",
                   "B12_mean")


sentinel2.efps$doy <- yday(sentinel2.efps$date_index)

# Imbalanced

# bot$sendMessage(chat_id = chat_id, text = paste("FFS starting for Imbalanced"))
# 
# cl <- makePSOCKcluster(40)
# registerDoParallel(cl)
# tic()
# ffsmodel_LLO <- ffs(sentinel2.efps[,predictors.v1], sentinel2.efps$GPP.day.smooth,
#                     metric="Rsquared",
#                     method="rf",
#                     verbose=TRUE, seed = 10,
#                     trControl=trainControl(method="cv",
#                                            index = indices.imbalanced$index, indexOut = indices.imbalanced$indexOut))
# toc()
# ffsmodel_LLO$selectedvars
# ffsmodel_LLO$selectedvars_perf
# 
# stopCluster(cl)
# toc()
# bot$sendMessage(chat_id = chat_id, text = paste("FFS done for Imbalanced"))

# saveRDS(ffsmodel_LLO, "~/dpabon/results_EFPs_ML/variable_selection_analysis/10_fold_meyer_var_selection_imbalanced.RDS")

# Undersampling 


# sentinel2.efps.undersampling$doy <- yday(sentinel2.efps.undersampling$date_index)
# 
# bot$sendMessage(chat_id = chat_id, text = paste("FFS starting for Undersampling"))
# cl <- makePSOCKcluster(40)
# registerDoParallel(cl)
# tic()
# 
# ffsmodel_LLO_undersampling <- ffs(sentinel2.efps.undersampling[,predictors.v1], sentinel2.efps.undersampling$GPP.day.smooth,
#                     metric="Rsquared",
#                     method="rf",
#                     verbose=TRUE, seed = 10,
#                     trControl=trainControl(method="cv",
#                                            index = indices.undersampling$index, indexOut = indices.undersampling$indexOut))
# 
# stopCluster(cl)
# toc()
# saveRDS(ffsmodel_LLO_undersampling, "~/dpabon/results_EFPs_ML/variable_selection_analysis/10_fold_meyer_var_selection_undersampling.RDS")
# 
# bot$sendMessage(chat_id = chat_id, text = paste("FFS ready for undersampling"))

# Oversampling

# bot$sendMessage(chat_id = chat_id, text = paste("FFS starting for oversampling"))
# 
# sentinel2.efps.oversampling$doy <- yday(sentinel2.efps.oversampling$date_index)
# 
# 
# cl <- makePSOCKcluster(40)
# registerDoParallel(cl)
# tic()
# ffsmodel_LLO_oversampling <- ffs(sentinel2.efps.oversampling[,predictors.v1], sentinel2.efps.oversampling$GPP.day.smooth,
#                                   metric="Rsquared",
#                                   method="rf",
#                                   verbose=TRUE, seed = 10,
#                                   trControl=trainControl(method="cv",
#                                                          index = indices.oversampling$index, indexOut = indices.oversampling$indexOut))
# 
# stopCluster(cl)
# toc()
# saveRDS(ffsmodel_LLO_oversampling, "~/dpabon/results_EFPs_ML/variable_selection_analysis/10_fold_meyer_var_selection_oversampling.RDS")

# bot$sendMessage(chat_id = chat_id, text = paste("FFS ready for oversampling"))

# smoter
# sentinel2.efps.smoter$doy <- yday(sentinel2.efps.smoter$date_index)
# 
# 
# bot$sendMessage(chat_id = chat_id, text = paste("FFS starting for smoter"))
# 
# cl <- makePSOCKcluster(40)
# registerDoParallel(cl)
# tic()
# ffsmodel_LLO_smoter <- ffs(sentinel2.efps.smoter[,predictors.v1], sentinel2.efps.smoter$GPP.day.smooth,
#                                  metric="Rsquared",
#                                  method="rf",
#                                  verbose=TRUE, seed = 10,
#                                  trControl=trainControl(method="cv",
#                                                         index = indices.smoter$index, indexOut = indices.smoter$indexOut))
# 
# stopCluster(cl)
# toc()
# saveRDS(ffsmodel_LLO_smoter, "~/dpabon/results_EFPs_ML/variable_selection_analysis/10_fold_meyer_var_selection_smoter.RDS")
# 
# bot$sendMessage(chat_id = chat_id, text = paste("FFS ready for smoter"))

# loading models

model_imbalanced <- readRDS("~/dpabon/results_EFPs_ML/variable_selection_analysis/10_fold_meyer_var_selection_imbalanced.RDS")

model_imbalanced$selectedvars
round(model_imbalanced$selectedvars_perf, 3)
round(model_imbalanced$selectedvars_perf_SE, 3)

plot(model_imbalanced)
model_imbalanced$selectedvars_perf_SE
str(model_imbalanced)
model_imbalanced$resample
model_imbalanced$results

plot(sentinel2.efps$GPP.day.smooth, model_imbalanced$finalModel$y)

sentinel2.efps$GPP.predict.imbalance <- model_imbalanced$finalModel$predicted

test_new_plot <- data.frame(observed = model_imbalanced$finalModel$y, predicted =  model_imbalanced$finalModel$predicted)

nrow(test_new_plot)
nrow(sentinel2.efps)

gpp_imbalance_plot <- ggplot(test_new_plot, aes(x = predicted, y = observed))+
  geom_hex()+
  # geom_smooth(method = "lm") +
  annotate("text", x = max(test_new_plot$predicted), y = min(test_new_plot$observed),
           label = paste0("paste(italic(R)[10-fold]^2,  \" =", 0.66,  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw() +
  scale_fill_viridis(direction = -1) +
  ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Observed)')) +
  xlab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Predicted)')) +
  ggtitle("Imbalanced")


model_undersampling <- readRDS("~/dpabon/results_EFPs_ML/variable_selection_analysis/10_fold_meyer_var_selection_undersampling.RDS")

round(model_undersampling$selectedvars_perf, 3)

model_undersampling$selectedvars
round(model_undersampling$selectedvars_perf_SE, 3)

model_undersampling$results

test_new_plot <- data.frame(observed = model_undersampling$finalModel$y, predicted =  model_undersampling$finalModel$predicted)

nrow(test_new_plot)
nrow(sentinel2.efps.undersampling)

gpp_undersampling_plot <- ggplot(test_new_plot, aes(x = predicted, y = observed))+
  geom_hex()+
  # geom_smooth(method = "lm") +
  annotate("text", x = max(test_new_plot$predicted), y = min(test_new_plot$observed),
           label = paste0("paste(italic(R)[10-fold]^2,  \" =", 0.68,  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw() +
  scale_fill_viridis(direction = -1) +
  ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Observed)')) +
  xlab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Predicted)')) +
  ggtitle("Undersampling")




model_oversampling <- readRDS("~/dpabon/results_EFPs_ML/variable_selection_analysis/10_fold_meyer_var_selection_oversampling.RDS")

model_oversampling$selectedvars

round(model_oversampling$selectedvars_perf, 2)

round(model_oversampling$selectedvars_perf_SE, 3)

model_oversampling$results

test_new_plot <- data.frame(observed = model_oversampling$finalModel$y, predicted =  model_oversampling$finalModel$predicted)

nrow(test_new_plot)
nrow(sentinel2.efps.oversampling)

gpp_oversampling_plot <- ggplot(test_new_plot, aes(x = predicted, y = observed))+
  geom_hex()+
  # geom_smooth(method = "lm") +
  annotate("text", x = max(test_new_plot$predicted), y = min(test_new_plot$observed),
           label = paste0("paste(italic(R)[10-fold]^2,  \" =", 0.67,  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw() +
  scale_fill_viridis(direction = -1) +
  ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Observed)')) +
  xlab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Predicted)')) +
  ggtitle("Oversampling")

gpp_oversampling_plot


model_smoter <- readRDS("~/dpabon/results_EFPs_ML/variable_selection_analysis/10_fold_meyer_var_selection_smoter.RDS")

round(model_smoter$selectedvars_perf, 3)

round(model_smoter$selectedvars_perf_SE, 3)

model_smoter$selectedvars

model_smoter$results
model_smoter$bestTune

sentinel2.efps.smoter$GPP.predict.smoter <- model_smoter$finalModel$predicted

test_new_plot <- data.frame(observed = model_smoter$finalModel$y, predicted =  model_smoter$finalModel$predicted)

nrow(test_new_plot)
nrow(sentinel2.efps.smoter)

gpp_smoter_plot <- ggplot(test_new_plot, aes(x = predicted, y = observed))+
  geom_hex()+
  # geom_smooth(method = "lm") +
  annotate("text", x = max(test_new_plot$predicted), y = min(test_new_plot$observed),
           label = paste0("paste(italic(R)[10-fold]^2,  \" =", 0.71,  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw() +
  scale_fill_viridis(direction = -1) +
  ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Observed)')) +
  xlab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Predicted)')) +
  ggtitle("SMOTER")


gpp_smoter_plot

ggarrange(gpp_imbalance_plot, gpp_undersampling_plot, gpp_oversampling_plot, gpp_smoter_plot, ncol = 1)

ggsave("~/dpabon/results_EFPs_ML/random_forest/global_datasets_fig2.png", width = 4, height = 13)
ggsave("~/dpabon/results_EFPs_ML/random_forest/global_datasets_fig2.pdf", width = 4, height = 13)


# Fixing plot 4 ----

plot4.data <- left_join(sentinel2.efps, sentinel2.efps.smoter, by = c("date_index", "site_name", "GPP.day.smooth", "igbp"))
nrow(plot4.data)
plot4.data$GPP.predict.imbalance
plots.cv.gpp <- list()
nrow(sentinel2.efps.smoter)
unique(sentinel2.efps.smoter$site_name)
# sentinel2.efps$gppsat.cv <- 0
# sentinel2.efps$gppsat.cv.smoter <- 0

for(i in 1:length(unique(plot4.data$site_name))){
    
    plots.cv.gpp[[as.character(unique(plot4.data$site_name)[i])]] <- plot4.data %>%
      filter(site_name == as.character(unique(plot4.data$site_name)[i])) %>% 
      ggplot() +
      geom_point(aes(x = date_index, y = GPP.day.smooth, col = "Observed")) +
      geom_point(aes(x = date_index, y = GPP.predict.imbalance, col = "Predicted Imbalanced")) +
      geom_point(aes(x = date_index, y = GPP.predict.smoter, col = "Predicted SMOTER")) +
      labs(color = "") +
      ggtitle(paste0(unique(plot4.data$site_name)[i], " (", plot4.data$igbp[which(plot4.data$site_name == unique(plot4.data$site_name)[i])[1]], ")")) +
      ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2} m^{-2} s^{-1} \\rbrack$')) +
      xlab("Date") +
      xlim(c(ymd("2015-04-01"), ymd("2018-12-31"))) +
      theme(axis.title.y=element_blank(), axis.title.x=element_blank()) +
      theme_bw()
}

plots.cv.gpp$`BE-Bra`


obs.predic.gpp <- plot4.data %>% 
  group_by(igbp) %>% 
  do(
    plots = ggplot(data = .) +
      geom_point(aes(x = date_index, y = GPP.day.smooth, col = "Observed")) +
      geom_point(aes(x = date_index, y = GPP.predict.imbalance, col = "Predicted Imbalanced")) +
      geom_point(aes(x = date_index, y = GPP.predict.smoter, col = "Predicted SMOTER")) +
      labs(color = "")+
      facet_wrap(~site_name, scales = "free_y") +
      ggtitle(paste(unique(.$igbp))) +
      ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2} m^{-2} s^{-1} \\rbrack$')) +
      xlab("Date")
  )

obs.predic.gpp$plots[[1]]
for (i in 1:nrow(obs.predic.gpp)) {
  ggsave(filename = paste0("~/dpabon/results_EFPs_ML/leave_one_out_cv/", obs.predic.gpp$igbp[i], ".png"), plot = obs.predic.gpp$plots[[i]], width = 14, height = 8)
}


# figure XXX

obs.predic.example.gpp <- ggarrange(plots.cv.gpp[["RU-Fy2"]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                                    plots.cv.gpp[["DE-Geb"]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                                    plots.cv.gpp[["DE-Hai"]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                                    plots.cv.gpp[["DE-Gri"]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                                    plots.cv.gpp[["CZ-wet"]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                                    plots.cv.gpp[["CH-Lae"]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                                    plots.cv.gpp[["ES-LM2"]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()),
                                    plots.cv.gpp[["IT-Lsn"]] + theme(axis.title.y = element_blank(), axis.title.x = element_blank()), ncol = 1, nrow = 8, common.legend = T, legend = "bottom", align = "h")

obs.predic.example.gpp

obs.predic.example.gpp <- annotate_figure(obs.predic.example.gpp,  left = "GPP umolCO2 m-2 s-1", bottom = "Date")

obs.predic.example.gpp

ggsave("~/dpabon/results_EFPs_ML/leave_one_out_cv/example_obs_vs_predict.png",  height = 10, width = 6)
ggsave("~/dpabon/results_EFPs_ML/leave_one_out_cv/example_obs_vs_predict.pdf",  height = 10, width = 6)


sentinel2.efps %>% 
  filter(GPP.day.smooth > 20)

# variables selected

vars.selected <- c(model_imbalanced$selectedvars,
model_undersampling$selectedvars,
model_oversampling$selectedvars,
model_smoter$selectedvars)
summary(as.factor(vars.selected))
hist(vars.selected)

vars.selected <- unique(vars.selected)

# variables selected ----

comp_vars_selected <- as.factor(c(model_imbalanced$selectedvars, model_oversampling$selectedvars, model_smoter$selectedvars, model_undersampling$selectedvars))

summary(comp_vars_selected)

# results folds ----

write_csv(model_imbalanced$resample, file = "~/dpabon/results_EFPs_ML/leave_one_out_cv/results_cross_validation_rf_imbalanced.csv")

write_csv(model_undersampling$resample, file = "~/dpabon/results_EFPs_ML/leave_one_out_cv/results_cross_validation_rf_undersampling.csv")

write_csv(model_oversampling$resample, file = "~/dpabon/results_EFPs_ML/leave_one_out_cv/results_cross_validation_rf_oversampling.csv")


write_csv(model_smoter$resample, file = "~/dpabon/results_EFPs_ML/leave_one_out_cv/results_cross_validation_rf_smoter.csv")


# dataset comparison ----
# comparison of the performance between the different datasets

summary_rf_datasets <- data.frame(metric_value = c(as.vector(as.matrix(model_imbalanced$resample[-4])), as.vector(as.matrix(model_undersampling$resample[-4])), as.vector(as.matrix(model_oversampling$resample[-4])), as.vector(as.matrix(model_smoter$resample[-4]))), metric_name = rep(c("RMSE", "R2", "MAE"), each = 500), dataset = rep(c("Imbalanced", "Undersampling", "Oversampling", "SMOTER"), each = 3*500), technique = "rf")

summary_rf_datasets <- as_tibble(summary_rf_datasets)

my_comparisons <- list(c("Imbalanced", "Oversampling"), c("Imbalanced", "SMOTER"), c("Imbalanced", "Undersampling"), c("Oversampling", "SMOTER"), c("Oversampling", "Undersampling"), c("SMOTER", "Undersampling"))

summary_rf_datasets %>% 
  ggplot(aes(x = dataset, y = metric_value, col = dataset)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  facet_wrap(~metric_name, scales = "free_y") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  xlab("Dataset") +
  ylab("Value") +
  theme(legend.position = "none") +
    stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label = "p.signif") +
  theme_bw()

ggsave("~/dpabon/results_EFPs_ML/random_forest/balanced_vs_imbalanced.png", height = 6, width = 8)

ggsave("~/dpabon/results_EFPs_ML/random_forest/balanced_vs_imbalanced.pdf", height = 6, width = 8)

model_imbalanced$bestTune$mtry
model_undersampling$bestTune$mtry
model_oversampling$bestTune$mtry
model_smoter$bestTune$mtry

# save(model_smoter, file = "~/dpabon/results_EFPs_ML/model_smoter.RData")

summary_rf_datasets %>% 
  group_by(dataset, metric_name) %>% 
  summarise(metric.v = mean(metric_value)) %>% 
  mutate(metric.v = round(metric.v, 2))

linear_4_comparison <- summary_cv %>% 
  filter(VI == "CIR") %>% 
  select(-VI)

lm_vs_rf <- rbind(summary_rf_datasets, linear_4_comparison)

my_comparisons <- list(c("rf", "lm"))

lm_vs_rf %>% 
  ggplot(aes(x = dataset, y = metric_value, col = technique)) +
  # geom_violin() +
  geom_boxplot(width = 0.1) +
  facet_wrap(~metric_name, scales = "free_y", nrow = 3) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  xlab("Dataset") +
  ylab("Value")


# predicted vs observed plots ----

# upscaling gpp example ----

source("dpabon/bin_EFPs_ML/EFPs_ML/sentinel-2_gpp_function.R")

example.image.gpp <- sentinel2.efps(sentinel2_product = "~/dpabon/data_EFPS_ML/example_entire_tile/S2A_MSIL2A_20200623T103031_N0214_R108_T32ULU_20200623T142851.SAFE/", model = model_smoter$finalModel, product = "gpp")

# example.image.gpp <- projectRaster(example.image.gpp, crs=CRS("+proj=longlat +datum=WGS84"))
library(sen2r)

example.image.rgb <- sen2r_rgb(sentinel2_product = "~/dpabon/data_EFPS_ML/example_entire_tile/S2A_MSIL2A_20200623T103031_N0214_R108_T32ULU_20200623T142851.SAFE/", resolution = "20m")

# example.image.rgb <- projectRaster(example.image.rgb, crs=CRS("+proj=longlat +datum=WGS84"))

test.rgb <- projectRaster(example.image.rgb, crs=CRS("+proj=longlat +datum=WGS84"))

test.gpp <-  projectRaster(example.image.gpp, crs=CRS("+proj=longlat +datum=WGS84"))
library(RStoolbox)
library(terrainr)

test1 <- ggRGB(test.rgb,r = 1, g = 2, b = 3, stretch = "lin")+
  xlab(expression(paste("Longitude (",degree,")"))) + ylab(expression(paste("Latitude (",degree,")"))) +
  coord_quickmap() +
  theme_bw()

test2 <- ggplot()+
  geom_raster(data = test.gpp, aes(x = x, y = y, fill = GPP)) +
  scale_fill_viridis(direction = -1) +
  coord_quickmap() +
  xlab(expression(paste("Longitude (",degree,")"))) + ylab(expression(paste("Latitude (",degree,")"))) +
  guides(fill = guide_colourbar(title = latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$') , title.position = "right"))+
  theme_bw() +
  theme(legend.title = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
  

gpp_product_example_plot <- ggarrange(test1, test2, ncol = 2, align = , label.x = "lon", label.y = "lat", common.legend = T, legend = "right")

ggsave(plot = gpp_product_example_plot,filename =  "~/dpabon/results_EFPs_ML/plots/example_gpp_tile.pdf", width = 7, height = 3)

# writeRaster(example.image, filename = "~/dpabon/results_EFPs_ML/example_product.tif", format="GTiff")
# writeRaster(example.image, filename = "~/dpabon/results_EFPs_ML/example_product.nc", format="CDF")
# writeRaster(example.image, filename = "~/dpabon/results_EFPs_ML/example_product.envi", format="ENVI")
writeRaster(example.image.gpp, filename = "~/dpabon/results_EFPs_ML/example_product.img", format="HFA", overwrite=T)
test.gpp
str(example.image)
# # Random Forest training ----
# 
# # per vegetation type
# 
gpp.rf <- list()
gppsat.rf <- list()

datasets <- c("global_imba", "global_under", "global_over", "global_smoter")



for(i in 1:length(datasets)){
  if(datasets[i] == "global_imba"){
    tempora <- sentinel2.efps

    gpp.rf[["global_imba"]] <- randomForest(x = tempora[,model_imbalanced$selectedvars], y = tempora$GPP.day.smooth, importance = T)
    # gppsat.rf[[i]] <- randomForest(x = tempora[,], y = tempora$gppsat, importance = T)
    sentinel2.efps$GPP.predict.imbalance <- predict(gpp.rf[[i]])
    # sentinel2.efps$GPPsat.predict.imbalance <- predict(gppsat.rf[[i]])

  }
  if(datasets[i] == "global_under"){
    tempora <- sentinel2.efps.undersampling
    gpp.rf[["global_under"]] <- randomForest(x = tempora[,model_undersampling$selectedvars], y = tempora$GPP.day.smooth, importance = T)
    # gppsat.rf[[i]] <- randomForest(x = tempora[,predictors.global], y = tempora$gppsat, importance = T)
    sentinel2.efps$GPP.predict.undersampling <- predict(gpp.rf[[i]], sentinel2.efps[,model_undersampling$selectedvars])
    # sentinel2.efps$GPPsat.predict.undersampling <- predict(gppsat.rf[[i]], sentinel2.efps[,predictors.global])

  }
  if(datasets[i] == "global_over"){
    tempora <- sentinel2.efps.oversampling
    gpp.rf[["global_over"]] <- randomForest(x = tempora[,model_oversampling$selectedvars], y = tempora$GPP.day.smooth, importance = T)
    # gppsat.rf[[i]] <- randomForest(x = tempora[,predictors.global], y = tempora$gppsat, importance = T)
    sentinel2.efps$GPP.predict.oversampling <- predict(gpp.rf[[i]], sentinel2.efps[,model_oversampling$selectedvars])
    # sentinel2.efps$GPPsat.predict.oversampling <- predict(gppsat.rf[[i]], sentinel2.efps[,predictors.global])
  }
  if(datasets[i] == "global_smoter") {
    tempora <- sentinel2.efps.smoter
    gpp.rf[["global_smoter"]] <- randomForest(x = tempora[,model_smoter$selectedvars], y = tempora$GPP.day.smooth, importance = T)
    # tempora <- sentinel2.efps.smoter.gppsat
    # gppsat.rf[[i]] <- randomForest(x = tempora[,predictors.global], y = tempora$gppsat, importance = T)
    sentinel2.efps$GPP.predict.smoter <- predict(gpp.rf[[i]], sentinel2.efps[,model_smoter$selectedvars])
    # sentinel2.efps$GPPsat.predict.smoter <- predict(gppsat.rf[[i]], sentinel2.efps[,predictors.global])
  }
  print(i)
}

# 
# 
# gpp.rf[["global_smoter"]]
# gpp.rf[["global_under"]]
# gpp.rf[["global_over"]]
# gpp.rf[["global_imba"]]

# Figure 2.

# GPP

# GPP Imbalanced

# sentinel2.efps$GPP.predict.imbalance <- predict(gpp.rf$global_imba, sentinel2.efps[,model_imbalanced$selectedvars])
# 
# gpp.imba.rf.plot <- sentinel2.efps %>% 
#   ggplot(aes(y = GPP.day.smooth, x = GPP.predict.imbalance)) +
#   geom_hex()+
#   # geom_smooth(method = "lm") +
#   annotate("text", x = max(sentinel2.efps$GPP.predict.imbalance), y = min(sentinel2.efps$GPP.day.smooth),
#            label = paste0("paste(italic(R)[10-fold]^2,  \" =", 0.69,  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
#   geom_abline(slope = 1, intercept = 0, col = "red") +
#   theme_bw() +
#   scale_fill_viridis() +
#   ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Observed)')) +
#   xlab(latex2exp::TeX('GPP $\\lbrack \\mu mol  CO_{2}  m^{-2}  s^{-1} \\rbrack$ (Predicted)')) +
#   ggtitle("Imbalanced")
# 
# gpp.imba.rf.plot
# 
# # GPP Undersampling
# 
# sentinel2.efps.undersampling$GPP.day.predicted <- predict(gpp.rf$global_under, sentinel2.efps.undersampling[,model_undersampling$selectedvars])
# 
# gpp.under.rf.plot <- sentinel2.efps.undersampling %>% 
#   ggplot(aes(y = GPP.day.smooth, x = GPP.day.predicted)) +
#   geom_hex()+
#   # geom_smooth(method = "lm") +
#   annotate("text", x = max(sentinel2.efps.undersampling$GPP.day.predicted), y = min(sentinel2.efps.undersampling$GPP.day.smooth),
#            label = paste0("paste(italic(R)[10-fold]^2,  \" =", 0.71,  "\")"), parse = TRUE, colour = "red", size = 3, vjust= -2, hjust=1) +
#   geom_abline(slope = 1, intercept = 0, col = "red") +
#   theme_bw() +
#   scale_fill_viridis() +
#   ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol CO_{2} m^{-2} s^{-1} \\rbrack$ (Observed)')) +
#   xlab(latex2exp::TeX('GPP $\\lbrack \\mu mol CO_{2} m^{-2} s^{-1} \\rbrack$ (Predicted)')) +
#   ggtitle("Udersampling")
# 
# gpp.under.rf.plot
#   
# # GPP Oversampling 
# 
# sentinel2.efps.oversampling$GPP.day.predicted <- predict(gpp.rf$global_over, sentinel2.efps.oversampling[,model_oversampling$selectedvars])
# 
# gpp.over.rf.plot <- sentinel2.efps.oversampling %>% 
#   ggplot(aes(y = GPP.day.smooth, x = GPP.day.predicted)) +
#   geom_hex()+
#   # geom_smooth(method = "lm") +
#   annotate("text", x = max(sentinel2.efps.oversampling$GPP.day.predicted), y = min(sentinel2.efps.oversampling$GPP.day.smooth),
#            label = paste0("paste(italic(R)[10-fold]^2,  \" =", 0.71,  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
#   geom_abline(slope = 1, intercept = 0, col = "red") +
#   theme_bw() +
#   scale_fill_viridis() +
#   ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol CO_{2} m^{-2} s^{-1} \\rbrack$ (Observed)')) +
#   xlab(latex2exp::TeX('GPP $\\lbrack \\mu mol CO_{2} m^{-2} s^{-1} \\rbrack$ (Predicted)')) +
#   ggtitle("Oversampling")
# 
# gpp.over.rf.plot
# 
# # GPP SMOTER
# 
# sentinel2.efps.smoter$GPP.day.predicted <- predict(gpp.rf$global_smoter, sentinel2.efps.smoter[,model_smoter$selectedvars])
# 
# gpp.smoter.rf.plot <- sentinel2.efps.smoter %>% 
#   ggplot(aes(y = GPP.day.smooth, x = GPP.day.predicted)) +
#   geom_hex()+
#   # geom_smooth(method = "lm") +
#   annotate("text", x = max(sentinel2.efps.smoter$GPP.day.predicted), y = min(sentinel2.efps.smoter$GPP.day.smooth),
#            label = paste0("paste(italic(R)[10-fold]^2,  \" =", 0.71,  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
#   geom_abline(slope = 1, intercept = 0, col = "red") +
#   theme_bw() +
#   scale_fill_viridis() +
#   ylab(latex2exp::TeX('GPP $\\lbrack \\mu mol CO_{2} m^{-2} s^{-1} \\rbrack$ (Observed)')) +
#   xlab(latex2exp::TeX('GPP $\\lbrack \\mu mol CO_{2} m^{-2} s^{-1} \\rbrack$ (Predicted)')) +
#   ggtitle("SMOTER")
# 
# gpp.smoter.rf.plot
# 
# ggarrange(gpp.imba.rf.plot, gpp.under.rf.plot, gpp.over.rf.plot, gpp.smoter.rf.plot, ncol = 1)
# 
# ggsave("~/dpabon/results_EFPs_ML/random_forest/global_datasets_fig2.png", width = 4, height = 12)
# ggsave("~/dpabon/results_EFPs_ML/random_forest/global_datasets_fig2.pdf", width = 4, height = 12)

# GPPsat

# GPPsat Imbalanced

sentinel2.efps$GPPsat.predict.imbalance <- predict(gppsat.rf$global_imba, sentinel2.efps[,predictors.global])

gppsat.imba.rf.plot <- sentinel2.efps %>% 
  ggplot(aes(y = gppsat, x = GPPsat.predict.imbalance)) +
  geom_hex()+
  geom_smooth(method = "lm") +
  annotate("text", x = max(sentinel2.efps$GPPsat.predict.imbalance), y = min(sentinel2.efps$gppsat),
           label = paste0("paste(italic(R)[10-fold]^2,  \" =", round(gppsat.rf$global_imba$results[which(gppsat.rf$global_imba$results[,1]== gppsat.rf$global_imba$bestTune$mtry),3], digits = 2),  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw() +
  scale_fill_viridis() +
  ylab(TeX('GPPsat $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$ (Observed)')) +
  xlab(TeX('GPPsat $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$ (Predicted)')) +
  ggtitle("Imbalanced")

gppsat.imba.rf.plot

# GPPsat Undersampling

sentinel2.efps.undersampling$GPPsat.day.predicted <- predict(gppsat.rf$global_under, sentinel2.efps.undersampling[,predictors.global])

gppsat.under.rf.plot <- sentinel2.efps.undersampling %>% 
  ggplot(aes(y = gppsat, x = GPPsat.day.predicted)) +
  geom_hex()+
  geom_smooth(method = "lm") +
  annotate("text", x = max(sentinel2.efps.undersampling$GPPsat.day.predicted), y = min(sentinel2.efps.undersampling$gppsat),
           label = paste0("paste(italic(R)[10-fold]^2,  \" =", round(gppsat.rf$global_under$results[which(gppsat.rf$global_under$results[,1]== gppsat.rf$global_under$bestTune$mtry),3], digits = 2),  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw() +
  scale_fill_viridis() +
  ylab(TeX('GPPsat $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$ (Observed)')) +
  xlab(TeX('GPPsat $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$ (Predicted)')) +
  ggtitle("Udersampling")

gppsat.under.rf.plot

# GPPsat Oversampling 

sentinel2.efps.oversampling$GPPsat.day.predicted <- predict(gppsat.rf$global_over, sentinel2.efps.oversampling[,predictors.global])

gppsat.over.rf.plot <- sentinel2.efps.oversampling %>% 
  ggplot(aes(y = gppsat, x = GPPsat.day.predicted)) +
  geom_hex()+
  geom_smooth(method = "lm") +
  annotate("text", x = max(sentinel2.efps.oversampling$GPPsat.day.predicted), y = min(sentinel2.efps.oversampling$gppsat),
           label = paste0("paste(italic(R)[10-fold]^2,  \" =", round(gppsat.rf$global_over$results[which(gppsat.rf$global_over$results[,1]== gppsat.rf$global_over$bestTune$mtry),3], digits = 2),  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw() +
  scale_fill_viridis() +
  ylab(TeX('GPPsat $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$ (Observed)')) +
  xlab(TeX('GPPsat $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$ (Predicted)')) +
  ggtitle("Oversampling")

gppsat.over.rf.plot

# GPPsat SMOTER

sentinel2.efps.smoter.gppsat$GPPsat.day.predicted <- predict(gppsat.rf$global_smoter, sentinel2.efps.smoter.gppsat[,predictors.global])

gppsat.smoter.rf.plot <- sentinel2.efps.smoter.gppsat %>% 
  ggplot(aes(y = gppsat, x = GPPsat.day.predicted)) +
  geom_hex()+
  geom_smooth(method = "lm") +
  annotate("text", x = max(sentinel2.efps.smoter.gppsat$GPPsat.day.predicted), y = min(sentinel2.efps.smoter.gppsat$gppsat),
           label = paste0("paste(italic(R)[10-fold]^2,  \" =", round(gppsat.rf$global_smoter$results[which(gppsat.rf$global_smoter$results[,1]== gppsat.rf$global_smoter$bestTune$mtry),3], digits = 2),  "\")"), parse = TRUE, colour = "red", size = 3, vjust=-2, hjust=1) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw() +
  scale_fill_viridis() +
  ylab(TeX('GPPsat $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$ (Observed)')) +
  xlab(TeX('GPPsat $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$ (Predicted)')) +
  ggtitle("SMOTER")

gppsat.smoter.rf.plot


ggarrange(gpp.imba.rf.plot, gppsat.imba.rf.plot,
          gpp.under.rf.plot, gppsat.under.rf.plot,
          gpp.over.rf.plot, gppsat.over.rf.plot,
          gpp.smoter.rf.plot, gppsat.smoter.rf.plot, ncol = 2, nrow = 4)


#ggsave("~/dpabon/results_EFPs_ML/random_forest/global_datasets_fig2.png", width = 7, height = 12)
#ggsave("~/dpabon/results_EFPs_ML/random_forest/global_datasets_fig2.pdf", width = 7, height = 12)

# leaved-one-site cross validation -----


# sites.names <- sentinel2.efps %>%
#   group_by(site_name) %>%
#   summarise(igbp = unique(igbp))
# 
# summary_efps <- data.frame(EFPs = rep(c("GPP.day.smooth"), each = length(sites.names$site_name)), EFPs.title = rep(c("GPP"), each = length(sites.names$site_name)), site.name = sites.names$site_name, igbp = sites.names$igbp)
# 
# metrics_names <- data.frame(RMSE = 0, r = 0, r2 = 0, MEF = 0, bias = 0, bias_sq = 0, var_err = 0, phase_err = 0)
# 
# summary_efps <- cbind(summary_efps, metrics_names)
# 
# summary_efps <- as_tibble(summary_efps)
# 
# summary_efps <- summary_efps %>%
#   mutate(model.rf = vector("list", length = nrow(summary_efps)),
#          model.rf.smoter = vector("list", length = nrow(summary_efps)))
# 
# 
# summary_efps
# 
# 
# for (i in 1:nrow(summary_efps)) {
#   training <- sentinel2.efps %>%
#     filter(site_name != summary_efps$site.name[i]) %>%
#     as.data.frame()
# 
#   testing <- sentinel2.efps %>%
#     filter(site_name == summary_efps$site.name[i]) %>%
#     as.data.frame()
# 
#   if(summary_efps$EFPs.title[i] == "GPP"){
#     mtry.local <- 23
#   }else{
#     mtry.local <- 26
#   }
# 
#    model.rf <- randomForest(y = training[,as.character(summary_efps$EFPs[i])], x = training[,as.character(model_imbalanced$selectedvars)] , importance = T, mtry = mtry.local)
# 
#    summary_efps$model.rf[[i]] <- model.rf
# 
#   rf_testing <- predict(model.rf, testing[,model_imbalanced$selectedvars])
# 
#   summary_efps[i,5:12] <- matrix(modelEval(obs = testing[,as.character(summary_efps$EFPs[i])], mod = rf_testing), nrow = 1, ncol = 8)
# 
# 
#   training <- sentinel2.efps.smoter %>%
#     filter(site_name != summary_efps$site.name[i]) %>%
#     as.data.frame()
# 
#   testing <- sentinel2.efps.smoter %>%
#     filter(site_name == summary_efps$site.name[i]) %>%
#     as.data.frame()
# 
#   model.rf.smoter <- randomForest(y = training[,as.character(summary_efps$EFPs[i])], x = training[,as.character(model_smoter$selectedvars)] , importance = T, mtry = mtry.local)
# 
#   summary_efps$model.rf.smoter[[i]] <- model.rf.smoter
# 
#   percent <- round((i * 100) / nrow(summary_efps), digits = 2)
#   print(paste(percent, "%"))
# }
# 
# 
# write_rds(summary_efps, "~/dpabon/results_EFPs_ML/leave_one_out_cv/results_leave_one_out_cv.RDS")

gpp.gppsat.cv <- read_rds("~/dpabon/results_EFPs_ML/leave_one_out_cv/results_leave_one_out_cv.RDS")
mean(gpp.gppsat.cv$r2)

plots.cv.gpp <- list()
plots.cv.gppsat <- list()

sentinel2.efps$gpp.cv <- 0
sentinel2.efps$gpp.cv.smoter <- 0
# sentinel2.efps$gppsat.cv <- 0
# sentinel2.efps$gppsat.cv.smoter <- 0

for(i in 1:nrow(gpp.gppsat.cv)){
  tempora <- sentinel2.efps %>% 
    filter(site_name == gpp.gppsat.cv$site.name[i])
  
  if(gpp.gppsat.cv$EFPs.title[i] == "GPP"){
    tempora$GPP.predict <- predict(gpp.gppsat.cv$model.rf[[i]], tempora[,model_imbalanced$selectedvars])
    
    tempora$GPP.predict.smoter <- predict(gpp.gppsat.cv$model.rf.smoter[[i]], tempora[,model_smoter$selectedvars])
    
    sentinel2.efps$gpp.cv[which(sentinel2.efps$site_name == gpp.gppsat.cv$site.name[i])] <- predict(gpp.gppsat.cv$model.rf[[i]], tempora[,model_imbalanced$selectedvars])
    
    sentinel2.efps$gpp.cv.smoter[which(sentinel2.efps$site_name == gpp.gppsat.cv$site.name[i])] <- predict(gpp.gppsat.cv$model.rf.smoter[[i]], tempora[,model_smoter$selectedvars])
    
    plots.cv.gpp[[as.character(gpp.gppsat.cv$site.name[i])]] <- ggplot(data = tempora) +
      geom_point(aes(x = date_index, y = GPP.day.smooth, col = "Observed")) +
      geom_point(aes(x = date_index, y = GPP.predict, col = "Predicted Imbalanced")) +
      geom_point(aes(x = date_index, y = GPP.predict.smoter, col = "Predicted SMOTER")) +
      labs(color = "") +
      ggtitle(paste0(gpp.gppsat.cv$site.name[i], " (", unique(tempora$igbp), ")")) +
      ylab(latex2exp::TeX('GPP $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$')) +
      xlab("Date") +
      xlim(c(ymd("2015-04-01"), ymd("2018-12-31"))) +
      theme(axis.title.y=element_blank(), axis.title.x=element_blank())
    
  }else{
    tempora$GPPsat.predict <- predict(gpp.gppsat.cv$model.rf[[i]], tempora[,predictors.global])
    
    tempora$GPPsat.predict.smoter <- predict(gpp.gppsat.cv$model.rf.smoter[[i]], tempora[,predictors.global])
    
    sentinel2.efps$gppsat.cv[which(sentinel2.efps$site_name == gpp.gppsat.cv$site.name[i])] <-  predict(gpp.gppsat.cv$model.rf[[i]], tempora[,predictors.global])
    
    sentinel2.efps$gppsat.cv.smoter[which(sentinel2.efps$site_name == gpp.gppsat.cv$site.name[i])] <-  predict(gpp.gppsat.cv$model.rf.smoter[[i]], tempora[,predictors.global])
    
    plots.cv.gppsat[[as.character(gpp.gppsat.cv$site.name[i])]] <- ggplot(data = tempora) +
      geom_point(aes(x = date_index, y = gppsat, col = "Observed")) +
      geom_point(aes(x = date_index, y = GPPsat.predict, col = "Predicted Imbalanced")) +
      geom_point(aes(x = date_index, y = GPPsat.predict.smoter, col = "Predicted SMOTER")) +
      labs(color = "") +
      ggtitle(paste0(gpp.gppsat.cv$site.name[i], " (", unique(tempora$igbp), ")")) +
      ylab(latex2exp::TeX('GPP $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$')) +
      xlab("Date") +
      theme(axis.title.y=element_blank(), axis.title.x=element_blank())
  }
}


obs.predic.gpp <- sentinel2.efps %>% 
  group_by(igbp) %>% 
  do(
    plots = ggplot(data = .) +
      geom_point(aes(x = date_index, y = GPP.day.smooth, col = "Observed")) +
      geom_point(aes(x = date_index, y = gpp.cv, col = "Predicted Imbalanced")) +
      geom_point(aes(x = date_index, y = gpp.cv.smoter, col = "Predicted SMOTER")) +
      labs(color = "")+
      facet_wrap(~site_name, scales = "free_y") +
      ggtitle(paste(unique(.$igbp))) +
      ylab(latex2exp::TeX('GPP $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$')) +
      xlab("Date")
  )


for (i in 1:nrow(obs.predic.gpp)) {
  ggsave(filename = paste0("~/dpabon/results_EFPs_ML/leave_one_out_cv/", obs.predic.gpp$igbp[i], ".png"), plot = obs.predic.gpp$plots[[i]], width = 14, height = 8)
}


obs.predic.gppsat <- sentinel2.efps %>% 
  group_by(igbp) %>% 
  do(
    plots = ggplot(data = .) +
      geom_point(aes(x = date_index, y = gppsat, col = "Observed")) +
      geom_point(aes(x = date_index, y = gppsat.cv, col = "Predicted Imbalanced")) +
      geom_point(aes(x = date_index, y = gppsat.cv.smoter, col = "Predicted SMOTER")) +
      labs(color = "") +
      facet_wrap(~site_name, scales = "free_y") +
      ggtitle(paste(unique(.$igbp))) +
      ylab(TeX('GPPsat $\\lbrack \\mu m^{-2} s^{-1} \\rbrack$')) +
      xlab("Date")
  )

for (i in 1:nrow(obs.predic.gppsat)) {
  ggsave(filename = paste0("~/dpabon/results_EFPs_ML/leave_one_out_cv/", obs.predic.gpp$igbp[i], "_gppsat", ".png"), plot = obs.predic.gppsat$plots[[i]], width = 14, height = 8)
}

# figure XXX

obs.predic.example.gpp <- ggarrange(plots.cv.gpp[["RU-Fy2"]],
                                    plots.cv.gpp[["DE-Geb"]],
                                    plots.cv.gpp[["DE-Hai"]],
                                    plots.cv.gpp[["DE-Gri"]],
                                    plots.cv.gpp[["CZ-wet"]],
                                    plots.cv.gpp[["CH-Lae"]],
                                    plots.cv.gpp[["ES-LM2"]],
                                    plots.cv.gpp[["IT-Lsn"]], ncol = 1, nrow = 8, common.legend = T, legend = "bottom", align = "h")

obs.predic.example.gpp

obs.predic.example.gpp <- annotate_figure(obs.predic.example.gpp,  left = "GPP umolCO2 m-2 s-1", bottom = "Date")

obs.predic.example.gpp
ggsave("~/dpabon/results_EFPs_ML/leave_one_out_cv/example_obs_vs_predict.png",  height = 10, width = 6)
ggsave("~/dpabon/results_EFPs_ML/leave_one_out_cv/example_obs_vs_predict.pdf",  height = 10, width = 6)

obs.predic.example.gppsat <- ggarrange(plots.cv.gppsat[["RU-Fy2"]],
                                       plots.cv.gppsat[["DE-Geb"]],
                                       plots.cv.gppsat[["DE-Hai"]],
                                       plots.cv.gppsat[["DE-Gri"]],
                                       plots.cv.gppsat[["CZ-wet"]],
                                       plots.cv.gppsat[["CH-Lae"]],
                                       plots.cv.gppsat[["ES-LM2"]],
                                       plots.cv.gppsat[["IT-Lsn"]], ncol = 1, nrow = 8, common.legend = T, legend = "bottom")

obs.predic.example.gppsat <- annotate_figure(obs.predic.example.gppsat,  left = "GPPsat umolCO2 m-2 s-1", bottom = "Date")

ggarrange(obs.predic.example.gpp, obs.predic.example.gppsat, ncol = 2, common.legend = T)

ggsave("~/dpabon/results_EFPs_ML/leave_one_out_cv/example_obs_vs_predict.png",  height = 12, width = 10)

##############

gppsat.comp <- sentinel2.efps %>% 
  ggplot(aes(y = gppsat, x = GPPsat.predict.imbalance)) +
  geom_hex()+
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  annotate("text", x = max(sentinel2.efps$GPPsat.predict.imbalance), y = min(sentinel2.efps$gppsat),
           label = paste0("paste(italic(R) ^ 2, cv10kfold,  \" =", round(gppsat.rf$global_imba$results[which(gppsat.rf$global_imba$results[,1]== gpp.rf$global_imba$bestTune$mtry),3], digits = 2),  "\")"), parse = TRUE, colour = "red", size = 5, vjust=-2, hjust=1) +
  theme_bw() +
  scale_fill_viridis() +
  ylab("GPPsat (observed)") +
  xlab("GPPsat (predicted)") +
  ggtitle("Model trained globally (Imbalance-Original dataset)")

gppsat.comp
ggarrange(gpp.comp, gppsat.comp)



## Comparison between vegetation types ----

# GPP
sentinel2.efps %>% 
  ggplot(aes(y = GPP.day.smooth, x= GPP.predict.pft)) +
  geom_point() +
  facet_wrap(~igbp) +
  scale_fill_viridis() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  #geom_smooth(method = "lm") +
  #stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  theme_bw() +
  ylab("GPP (observed)") +
  xlab("GPP (predicted)") +
  ggtitle("Model trained per PFT (Imbalance-Original dataset)")


# GPPsat
sentinel2.efps %>% 
  ggplot(aes(y = gppsat, x= GPPsat.predict.pft)) +
  geom_point() +
  facet_wrap(~igbp) +
  scale_fill_viridis() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  #geom_smooth(method = "lm") +
  #stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  theme_bw() +
  ylab("GPP (observed)") +
  xlab("GPP (predicted)") +
  ggtitle("Model trained per PFT (Imbalance-Original dataset)")


#### 
imbalance <- sentinel2.efps %>% 
  ggplot(aes(y = GPP.day.smooth, x = GPP.predict.imbalance)) +
  geom_hex()+
  #geom_smooth(method = "lm")+
  #stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  scale_fill_viridis() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw() +
  ylab("GPP (observed)") +
  xlab("GPP (predicted)") +
  ggtitle("Model trained globally (Imbalance-Original dataset)")


undersampling <- sentinel2.efps %>% 
  ggplot(aes(y = GPP.day.smooth, x = GPP.predict.undersampling)) +
  geom_hex()+
  #geom_smooth(method = "lm") +
  #stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  scale_fill_viridis() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw()+
  ylab("GPP (observed)") +
  xlab("GPP (predicted)") +
  ggtitle("Model trained globally (Balanced, Undersampling)")

undersampling

oversampling <- sentinel2.efps %>% 
  ggplot(aes(y = GPP.day.smooth, x = GPP.predict.oversampling)) +
  geom_hex()+
  #geom_smooth(method = "lm") +
  #stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  scale_fill_viridis() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw()+
  ylab("GPP (observed)") +
  xlab("GPP (predicted)") +
  ggtitle("Model trained globally (Balanced, Oversampling)")

oversampling

smoter <- sentinel2.efps %>% 
  ggplot(aes(y = GPP.day.smooth, x = GPP.predict.smoter)) +
  geom_hex()+
  #geom_smooth(method = "lm") +
  #stat_cor(aes(label = ..rr.label..), color = "red", geom = "label") +
  scale_fill_viridis() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw()+
  ylab("GPP (observed)") +
  xlab("GPP (predicted)") +
  ggtitle("Model trained globally (Balanced, SMOTER)")


ggarrange(imbalance, undersampling, oversampling, smoter, ncol = 1)


# GPP and GPPsat prediction global without balancing ----

# Metrics:
# Regression using entire dataset ()
# leaved-one-site-out cross validation
# k-fold cross validation (5 or 10 k)

predictors.v2 <- predictors.v1[gpp.vsurf$varselect.pred]

# using the entire data set 
gpp.rf <- randomForest(x = sentinel2.efps[,predictors.v2], y = sentinel2.efps$GPP.day.smooth, importance = T)

sentinel2.efps$gpp.predict <- predict(gpp.rf, sentinel2.efps[,predictors.v2])

modelEval(obs = sentinel2.efps$GPP.day.smooth, mod = sentinel2.efps$gpp.predict)

sentinel2.efps %>% 
  ggplot(aes(x = GPP.day.smooth, y = gpp.predict)) +
  geom_hex() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw()

# leave-one-side-out cross validation

sites.names <- sentinel2.efps %>%
  group_by(site_name) %>%
  summarise(igbp = unique(igbp))

# sentinel2.efps$gpp.predict.cv.one <- NA

sites.names$site_name

for (i in 1:length(sites.names$site_name)) {
  training <- sentinel2.efps %>%
    filter(site_name != sites.names$site_name[i]) %>%
    as.data.frame()

  testing <- sentinel2.efps %>%
    filter(site_name == sites.names$site_name[i]) %>%
    as.data.frame()

  regressor <- randomForest(y = training$GPP.day.smooth, x = training[,predictors.v2] , importance = T)

  rf_testing <- predict(regressor, testing[,predictors.v2])
  
  sentinel2.efps$gpp.predict.cv.one[which(sentinel2.efps$site_name == sites.names$site_name[i])] <- rf_testing
}

modelEval(obs = sentinel2.efps$GPP.day.smooth, mod = sentinel2.efps$gpp.predict.cv.one)

sentinel2.efps %>% 
  ggplot(aes(x = GPP.day.smooth, y = gpp.predict.cv.one)) +
  geom_hex() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  theme_bw()

# k-fold cross validation

fitcontrol <- trainControl(method = "cv", number = 10, search = "random", savePredictions = T)

gpp.kfold <- caret::train(x = sentinel2.efps[,predictors.v2], y =  sentinel2.efps$GPP.day.smooth, method = "rf", trControl = fitcontrol, ntree = 500, tuneLength = 10)

gpp.kfold


