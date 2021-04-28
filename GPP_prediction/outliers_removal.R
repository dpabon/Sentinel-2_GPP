# the outliers removal has three steps
# 1. Automatic outlier detection using time series decomposition over the spectral bands.
# 2. Clouds filter at 900 meters radious for the buffer area with a threshold of 0.7.
# 3. Manual outlier removal (For inconsistences observed manually per site).

library(tidyverse)
library(anomalize)
library(plotly)
library(lubridate)
library(outliers)

method.out <- "iqr"

# original dataset  

sentinel2.efps <- read_csv("~/dpabon/data_EFPS_ML/sentinel2.efpsv4.csv")
sentinel2.efps %>% nrow()


sentinel2.efps <- sentinel2.efps %>% 
  drop_na()

sentinel2.efps %>% nrow()


# 1. Automatic outlier detection using time series decomposition over the spectral bands. ----

sites_detection <- unique(sentinel2.efps$site_name)

sites_detection <- sites_detection[-c(53, 55, 59)]

test <- list()

for(i in 1:length(sites_detection)){
  test[[i]] <- sentinel2.efps %>% 
    filter(site_name == sites_detection[i]) %>% 
    dplyr::select(date_index, B1_mean, B2_mean, B3_mean, B4_mean, B5_mean, B11_mean, B12_mean) %>%
    gather("variable", "value", -date_index) %>% 
    group_by(variable) %>% 
    time_decompose(value, method = "stl") %>% 
    anomalize(remainder, method = method.out, alpha = 0.05, max_anoms = 0.2) %>% 
    mutate(site_name = sites_detection[i]) %>% 
    group_by(date_index) %>% 
    summarise(anomaly = length(which(anomaly == "Yes")),
              site_name = unique(site_name))
}


for(i in 1:length(sites_detection)){
  tets <- sentinel2.efps %>% 
    filter(site_name == sites_detection[i]) %>% 
    dplyr::select(date_index, B1_mean, B2_mean, B3_mean, B4_mean, B5_mean, B11_mean, B12_mean) %>%
    gather("variable", "value", -date_index) %>% 
    group_by(variable) %>% 
    time_decompose(value, method = "stl") %>% 
    anomalize(remainder, method = method.out, alpha = 0.05, max_anoms = 0.2) %>%
    time_recompose() %>%
    plot_anomalies(time_recomposed = TRUE, ncol = 1, alpha_dots = 0.25) +
    labs(title = sites_detection[i], subtitle = paste("STL +", method.out, "Methods"))
  
  ggsave(filename = paste0("~/dpabon/results_EFPs_ML/bands_annomalies_detection_", method.out, "/", sites_detection[i], ".png"), plot = tets, height = 12, width = 6)
  ggsave(filename = paste0("~/dpabon/results_EFPs_ML/bands_annomalies_detection_", method.out, "/", sites_detection[i], ".pdf"), plot = tets,  width = 210, height = 297, units = "mm")
}



to.merge <- bind_rows(test)

# some statistics in this step

statistics.step.1 <- to.merge %>% 
  group_by(site_name) %>% 
  summarise(days = length(site_name), number.anomalies = length(which(anomaly != 0))) %>% 
  mutate(percent.anomalies = (number.anomalies * 100)/days)


statistics.step.1 %>% 
  ggplot(aes(x=percent.anomalies)) +
  geom_histogram()

statistics.step.1 %>% 
  summarise(mean(percent.anomalies))

#####

unique(to.merge$site_name)

sentinel2.efps <- full_join(to.merge, sentinel2.efps)

sentinel2.efps$anomaly

# sentinel2.efps <- sentinel2.efps %>%
  # filter(anomaly == 0 | is.na(anomaly) == T)

# 2. Clouds filter at 900 meters radious for the buffer area with a threshold of 0.7. ----

# incoporating cloud percentage at 900 meters 

clouds <- read_csv("sentinel2_900m_info.csv")

clouds <- clouds %>%
  group_by(site_name, date_index) %>%
  summarise_all(mean) %>%
  ungroup()

# merging datasets 

nrow(sentinel2.efps)
sentinel2.efps <- inner_join(sentinel2.efps, clouds)

nrow(sentinel2.efps)

sentinel2.efps %>% 
  dplyr::select(qc_clouds, site_name) %>% 
  ggplot(data = ., aes(x = site_name, y = qc_clouds)) +
  geom_boxplot() +
  theme_bw() +
  scale_x_discrete(guide = guide_axis(check.overlap = TRUE))

ggplotly(sentinel2.efps %>% 
           dplyr::select(qc_clouds, site_name) %>% 
           ggplot(data = ., aes(x = site_name, y = qc_clouds)) +
           geom_boxplot() +
           theme_bw() +
           scale_x_discrete(guide = guide_axis(check.overlap = TRUE)))



sentinel2.efps <- sentinel2.efps %>% 
  filter(qc_clouds < 0.7)

# 3. Manual outlier removal (For clear inconsistences observed manually per site). -----

plots.time.series <- sentinel2.efps %>% 
  group_by(site_name) %>% 
  do(
    plots = dplyr::select(.data=., GPP.day.smooth, GPP.day, ireci_mean, anomaly, date_index, igbp, qc_clouds) %>% 
      gather("variable", "value", GPP.day.smooth, GPP.day, ireci_mean) %>% 
      ggplot(data = ., aes(x = date_index, y = value)) +
      geom_point() +
      geom_line() +
      facet_wrap(ncol = 1, ~variable, scales = "free_y") +
      theme_bw() +
      ggtitle(paste(unique(.$site_name), unique(.$igbp)))
  )


for(i in 1:nrow(plots.time.series)){
  try(ggsave(filename = paste0("~/dpabon/results_EFPs_ML/diagnostic_GPP_GPP.smoothed_IRECI_v4_", method.out, "/", plots.time.series$site_name[i], ".png"), plot = plots.time.series$plots[[i]]))
}

ggplotly(plots.time.series$plots[[which(plots.time.series$site_name == "BE-Bra")]])


# loading manual dataset

manual.outlier <- read_csv("~/dpabon/data_EFPS_ML/manual_filter.csv")

manual.outlier$date_index <- ymd(manual.outlier$date_index)

nrow(sentinel2.efps)
sentinel2.efps <- full_join(sentinel2.efps, manual.outlier)

sentinel2.efps <- sentinel2.efps %>% 
  filter(is.na(manual_outlier) == T) 


## Number of valid images per vegetation type

test.externo <-  sentinel2.efps %>% 
  mutate(go = if_else((anomaly == 0 | is.na(anomaly) == T) & (qc_clouds < 0.7) & (is.na(manual_outlier) == T), TRUE, FALSE)) 


test.externo %>%
  group_by(igbp, go) %>% 
  summarise(number = n()) %>% 
  ggplot(aes(x = igbp, y = number, fill = go)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_discrete(name = "Sentinel-2 images", labels = c("Not Valid", "Valid")) +
  xlab("Vegetation type") +
  ylab("Number of observations") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 650))
  
ggsave("~/dpabon/results_EFPs_ML/plots/gpp_sentinel2_imbalanced_vt.png", width = 8.71, height = 4.85)  


## some statistics for annomaly -----



statistics.step.1 <- sentinel2.efps %>% 
  group_by(site_name) %>% 
  summarise(days = length(site_name), number.anomalies = length(which(anomaly != 0))) %>% 
  mutate(percent.anomalies = (number.anomalies * 100)/days)


statistics.step.1 %>% 
  ggplot(aes(x=percent.anomalies)) +
  geom_histogram() +
  theme_bw() +
  xlab("Percent of images detected as outliers") +
  ylab("Number of sites") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))

ggsave(filename = "~/dpabon/results_EFPs_ML/bands_annomalies_detection_iqr/percent_sites_outliers.png", height = 6, width = 8)

statistics.step.1 %>% 
  summarise(mean(percent.anomalies))

statistics.step.1 %>% 
  filter(number.anomalies > 0) %>% 
  summarise(min(number.anomalies))

statistics.step.1 %>% 
  filter(number.anomalies > 0) %>% 
  summarise(max(number.anomalies))

statistics.step.1 %>% 
  filter(number.anomalies > 0) %>% 
  summarise(mean(number.anomalies))

######


write_csv(sentinel2.efps, paste0("~/dpabon/data_EFPS_ML/sentinel2.efpsv4_wo_outliers_",method.out,".csv"))

sentinel2.efps %>% nrow()
