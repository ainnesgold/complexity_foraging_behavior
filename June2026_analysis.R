library(tidyverse)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(gridExtra)

biodiversity_raw <- read.csv("Biodiversity Surveys.csv")
diet_raw <- read.csv("Diet Categories.csv")
foraging_raw <- read.csv("Foraging Surveys.csv")
video_raw <- read.csv("Video Surveys.csv")

#counting how many video surveys
video_raw %>%
  count(Date, Site, Complexity_Level)


# Merge the biodiversity data set and diet categories by "Genus" and "Species" columns
bd_merged_data <- merge(biodiversity_raw, diet_raw, by = c("Genus", "Species"), all.x = TRUE)

# Remove data points with genus = "chaetodontidae" and species = "unknown"
bd_merged_data <- bd_merged_data %>%
  filter(Genus != "chaetodontidae")

# Fixing the complexity levels to go from lowest to highest
bd_merged_data$Complexity_Level <- factor(bd_merged_data$Complexity_Level, levels = c("Low", "Medium", "High"))

# Remove data points with Site equal to 3
bd_merged_data <- bd_merged_data %>%
  filter(Site != 3)

# Remove "Adult" and "Juvenile" columns
bd_merged_data <- bd_merged_data %>%
  select(-Adult, -Juvenile)

# Getting a total fish count by survey
bd_totals <- bd_merged_data %>%
  group_by(Date, Site, Complexity_Level) %>%
  summarise(total_abundance = sum(Fish_Quantity))








########## Foraging Behavior Survey Data ###########
# Merge the foraging data sets and diet categories by the 'Genus' and 'Species' columns
fd_merged_data <- merge(foraging_raw, diet_raw, by = c("Genus", "Species"), all.x = TRUE)

# Remove rows where Site is equal to 3
fd_merged_data <- fd_merged_data %>%
  filter(Site != 3)

# Merge foraging and biodiversity data by Date, Site, and Complexity_Level
fd_merged_data <- merge(fd_merged_data, bd_totals, by = c("Date", "Site", "Complexity_Level"),
                        all = TRUE)

# Remove rows with NA values
fd_merged_data <- fd_merged_data[complete.cases(fd_merged_data), ]








########## Video Survey Data ###########
# Remove specified columns from video_raw
video_raw <- subset(video_raw, select = -c(Notes, X, X.1))

# Changing the complexity levels from numeric to character variables
video_raw$Complexity_Level <- factor(video_raw$Complexity_Level,
                                     levels = c(3, 2, 1),
                                     labels = c("Low", "Medium", "High"))

# Merge the data sets by "Genus" and "Species" columns
vd_merged_data <- merge(video_raw, diet_raw, by = c("Genus", "Species"), all.x = TRUE)

# Combining the video data set by fixed columns
vd_merged_data <- merge(vd_merged_data, bd_totals, by = c("Date", "Site", "Complexity_Level"),
                        all = TRUE)

# Remove rows with NA values
vd_merged_data <- vd_merged_data[complete.cases(vd_merged_data), ]




########### Merging In-situ & Video Data Sets ###########
# Creating new columns to accommodate extra observation between video and foraging data sets
fd_merged_data$Type <- "Observer"
vd_merged_data$Type <- "Video"

# Ensure df1 and df2 have the same column names
common_columns <- intersect(names(fd_merged_data), names(vd_merged_data))

# Add missing columns to df1
missing_columns_df1 <- setdiff(names(vd_merged_data), common_columns)
fd_merged_data[missing_columns_df1] <- NA

# Add missing columns to df2
missing_columns_df2 <- setdiff(names(fd_merged_data), common_columns)
vd_merged_data[missing_columns_df2] <- NA

# Reorder columns to match
fd_merged_data <- fd_merged_data[names(vd_merged_data)]
vd_merged_data <- vd_merged_data[names(fd_merged_data)]  

# Combining data sets
full_data <- rbind(fd_merged_data, vd_merged_data)

# mutation the RTS_Location column to reflect distance of RTS_location from reef in meters
full_data <- full_data %>%
  mutate(RTS_Location = (RTS_Location) * 0.5)

# Removing columns with missing data between foraging and video data sets
clean_data <- full_data %>%
  select(-RTS, -Camera_Position, -Limu_Bait, -TimeIn_Frame, -TimeOut_Frame, -Time_Seconds)

# Fixing the complexity levels to go from lowest to highest
clean_data$Complexity_Level <- factor(clean_data$Complexity_Level, levels = c("Low", "Medium", "High"))


add_trailing_s <- function(x) {
  ifelse(grepl("s$", x), x, paste0(x, "s"))
}

clean_data$Diet_Category <- add_trailing_s(clean_data$Diet_Category)




########### Species Richness & Abundance ###########

# Calculate species richness and abundance
richness_abundance_sum <- clean_data %>%
  group_by(Complexity_Level, Site, Date) %>%
  summarise(richness = n_distinct(paste(Genus, Species)),
            abundance = sum(Fish_Quantity),
            shannon_index = -sum(Fish_Quantity/sum(Fish_Quantity) * log(Fish_Quantity/sum(Fish_Quantity))))


# Merge the richness and abundance data with si_clean_data
clean_data <- clean_data %>%
  left_join(richness_abundance_sum, by = c("Complexity_Level", "Site", "Date")) 

# Remove rows with NA values
clean_data <- na.omit(clean_data)


####### Sub-setting abundance for only HMD (herbivore diet category)
HMD_abundance <- clean_data[clean_data$Diet_Category == "Herbivores, microvores, and detritivores",] %>%
  group_by(Complexity_Level, Site, Date) %>%
  summarise(HMD_richness = n_distinct(paste(Genus, Species)),
            HMD_abundance = sum(Fish_Quantity))

# Looking at the effect of complexity level on abundance, richness, shannon index, HMD abundance
richness_abundance_sum
HMD_abundance
richness_abundance_sum <- left_join(richness_abundance_sum, HMD_abundance, by = c("Complexity_Level", "Site", "Date"))

# Replace all NA values with 0, since those were dates when there were no HMD observed
richness_abundance_sum[is.na(richness_abundance_sum)] <- 0






################################# continuous complexity data ###############################################################
source("complexity_metrics.R")

complexity_sum <- RDH_allreef_sum
#complexity_sum <- RDH_allplot_sum



## Summary figure of complexity and fish community data


#log abundances for visualizations
richness_abundance_sum$log_abundance <- log(richness_abundance_sum$abundance)
richness_abundance_sum$log_HMD_abundance <- log(richness_abundance_sum$HMD_abundance)

community_long <- richness_abundance_sum %>%
  pivot_longer(
    cols = c(richness, log_abundance, log_HMD_abundance),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      "richness" = "Richness",
      "log_abundance" = "log(Abundance)",
      "log_HMD_abundance" = "log(Herbivore Abundance)"
    ))

fish_community_plot <- ggplot(community_long, aes(x = as.factor(Site), y = Value, fill = as.factor(Site))) +
  geom_boxplot()+
  geom_jitter(
    width = 0.15,
    alpha = 0.4,
    size = 2
  ) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(
    values = c("#0072B2", "#D55E00"),  # Okabe-Ito palette
    name = ""
  ) +
  labs(
    x = "Site",
    y = "Value"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  )


summary_stats_plot <- ggarrange(fish_community_plot, allreefplot_v2, nrow=2, ncol=1)


ggsave(
  filename = "outputs/summary_stats_plot.pdf",
  plot = summary_stats_plot,
  device = "pdf",
  width = 10,
  height = 10
)





################################# complexity + fish community ###############################################################

## merging two datasets
merged_df <- richness_abundance_sum %>%
  left_join(
    complexity_sum,
    by = c("Site", "Complexity_Level")
  )

######## plots ###########
## richness
p1<-ggplot(merged_df, aes(x = HeightRange, y = richness)) +
  geom_point(color = "#0072B2", shape = 16) +
  labs(
    title = "A.",
    x = "Height Range",
    y = "Species Richness"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

p2<-ggplot(merged_df, aes(x = Rugosity, y = richness)) +
  geom_point(color = "#E69F00", shape = 16) +
  labs(
    title = "B.",
    x = "Rugosity",
    y = "Species Richness"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


p3<-ggplot(merged_df, aes(x = FractalDimension, y = richness)) +
  geom_point(color = "#009E73", shape = 16) +
  labs(
    title = "C.",
    x = "Fractal Dimension",
    y = "Species Richness"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

## shannon

p4<-ggplot(merged_df, aes(x = HeightRange, y = shannon_index)) +
  geom_point(color = "#0072B2", shape = 17) +
  labs(
    title = "D.",
    x = "Height Range",
    y = "Shannon Index"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

p5<-ggplot(merged_df, aes(x = Rugosity, y = shannon_index)) +
  geom_point(color = "#E69F00", shape = 17) +
  labs(
    title = "E.",
    x = "Rugosity",
    y = "Shannon Index"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


p6<-ggplot(merged_df, aes(x = FractalDimension, y = shannon_index)) +
  geom_point(color = "#009E73", shape = 17) +
  labs(
    title = "F.",
    x = "Fractal Dimension",
    y = "Shannon Index"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


## abundance
p7<-ggplot(merged_df, aes(x = HeightRange, y = log(abundance))) +
  geom_point(color = "#0072B2", shape = 15) +
  labs(
    title = "G.",
    x = "Height Range",
    y = "log(Overall Abundance)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

p8<-ggplot(merged_df, aes(x = Rugosity, y = log(abundance))) +
  geom_point(color = "#E69F00", shape = 15) +
  labs(
    title = "H.",
    x = "Rugosity",
    y = "log(Overall Abundance)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

p9<-ggplot(merged_df, aes(x = FractalDimension, y = log(abundance))) +
  geom_point(color = "#009E73", shape = 15) +
  labs(
    title = "I.",
    x = "Fractal Dimension",
    y = "log(Overall Abundance)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

## herbivore abundance

p10<-ggplot(merged_df, aes(x = HeightRange, y = log1p(HMD_abundance))) +
  geom_point(color = "#0072B2", shape = 18) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  )+
  labs(
    title = "J.",
    x = "Height Range",
    y = "log(Herbivore Abundance) + 1"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


p11<-ggplot(merged_df, aes(x = Rugosity, y = log1p(HMD_abundance))) +
  geom_point(color = "#E69F00", shape = 18) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  )+
  labs(
    title = "K.",
    x = "Rugosity",
    y = "log(Herbivore Abundance) + 1"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


p12<-ggplot(merged_df, aes(x = FractalDimension, y = log1p(HMD_abundance))) +
  geom_point(color = "#009E73", shape = 18) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.05))
  )+
  labs(
    title = "L.",
    x = "Fractal Dimension",
    y = "log(Herbivore Abundance) + 1"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )



complexity_community_plot <- ggarrange(p1 + rremove("xlab"), p2 + rremove("xlab") + rremove("ylab"), p3 + rremove("xlab") + rremove("ylab"),
          p4 + rremove("xlab"), p5 + rremove("xlab") + rremove("ylab"), p6 + rremove("xlab") + rremove("ylab"),
          p7 + rremove("xlab"), p8 + rremove("xlab") + rremove("ylab"), p9 + rremove("xlab") + rremove("ylab"),
          p10, p11 + rremove("ylab"), p12 + rremove("ylab"), nrow = 4, ncol = 3)


# ggsave(
#   filename = "outputs/complexity_community_plot.pdf",
#   plot = complexity_community_plot,
#   device = "pdf",
#   width = 15,
#   height = 15
# )


######## models ###########

#using glmer for count data: richness, abundance, HMD abundance
  #poisson if there is no overdispersion, negative binomial if there is
#use lme with normal distribution for continuous, normal data: shannon index

hist(merged_df$richness) #skewed
hist(merged_df$shannon_index) #normal
hist(merged_df$abundance) #skewed
hist(merged_df$HMD_abundance) #skewed


# Effect of complexity on richness

#overdispersion if var >> mean
mean(merged_df$richness)
var(merged_df$richness)

#they are very similar - probably use poisson distribution 

#Site is removed, no effect and lower AIC without it

richness_model_1 <- glmer(
  richness ~ HeightRange + (1 | Date),
  family = poisson,
  data = merged_df
)

summary(richness_model_1)

AIC(richness_model_1)


#overdispersion test
sim_pois <- simulateResiduals(richness_model_1)
plot(sim_pois)
testDispersion(sim_pois)

#large p value so no overdispersion, poisson model is good


richness_model_2 <- glmer(
  richness ~ Rugosity + (1 | Date),
  family = poisson,
  data = merged_df
)

summary(richness_model_2)
AIC(richness_model_2)


richness_model_3 <- glmer(
  richness ~ HeightRange*Rugosity + (1 | Date),
  family = poisson,
  data = merged_df
)

summary(richness_model_3)
AIC(richness_model_3)


richness_model_4 <- glmer(
  richness ~ FractalDimension + (1 | Date),
  family = poisson,
  data = merged_df
)

summary(richness_model_4)
AIC(richness_model_4)






# Effect of complexity on shannon: lmer
#no significant effects so not included in manuscript since no differences from richness in results

shannon_model_1 <- lmer(shannon_index ~ HeightRange + 
                           (1 | Date),
                         data = merged_df)
summary(shannon_model_1)


shannon_model_2 <- lmer(shannon_index ~ Rugosity + 
                          (1 | Date),
                        data = merged_df)
summary(shannon_model_2)

shannon_model_3 <- lmer(shannon_index ~ HeightRange*Rugosity + 
                          (1 | Date),
                        data = merged_df)
summary(shannon_model_3)

shannon_model_4 <- lmer(shannon_index ~ FractalDimension + 
                          (1 | Date),
                        data = merged_df)
summary(shannon_model_4)



# Effect of complexity on abundance

## if var >> mean, means it is overdispersed and need negative binomial distribution

mean(merged_df$abundance)
var(merged_df$abundance)
##var is >> mean, probably negative binomial distribution


#removed site as a random effect, lower model AIC without it and it did not account for much variance. 
#No change in significant results once site was removed

abundance_model_1 <- glmmTMB(
  abundance ~ HeightRange + (1 | Date),
  family = nbinom2,
  data = merged_df
)

summary(abundance_model_1)
AIC(abundance_model_1)

#test overdispersion
sim <- simulateResiduals(abundance_model_1)
testDispersion(sim)


abundance_model_2 <- glmmTMB(
  abundance ~ Rugosity + 
    (1 | Date),
  family = nbinom2,
  data = merged_df
)

summary(abundance_model_2)
AIC(abundance_model_2)


abundance_model_3 <- glmmTMB(
  abundance ~ HeightRange*Rugosity + (1 | Date),
  family = nbinom2,
  data = merged_df
)

summary(abundance_model_3)
AIC(abundance_model_3)


abundance_model_4 <- glmmTMB(
  abundance ~ FractalDimension + (1 | Date),
  family = nbinom2,
  data = merged_df
)

summary(abundance_model_4)
AIC(abundance_model_4)




# Effect of complexity on herbivore abundance


mean(merged_df$HMD_abundance)
var(merged_df$HMD_abundance)
##var is >> mean, probably negative binomial distribution


#removed site as a random effect, lower model AIC without it and it did not account for much variance. 
#No change in significant results once site was removed

HMD_abundance_model_1 <- glmmTMB(
  HMD_abundance ~ HeightRange + (1 | Date),
  family = nbinom2,
  data = merged_df
)

summary(HMD_abundance_model_1)
AIC(HMD_abundance_model_1)

#test overdispersion
sim <- simulateResiduals(HMD_abundance_model_1)
testDispersion(sim)

HMD_abundance_model_2 <- glmmTMB(
  HMD_abundance ~ Rugosity + (1 | Date),
  family = nbinom2,
  data = merged_df
)

summary(HMD_abundance_model_2)
AIC(HMD_abundance_model_2)



HMD_abundance_model_3 <- glmmTMB(
  HMD_abundance ~ HeightRange*Rugosity + (1 | Date),
  family = nbinom2,
  data = merged_df
)

summary(HMD_abundance_model_3)
AIC(HMD_abundance_model_3)

HMD_abundance_model_4 <- glmmTMB(
  HMD_abundance ~ FractalDimension + (1 | Date),
  family = nbinom2,
  data = merged_df
)

summary(HMD_abundance_model_4)
AIC(HMD_abundance_model_4)



### Figures of significant predictors

#Abundance ~ Rugosity
pred_abundance <- ggpredict(
  abundance_model_2,
  terms = "Rugosity [all]"
)

abundance_p1 <- ggplot() +
  # raw data
  geom_jitter(
    data = merged_df, aes(x = Rugosity, y = abundance),
    height = 0.02, width = 0, alpha = 0.2) +
  # confidence interval
  geom_ribbon(
    data = pred_abundance,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.2) +
  # model prediction
  geom_line(
    data = pred_abundance, aes(x = x, y = predicted), linewidth = 1.2) +
  labs(
    title = "",
    x = "Rugosity",
    y = "Abundance"
  ) +
  theme_minimal(base_size = 16)

#Abundance ~ Height Range * Rugosity
newdat <- expand.grid(
  HeightRange = seq(min(merged_df$HeightRange, na.rm = TRUE),
                    max(merged_df$HeightRange, na.rm = TRUE),
                    length.out = 100),
  Rugosity = seq(min(merged_df$Rugosity, na.rm = TRUE),
                 max(merged_df$Rugosity, na.rm = TRUE),
                 length.out = 100),
  Date = merged_df$Date[1]  # pick one level
)

newdat$pred <- predict(abundance_model_3, newdata = newdat, type = "response")

abundance_p2 <- ggplot(newdat, aes(x = Rugosity, y = HeightRange, fill = pred)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +  # colorblind-friendly
  labs(
    x = "Rugosity",
    y = "Height Range",
    fill = "Abundance",
    title = ""
  ) +
  theme_minimal(base_size = 16)



#Herbivore abundance ~ height range
pred_HMD_abundance <- ggpredict(
  HMD_abundance_model_3,
  terms = "HeightRange [all]"
)

HMD_abundance_p1 <- ggplot() +
  # raw data
  geom_jitter(
    data = merged_df, aes(x = HeightRange, y = HMD_abundance),
    height = 0.02, width = 0, alpha = 0.2) +
  # confidence interval
  geom_ribbon(
    data = pred_HMD_abundance,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.2) +
  # model prediction
  geom_line(
    data = pred_HMD_abundance, aes(x = x, y = predicted), linewidth = 1.2) +
  labs(
    title = "",
    x = "Height Range",
    y = "Herbivore Abundance"
  ) +
  theme_minimal(base_size = 16)


#Herbivore abundance ~ rugosity
pred_HMD_abundance <- ggpredict(
  HMD_abundance_model_3,
  terms = "Rugosity [all]"
)

HMD_abundance_p2 <- ggplot() +
  # raw data
  geom_jitter(
    data = merged_df, aes(x = Rugosity, y = HMD_abundance),
    height = 0.02, width = 0, alpha = 0.2) +
  # confidence interval
  geom_ribbon(
    data = pred_HMD_abundance,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    alpha = 0.2) +
  # model prediction
  geom_line(
    data = pred_HMD_abundance, aes(x = x, y = predicted), linewidth = 1.2) +
  labs(
    title = "",
    x = "Rugosity",
    y = "Herbivore Abundance"
  ) +
  theme_minimal(base_size = 16)


#Herbivore abundance ~ height range * rugosity
newdat <- expand.grid(
  HeightRange = seq(min(merged_df$HeightRange, na.rm = TRUE),
                    max(merged_df$HeightRange, na.rm = TRUE),
                    length.out = 100),
  Rugosity = seq(min(merged_df$Rugosity, na.rm = TRUE),
                 max(merged_df$Rugosity, na.rm = TRUE),
                 length.out = 100),
  Date = merged_df$Date[1]  # pick one level
)

newdat$pred <- predict(HMD_abundance_model_3, newdata = newdat, type = "response")

HMD_abundance_p3 <- ggplot(newdat, aes(x = Rugosity, y = HeightRange, fill = pred)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C") +  # colorblind-friendly
  labs(
    x = "Rugosity",
    y = "Height Range",
    fill = "Herbivore\nAbundance",
    title = ""
  ) +
  theme_minimal(base_size = 16)

blank <- ggplot() + theme_void()

fish_community_plot <- ggarrange(abundance_p2, HMD_abundance_p3, blank, abundance_p1, HMD_abundance_p1, HMD_abundance_p2, nrow = 2, ncol = 3)

ggsave(
  filename = "outputs/fish_community_plot.pdf",
  plot = fish_community_plot,
  device = "pdf",
  width = 14,
  height = 6
)












############################################# complexity and foraging behavior #######################################################


################################# foraging freq ###############################################################

behavior_frequencies <- clean_data %>%
  group_by(RTS_Location, Complexity_Level, Diet_Category, Site, Date) %>%
  summarise(foraging_frequency = sum(Behavior == "foraging") / n(),
            swimming_frequency = sum(Behavior == "swim") / n(),
            .groups = 'drop')


## merging two datasets

behavior_merged <- behavior_frequencies %>%
  left_join(
    complexity_sum,
    by = c("Site", "Complexity_Level")
  )




#################### models ####################

#Foraging freq

#response variable is proportion, with many zeros. Use beta family.

hist(behavior_merged$foraging_frequency)

# transform data since beta distribution cannot handle 0 and 1s. scale everything so it fits in the distribution.
n <- nrow(behavior_merged)
behavior_merged$ff_beta <- (behavior_merged$foraging_frequency * (n - 1) + 0.5) / n

#i have removed random effects because they do not explain variance and lower AIC without them. No change in significance with or without them
freq_model_1 <- glmmTMB(
  ff_beta ~ HeightRange + RTS_Location + Diet_Category,
  family = beta_family(),
  data = behavior_merged
)

summary(freq_model_1)
AIC(freq_model_1)



#rugosity
freq_model_2 <- glmmTMB(
  ff_beta ~ Rugosity + RTS_Location + Diet_Category,
  family = beta_family(),
  data = behavior_merged
)

summary(freq_model_2)
AIC(freq_model_2)


#interaction model?
freq_model_3 <- glmmTMB(
  ff_beta ~ HeightRange * Rugosity + RTS_Location + Diet_Category,
  family = beta_family(),
  data = behavior_merged
)

summary(freq_model_3)
AIC(freq_model_3)


#visualize rugosity and height range x foraging freq
pred_height <- ggpredict(
  freq_model_3,
  terms = "HeightRange [all]"
)

freq_p1 <- ggplot() +
  # raw data
  geom_jitter(
    data = behavior_merged,
    aes(x = HeightRange,
        y = ff_beta),
    height = 0.02, width = 0, alpha = 0.2) +
  # confidence interval
  geom_ribbon(data = pred_height, aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  # model prediction
  geom_line(
    data = pred_height,
    aes(x = x, y = predicted), linewidth = 1.2) +
  labs(title = "A.", x = "Height Range", y = "Foraging Frequency") +
  theme_minimal(base_size = 16)



pred_rug <- ggpredict(
  freq_model_3,
  terms = "Rugosity [all]"
)

freq_p2 <- ggplot() +
  # raw data
  geom_jitter(
    data = behavior_merged,
    aes(x = Rugosity,
        y = ff_beta),
    height = 0.02, width = 0, alpha = 0.2) +
  # confidence interval
  geom_ribbon(data = pred_rug, aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.2) +
  # model prediction
  geom_line(
    data = pred_rug,
    aes(x = x, y = predicted), linewidth = 1.2) +
  labs(title = "B.", x = "Rugosity", y = "Foraging Frequency") +
  theme_minimal(base_size = 16)



#visualize interaction
newdat <- expand.grid(
  HeightRange = seq(min(behavior_merged$HeightRange, na.rm = TRUE),
                    max(behavior_merged$HeightRange, na.rm = TRUE),
                    length.out = 100),
  Rugosity = seq(min(behavior_merged$Rugosity, na.rm = TRUE),
                 max(behavior_merged$Rugosity, na.rm = TRUE),
                 length.out = 100),
  RTS_Location = mean(behavior_merged$RTS_Location, na.rm = TRUE),
  Diet_Category = behavior_merged$Diet_Category[1]  # pick reference level
)

newdat$pred <- predict(freq_model_3, newdata = newdat, type = "response")

freq_p3 <- ggplot(newdat, aes(x = Rugosity, y = HeightRange, fill = pred)) +
  geom_tile() +
  # optional: raw data on top
  # geom_point(
  #   data = behavior_merged,
  #   aes(x = HeightRange, y = Rugosity,
  #   inherit.aes = FALSE,
  #   alpha = 0.3,
  #   size = 1
  # ) +
  scale_fill_viridis_c(option = "C") +  # colorblind-friendly
  labs(
    x = "Rugosity",
    y = "Height Range",
    fill = "Foraging Frequency",
    title = "C."
  ) +
  theme_minimal(base_size = 16)








freq_model_4 <- glmmTMB(
  ff_beta ~ FractalDimension + RTS_Location + Diet_Category,
  family = beta_family(),
  data = behavior_merged
)

summary(freq_model_4)
AIC(freq_model_4)

pred_fd <- ggpredict(
  freq_model_4,
  terms = "FractalDimension [all]"
)

freq_p4 <- ggplot() +
  # raw data
  geom_jitter(
    data = behavior_merged,
    aes(x = FractalDimension,
        y = ff_beta),
    height = 0.02,
    width = 0,
    alpha = 0.2
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_fd,
    aes(x = x,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2
  ) +
  
  # model prediction
  geom_line(
    data = pred_fd,
    aes(x = x,
        y = predicted),
    linewidth = 1.2
  ) +
  
  labs(
    title = "",
    x = "Fractal Dimension",
    y = "Foraging Frequency"
  ) +
  
  theme_minimal(base_size = 16)







frequencyplot <- ggarrange(freq_p1,freq_p2, freq_p3, nrow=1, ncol=3)


ggsave(
  filename = "outputs/frequencyplot.pdf",
  plot = frequencyplot,
  device = "pdf",
  width = 18,
  height = 5
)


ggsave(
  filename = "outputs/freq_fd.pdf",
  plot = freq_p4,
  device = "pdf",
  width = 5,
  height = 5
)









################################# distance from reef ###############################################################


foraging_merged <- clean_data %>%
  left_join(
    complexity_sum,
    by = c("Site", "Complexity_Level")
  )

foraging_merged <- foraging_merged %>%
  mutate(
    Diet_Category = recode(
      Diet_Category,
      "Herbivores, microvores, and detritivores" = "Herbivores"
    )
  )

hist(foraging_merged$RTS_Location)

#glmm, gaussian family

#took out Site random effects - AIC lower without it and no changes to significance

distance_model_1 <- glmmTMB(
  RTS_Location ~ HeightRange + Diet_Category + (1|Date),
  family = gaussian(),
  data = foraging_merged
)

summary(distance_model_1)
AIC(distance_model_1)


sim <- simulateResiduals(distance_model_1)
plot(sim)
testUniformity(sim)
testDispersion(sim)
testZeroInflation(sim)


#do rest of distance models
distance_model_2 <- glmmTMB(
  RTS_Location ~ Rugosity + Diet_Category + (1|Date),
  family = gaussian(),
  data = foraging_merged
)

summary(distance_model_2)
AIC(distance_model_2)


distance_model_3 <- glmmTMB(
  RTS_Location ~ HeightRange*Rugosity + Diet_Category + (1|Date),
  family = gaussian(),
  data = foraging_merged
)
summary(distance_model_3)
AIC(distance_model_3)


distance_model_4 <- glmmTMB(
  RTS_Location ~ FractalDimension + Diet_Category + (1|Date),
  family = gaussian(),
  data = foraging_merged
)

summary(distance_model_4)
AIC(distance_model_4)





#visualization from model 3
#visualization of height range effect
pred_height <- ggpredict(
  distance_model_3,
  terms = "HeightRange [all]"
)

dist_p1 <- ggplot() +
  # raw data
  geom_jitter(
    data = foraging_merged,
    aes(x = HeightRange,
        y = RTS_Location),
    height = 0.02,
    width = 0,
    alpha = 0.2
  ) +
  # confidence interval
  geom_ribbon(
    data = pred_height,
    aes(x = x,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2
  ) +
  # model prediction
  geom_line(
    data = pred_height,
    aes(x = x,
        y = predicted),
    linewidth = 1.2
  ) +
  labs(
    title = "A.",
    x = "Height Range",
    y = "Foraging Distance"
  ) +
  theme_minimal(base_size = 16)


pred_rugosity <- ggpredict(
  distance_model_3,
  terms = "Rugosity [all]"
)

dist_p2 <- ggplot() +
  # raw data
  geom_jitter(
    data = foraging_merged,
    aes(x = Rugosity,
        y = RTS_Location),
    height = 0.02,
    width = 0,
    alpha = 0.2
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_rugosity,
    aes(x = x,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2
  ) +
  
  # model prediction
  geom_line(
    data = pred_rugosity,
    aes(x = x,
        y = predicted),
    linewidth = 1.2
  ) +
  
  labs(
    title = "B.",
    x = "Rugosity",
    y = "Foraging Distance"
  ) +
  
  theme_minimal(base_size = 16)



#visualize interaction
newdat <- expand.grid(
  HeightRange = seq(min(foraging_merged$HeightRange, na.rm = TRUE),
                    max(foraging_merged$HeightRange, na.rm = TRUE),
                    length.out = 100),
  Rugosity = seq(min(foraging_merged$Rugosity, na.rm = TRUE),
                 max(foraging_merged$Rugosity, na.rm = TRUE),
                 length.out = 100),
  Diet_Category = foraging_merged$Diet_Category[1],  # pick reference level
  Date = foraging_merged$Date[1] 
)

newdat$pred <- predict(distance_model_3, newdata = newdat, type = "response")

dist_interaction_plot <- ggplot(newdat, aes(x = Rugosity, y = HeightRange, fill = pred)) +
  geom_tile() +
  
  # optional: raw data on top
  # geom_point(
  #   data = foraging_merged,
  #   aes(x = HeightRange, y = Rugosity,
  #   inherit.aes = FALSE,
  #   alpha = 0.3,
  #   size = 1
  # ) +
  
  scale_fill_viridis_c(option = "C") +  # colorblind-friendly
  labs(
    x = "Rugosity",
    y = "Height Range",
    fill = "Foraging \nDistance",
    title = "C."
  ) +
  theme_minimal(base_size = 16)



combined_dist_plot <- ggarrange(dist_p1, dist_p2, dist_interaction_plot, nrow=1, ncol=3)


ggsave(
  filename = "outputs/combined_dist_plot.pdf",
  plot = combined_dist_plot,
  device = "pdf",
  width = 15,
  height = 5
)





ggsave(
  filename = "outputs/dist_diet_plot_m1.pdf",
  plot = dist_diet_plot_m1,
  device = "pdf",
  width = 10,
  height = 8
)


#post hoc diet category



emm_m3 <- emmeans(distance_model_3, ~ Diet_Category, type = "response")
pairs(emm_m3, adjust = "tukey", type = "response")
emm_df <- as.data.frame(emm_m3)

# Table of significant postdoc Tukey tests
pairs_out <- pairs(emmeans(distance_model_3, ~ Diet_Category),
                   adjust = "tukey")
pairs_df <- as.data.frame(pairs_out)
sig_pairs <- subset(pairs_df, p.value < 0.05)
sig_pairs

sig_table <- sig_pairs %>%
  mutate(
    estimate = round(estimate, 2),
    SE = round(SE, 2),
    z.ratio = round(z.ratio, 2),
    p.value = signif(p.value, 3)
  ) %>%
  select(
    Contrast = contrast,
    Estimate = estimate,
    SE,
    df,
    `z value` = z.ratio,
    `Adjusted p` = p.value
  )


table_plot <- tableGrob(sig_table, rows = NULL)
grid.newpage()
grid.draw(table_plot)
pdf("outputs/distance_m3_diet_table.pdf", width = 10, height = 6)
grid.newpage()
grid.draw(table_plot)
dev.off()




dist_diet_plot_m3 <- ggplot() +
  # RAW DATA
  geom_jitter(
    data = foraging_merged,
    aes(x = Diet_Category, y = RTS_Location),
    alpha = 0.15
  ) +
  # MODEL ESTIMATES (now follow factor order automatically)
  geom_point(
    data = emm_df,
    aes(x = Diet_Category, y = emmean),
    size = 2
  ) +
  geom_errorbar(
    data = emm_df,
    aes(x = Diet_Category,
        ymin = asymp.LCL,
        ymax = asymp.UCL),
    width = 0.2
  ) +
  coord_flip() +
  labs(title = "",
       x = "Diet Category",
       y = "Foraging Distance (m)"
  ) +
  theme_minimal(base_size = 16)




ggsave(
  filename = "outputs/dist_diet_plot_m3.pdf",
  plot = dist_diet_plot_m3,
  device = "pdf",
  width = 6,
  height = 6
)



################################# foraging distance ###############################################################

foraging_only <- foraging_merged %>%
  filter(Behavior == "foraging")

hist(foraging_only$RTS_Location)

#glmm, gaussian family
#took out random effects - AIC lower without them and no changes to significance

forage_dist_model_1 <- glmmTMB(
  RTS_Location ~ HeightRange + Diet_Category + (1|Date),
  family = gaussian(),
  data = foraging_only
)

summary(forage_dist_model_1)
AIC(forage_dist_model_1)

forage_dist_model_2 <- glmmTMB(
  RTS_Location ~ Rugosity + Diet_Category + (1|Date),
  family = gaussian(),
  data = foraging_only
)

summary(forage_dist_model_2)
AIC(forage_dist_model_2)


forage_dist_model_3 <- glmmTMB(
  RTS_Location ~ HeightRange*Rugosity + Diet_Category + (1|Date),
  family = gaussian(),
  data = foraging_only
)

summary(forage_dist_model_3)
AIC(forage_dist_model_3)


forage_dist_model_4 <- glmmTMB(
  RTS_Location ~ FractalDimension + Diet_Category + (1|Date),
  family = gaussian(),
  data = foraging_only
)

summary(forage_dist_model_4)
AIC(forage_dist_model_4)


## model 3 visualizations

#foraging distance ~ height range

pred_height <- ggpredict(
  forage_dist_model_3,
  terms = "HeightRange [all]"
)

foraging_dist_p1 <- ggplot() +
  # raw data
  geom_jitter(
    data = foraging_only,
    aes(x = HeightRange,
        y = RTS_Location),
    height = 0.02,
    width = 0,
    alpha = 0.2
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_height,
    aes(x = x,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2
  ) +
  
  # model prediction
  geom_line(
    data = pred_height,
    aes(x = x,
        y = predicted),
    linewidth = 1.2
  ) +
  
  labs(
    title = "A.",
    x = "Height Range",
    y = "Foraging Distance"
  ) +
  
  theme_minimal(base_size = 16)



#foraging distance ~ rugosity
pred_rugosity <- ggpredict(
  forage_dist_model_3,
  terms = "Rugosity [all]"
)

foraging_dist_p2 <- ggplot() +
  # raw data
  geom_jitter(
    data = foraging_only,
    aes(x = Rugosity,
        y = RTS_Location),
    height = 0.02,
    width = 0,
    alpha = 0.2
  ) +
  
  # confidence interval
  geom_ribbon(
    data = pred_rugosity,
    aes(x = x,
        ymin = conf.low,
        ymax = conf.high),
    alpha = 0.2
  ) +
  
  # model prediction
  geom_line(
    data = pred_rugosity,
    aes(x = x,
        y = predicted),
    linewidth = 1.2
  ) +
  
  labs(
    title = "B.",
    x = "Rugosity",
    y = "Foraging Distance"
  ) +
  
  theme_minimal(base_size = 16)


## combined plot

combined_forage_dist_plot <- ggarrange(foraging_dist_p1, foraging_dist_p2, nrow = 1, ncol=2)

ggsave(
  filename = "outputs/foraging_dist_plot_m3.pdf",
  plot = combined_forage_dist_plot,
  device = "pdf",
  width = 10,
  height = 6
)



#diet post hoc tests
emm_m3 <- emmeans(forage_dist_model_3, ~ Diet_Category, type = "response")
pairs(emm_m3, adjust = "tukey", type = "response")
emm_df <- as.data.frame(emm_m3)

# Table of significant postdoc Tukey tests
pairs_out <- pairs(emmeans(forage_dist_model_3, ~ Diet_Category),
                   adjust = "tukey")
pairs_df <- as.data.frame(pairs_out)
sig_pairs <- subset(pairs_df, p.value < 0.05)
sig_pairs

sig_table <- sig_pairs %>%
  mutate(
    estimate = round(estimate, 2),
    SE = round(SE, 2),
    z.ratio = round(z.ratio, 2),
    p.value = signif(p.value, 3)
  ) %>%
  select(
    Contrast = contrast,
    Estimate = estimate,
    SE,
    df,
    `z value` = z.ratio,
    `Adjusted p` = p.value
  )


table_plot <- tableGrob(sig_table, rows = NULL)
grid.newpage()
grid.draw(table_plot)
pdf("outputs/foraging_distance_m3_diet_table.pdf", width = 10, height = 6)
grid.newpage()
grid.draw(table_plot)
dev.off()




foraging_dist_diet_plot_m3 <- ggplot() +
  # RAW DATA
  geom_jitter(
    data = foraging_only,
    aes(x = Diet_Category, y = RTS_Location),
    alpha = 0.15
  ) +
  # MODEL ESTIMATES (now follow factor order automatically)
  geom_point(
    data = emm_df,
    aes(x = Diet_Category, y = emmean),
    size = 2
  ) +
  geom_errorbar(
    data = emm_df,
    aes(x = Diet_Category,
        ymin = asymp.LCL,
        ymax = asymp.UCL),
    width = 0.2
  ) +
  coord_flip() +
  labs(title = "",
       x = "Diet Category",
       y = "Foraging Distance (m)"
  ) +
  theme_minimal(base_size = 16)




ggsave(
  filename = "outputs/foraging_dist_diet_plot_m3.pdf",
  plot = foraging_dist_diet_plot_m3,
  device = "pdf",
  width = 6,
  height = 6
)











################################# number of bites ###############################################################


bites_df <- clean_data %>%
  filter(Bites > 0)

bites_df <- bites_df %>%
  left_join(
    complexity_sum,
    by = c("Site", "Complexity_Level")
  )
hist(bites_df$Bites)

# var >> mean, so nbinom
mean(bites_df$Bites)
var(bites_df$Bites)

bites_model_1 <- glmmTMB(
  Bites ~ HeightRange + RTS_Location + Diet_Category + (1 | Date),
  family = nbinom2,
  data = bites_df
)

AIC(bites_model_1)
summary(bites_model_1)

#colinearity test, <2 means there isnʻt any
car::vif(glm(Bites ~ HeightRange + Diet_Category + RTS_Location, family=poisson, data=bites_df))


bites_model_2 <- glmmTMB(
  Bites ~ Rugosity + RTS_Location + Diet_Category + (1 | Date),
  family = nbinom2,
  data = bites_df
)

summary(bites_model_2)
AIC(bites_model_2)

bites_model_3 <- glmmTMB(
  Bites ~ HeightRange*Rugosity + RTS_Location + Diet_Category + (1 | Date),
  family = nbinom2,
  data = bites_df
)

summary(bites_model_3)
AIC(bites_model_3)


bites_model_4 <- glmmTMB(
  Bites ~ FractalDimension + RTS_Location + Diet_Category + (1 | Date),
  family = nbinom2,
  data = bites_df
)

summary(bites_model_4)
AIC(bites_model_4)



## visualization

pred_rts <- ggpredict(bites_model_3, terms = "RTS_Location")

bites_rts_plot <- ggplot(pred_rts, aes(x, predicted)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_jitter(data = bites_df,
              aes(x = RTS_Location, y = Bites),
              alpha = 0.15) +
  labs(x = "Distance From Reef (m)", y = "Predicted Bites") +
  theme_minimal(base_size = 16)


ggsave(
  filename = "outputs/bites_rts_plot.pdf",
  plot = bites_rts_plot,
  device = "pdf",
  width = 5,
  height = 5
)



