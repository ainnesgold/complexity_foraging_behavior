library(tidyverse)
library(gridExtra)
library(lmerTest)
library(emmeans)
library(multcompView)
library(ggpubr)
library(grid)

biodiversity_raw <- read.csv("Biodiversity Surveys.csv")
diet_raw <- read.csv("Diet Categories.csv")
foraging_raw <- read.csv("Foraging Surveys.csv")
video_raw <- read.csv("Video Surveys.csv")


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



# Boxplot of Species Richness by Complexity Level
p1 <- ggplot(richness_abundance_sum, aes(x = Complexity_Level, y = richness, fill = Complexity_Level)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.4,
    size = 2
  ) +
  scale_fill_brewer(palette = "Set2") + 
  labs(
    title = "A.",
    x = "Complexity Level",
    y = "Species Richness"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

# Boxplot of Shannon Index by Complexity Level
p2 <- ggplot(richness_abundance_sum, aes(x = Complexity_Level, y = shannon_index, fill = Complexity_Level)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.4,
    size = 2
  ) +
  scale_fill_brewer(palette = "Set2") + 
  labs(
    title = "B.",
    x = "Complexity Level",
    y = "Shannon Index"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

# Boxplot of Abundance by Complexity Level 
#logged because of outliers making it hard to see

p3 <- ggplot(richness_abundance_sum, aes(x = Complexity_Level, y = log(abundance), fill = Complexity_Level)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.4,
    size = 2
  ) +
  scale_fill_brewer(palette = "Set2") + 
  labs(
    title = "C.",
    x = "Complexity Level",
    y = "log(Species Abundance)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )


# Boxplot of Herbivore Abundance by Complexity Level
p4 <- ggplot(richness_abundance_sum, 
             aes(x = Complexity_Level, 
                 y = log(HMD_abundance), 
                 fill = Complexity_Level)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.4,
    size = 2
  ) +
  scale_fill_brewer(palette = "Set2") +
  coord_cartesian(clip = "off") +
  labs(
    title = "D.",
    x = "Complexity Level",
    y = "log(HMD Abundance)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

# Arrange the plots in a 2x2 grid
combined_plot <- grid.arrange(p1 + rremove("xlab"), p2 + rremove("xlab"), p3 + rremove("xlab"), p4 + rremove("xlab"), ncol = 2,
                              bottom = textGrob("Complexity Level", gp = gpar(fontsize = 16)))

ggsave(
  filename = "outputs/combined_plot.pdf",
  plot = combined_plot,
  device = "pdf",
  width = 10,
  height = 8
)






############################### Statistical Analysis ####################################################





###### MODELS #########

# Effect of complexity on richness
richness_model <- aov(richness ~ Complexity_Level, data = richness_abundance_sum)
summary(richness_model)

#LME for including mixed effects
richness_model <- lmer(richness ~ Complexity_Level + 
                (1 | Site) + 
                (1 | Date),
              data = richness_abundance_sum)

anova(richness_model)


# Effect of complexity on shannon
shannon_model <- aov(shannon_index ~ Complexity_Level, data = richness_abundance_sum)
summary(shannon_model)

#LME for including mixed effects
shannon_model <- lmer(shannon_index ~ Complexity_Level + 
                         (1 | Site) + 
                         (1 | Date),
                       data = richness_abundance_sum)

anova(shannon_model)


# Effect of complexity on abundance
abundance_model <- aov(abundance ~ Complexity_Level, data = richness_abundance_sum)
summary(abundance_model)

#LME for including mixed effects
abundance_model <- lmer(abundance ~ Complexity_Level + 
                         (1 | Site) + 
                         (1 | Date),
                       data = richness_abundance_sum)

anova(abundance_model)



# Effect of complexity on HMD abundance
HMD_abundance_model <- aov(HMD_abundance ~ Complexity_Level, data = richness_abundance_sum)
summary(HMD_abundance_model)

#LME for including mixed effects
HMD_abundance_model <- lmer(HMD_abundance ~ Complexity_Level + 
                          (1 | Site) + 
                          (1 | Date),
                        data = richness_abundance_sum)

anova(HMD_abundance_model)










### Remake that first figure based on model estimates

emm_richness <- as.data.frame(emmeans(richness_model, ~ Complexity_Level))
emm_shannon  <- as.data.frame(emmeans(shannon_model, ~ Complexity_Level))
emm_abund    <- as.data.frame(emmeans(abundance_model, ~ Complexity_Level))
emm_hmd      <- as.data.frame(emmeans(HMD_abundance_model, ~ Complexity_Level))


p1 <- ggplot() +
  # RAW DATA (grey so model stands out)
  geom_jitter(
    data = richness_abundance_sum,
    aes(Complexity_Level, richness),
    width = 0.15,
    alpha = 0.4,
    size = 2,
    color = "grey60"
  ) +
  # MODEL MEAN (colored)
  geom_point(
    data = emm_richness,
    aes(Complexity_Level, emmean, color = Complexity_Level),
    size = 4
  ) +
  # MODEL CI (colored)
  geom_errorbar(
    data = emm_richness,
    aes(Complexity_Level,
        ymin = lower.CL,
        ymax = upper.CL,
        color = Complexity_Level),
    width = 0.15
  ) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "A.",
    x = "Complexity Level",
    y = "Species Richness"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")


p2 <- ggplot() +
  geom_jitter(
    data = richness_abundance_sum,
    aes(Complexity_Level, shannon_index),
    width = 0.15,
    alpha = 0.25,
    size = 2
  ) +
  geom_point(
    data = emm_shannon,
    aes(Complexity_Level, emmean, color = Complexity_Level),
    size = 4
  ) +
  geom_errorbar(
    data = emm_shannon,
    aes(Complexity_Level, ymin = lower.CL, ymax = upper.CL, color = Complexity_Level),
    width = 0.15
  ) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "B.",
       x = "Complexity Level",
       y = "Shannon Index") +
  theme_minimal(base_size = 16)  +
  theme(legend.position = "none")


p3 <- ggplot() +
  geom_jitter(
    data = richness_abundance_sum,
    aes(Complexity_Level, abundance),
    width = 0.15,
    alpha = 0.25,
    size = 2
  ) +
  geom_point(
    data = emm_abund,
    aes(Complexity_Level, emmean, color = Complexity_Level),
    size = 3
  ) +
  geom_errorbar(
    data = emm_abund,
    aes(Complexity_Level, ymin = lower.CL, ymax = upper.CL, color = Complexity_Level),
    width = 0.15
  ) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "C.",
       x = "Complexity Level",
       y = "Species Abundance") +
  theme_minimal(base_size = 16)  +
  theme(legend.position = "none")


p4 <- ggplot() +
  geom_jitter(
    data = richness_abundance_sum,
    aes(Complexity_Level, HMD_abundance),
    width = 0.15,
    alpha = 0.25,
    size = 2
  ) +
  geom_point(
    data = emm_hmd,
    aes(Complexity_Level, emmean, color = Complexity_Level),
    size = 3
  ) +
  geom_errorbar(
    data = emm_hmd,
    aes(Complexity_Level, ymin = lower.CL, ymax = upper.CL, color = Complexity_Level),
    width = 0.15
  ) +
  scale_color_brewer(palette = "Set2") +
  labs(title = "D.",
       x = "Complexity Level",
       y = "HMD Abundance") +
  theme_minimal(base_size = 16)  +
  theme(legend.position = "none")


combined_plot <- grid.arrange(
  p1 + rremove("xlab"), p2 + rremove("xlab"), p3+ rremove("xlab"),  p4  + rremove("xlab"),
  ncol = 2,
  bottom = textGrob("Complexity Level",
                    gp = gpar(fontsize = 16))
)

ggsave(
  filename = "outputs/model_based_combined_plot.pdf",
  plot = combined_plot,
  width = 10,
  height = 8
)






####### Foraging results


# Frequency of foraging

# Calculate behavior frequencies

behavior_frequencies <- clean_data %>%
  group_by(RTS_Location, Complexity_Level, Diet_Category, Site, Date) %>%
  summarise(foraging_frequency = sum(Behavior == "foraging") / n(),
            swimming_frequency = sum(Behavior == "swim") / n(),
            .groups = 'drop')


#Foraging freq with Diet
foraging_freq_model <- lmer(foraging_frequency ~ Complexity_Level + RTS_Location + Diet_Category +
                              (1 | Site) + 
                              (1 | Date),
                            data = behavior_frequencies)
anova(foraging_freq_model)

# We see that complexity level and diet category are significant predictors of foraging frequency

# Post hoc tests
emmeans(foraging_freq_model, pairwise ~ Complexity_Level, adjust = "tukey")
# We see that there low and medium complexity are not significantly different, but low and high, and medium and high, are significantly different

emmeans(foraging_freq_model, pairwise ~ Diet_Category, adjust = "tukey")
#We see that most diet groups are not significantly different from one another. Exception is significant difference between crustacivores and planktivores.


## Visualization
emm <- emmeans(foraging_freq_model, ~ Complexity_Level)
pairs(emm, adjust = "tukey")

# Manually putting p values to get significance letters for plot
pvals <- c(
  "Low-Medium" = 0.3760,
  "Low-High" = 0.0002,
  "Medium-High" = 0.0288
)

letters <- multcompLetters(pvals, threshold = 0.05)
letters

cld_df <- data.frame(
  Complexity_Level = c("Low", "Medium", "High"),
  letter = c(
    "a",   # Low (not different from Medium, but different from High)
    "a",   # Medium (differs from High only)
    "b"    # High
  )
)

label_df <- merge(cld_df,
                  aggregate(foraging_frequency ~ Complexity_Level,
                            behavior_frequencies, max),
                  by = "Complexity_Level")



freq_plot <- ggplot(behavior_frequencies, 
       aes(Complexity_Level, y = foraging_frequency, fill = Complexity_Level)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(
    data = label_df,
    aes(x = Complexity_Level,
        y = foraging_frequency * 1.05,
        label = letter),
    size = 6
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "A.",
    x = "Complexity Level",
    y = "Foraging Frequency"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )




#version with model estimate


emm <- emmeans(foraging_freq_model, ~ Complexity_Level)
emm_df <- as.data.frame(emm)

label_df <- merge(cld_df,
                  aggregate(foraging_frequency ~ Complexity_Level,
                            behavior_frequencies, max),
                  by = "Complexity_Level")

freq_plot_v2 <- ggplot() +
  # RAW DATA 
  geom_jitter(
    data = behavior_frequencies,
    aes(x = Complexity_Level, y = foraging_frequency),
    width = 0.15,
    alpha = 0.15,
    size = 2
  ) +
  # MODEL ESTIMATES (colored)
  geom_point(
    data = emm_df,
    aes(x = Complexity_Level, y = emmean, color = Complexity_Level),
    size = 4
  ) +
  # MODEL CI (colored)
  geom_errorbar(
    data = emm_df,
    aes(x = Complexity_Level,
        ymin = lower.CL,
        ymax = upper.CL,
        color = Complexity_Level),
    width = 0.15
  ) +
  # COLOR BLIND FRIENDLY PALETTE
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  geom_text(
    data = label_df,
    aes(x = Complexity_Level,
        y = foraging_frequency * 1.1,
        label = letter),
    size = 6
  ) +
  labs(
    title = "A.",
    x = "Complexity Level",
    y = "Foraging Frequency"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none"
  )















 ###### Foraging distance ######

foraging_dist_model <- lmer(RTS_Location ~ Complexity_Level + Diet_Category +
                              (1 | Site) + 
                              (1 | Date),
                            data = clean_data)
anova(foraging_dist_model)


# Post hoc tests
emmeans(foraging_dist_model, pairwise ~ Complexity_Level, adjust = "tukey")
#We see that there is a significant difference between low and high, and medium and high, but not low and medium.

emmeans(foraging_dist_model, pairwise ~ Diet_Category, adjust = "tukey")
pairs_out <- pairs(emmeans(foraging_dist_model, ~ Diet_Category),
                   adjust = "tukey")
pairs_df <- as.data.frame(pairs_out)
sig_pairs <- subset(pairs_df, p.value < 0.05)
sig_pairs
#We see significant differences in several diet groups, including corallivore-macroinvertivore, corallivore-planktivore, corallivore-sessile invertivore, 
# crustacivore - macroinvertivore, crustacivore - sessile invertivore, HMD - macroinvertivore, HMD - sessile invertivore, macroinvertivore - microinvertivore, 
#macroinvertivore - planktivore, micro-invertivore - sessile invertivore, .


## Visualization
emm <- emmeans(foraging_dist_model, ~ Complexity_Level)
pairs(emm, adjust = "tukey")

# Manually putting p values to get significance letters for plot
pvals <- c(
  "Low-Medium" = 0.8386,
  "Low-High" = 0.0001,
  "Medium-High" = 0.0001
)

letters <- multcompLetters(pvals, threshold = 0.05)
letters

cld_df <- data.frame(
  Complexity_Level = c("Low", "Medium", "High"),
  letter = c(
    "a",   # Low (not different from Medium, but different from High)
    "a",   # Medium (differs from High only)
    "b"    # High
  )
)

label_df <- merge(cld_df,
                  aggregate(RTS_Location ~ Complexity_Level,
                            clean_data, max),
                  by = "Complexity_Level")

dist_plot <- ggplot(clean_data, 
                    aes(Complexity_Level, y = RTS_Location, fill = Complexity_Level)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 2) +
  scale_fill_brewer(palette = "Set2") +
  geom_text(
    data = label_df,
    aes(x = Complexity_Level,
        y = RTS_Location * 1.1,
        label = letter),
    size = 6
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "B.",
    x = "Complexity Level",
    y = "Foraging Distance"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

#version with model estimate


emm <- emmeans(foraging_dist_model, ~ Complexity_Level)
emm_df <- as.data.frame(emm)

dist_plot_v2 <- ggplot() +
  # RAW DATA (keep light / neutral)
  geom_jitter(
    data = clean_data,
    aes(x = Complexity_Level, y = RTS_Location),
    width = 0.15,
    alpha = 0.15,
    size = 2
  ) +
  # MODEL ESTIMATES (colored)
  geom_point(
    data = emm_df,
    aes(x = Complexity_Level, y = emmean, color = Complexity_Level),
    size = 4
  ) +
  # MODEL CI (colored)
  geom_errorbar(
    data = emm_df,
    aes(x = Complexity_Level,
        ymin = lower.CL,
        ymax = upper.CL,
        color = Complexity_Level),
    width = 0.15
  ) +
  # COLOR BLIND FRIENDLY PALETTE
  scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  # LETTERS
  geom_text(
    data = label_df,
    aes(x = Complexity_Level,
        y = RTS_Location * 1.3,
        label = letter),
    size = 6
  ) +
  labs(
    title = "B.",
    x = "Complexity Level",
    y = "Foraging Distance"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none"
  )




#foraging_plot <- ggarrange(freq_plot, dist_plot, nrow = 1, ncol = 2)

foraging_plot <- grid.arrange(freq_plot + rremove("xlab"), dist_plot + rremove("xlab"), ncol = 2, nrow = 1,
                              bottom = textGrob("Complexity Level", gp = gpar(fontsize = 16)))

foraging_plot_v2 <- grid.arrange(freq_plot_v2 + rremove("xlab"), dist_plot_v2 + rremove("xlab"), ncol = 2, nrow = 1,
                              bottom = textGrob("Complexity Level", gp = gpar(fontsize = 16)))


ggsave(
  filename = "outputs/foraging_plot.pdf",
  plot = foraging_plot,
  device = "pdf",
  width = 10,
  height = 6
)

ggsave(
  filename = "outputs/foraging_plot_v2.pdf",
  plot = foraging_plot_v2,
  device = "pdf",
  width = 10,
  height = 6
)



#Diet groups visualization??


emm_diet <- emmeans(foraging_dist_model, ~ Diet_Category)
emm_df <- as.data.frame(emm_diet)

distance_diet_plot <- ggplot() +
  # RAW DATA
  geom_jitter(
    data = clean_data,
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
        ymin = lower.CL,
        ymax = upper.CL),
    width = 0.2
  ) +
  coord_flip() +
  labs(title = "A.",
    x = "Diet Category",
    y = "Foraging Distance"
  ) +
  theme_minimal(base_size = 16)






###### Number of bites, for those that were foraging ######



bites_model <- lmer(Bites ~ Complexity_Level + Diet_Category + RTS_Location +
                              (1 | Site) + 
                              (1 | Date),
                            data = clean_data)
anova(bites_model)

emmeans(bites_model, pairwise ~ Diet_Category, adjust = "tukey")

#Extracting just the significant pairs
pairs_out <- pairs(emmeans(bites_model, ~ Diet_Category),
                   adjust = "tukey")
pairs_df <- as.data.frame(pairs_out)
sig_pairs <- subset(pairs_df, p.value < 0.05)
sig_pairs

# visualization of bites
emm_bites <- emmeans(bites_model, ~ Diet_Category)
emm_df_bites <- as.data.frame(emm_bites)

bites_diet_plot <- ggplot() +
  # RAW DATA
  geom_jitter(
    data = clean_data,
    aes(x = Diet_Category, y = Bites),
    alpha = 0.15
  ) +
  # MODEL ESTIMATES
  geom_point(
    data = emm_df_bites,
    aes(x = Diet_Category, y = emmean),
    size = 2
  ) +
  geom_errorbar(
    data = emm_df_bites,
    aes(x = Diet_Category,
        ymin = lower.CL,
        ymax = upper.CL),
    width = 0.2
  ) +
  coord_flip() +
  labs(title = "B.",
    x = "Diet Category",
    y = "Number of Bites"
  ) +
  theme_minimal(base_size = 16)


diet_plot <- grid.arrange(distance_diet_plot + rremove("ylab"), bites_diet_plot + rremove("ylab"), ncol = 1, nrow = 2)

ggsave(
  filename = "outputs/diet_plot.pdf",
  plot = diet_plot,
  device = "pdf",
  width = 10,
  height = 10
)
