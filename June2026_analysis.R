library(tidyverse)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(glmmTMB)



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
  geom_point() +
  geom_smooth(method = "lm") +
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
  geom_point() +
  geom_smooth(method = "lm") +
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
  geom_point() +
  geom_smooth(method = "lm") +
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
  geom_point() +
  geom_smooth(method = "lm") +
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
  geom_point() +
  geom_smooth(method = "lm") +
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
  geom_point() +
  geom_smooth(method = "lm") +
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
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "G.",
    x = "Height Range",
    y = "log(Species Abundance)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

p8<-ggplot(merged_df, aes(x = Rugosity, y = log(abundance))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "H.",
    x = "Rugosity",
    y = "log(Species Abundance)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )

p9<-ggplot(merged_df, aes(x = FractalDimension, y = log(abundance))) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "I.",
    x = "Fractal Dimension",
    y = "log(Species Abundance)"
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
  geom_point() +
  geom_smooth(method = "lm") +
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
  geom_point() +
  geom_smooth(method = "lm") +
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
  geom_point() +
  geom_smooth(method = "lm") +
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


ggsave(
  filename = "outputs/complexity_community_plot.pdf",
  plot = complexity_community_plot,
  device = "pdf",
  width = 15,
  height = 15
)


######## models ###########

# Effect of complexity on richness
# richness_model_1 <- lmer(richness ~ HeightRange + 
#                          (1 | Site) + 
#                          (1 | Date),
#                        data = merged_df)
# anova(richness_model_1)

hist(merged_df$richness)


richness_model_1 <- glmer(
  richness ~ HeightRange + (1 | Site) + (1 | Date),
  family = poisson,
  data = merged_df
)

summary(richness_model_1)

# library(DHARMa)
# sim_pois <- simulateResiduals(m_pois)
# plot(sim_pois)
# testDispersion(sim_pois)


richness_model_2 <- glmer(
  richness ~ Rugosity + (1 | Site) + (1 | Date),
  family = poisson,
  data = merged_df
)

summary(richness_model_2)


richness_model_3 <- glmer(
  richness ~ FractalDimension + (1 | Site) + (1 | Date),
  family = poisson,
  data = merged_df
)

summary(richness_model_3)















######## EDIT MODELS STARTING HERE - change to poisson models


# Effect of complexity on shannon
shannon_model_1 <- lmer(shannon_index ~ HeightRange + 
                           (1 | Site) + 
                           (1 | Date),
                         data = merged_df)
anova(shannon_model_1)

shannon_model_2 <- lmer(shannon_index ~ Rugosity + 
                          (1 | Site) + 
                          (1 | Date),
                        data = merged_df)
anova(shannon_model_2)

shannon_model_3 <- lmer(shannon_index ~ FractalDimension + 
                          (1 | Site) + 
                          (1 | Date),
                        data = merged_df)
anova(shannon_model_3)

shannon_model_4 <- lmer(shannon_index ~ Rugosity * HeightRange + 
                          (1 | Site) + 
                          (1 | Date),
                        data = merged_df)
anova(shannon_model_4)


# Effect of complexity on abundance
abundance_model_1 <- lmer(abundance ~ HeightRange + 
                          (1 | Site) + 
                          (1 | Date),
                        data = merged_df)
anova(abundance_model_1)

abundance_model_2 <- lmer(abundance ~ Rugosity + 
                            (1 | Site) + 
                            (1 | Date),
                          data = merged_df)
anova(abundance_model_2)

abundance_model_3 <- lmer(abundance ~ FractalDimension + 
                            (1 | Site) + 
                            (1 | Date),
                          data = merged_df)
anova(abundance_model_3)

abundance_model_4 <- lmer(abundance ~ Rugosity * HeightRange + 
                            (1 | Site) + 
                            (1 | Date),
                          data = merged_df)
anova(abundance_model_4)

# Effect of complexity on herbivore abundance
HMD_abundance_model_1 <- lmer(HMD_abundance ~ HeightRange + 
                            (1 | Site) + 
                            (1 | Date),
                          data = merged_df)
anova(HMD_abundance_model_1)

HMD_abundance_model_2 <- lmer(HMD_abundance ~ Rugosity + 
                                (1 | Site) + 
                                (1 | Date),
                              data = merged_df)
anova(HMD_abundance_model_2)

HMD_abundance_model_3 <- lmer(HMD_abundance ~ FractalDimension + 
                                (1 | Site) + 
                                (1 | Date),
                              data = merged_df)
anova(HMD_abundance_model_3)

HMD_abundance_model_4 <- lmer(HMD_abundance ~ Rugosity * HeightRange + 
                                (1 | Site) + 
                                (1 | Date),
                              data = merged_df)
anova(HMD_abundance_model_4)





############################################# complexity and foraging behavior #######################################################

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





#Foraging freq - edit these models

hist(behavior_merged$foraging_frequency)

install.packages("glmmTMB")
library(glmmTMB)

model_zi <- glmmTMB(
  foraging_frequency ~ HeightRange + RTS_Location + Diet_Category +
    (1 | Site) + (1 | Date),
  ziformula = ~1,
  family = nbinom2,
  data = behavior_merged
)

summary(model_zi)

#tukey for categories
emmeans(foraging_freq_model, pairwise ~ Diet_Category, adjust = "tukey")

#plot of it
ggplot(behavior_merged, aes(x=HeightRange, y=foraging_frequency)) + geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "",
    x = "Height Range",
    y = "Foraging Frequency"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.position = "none"
  )
















###### Foraging distance ######
foraging_merged <- clean_data %>%
  left_join(
    complexity_sum,
    by = c("Site", "Complexity_Level")
  )

foraging_dist_model <- lmer(RTS_Location ~ range_height + Diet_Category +
                              (1 | Site) + 
                              (1 | Date),
                            data = foraging_merged)
anova(foraging_dist_model)
