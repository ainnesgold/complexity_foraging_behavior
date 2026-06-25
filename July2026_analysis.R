library(tidyverse)
library(lmerTest)
library(emmeans)
library(ggpubr)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(gridExtra)
library(grid)
library(multcomp)

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



summary_stats_plot <- ggarrange(allreefplot_v2 + rremove("xlab"), fish_community_plot, nrow=2, ncol=1)


ggsave(
  filename = "outputs/summary_stats_plot.pdf",
  plot = summary_stats_plot,
  device = "pdf",
  width = 10.5,
  height = 10
)





################################# complexity + fish community ###############################################################

## merging two datasets
merged_df <- richness_abundance_sum %>%
  left_join(
    complexity_sum,
    by = c("Site", "Complexity_Level")
  )




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


## Model 3 visualization


HR_mean <- mean(merged_df$HeightRange, na.rm = TRUE)
HR_sd   <- sd(merged_df$HeightRange, na.rm = TRUE)

# New prediction grid
newdat <- expand.grid(
  Rugosity = seq(min(merged_df$Rugosity),
                 max(merged_df$Rugosity),
                 length.out = 100),
  HeightRange = c(HR_mean - HR_sd,
                  HR_mean,
                  HR_mean + HR_sd)
)



# Predict (response scale)
newdat$pred <- predict(
  abundance_model_3,
  newdata = newdat,
  type = "response",
  re.form = NA
)

# Label SD groups
newdat$HR_group <- factor(
  newdat$HeightRange,
  labels = c("-1 SD", "Mean", "+1 SD")
)


abundance_plot <- ggplot() +

  # raw data
  geom_point(
    data = merged_df,
    aes(x = Rugosity, y = abundance),
    alpha = 0.25
  ) +

  # model lines
  geom_line(
    data = newdat,
    aes(x = Rugosity,
        y = pred,
        color = HR_group,
        linetype = HR_group),
    linewidth = 1.2
  ) +
  scale_y_log10()+
  #scale_y_continuous(trans = "sqrt") +
  scale_color_manual(values = c("#0072B2", "black", "#D55E00")) +
  scale_linetype_manual(values = c("dashed", "solid", "dotted")) +

  labs(
    x = "Rugosity",
    y = "Abundance",
    color = "Height Range",
    linetype = "Height Range",
    title = "A"
  ) +

  theme_minimal(base_size = 16)
















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



### visualization for herbivore model 3


# Predict (response scale)
newdat$HMD_pred <- predict(
  HMD_abundance_model_3,
  newdata = newdat,
  type = "response",
  re.form = NA
)

# Label SD groups
newdat$HR_group <- factor(
  newdat$HeightRange,
  labels = c("-1 SD", "Mean", "+1 SD")
)


HMD_abundance_plot <- ggplot() +
  
  # raw data
  geom_point(
    data = merged_df,
    aes(x = Rugosity, y = HMD_abundance),
    alpha = 0.25
  ) +
  
  # model lines
  geom_line(
    data = newdat,
    aes(x = Rugosity,
        y = HMD_pred,
        color = HR_group,
        linetype = HR_group),
    linewidth = 1.2
  ) +
  scale_y_log10()+
  #scale_y_continuous(trans = "sqrt") +
  scale_color_manual(values = c("#0072B2", "black", "#D55E00")) +
  scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
  
  labs(
    x = "Rugosity",
    y = "Herbivore Abundance",
    color = "Height Range",
    linetype = "Height Range",
    title = "B"
  ) +
  
  theme_minimal(base_size = 16)



combined_abundance_plot <- ggarrange(abundance_plot, HMD_abundance_plot, nrow = 1, ncol = 2, common.legend = TRUE, legend = "top")

ggsave(
  filename = "outputs/combined_abundance_plot.pdf",
  plot = combined_abundance_plot,
  device = "pdf",
  width = 10,
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


behavior_merged$Diet_Category <- factor(behavior_merged$Diet_Category)


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


#interaction model
freq_model_3 <- glmmTMB(
  ff_beta ~ HeightRange * Rugosity + RTS_Location + Diet_Category,
  family = beta_family(),
  data = behavior_merged
)

summary(freq_model_3)
AIC(freq_model_3)




# model 3 results visualization
HR_mean <- mean(behavior_merged$HeightRange, na.rm = TRUE)
HR_sd   <- sd(behavior_merged$HeightRange, na.rm = TRUE)


newdat <- expand.grid(
  Rugosity = seq(min(behavior_merged$Rugosity, na.rm=TRUE),
                 max(behavior_merged$Rugosity, na.rm=TRUE),
                 length.out=100),
  
  HeightRange = c(HR_mean - HR_sd, HR_mean, HR_mean + HR_sd),
  
  RTS_Location = mean(behavior_merged$RTS_Location, na.rm=TRUE),
  
  Diet_Category = levels(behavior_merged$Diet_Category)
)

newdat$pred <- predict(freq_model_3, newdata=newdat, type="response", re.form=NA)

# average over diet
newdat <- newdat %>%
  group_by(Rugosity, HeightRange) %>%
  summarise(pred = mean(pred), .groups="drop")




# Label SD groups
newdat$HR_group <- factor(
  newdat$HeightRange,
  labels = c("-1 SD", "Mean", "+1 SD")
)


freq_plot <- ggplot() +
  
  # raw data
  geom_point(
    data = behavior_merged,
    aes(x = Rugosity, y = ff_beta),
    alpha = 0.25
  ) +
  
  # model lines
  geom_line(
    data = newdat,
    aes(x = Rugosity,
        y = pred,
        color = HR_group,
        linetype = HR_group),
    linewidth = 1.2
  ) +
  #scale_y_continuous(trans = "sqrt") +
  scale_color_manual(values = c("#0072B2", "black", "#D55E00")) +
  scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
  
  labs(
    x = "Rugosity",
    y = "Foraging Frequency",
    color = "Height Range",
    linetype = "Height Range",
    title = "A"
  ) +
  
  theme_minimal(base_size = 16)












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

foraging_merged$Diet_Category <- factor(foraging_merged$Diet_Category)


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





#visualization from model 3

HR_mean <- mean(foraging_merged$HeightRange, na.rm = TRUE)
HR_sd   <- sd(foraging_merged$HeightRange, na.rm = TRUE)


newdat <- expand.grid(
  Rugosity = seq(min(foraging_merged$Rugosity, na.rm=TRUE),
                 max(foraging_merged$Rugosity, na.rm=TRUE),
                 length.out=100),
  
  HeightRange = c(HR_mean - HR_sd, HR_mean, HR_mean + HR_sd),
  
  Diet_Category = levels(foraging_merged$Diet_Category)
)

newdat$pred <- predict(distance_model_3, newdata=newdat, type="response", re.form=NA)

# average over diet
newdat <- newdat %>%
  group_by(Rugosity, HeightRange) %>%
  summarise(pred = mean(pred), .groups="drop")




# Label SD groups
newdat$HR_group <- factor(
  newdat$HeightRange,
  labels = c("-1 SD", "Mean", "+1 SD")
)


dist_plot <- ggplot() +
  
  # raw data
  geom_point(
    data = foraging_merged,
    aes(x = Rugosity, y = RTS_Location),
    alpha = 0.25
  ) +
  
  # model lines
  geom_line(
    data = newdat,
    aes(x = Rugosity,
        y = pred,
        color = HR_group,
        linetype = HR_group),
    linewidth = 1.2
  ) +
  #scale_y_continuous(trans = "sqrt") +
  scale_color_manual(values = c("#0072B2", "black", "#D55E00")) +
  scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
  
  labs(
    x = "Rugosity",
    y = "Distance from Reef",
    color = "Height Range",
    linetype = "Height Range",
    title = "B"
  ) +
  theme_minimal(base_size = 16)







#post hoc diet category
emm_m3 <- emmeans(distance_model_3, ~ Diet_Category, type = "response")
pairs(emm_m3, adjust = "tukey", type = "response")
emm_df <- as.data.frame(emm_m3)


# Figure showing significant diet results


letters_df <- cld(
  emm_m3,
  adjust = "sidak",
  Letters = letters
)

dist_diet_plot <- ggplot() +
  geom_jitter(
    data = foraging_merged,
    aes(x = Diet_Category, y = RTS_Location),
    alpha = 0.15
  ) +
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
  #significance letters
  geom_text(
    data = letters_df,
    aes(x = Diet_Category,
        y = emmean,
        label = .group),
    vjust = -0.5, 
    size = 5
  ) +
  coord_flip() +
  labs(
    x = "Diet Category",
    y = "Distance from Reef (m)",
    title = "A"
  ) +
  theme_minimal(base_size = 16)
















################################# foraging distance ###############################################################

foraging_only <- foraging_merged %>%
  filter(Behavior == "foraging")

hist(foraging_only$RTS_Location)

foraging_only$Diet_Category <- factor(foraging_only$Diet_Category)


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




## model 3 visualization


HR_mean <- mean(foraging_only$HeightRange, na.rm = TRUE)
HR_sd   <- sd(foraging_only$HeightRange, na.rm = TRUE)

newdat <- expand.grid(
  Rugosity = seq(min(foraging_only$Rugosity, na.rm=TRUE),
                 max(foraging_only$Rugosity, na.rm=TRUE),
                 length.out=100),
  
  HeightRange = c(HR_mean - HR_sd, HR_mean, HR_mean + HR_sd),
  
  Diet_Category = levels(foraging_only$Diet_Category)
)

newdat$pred <- predict(forage_dist_model_3, newdata=newdat, type="response", re.form=NA)

# average over diet
newdat <- newdat %>%
  group_by(Rugosity, HeightRange) %>%
  summarise(pred = mean(pred), .groups="drop")



# Label SD groups
newdat$HR_group <- factor(
  newdat$HeightRange,
  labels = c("-1 SD", "Mean", "+1 SD")
)


foraging_dist_plot <- ggplot() +
  
  # raw data
  geom_point(
    data = foraging_only,
    aes(x = Rugosity, y = RTS_Location),
    alpha = 0.25
  ) +
  
  # model lines
  geom_line(
    data = newdat,
    aes(x = Rugosity,
        y = pred,
        color = HR_group,
        linetype = HR_group),
    linewidth = 1.2
  ) +
  #scale_y_continuous(trans = "sqrt") +
  scale_color_manual(values = c("#0072B2", "black", "#D55E00")) +
  scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
  
  labs(
    x = "Rugosity",
    y = "Foraging Distance",
    color = "Height Range",
    linetype = "Height Range",
    title = "C"
  ) +
  theme_minimal(base_size = 16) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "*Marginally significant interaction (p = 0.057)",
    hjust = -0.1,
    vjust = 1.5,
    size = 5
  )




#diet post hoc tests


#post hoc diet category
forage_emm_m3 <- emmeans(forage_dist_model_3, ~ Diet_Category, type = "response")
pairs(forage_emm_m3, adjust = "tukey", type = "response")
forage_emm_m3_df <- as.data.frame(forage_emm_m3)


# Figure showing significant diet results


forage_letters_df <- cld(
  forage_emm_m3,
  adjust = "sidak",
  Letters = letters
)

forage_dist_diet_plot <- ggplot() +
  geom_jitter(
    data = foraging_only,
    aes(x = Diet_Category, y = RTS_Location),
    alpha = 0.15
  ) +
  geom_point(
    data = forage_emm_m3_df,
    aes(x = Diet_Category, y = emmean),
    size = 2
  ) +
  geom_errorbar(
    data = forage_emm_m3_df,
    aes(x = Diet_Category,
        ymin = asymp.LCL,
        ymax = asymp.UCL),
    width = 0.2
  ) +
  #significance letters
  geom_text(
    data = forage_letters_df,
    aes(x = Diet_Category,
        y = emmean,
        label = .group),
    vjust = -0.5, 
    size = 5
  ) +
  coord_flip() +
  labs(
    x = "Diet Category",
    y = "Foraging Distance (m)",
    title = "B"
  ) +
  theme_minimal(base_size = 16)



### All Foraging complexity plots together



foraging_complexity_plots <- ggarrange(freq_plot, dist_plot, foraging_dist_plot, nrow = 1, ncol = 3, common.legend = TRUE, legend = "top")

ggsave(
  filename = "outputs/foraging_complexity_plots.pdf",
  plot = foraging_complexity_plots,
  device = "pdf",
  width = 15,
  height = 7
)



## Diet plots together

diet_plots <- ggarrange(dist_diet_plot, forage_dist_diet_plot, nrow = 1, ncol = 2)

ggsave(
  filename = "outputs/diet_plots.pdf",
  plot = diet_plots,
  device = "pdf",
  width = 14,
  height = 8
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





## visualization

pred_rts <- ggpredict(bites_model_3, terms = "RTS_Location")

bites_rts_plot <- ggplot(pred_rts, aes(x, predicted)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_jitter(data = bites_df,
              aes(x = RTS_Location, y = Bites),
              alpha = 0.15) +
  labs(x = "Distance From Reef (m)", y = "Number of Bites") +
  theme_minimal(base_size = 16)


ggsave(
  filename = "outputs/bites_rts_plot.pdf",
  plot = bites_rts_plot,
  device = "pdf",
  width = 5,
  height = 5
)



