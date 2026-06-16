
library(tidyverse)

##### first version of data: all reef
RDH_allreef <- read.csv("RDH_AllReef.csv")

RDH_allreef_sum <- RDH_allreef %>%
  group_by(Site, Complexity_Level) %>%
  summarise(Rugosity = mean(Rugosity),
            FractalDimension = mean(FractalDimension),
            HeightRange = mean(HeightRange))

RDH_long <- RDH_allreef_sum %>%
  mutate(
    Complexity_Level = factor(
      Complexity_Level,
      levels = c("Low", "Medium", "High")
    )
  ) %>%
  pivot_longer(
    cols = c(Rugosity, FractalDimension, HeightRange),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      "Rugosity" = "Rugosity",
      "FractalDimension" = "Fractal Dimension",
      "HeightRange" = "Height Range"
    )
  )



allreefplot <- ggplot(
  RDH_long,
  aes(
    x = Complexity_Level,
    y = Value,
    fill = factor(Site)
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(
    values = c("#0072B2", "#D55E00"),  # Okabe-Ito palette
    name = "Site"
  ) +
  labs(
    x = "Complexity Level",
    y = "Value"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )



# ggsave(
#   filename = "outputs/allreefplot.pdf",
#   plot = allreefplot,
#   device = "pdf",
#   width = 10,
#   height = 5
# )




## trying out a different plot
RDH_allreef <- read.csv("RDH_AllReef_V2.csv")

# RDH_allreef_sum <- RDH_allreef %>%
#   group_by(Site) %>%
#   summarise(Rugosity = mean(Rugosity),
#             FractalDimension = mean(FractalDimension),
#             HeightRange = mean(HeightRange))

RDH_long <- RDH_allreef %>%
  pivot_longer(
    cols = c(Rugosity, FractalDimension, HeightRange),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      "Rugosity" = "Rugosity",
      "FractalDimension" = "Fractal Dimension",
      "HeightRange" = "Height Range"
    ))



allreefplot_v2 <- ggplot(RDH_long, aes(x = as.factor(Site), y = Value, fill = as.factor(Site))) +
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
  

# ggsave(
#   filename = "outputs/allreefplot_v2.pdf",
#   plot = allreefplot_v2,
#   device = "pdf",
#   width = 10,
#   height = 5
# )
















 #### second version of data: all plot


RDH_allplot <- read.csv("RDH_AllPlot.csv")

RDH_allplot_sum <- RDH_allplot %>%
  group_by(Site, Complexity_Level) %>%
  summarise(Rugosity = mean(Rugosity),
            FractalDimension = mean(FractalDimension),
            HeightRange = mean(HeightRange))

RDH_long_allplot <- RDH_allplot_sum %>%
  mutate(
    Complexity_Level = factor(
      Complexity_Level,
      levels = c("Low", "Medium", "High")
    )
  ) %>%
  pivot_longer(
    cols = c(Rugosity, FractalDimension, HeightRange),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      "Rugosity" = "Rugosity",
      "FractalDimension" = "Fractal Dimension",
      "HeightRange" = "Height Range"
    )
  )

allplotplot <- ggplot(
  RDH_long_allplot,
  aes(
    x = Complexity_Level,
    y = Value,
    fill = factor(Site)
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(
    values = c("#0072B2", "#D55E00"),  # Okabe-Ito palette
    name = "Site"
  ) +
  labs(
    x = "Complexity Level",
    y = "Value"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    panel.grid.major.x = element_blank()
  )



ggsave(
  filename = "outputs/allplotplot.pdf",
  plot = allplotplot,
  device = "pdf",
  width = 10,
  height = 5
)















# #### OLD
# ### REDO WITH FINAL COMPLEXITY DATA ####
# complexity_raw <- read.csv("Complexity_Metrics.csv")
# 
# df_long <- complexity_raw %>%
#   pivot_longer(
#     cols = -Complexity_Level,
#     names_to = c("metric", "stat"),
#     names_pattern = "(.*)_(Mean|SD)",
#     values_to = "value"
#   ) %>%
#   pivot_wider(names_from = stat, values_from = value)
# 
# df_long$Complexity_Level <- factor(
#   df_long$Complexity_Level,
#   levels = c("Low", "Medium", "High")
# )
# 
# df_long$metric <- recode(df_long$metric,
#                          "Fractal_Dimension" = "Fractal Dimension")
# 
# complexity_plot <- ggplot(df_long, aes(x = Complexity_Level, y = Mean, color = Complexity_Level)) +
#   geom_point(size = 3) +
#   geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15) +
#   facet_wrap(~ metric, scales = "free_y") +
#   scale_color_brewer(palette = "Set2") +
#   labs(
#     x = "Complexity Level",
#     y = "Mean ± SD"
#   ) +
#   theme_minimal(base_size = 16) +
#   theme(legend.position = "none")
# 
# 
# ggsave(
#   filename = "outputs/complexity_plot.pdf",
#   plot = complexity_plot,
#   device = "pdf",
#   width = 10,
#   height = 5
# )
# 
# 
# ##### June 9 2026, starting with raw data #####
# #dont think i have the right raw data file, need to update
# 
# complexity_raw <- read.csv("Complexity_Raw.csv")
# 
# unique(complexity_raw$complexity)
# 
# complexity_sum <- complexity_raw %>%
#   group_by(site, complexity) %>%
#   summarise(height_mean = mean(height), height_range = diff(range(height)),
#             rugosity_mean = mean(rugosity), rugosity_range = diff(range(rugosity)),
#             fractal_mean_mean = mean(fractal_mean), fractal_mean_range = diff(range(fractal_mean)),
#             fractal_raw_mean = mean(fractal_raw), fractal_raw_range = diff(range(fractal_raw)))
# 
# 
# 
# 
# 
# 
# 
#   