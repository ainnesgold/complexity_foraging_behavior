
library(tidyverse)

complexity_raw <- read.csv("Complexity_Metrics.csv")

df_long <- complexity_raw %>%
  pivot_longer(
    cols = -Complexity_Level,
    names_to = c("metric", "stat"),
    names_pattern = "(.*)_(Mean|SD)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = stat, values_from = value)

df_long$Complexity_Level <- factor(
  df_long$Complexity_Level,
  levels = c("Low", "Medium", "High")
)

df_long$metric <- recode(df_long$metric,
                         "Fractal_Dimension" = "Fractal Dimension")

complexity_plot <- ggplot(df_long, aes(x = Complexity_Level, y = Mean, color = Complexity_Level)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.15) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = "Complexity Level",
    y = "Mean ± SD"
  ) +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none")


ggsave(
  filename = "outputs/complexity_plot.pdf",
  plot = complexity_plot,
  device = "pdf",
  width = 10,
  height = 5
)


##### June 9 2026, starting with raw data #####
#dont think i have the right raw data file, need to update

complexity_raw <- read.csv("Complexity_Raw.csv")

unique(complexity_raw$complexity)

complexity_sum <- complexity_raw %>%
  group_by(site, complexity) %>%
  summarise(height_mean = mean(height), height_range = diff(range(height)),
            rugosity_mean = mean(rugosity), rugosity_range = diff(range(rugosity)),
            fractal_mean_mean = mean(fractal_mean), fractal_mean_range = diff(range(fractal_mean)),
            fractal_raw_mean = mean(fractal_raw), fractal_raw_range = diff(range(fractal_raw)))







  