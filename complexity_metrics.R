

complexity_raw <- read.csv("Complexity_Metrics.csv")


library(tidyr)
library(ggplot2)
library(tidyr)
library(dplyr)



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

ggplot(df_long, aes(x = Complexity_Level, y = Mean, color = Complexity_Level)) +
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
