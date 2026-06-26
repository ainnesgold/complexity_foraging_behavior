
library(tidyverse)


RDH_allreef <- read.csv("RDH_AllReef_V2.csv")


RDH_long <- RDH_allreef %>%
  pivot_longer(
    cols = c(Rugosity, HeightRange),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      "Rugosity" = "Rugosity",
      "HeightRange" = "Height Range"
    ))

RDH_allreef %>%
  filter(Site == 1) %>%
  summarise(min(Rugosity, na.rm = TRUE))

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
  


