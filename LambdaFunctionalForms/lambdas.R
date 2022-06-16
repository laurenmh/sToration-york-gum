library(tidyverse)

read_csv("water_full_env.csv", col_types = cols()) %>%
  select(-`...1`) %>%
  rename_with(function(colname) str_remove(colname, ".x")) %>%
  rename(Phosphorus = Colwell.P) %>%
  mutate(Invasive = if_else(Focal.sp %in% c("A", "H"), "exotic", "native")) %>%
  mutate(Focal.sp = recode(Focal.sp,
                           "A" = "Arctotheca calendula", "T" = "Trachymene cyanopetala",
                           "W" = "Waitzia acuminata", "H" = "Hypochaeris glabra")) %>%
  mutate(Focal.sp = str_c(Focal.sp, if_else(Focal.sp %in% c("Arctotheca calendula",
                                                            "Hypochaeris glabra"),
                                            " (exotic)", " (native)"))) %>%
  filter(Trt.comp == "S") %>%
  group_by(Reserve, Plot.ID, Quadrant, Focal.sp, Canopy, Invasive) %>%
  summarise(`Mean flowers` = mean(Number.flowers.total), .groups = "drop") %>%
  drop_na() %>%
  mutate(Canopy = scale(Canopy)) %>%
  ggplot(aes(x = Canopy, y = `Mean flowers`, colour = Reserve)) +
  geom_point(position = "jitter") +
  facet_wrap(~ Focal.sp, scale = "free") +
  theme_bw()
