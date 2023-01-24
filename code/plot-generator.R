library(tidyverse)
library(patchwork)

#### Data files ----------------------------------------------------------------
fiber_data <- read_rds(here::here('data/processed/fiber-data.rds'))
concentrations <- read_rds(here::here('data/processed/concentrations.rds'))

gm <- function(x) { exp(mean(log(x))) }
gsd <- function(x) { exp(sd(log(x[x>0]))) }
sum_na <- function(v) { 
  return(ifelse(all(is.na(v)), NA, sum(v, na.rm = T)))
}
fiber_smry <- fiber_data %>% 
  group_by(chem_risk_id, distance_ft, fiber_type, direction, degrees, orientation) %>% 
  summarize(n = n(),
            length = gm(length), 
            width = gm(width), 
            aspect_ratio = gm(aspect_ratio[!is.infinite(aspect_ratio)]), 
            debris = sum_na(debris)/n
            # debris = case_when(debris != 0 ~ debris), 
  ) %>% 
  arrange(distance_ft, degrees)
fiber_conc <- fiber_smry %>% 
  ungroup() %>% 
  left_join(concentrations, by = c('chem_risk_id', 'distance_ft', 'direction', 'degrees', 'orientation')) %>% 
  mutate(fiber_cc = n/total_fibers*total_fiber_cc) %>% 
  ungroup() %>% 
  group_by(distance_ft, fiber_type, direction, degrees) %>% 
  summarize(length = mean(length), 
            width = mean(width), 
            aspect_ratio = mean(aspect_ratio),  
            fiber_cc = mean(fiber_cc),  
            aed = (fiber_cc*log(2*aspect_ratio))^1/2 * width,
            total_fiber_cc = mean(total_fiber_cc)) 


plume_window <- c("E", "NE", "SE")
fiber_smry2 <- fiber_smry %>% filter(direction %in% plume_window)
fiber_conc2 <- fiber_conc %>% 
  filter(direction %in% plume_window) %>% 
  mutate(xx = fiber_cc*case_when(distance_ft == 8 ~ 1)) %>% 
  group_by(fiber_type) %>% 
  fill(xx, .direction = 'down') %>% 
  mutate(norm_conc = fiber_cc/xx, 
         norm_dist = distance_ft - 8) %>% 
  select(-xx)

debris_smry <- fiber_data %>% 
  filter(!is.na(debris))
debris_conc <- fiber_data %>% 
  # filter(!is.na(debris)) %>% 
  group_by(chem_risk_id, distance_ft, fiber_type, direction, degrees, orientation, debris) %>% 
  summarize(n = n(),
            length = gm(length), 
            width = gm(width), 
            aspect_ratio = gm(aspect_ratio)) %>% 
  arrange(distance_ft, degrees) %>% 
  left_join(concentrations, by = c('chem_risk_id', 'distance_ft', 'direction', 'degrees', 'orientation')) %>% 
  mutate(fiber_cc = n/total_fibers*total_fiber_cc) 
debris_conc2 <- debris_conc %>% 
  filter(!is.na(debris), direction %in% plume_window) %>% 
  mutate(xx = fiber_cc*case_when(distance_ft == 8 ~ 1)) %>% 
  ungroup() %>% 
  group_by(fiber_type, debris) %>% 
  fill(xx, .direction = 'down') %>% 
  mutate(norm_conc = fiber_cc/xx, 
         norm_dist = distance_ft - 8) %>% 
  select(-xx)

#### Plots ---------------------------------------------------------------------
## Sampling map
fiber_data %>% 
  distinct(distance_ft, degrees, direction) %>% 
  mutate(loc = paste0(distance_ft, ' ft ', direction), 
         source = F) %>%
  add_row(distance_ft = 0, degrees = 0, direction = NA, loc = NA, source = T) %>% 
  ggplot(aes(x = distance_ft, y = degrees, label = loc, shape = source)) + 
  coord_polar(theta = 'y', start = -1.57, direction = -1) + 
  geom_point(size = 4) +
  # ggrepel::geom_text_repel(size = 3) +
  geom_text(hjust = 0.5, vjust = -1, size = 3) + 
  labs(x = 'Distance (ft)') + 
  scale_shape_manual(values = c(10, 18)) +
  scale_size_continuous(range = c(1, 12)) + 
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_y_continuous(limits = c(0, 360), breaks = seq(0, 360, by = 45)) + 
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size = 14), 
        axis.title.x = element_blank())
ggsave('output/sampling-locations.jpg', height = 8, width = 8, units = 'in')

## AED and aspect ratio cumulative probability by fiber type
p1 <- ggplot(fiber_data %>% filter(width>0), aes(x = aed, linetype = fiber_type)) + 
  stat_ecdf(size = 0.8) + 
  scale_x_log10() +
  annotation_logticks(sides = 'b') + 
  labs(x = 'AED (um), log scale', y = 'Cumulative probability') + 
  scale_linetype_discrete(name = 'Fiber Type') + 
  theme_bw() + 
  theme(text = element_text(size = 12.5))
p2 <- ggplot(fiber_data %>% filter(width>0), aes(x = aspect_ratio, linetype = fiber_type)) + 
  stat_ecdf(size = 0.8) + 
  scale_x_log10() +
  annotation_logticks(sides = 'b') + 
  labs(x = 'Aspect ratio, log scale', y = 'Cumulative probability') + 
  scale_linetype_discrete(name = 'Fiber Type') + 
  theme_bw() + 
  theme(text = element_text(size = 12.5))
p2 + p1 + plot_layout(guides = 'collect')
ggsave('output/aspect-ratio-and-aed-cum-prob.jpg', height = 4, width = 8, units = 'in')

## Fiber deposition (polar coordinates)
fiber_conc %>% ungroup() %>% add_row(distance_ft = 0, degrees = 0, fiber_cc = 0, fiber_type ='Chrysotile') %>% 
  ggplot(aes(x = distance_ft, y = degrees, size = fiber_cc, shape = fiber_type)) + 
  coord_polar(theta = 'y', start = -1.57, direction = -1) + 
  geom_point(alpha = 0.8) + 
  facet_wrap(~ fiber_type) + 
  labs(x = 'Distance (ft)') +
  scale_shape_manual(values = c(19, 1)) +
  scale_size_continuous(range = c(1, 12)) + 
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) + 
  scale_y_continuous(limits = c(0, 360), breaks = seq(0, 360, by = 45)) + 
  theme_bw() + 
  theme(legend.position = "none",
        text = element_text(size = 12), 
        axis.title.x = element_blank())
ggsave('output/deposition-polar.jpg', height = 4, width = 8, units = 'in')

## Fiber deposition (cartesian coordinates)
ggplot(fiber_conc, aes(x = distance_ft, y = fiber_cc, shape = fiber_type, linetype = fiber_type)) + 
  facet_wrap(~ reorder(direction, degrees)) + 
  geom_point(size = 2.5) + 
  geom_line(size = 0.6) + 
  scale_shape_manual(values = c(19, 1)) + 
  labs(x = 'Distance (ft)', y = 'Concentrations (f/cc)') + 
  theme_bw() +    
  theme(legend.title = element_blank(), 
        text = element_text(size = 12.5))
ggsave('output/deposition-cartesian.jpg', height = 5, width = 8, units = 'in')


## Fiber deposition rate
ggplot(fiber_conc2, aes(x = norm_dist, y = norm_conc, shape = fiber_type, linetype = fiber_type)) + 
  geom_point(size = 2.5) + 
  # facet_wrap(~ fiber_type) +
  stat_smooth(method = 'lm', formula = y ~ x - 1, se = F, color = 'black', size = 0.8) + 
  labs(x = 'Distance from the 8 ft collector (ft)', y = 'Normalized concentration, log scale') +
  scale_x_continuous(limits = c(0, 100)) + 
  scale_y_log10() +
  scale_shape_manual(name = 'Fiber type', values = c(19, 1)) + 
  scale_linetype_discrete(name = 'Fiber type') + 
  annotation_logticks(sides = 'l') +
  theme_bw() +
  theme(text = element_text(size = 12.5))
ggsave('output/deposition-rate-fit.jpg', height = 4, width = 8, units = 'in')


## Effect of fiber dimension on debris binding
ggplot(debris_smry, aes(x = aed, linetype = debris)) + 
  facet_wrap(~fiber_type) + 
  stat_ecdf(size = 0.8) + 
  scale_x_log10() +
  annotation_logticks(sides = 'b') + 
  labs(x = 'AED (um)', y = 'Cumulative probability') + 
  scale_linetype_discrete(name = 'Fiber bound\nto debris') +
  theme_bw() +
  theme(text = element_text(size = 12.5))
ggsave('output/debris-binding-width-cum-prob.jpg', height = 4, width = 8, units = 'in')

ggplot(debris_smry, aes(x = length, linetype = debris)) +
  facet_wrap(~fiber_type) + 
  stat_ecdf(size = 0.8) +
  scale_x_log10() +
  annotation_logticks(sides = 'b') +
  labs(x = 'Length (um)', y = 'Cumulative probability') +
  scale_linetype_discrete(name = 'Fiber bound\nto debris') +
  theme_bw() +
  theme(text = element_text(size = 12.5))
ggsave('output/debris-binding-length-cum-prob.jpg', height = 4, width = 8, units = 'in')

debris_conc2 %>% 
  group_by(fiber_type, debris) %>% 
  summarize(wt_dist = sum(fiber_cc*distance_ft)/sum(fiber_cc)) %>% 
  ggplot(aes(x = fiber_type, y = wt_dist, fill = debris)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  geom_text(aes(label = round(wt_dist, 1)), vjust = -0.4, position = position_dodge(0.9)) +
  labs(x = 'Fiber type', y = 'Distance (ft)') + 
  scale_fill_grey(name = 'Fiber bound\nto debris') + 
  theme_bw() +
  theme(text = element_text(size = 12.5))
ggsave('output/dist-traveled-by-fiber-type-and-debris.jpg', height = 4, width = 8, units = 'in')
