library(tidyverse)
library(readxl)

polar_coords <- function() { 
  data.frame( 
    direction = c("N","E","S","W","SW","NE","NW","SE"), 
    degrees = c(90, 0, 270, 180, 225, 45, 135, 315)
  ) %>% 
    mutate(radians = 2*pi*degrees/360) 
}



fiber_data <- read_csv(here::here('data/Trial 1 All Fiber Data.csv'), col_types = cols(.default = 'c')) %>% 
  janitor::clean_names() %>% 
  rename(length = length_um, width = width_um) %>%
  mutate_at(vars(rjlg_sample_id, fiber_no, length, width, aspect_ratio, distance_ft), as.numeric) %>% 
  mutate(units = 'um', 
         field = as.numeric(str_replace(field, 'Field ', '')), 
         fiber_type = ifelse(fiber_type == 'A', 'Crocidolite', 'Chrysotile'), 
         fiber_density = ifelse(fiber_type == 'Crocidolite', 3.2, 2.4), 
         aed_henn = sqrt(fiber_density*log(2*aspect_ratio))*width, # Henn (1996) 
         aed = 66*width*(aspect_ratio/(2 + 4*aspect_ratio))^2.2, # Timbrell (1965)
         aed_hinds = width * fiber_density^0.5
  ) %>%
  filter(!is.na(length)) %>% 
  left_join(polar_coords(), by = 'direction') 
write_rds(fiber_data, 'data/processed/fiber-data.rds')

debris <- excel_sheets('data/TOTAL Fiber Dispersion Compilation_Debris.xlsx') %>% 
  map_dfr(~ read_xlsx('data/TOTAL Fiber Dispersion Compilation_Debris.xlsx', sheet = .x, col_types = 'text') %>% 
            mutate(id = .x)) %>% 
  janitor::clean_names() %>% 
  select(id, fiber_no = fiber, fiber_type, length = length_um, width = width_um, debris) %>% 
  separate(id, c("rjlg_sample_id", 'chem_risk_id'), sep = ' ') %>% 
  mutate_at(vars(rjlg_sample_id, fiber_no, length, width), as.numeric) %>% 
  mutate(debris = debris == 'yes', 
         aspect_ratio = length/width,
         fiber_type = ifelse(fiber_type == 'A', 'Crocidolite', 'Chrysotile')) %>% 
  filter(!is.na(debris)) |>
  # group_by(chem_risk_id) |> 
  # fill(field, .direction = 'down') |> 
  left_join(fiber_data |> select(rjlg_sample_id, chem_risk_id, distance_ft, direction, orientation) |> distinct(), 
            by = c('rjlg_sample_id', 'chem_risk_id')) |> 
  left_join(polar_coords(), by = 'direction') 
write_rds(debris, 'data/processed/debris-data.rds')

concentrations <- read_xlsx('data/concentration-data.xlsx') %>% 
  janitor::clean_names() %>% 
  filter(str_detect(chem_risk_id, 'CH_1')) %>% 
  rename(total_fiber_cc = fiber_cc, total_fibers = number_fibers) %>%
  left_join(fiber_data %>% select(chem_risk_id, distance_ft, direction, orientation, degrees, radians) %>% distinct(), by = 'chem_risk_id')
write_rds(concentrations, 'data/processed/concentrations.rds')
