## Depends upon

### 1- Housekeeping ----
# Load Packages
library(tidyverse)
library(survival)
library(ggplot2)
library(janitor)
library(lubridate)
library(hrbrthemes)

# Run script 01 to create data, note will just read the created file if it has already been created
source(here::here("code/full_gtcs", "01_create_basefile.R"))

# Read in teachers data as denominator
teachers <- readRDS("Scratch_data/cleaned_teachers.Rds")

# Read in denominator data for whole pop
denom <- readRDS(here::here("Data", "HB2019_pop_est_1981_2019.rds")) %>% 
  filter(year == 2019) %>% 
  # Limit to working age (21 - 65)
  filter(age >= 21 & age <=65)
denom <- read_csv("Data/simd-19-tab3.csv")
denom <- denom %>% 
  select(-`90+`) %>% 
  gather(key = "age", "n", `0`:`89`) %>% 
  mutate(age =  as.integer(age),
         male = if_else(Sex == "Males", 1L, 0L),
         simd = case_when(
           `SIMD decile` %in% 1:2 ~ 1,
           `SIMD decile` %in% 3:4 ~ 2,
           `SIMD decile` %in% 5:6 ~ 3,
           `SIMD decile` %in% 7:8 ~ 4,
           `SIMD decile` %in% 9:10 ~ 5
         )) %>% 
  filter(age >= 21 & age <=65) %>% 
  group_by(simd, male, age) %>% 
  summarise(n = sum(n)) %>% 
  ungroup() 


# Read in HCW data to create denom
hcw <- readRDS("TOPLEVEL/HCW/hcw2021/processed_data/hcw_for_tch_analysis20210728.Rds")
mar20 <- readRDS("TOPLEVEL/Case_control/data/2021-07-28/hcw_mar20_linked_ANON_2021-07-28.rds") %>% 
  select(anon_id, age_years, sex, simd2020_sc_quintile)
nov20 <- readRDS("TOPLEVEL/Case_control/data/2021-07-28/hcw_nov20_linked_ANON_2021-07-28.rds") %>% 
  select(anon_id, age_years, sex, simd2020_sc_quintile)
apr21 <- readRDS("TOPLEVEL/Case_control/data/2021-07-28/hcw_apr21_linked_ANON_2021-07-28.rds") %>% 
  select(anon_id, age_years, sex, simd2020_sc_quintile)
hcw2 <- bind_rows(mar20,
                  nov20,
                  apr21) %>% 
  distinct(anon_id, .keep_all = TRUE)
rm(mar20, nov20, apr21)
hcw <- map(hcw, ~ inner_join(.x, hcw2))
rm(hcw2)

## Summarise all denominator data
hcw <- map(hcw, ~ .x %>% 
             rename(simd = simd2020_sc_quintile,
                    age_year = age_years) %>% 
             filter(age_year >=21 & age_year <= 65) %>% 
             mutate(age_grp = case_when(age_year >= 21 & age_year <= 30 ~ "21 - 30",
                                        age_year >= 31 & age_year <= 40 ~ "31 - 40",
                                        age_year >= 41 & age_year <= 50 ~ "41 - 50",
                                        age_year >= 51 & age_year <= 65 ~ "51 - 65",
                                        TRUE ~ "Unknown"),
                    male = if_else(sex == "1", 1L, 0L)) %>% 
             count(age_grp, male, simd))

# Pull out healthcare workers data for case control
hcw_denom <- hcw$pf 
hhd_denom <- hcw$pf_hhld 
npf_denom <- hcw$npf 
excld_denom <-  hcw$drop_hcw 
hcw_smry <- list(
  hcw = sum(hcw_denom$n),
  hhd = sum(hhd_denom$n),
  npf = sum(npf_denom$n),
  excld = sum(excld_denom$n))
saveRDS(hcw_smry, "output/number_hcws_hhds.Rds")

### Denominator Wrangling (Teachers) ----
teachers_denom <- teachers %>% 
  rename(age_year = age) %>% 
  # Limit to working age (21 - 65)
  filter(age_year >= 21 & age_year <=65) %>% 
  
  mutate(age_grp = case_when(age_year >= 21 & age_year <= 30 ~ "21 - 30",
                             age_year >= 31 & age_year <= 40 ~ "31 - 40",
                             age_year >= 41 & age_year <= 50 ~ "41 - 50",
                             age_year >= 51 & age_year <= 65 ~ "51 - 65",
                             TRUE ~ "Unknown")) %>% 
  mutate(male = if_else(sex == "1", 1L, 0L)) %>% 
  group_by(age_grp, male, simd) %>% 
  summarise(teacher_pop = n()) %>% 
  ungroup() 

### Denominator Wrangling (HCW) ----
hcw_denom <- hcw_denom %>% 
  rename(hcw_pop = n) 
npf_denom <- npf_denom %>% 
  rename(npf_pop = n)
excld_denom <- excld_denom %>% 
  rename(exclude_pop = n) 

### Denominator Wrangling (Whole Pop) ----
denom <- denom %>% 
  mutate(age_grp = case_when(age >= 21 & age <= 30 ~ "21 - 30",
                             age >= 31 & age <= 40 ~ "31 - 40",
                             age >= 41 & age <= 50 ~ "41 - 50",
                             age >= 51 & age <= 65 ~ "51 - 65",
                             TRUE ~ "Unknown")) %>% 
  group_by(age_grp,male, simd) %>% 
  summarise(whole_pop = sum(n)) %>% 
  ungroup() 

### Create whole pop file with teacher/hcw split ----
pop_teacher_hcw_split <- denom %>% 
  left_join(teachers_denom, by = c("age_grp", "male", "simd")) %>% 
  left_join(hcw_denom, by = c("age_grp", "male", "simd")) %>% 
  left_join(npf_denom, by = c("age_grp", "male", "simd")) %>% 
  left_join(excld_denom, by = c("age_grp", "male", "simd")) %>% 
  mutate_at(vars(teacher_pop, hcw_pop, npf_pop, exclude_pop), function(x) if_else(is.na(x), 0L, x)) %>% 
  mutate(pop = whole_pop - teacher_pop - hcw_pop -npf_pop - exclude_pop) %>% 
  select(-teacher_pop, -whole_pop, -hcw_pop, -npf_pop, -exclude_pop) %>% 
  mutate(hcw_teacher = "Neither") %>% 
  bind_rows(teachers_denom %>% rename(pop = teacher_pop) %>% 
              mutate(hcw_teacher = "Teacher")) %>%
  bind_rows(hcw_denom %>% rename(pop = hcw_pop) %>% 
              mutate(hcw_teacher = "Healthcare worker, patient facing")) %>% 
  bind_rows(npf_denom %>% rename(pop = npf_pop) %>% 
              mutate(hcw_teacher = "Healthcare worker, non-patient facing")) 
saveRDS(pop_teacher_hcw_split, "Scratch_data/denom_pop_tch_hcw.Rds")

### CC Data Wrangling ----
# Create a time variable of time since 24/02/2020
# Choosing 1 week before 1st Mar to avoid negative times
cc_data <- cc_data %>% 
  mutate(time = as.double(specimendate - as.Date("2020-02-24"))) %>% 
  mutate(male = if_else(sex == 1, 1, 0))

### Calculate Cumulative Incidence for any cases ----
# Run survfit to get proportion of cases for each combo of age/sex/simd/hcw/teacher. Note
# no censoring so just used for convenience
cc_data$hcw_teacher[cc_data$hcw_teacher == "Healthcare Worker, patient facing"] <- "Healthcare worker, patient facing"
cc_data <- cc_data %>% 
  mutate(hcw_teacher2 = str_replace(hcw_teacher, "\\,", "_"))
km1 <- survfit(Surv(time, is_case) ~ age_grp + male + hcw_teacher2, 
               data = cc_data %>% filter(is_case == TRUE))
## The surv_summary function has less clean output when a string with a comma is included so modify text
MakeCiPlotData <- function(km_choose){
  a <- survminer::surv_summary(km_choose) %>% 
    as_tibble() 
  a <- a %>% 
    mutate(hcw_teacher = str_replace(hcw_teacher2, "_", ","))
  a <- a %>% 
    select(-hcw_teacher2)
  # Get no of cases for each combo of age/sex/teacher
  sum_data <- a %>% 
    group_by(age_grp, male, hcw_teacher) %>% 
    summarise(n_cases = sum(n.event))
  
  # Get proportions
  prop_data <- a %>% 
    select(age_grp, male, hcw_teacher, time, n.risk, n.event, surv, upper, lower) 
  
  # Get cumulative numbers of cases over time by multiplying cases by the proportion of case
  cum_data <- prop_data %>% 
    left_join(sum_data, by = c("age_grp", "male", "hcw_teacher")) %>% 
    mutate(cum_cases = (1-surv)*n_cases, 
           upper_cases = (1-lower)*n_cases,
           lower_cases = (1-upper)*n_cases) 
  cum_data <- cum_data %>% 
    mutate_at(vars(age_grp, hcw_teacher), as.character) %>% 
    mutate(male = if_else(male =="1", 1L, 0L))
  
  cum_data <- cum_data %>% 
    
    # match on population data (split by teachers)
    left_join(pop_teacher_hcw_split %>% 
                group_by(age_grp, male, hcw_teacher) %>% 
                summarise(pop = sum(pop)), 
              by = c("age_grp", "male", "hcw_teacher")) %>% 
    
    mutate(cum_inc = cum_cases/pop,
           cum_inc_upper = upper_cases/pop,
           cum_inc_lower = lower_cases/pop) %>% 
    
    mutate(date = time + as.Date("2020-02-24"))
  
  
  cum_data_strt <- cum_data %>% 
    group_by(hcw_teacher, age_grp, male) %>% 
    slice(1) %>% 
    mutate(date = min(cum_data$date, na.rm = TRUE),
           cum_inc = 0, 
           cum_inc_lower =0,
           cum_inc_upper = 0) %>% 
    ungroup()
  
  cum_data_end <- cum_data %>% 
    arrange(desc(cum_inc)) %>% 
    group_by(hcw_teacher, age_grp, male) %>% 
    slice(1) %>% 
    mutate(date = as.Date("2021-06-30")) %>% 
    ungroup()
  
  bind_rows(cum_data_strt, cum_data, cum_data_end)
  
}
cum_data2 <- MakeCiPlotData(km1)

male_lookup <- c(
  `0` = "Women", 
  `1` = "Men")

plot1 <- ggplot() +
  geom_step(data = cum_data2, mapping=aes(x=date, y=100*cum_inc, group = hcw_teacher,
                                          colour = hcw_teacher)) +
  facet_grid(age_grp ~ male, labeller = labeller(.cols = male_lookup)) +
  xlab("Date") +
  ylab("% Risk - Any Case") +
  scale_x_date("Week beginning", date_breaks = "months", date_labels = "%b-%y") +
  theme_ipsum() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_ipsum()
plot1

# save plot
ggsave(plot = plot1, filename = here::here("output", "cum_inc_any_test.png"))
saveRDS(plot1, here::here("output", "cum_inc_case.RDS"))

### Calculate Cumulative Incidence for hospitalisation ----
km2 <- survfit(Surv(time, hosp) ~ age_grp + male + hcw_teacher2, 
               data = cc_data %>% filter(hosp == 1))
cum_data2 <- MakeCiPlotData(km2)
plot2 <- ggplot() +
  geom_step(data = cum_data2 %>% 
              mutate(age_grp = paste(age_grp, "yrs"),
                     hcw_teacher = str_to_sentence(hcw_teacher)), mapping=aes(x=date, y=100*cum_inc, group = hcw_teacher,
                                                                              colour = hcw_teacher)) +
  facet_grid(age_grp ~ male, labeller = labeller(.cols = male_lookup)) +
  geom_vline(xintercept = as.Date("2020-08-24"), linetype = 2, colour = "gray45") +
  geom_vline(xintercept = as.Date("2020-12-23"), linetype = 2, colour = "gray45") +
  geom_vline(xintercept = as.Date("2021-02-26"), linetype = 2, colour = "gray45") +
  geom_vline(xintercept = as.Date("2021-04-23"), linetype = 2, colour = "gray45") +
  xlab("Date") +
  ylab("Risk - Hospitalisations (%)") +
  scale_x_date("Week beginning", date_breaks = "2 months", date_labels = "%b-%y") +
  theme_ipsum() +
  scale_colour_manual("",
                      breaks = c("Healthcare worker, patient facing",
                                 "Healthcare worker, non-patient facing",
                                 "Teacher",
                                 "Neither"),
                      values = c("#422d54",
                                 "#758bd1",
                                 "#2d543d",
                                 "#d175b8")) 
plot2

# save plot
ggsave(plot = plot2, filename = here::here("output", "cum_inc_hosp.png"))
saveRDS(plot2, here::here("output", "cum_inc_hosp.RDS"))

tiff("output/plot_ci.tiff", res = 300, width = 30, height = 16, unit = "cm", compression = "lzw")
plot2 + theme_ipsum(axis_text_size = 8)
dev.off()

