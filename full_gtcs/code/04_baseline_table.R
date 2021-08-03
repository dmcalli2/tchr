## Generates tables by re-weighting the control arm to the scottish population
## Note this is quite slow

# Depends
# Runs code/full_gtcs/01_create_basefile.R for case control study
# Runs also code/full_gtcs/02_clean_ecoss.R for tsting results 
# Read in denominator data for age, sex, simd and healthboard for Scotland "Data/simd-19-tab3.csv"

## Can use practice number or health board and SIMD
## Gives very similar results, however (using 5-year bands) about 5% of the population is 
## Missed off based on practice number. Can get the practice number analysis from version control history. Dropped as more complex and gives same answer

library(tableone)
library(survey)

# read case control data ----
source("code/full_gtcs/01_create_basefile.R")

# read ecoss data ----
source("code/full_gtcs/02_process_ecoss.R")

## limit to controls and to cases ----
cases <- cc_data %>% 
  filter(is_case==1)
cc_data <- cc_data %>% 
  filter(is_case ==0) 

## Limit to first occurrence of control in data ----
cc_data <- cc_data %>% 
  arrange(specimendate) %>% 
  distinct(anon_id, .keep_all = TRUE)

# Calculate inclusion probabilities by age, sex, simd and hb (hb/simd combo as a proxy for gp practice), denom all of Scotland ----
denom <- read_csv("Data/simd-19-tab3.csv")
denom <- denom %>% 
  gather("age", "n", -`NHS health board`, -`SIMD decile`, -Sex)
denom <- denom %>% 
  rename(simd = `SIMD decile`) %>% 
  mutate(age = if_else(age == "90+", "92", age),
         age = as.integer(age),
         hb2019name = paste0("NHS ", `NHS health board`),
         sex = if_else(Sex == "Males", 1L, 2L),
         simd = case_when(
           simd %in% 1:2 ~ 1,
           simd %in% 3:4 ~ 2,
           simd %in% 5:6 ~ 3,
           simd %in% 7:8 ~ 4,
           simd %in% 9:10 ~ 5)
  ) %>% 
  select(age_year = age, sex, simd, hb2019name, n)
denom <- denom %>% 
  group_by(age_year, sex, simd, hb2019name) %>%  
  summarise(n = sum(n)) %>% 
  ungroup()
numer <- cc_data %>% 
  group_by(age_year, sex, simd, hb2019name) %>% 
  count() %>% 
  ungroup()

# note one individual in a small board and young who is not in denominator. Therefore set probability to one for this individual
hb_simd <- denom %>%
  mutate(sex = as.character(sex)) %>% 
  inner_join(numer %>% rename(x = n)) %>% 
  mutate(prbs = x/n,
         prbs = if_else(prbs >1, 1, prbs))  %>% 
  select(-x, -n)

## impute probability for missing SIMD, 0.13%
cc_data <- cc_data %>% 
  left_join(hb_simd %>% rename(hb_simd_prb = prbs))
cc_data <- cc_data %>% 
  group_by(age_year, sex, simd, hb2019name) %>% 
  mutate(hb_simd_prb = if_else(is.na(hb_simd_prb), mean(hb_simd_prb, na.rm = TRUE), 
                               hb_simd_prb)) %>% 
  ungroup()

## impute probability for missing hb and SIMD
cc_data <- cc_data %>% 
  group_by(age_year, sex) %>% 
  mutate(hb_simd_prb = if_else(is.na(hb_simd_prb), mean(hb_simd_prb, na.rm = TRUE), 
                               hb_simd_prb)) %>% 
  ungroup()

# rename variables for tabulation
cc_data <- cc_data %>%
  mutate(female = if_else(sex == 2, 1, 0),
         ethnicity = case_when(white == 1 ~ "White",
                               s_asian == 1 ~ "Asian",
                               chinese == 1 ~ "Chinese",
                               black == 1 ~ "Black",
                               other == 1 ~ "Other",
                               unknown == 1 ~ "Unknown"),
         como_count = case_when(como_count == 0 ~ "0",
                                como_count == 1 ~ "1",
                                como_count >1 ~ "2+"),
         como_any = case_when(ihd_any ==1 | heart_other_any == 1 | circulatory_other == 1 |
                                ckd_any == 1 | oad_any == 1 | neuro_any == 1 | liver_any == 1
                              | immune_any == 1 | neoplasm_any == 1 
                              | neoplasm_last_year ==1 |esoph_stomach_duod == 1 ~ 1, TRUE ~ 0),
         sex.male = case_when(sex == "Male" ~ 1, TRUE ~ 0),
         sex.female = case_when(sex == "Female" ~ 1, TRUE ~ 0))
## vaccination data
cc_data <- cc_data %>% 
  mutate(dose_ever = case_when(
    !is.na(vt2_dt) ~ "c_post2nd",
    !is.na(vt1_dt) ~ "b_post1st",
    is.na(vt1_dt) ~ "a_unvacc"))

# remove duplicates caused because anonymous ID can appear twice
cc_data <- cc_data %>% 
  distinct(anon_id, .keep_all = TRUE)

## Estimate vaccination CI ----
cc_data_nst <- cc_data %>% 
  filter(!is_case ==1) %>% 
  select(anon_id, vt1_dt, vt2_dt, hcw_teacher, age_grp, sex, hb_simd_prb) %>% 
  mutate_at(vars(vt1_dt, vt2_dt), as.integer) %>% 
  mutate_at(vars(vt1_dt, vt2_dt), ~ if_else(is.na(.x), as.integer(as.Date("2999-01-01")), .x )) 

spec_list <- unique(c(cc_data$vt1_dt, cc_data$vt2_dt)) %>% sort()
spec_list_int <- as.integer(spec_list)

## By sex
cc_data_nst_grp <- cc_data_nst %>% 
  select(hcw_teacher, vt1_dt, vt2_dt, hb_simd_prb) %>% 
  group_by(hcw_teacher) %>% 
  nest() %>% 
  ungroup()
cc_data_nst_grp$ci <- map(cc_data_nst_grp$data, function(mydf){
  tibble(vax_dt = spec_list,
         v1_ci  = map_dbl(spec_list_int, ~ weighted.mean(mydf$vt1_dt <= .x, 1/mydf$hb_simd_prb)),
         v2_ci  = map_dbl(spec_list_int, ~ weighted.mean(mydf$vt2_dt <= .x, 1/mydf$hb_simd_prb)))
})
cc_data_nst_grp <- cc_data_nst_grp %>% 
  select(-data) %>% 
  unnest(ci)

cc_data_nst_grp2 <- cc_data_nst_grp %>% 
   gather("dose", "ci", v1_ci, v2_ci)
plot1 <- ggplot(cc_data_nst_grp2,
                aes(x = vax_dt, y = ci, colour = hcw_teacher)) +  
  geom_step() +
  facet_grid(~ dose)
plot1 +
  facet_grid(~dose) +
  scale_y_continuous("Proportion vaccinated") +
  scale_x_date("Date of dose") +
  scale_color_discrete("") +
  scale_linetype("")

## By age and sex
cc_data_nst_grp <- cc_data_nst %>% 
  select(hcw_teacher, age_grp, sex, vt1_dt, vt2_dt, hb_simd_prb) %>% 
  group_by(hcw_teacher, age_grp, sex) %>% 
  nest() %>% 
  ungroup()
cc_data_nst_grp$ci <- map(cc_data_nst_grp$data, function(mydf){
  tibble(vax_dt = spec_list,
         v1_ci  = map_dbl(spec_list_int, ~ weighted.mean(mydf$vt1_dt <= .x,1/ mydf$hb_simd_prb)),
         v2_ci  = map_dbl(spec_list_int, ~ weighted.mean(mydf$vt2_dt <= .x, 1/mydf$hb_simd_prb)))
})
cc_data_nst_grp <- cc_data_nst_grp %>% 
  select(-data) %>% 
  unnest(ci)

cc_data_nst_grp <- cc_data_nst_grp %>% 
  gather("dose", "ci", v1_ci, v2_ci)
plot2 <- ggplot(cc_data_nst_grp %>% 
                  mutate(age_grp = paste(age_grp, "yrs"),
                         sex = factor(sex, levels = c("2", "1"), labels = c("Women", "Men")),
                         dose = factor(dose, levels = c("v1_ci", "v2_ci"), labels = c("1st dose", "2nd dose")),
                         hcw_teacher = str_to_sentence(hcw_teacher)
                         ),
                aes(x = vax_dt, y = 100*ci, colour = hcw_teacher, linetype = sex)) +
  geom_step() +
  facet_grid(age_grp~dose) +
  scale_y_continuous("Vaccinated dose (%)") +
  scale_x_date("Week beginning", date_breaks = "2 months", date_labels = "%b-%y", limits = as.Date(c("2020-12-08", "2021-06-30"))) +
  theme_ipsum() +
  scale_colour_manual("",
                      breaks = c("Healthcare worker, patient facing",
                                 "Healthcare worker, non-patient facing",
                                 "Teacher",
                                 "Neither"),
                      values = c("#422d54",
                                 "#758bd1",
                                 "#2d543d",
                                 "#d175b8")) +
  geom_vline(xintercept = as.Date("2020-08-24"), linetype = 2, colour = "gray45") +
  geom_vline(xintercept = as.Date("2020-12-23"), linetype = 2, colour = "gray45") +
  geom_vline(xintercept = as.Date("2021-02-26"), linetype = 2, colour = "gray45") +
  geom_vline(xintercept = as.Date("2021-04-23"), linetype = 2, colour = "gray45") +
  scale_linetype("")
plot2
saveRDS(plot2, "output/Vaccine_update_in_teachers.Rds")
tiff("output/Vaccine_update_in_teachers.tiff", res = 300, compression = "lzw", width = 10, height = 10, unit = "in")
plot2
dev.off()

### Create Table ----
# List of variables to summarise
myVars <- c("age_year", "female", "simd", "ethnicity", "shielding",  "como_count",  "hhd", "hhd_tch", "dose_ever")

# List of categorical variables
catVars <- c("female", "simd", "ethnicity", "shielding", "como_count",  "hhd", "hhd_tch", "dose_ever") 

# Create survey object
cc_data_svy_hb <- svydesign(~anon_id, probs = ~hb_simd_prb, data = cc_data)

baseline_chars <- CreateTableOne(vars = myVars, strata = "hcw_teacher",
                                 data = cc_data, factorVars = catVars, test = FALSE, includeNA = TRUE) 
baseline_chars_hb <- svyCreateTableOne(vars = myVars, strata = "hcw_teacher",
                                       data = cc_data_svy_hb, factorVars = catVars, test = FALSE, includeNA = TRUE)  

## Create object with type of teacher
baseline_sect <- CreateTableOne(vars = myVars, strata = "sector_grp",
                                data = cc_data, factorVars = catVars, test = FALSE, includeNA = TRUE)
baseline_sect_hb <- svyCreateTableOne(vars = myVars, strata = "sector_grp",
                                      data = cc_data_svy_hb, factorVars = catVars, test = FALSE, includeNA = TRUE)  


saveRDS(list(baseline_chars = baseline_chars, baseline_chars_hb = baseline_chars_hb, 
             baseline_sect = baseline_sect, baseline_sect_hb= baseline_sect_hb),
        "output/table1_components.Rds")

## check totals from weighted mean versus survey package
baseline_totl <- svyCreateTableOne(vars = "dose_ever",
                                data = cc_data_svy_hb, factorVars = "dose_ever", test = FALSE)
dt_lst <- as.Date(c("2020-12-08", "2021-07-14"))
dt_lst <- dt_lst[1]:dt_lst[2]


res <- map_dbl(dt_lst, function(x){
  weighted.mean( (!is.na(cc_data$vt1_dt) & ( as.integer(cc_data$vt1_dt) <= x))[], 
                 1/cc_data$hb_simd_prb[])
})
points(dt_lst, res, col = "red")
tail(res, 1)
