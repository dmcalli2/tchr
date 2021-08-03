library(tidyverse)

source("code/full_gtcs/01_create_basefile.R")


s1s <- readRDS("output/main_analysis/hazard_ratio_supplementary_appendix.Rds")



## results Tables 1 and 2
t1_2 <- bind_rows(s1s$s1_over %>% 
                    filter(Individual == "Teachers overall", Household == "Teachers overall"),
                  s1s$s1_rest %>% 
                    filter(Individual == "Teachers overall", Household == "Teachers overall")
) 

## Tables
# Table 2 manuscript, teachers and HCWs and households
t2 <- t1_2 %>% 
  gather("Outcome", "value", `Any case`, Hospitalisation, Severe) %>% 
  spread(Parameter, value) 
t2 <- t2 %>% 
  mutate(Neither = "1") %>% 
  select(-Individual, - Household) %>% 
  filter(!is.na(Stratification))

t2 <- t2[ , c("Adjustment", "Time period", "Sex", "Age group (years)", "Stratification", 
              "Outcome", 
              "Neither",
              "Healthcare worker, patient facing", "Household member of healthcare worker", 
              "Household member of teacher", "Teacher",
              "Healthcare worker, non-patient facing")]

t3 <- bind_rows(s1s$s1_over %>% 
                  filter(Individual == "By sector", Household == "By sector"),
                s1s$s1_rest %>% 
                  filter(Individual == "By sector", Household == "By sector")
) 
t3 <- t3 %>% 
  filter(Parameter %in% c(  hh_tch_detPrimary = "Household member of primary teacher", hh_tch_detSecondary = "Household member of secondary teacher", 
                            hh_tch_detOther = "Household member other teacher", `sector_grp_tvNursery/Primary or Nursery` = "Nursery/Primary or Nursery", 
                            sector_grp_tvPrimary = "Primary", sector_grp_tvSecondary = "Secondary", 
                            sector_grp_tvOther = "Teacher in other sector")) %>% 
  gather("Outcome", "value", `Any case`, Hospitalisation) %>% 
  spread(Parameter, value) 

t3 <- t3 %>%
  filter(!is.na(Stratification)) %>% 
  mutate(Neither = "1") 
t3 <- t3 [, c("Adjustment", "Time period", "Sex", 
              "Age group (years)", "Stratification", "Outcome", "Household member of primary teacher", 
              "Neither",
              "Household member of secondary teacher", "Household member other teacher", 
              "Nursery/Primary or Nursery", "Primary", "Secondary", "Teacher in other sector")]


t2 <- t2 %>% 
  arrange(`Time period`, Sex, `Age group (years)`, Outcome, desc(Adjustment))
t3 <- t3 %>% 
  arrange(`Time period`, Sex, `Age group (years)`, Outcome, desc(Adjustment))


bystrata <- unique(t2$Stratification)
names(bystrata) <- bystrata
bystrata <- map(bystrata, ~ list(t2 = t2 %>% filter(Stratification == .x) %>% select(-Stratification),
                                 t3 = t3 %>% filter(Stratification == .x) %>% select(-Stratification)))

saveRDS(bystrata, "output/main_analysis/table_reg_appendix.Rds")



## Counts for Stable 2

## count individuals and household members overall
cc_data <- cc_data %>% 
  mutate(hcw_teacher2 = case_when(
    hcw_teacher == "Teacher" & tch ==0 ~ "Neither",
    hcw_teacher == "Teacher" ~ "Teacher",
    hcw_teacher == "Healthcare Worker, patient facing"     ~ "Healthcare worker, patient facing",
    hcw_teacher == "Healthcare worker, non-patient facing" ~ "Healthcare worker, non-patient facing",
    hhd_tch ==1 ~ "Household member of teacher",
    hhd ==1 ~ "Household member of healthcare worker",
    TRUE ~ "Neither"),
    factor(closure_status,
           levels = c("a_pre_reop", "b_reop2020", "c_closed", "d_mixed", "e_reop2021"),
           labels = c("Spring/Summer 2020 closed", "Autumn 2020 open", "Winter 2020 closed", "Winter 2021 mixed", "Spring 2021 reopened") ))


cc_data_smry <- cc_data %>% 
  arrange(desc(severe), desc(hosp), desc(is_case)) %>% 
  select(hcw_teacher2, closure_status, age_grp, sex, severe, hosp, case = is_case, anon_id) %>% 
  distinct(anon_id, .keep_all = TRUE)

cc_data_smry2 <- cc_data_smry %>% 
  group_by(age_grp, sex, closure_status, hcw_teacher2) %>% 
  summarise(
    cases   = sum(case ==1 & hosp ==0 & severe ==0),
    control = sum(as.integer(case ==0)),
    hosp_not_severe   = sum(hosp == 1 & severe ==0),
    severe_not_hospitalised = sum(severe ==1 & hosp ==0),
    hosp_severe = sum(severe ==1 & hosp ==1),
    tot = length(anon_id)) %>% 
  ungroup() %>% 
  mutate(sex = if_else(sex ==1, "Men", "Women")) %>% 
  select(Sex = sex, `Time period` = closure_status,
         `Age group (years)` = age_grp,
         `Group` = hcw_teacher2,
         Controls = control,
         `Cases only` = cases,
         `Hospitalised not severe` = hosp_not_severe,
         `Severe not hospitalised` = severe_not_hospitalised,
         `Hospitalised and severe` = `hosp_severe`)

saveRDS(cc_data_smry2, "output/main_analysis/table_count_appendix.Rds")
