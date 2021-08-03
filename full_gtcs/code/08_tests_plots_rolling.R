# 11_tests_as_denom

## Note that this runs on the 31/1/2021 data, not the more recent data.

## The following code attempts to account for the increased testing rates in teachers and HCWs relative to the general Scottish population.
## We generate summary statistics for the number tested and number tested positive.
## We also analyse the test positivity (number testing positive/number tested) within the case-control data.
## For the latter we generate new strata

## There are three datasets we use for this anlaysis
## The denominator counts of teachers, healthcare workers and the remaining population (by age, sex and SIMD - obtained from the teacher and HCW cohorts and from the mid-year estimates)
## The case-control study (1 case and 10 gp-age-sex matched controls)
## The ECOSS test database - numers of positive and negative tests for each person

## The analysis converts dates to weeks using the Week2020() function. 
## it follows the convention of lubridate::week() "returns the number of complete seven day periods that have occurred between the date and January 1st, plus one"
## BUT SO that it also WORKS after DECEMBER 31ST 2020 it uses January 1st 2020 as the reference so that weeks can be > 52


# Depends
# "code/full_gtcs/01_create_basefile.R" - case-control data
# "code/full_gtcs/02_process_ecoss.R" - lab data
# "Scratch_data/denom_pop_tch_hcw.Rds" - population denominator written in 03_cumulative_incidence_plots.R
# "output/model_rolling_rate_ratio_cleaned.Rds" - rolling models for standard analysis fo case control data (ie not test positivity)

# Outputs
# output/plot_test_outcomes.pdf panel plot of testing, odds ratios and hazard ratios for nay case and test positivity
# output/plot_test_outcomes.tiff panel plot of testing, odds ratios and hazard ratios for nay case and test positivity
# output/plot_test_outcomes.Rds panel plot of testing, odds ratios and hazard ratios for nay case and test positivity
# output/test_pos_models_periods.Rds -         model results and counts for test positivity by closure status for occupational groups 
# output/test_pos_models_periods_summary.RDS - model results and counts for test positivity by closure status for occupational groups - in summary tables
# output/test_pos_models_reopen_agesex.Rds   - model results and counts for test positivity by closure status age and sex for occupational groups
# output/model_rolling.Rds                   - rolling effect estimates with 3-week window for each week

## NOTE
# Some of the models take a long time to run, so there are conditional statements not to run the models if the output already exists.
# delete these files to run a fresh analysis
# if(!file.exists("Outputs/test_pos_models_periods.Rds"))
# if(!file.exists("Outputs/test_pos_models_reopen_agesex.Rds"))
# if(!file.exists("Outputs/model_rolling.Rds")){



TESTWEEKCUT <- 68

library(tidyverse)
library(survival)
library(lubridate)
library(janitor)

## Functions ----
source("code/full_gtcs/01_create_basefile.R")
source("code/full_gtcs/02_process_ecoss.R")

## Summary statistics for testing and testing positive each week for controls ----
## First identify all controls regardless of whether or not they are tested by occupation, age and sex
## Note total number of individuals who are "controls" and total number of controls (Allowing double counting of individuals) are very similar
cntrls_all <- cc_data %>% 
  filter(!is_case ==1) %>% 
  group_by(hcw, tch, age_grp, sex) %>% 
  summarise(cntrls = sum(!duplicated(anon_id))) %>% 
  ungroup()

## Limit case control set to tested individuals for subsequent analyses ----
## This leads to multiple rows per case and control
cc_data <- cc_data %>% 
  filter(!is.na(simd)) %>% 
  semi_join(ecoss %>% select(anon_id)) %>% 
  left_join(ecoss) %>% 
  as_tibble()

## Also limit to vaccinated individuals
cc_data <- cc_data %>% 
  filter(vt == "a_unvacc")

## For this "test positivity analysis" exclude any case who was a case despite never testing positive (eg based on diagnosis) even though had a test, should
## be fine if same for all groups
## very few cases are lost by excluding those with negative tests after their positive ones
cc_data %>% 
  count(is_case, ncov_result)
sum(!duplicated(cc_data$anon_id[cc_data$is_case]))
cc_data <- cc_data %>% 
  filter(! (is_case & ncov_result == "Negative"))
sum(!duplicated(cc_data$anon_id[cc_data$is_case]))

## Create new strata by health board and test week ----
## 4 implausible dates, drop these as too early. Also drop ones before 1st of February
## also drop dates before
cc_data <- cc_data %>% 
  filter(!test_date %in% as.Date(c("2022-02-17", "2022-04-29", "2022-11-28", "2032-01-25")),
         !test_date < as.Date("2020-02-01"))

cc_data <- cc_data %>% 
  mutate(test_week = Week2020(test_date)) %>% 
  arrange(hb2019name, test_week) %>% 
  mutate(stratum_new = cumsum(!duplicated(paste(hb2019name, test_week, sep = "_")))) 
nrow(cc_data)

## Individuals cannot be controls more than once in a stratum, or a control to themself. This code removes these individuals
cc_data <- cc_data %>% 
  arrange(desc(is_case)) %>% 
  distinct(stratum_new, anon_id, .keep_all = TRUE)
person <- readRDS("Data/person_anon.Rds")
nrow(cc_data)
cc_data <- cc_data %>%
  left_join(person %>% mutate(anon_id = as.character(anon_id))) %>%
  distinct(stratum_new, patid, .keep_all = TRUE)
nrow(cc_data)

## Limit data to where have one of the exposures of interest within a stratum for modelling only, not for event rates
cc_data_limit <- cc_data %>% 
  mutate(anyexp = (tch + hcw + hhd + hhd_tch + npf) >=1) %>% 
  group_by(stratum_new) %>% 
  mutate(any_tch_hcw = any(anyexp)) %>% 
  ungroup() %>% 
  filter(any_tch_hcw) %>% 
  distinct(anon_id) %>% 
  pull()
nrow(cc_data)

# Create lookup table from test date to start of week
# the first week is the 29th of January
test_week_st <- cc_data %>% 
  distinct(test_week) %>% 
  arrange(test_week) %>% 
  mutate(week_start = StartWeek(test_week))

## Count the number of tests among the CONTROLS. This includes all controls whether or not 
## there is a case within the same test_week and healthboard
roll <- cc_data$test_week %>% unique() %>% sort()
cntrl_tst <- map(roll, ~ cc_data %>% 
                     filter(test_week < TESTWEEKCUT + 1, test_week == .x, is_case == 0) %>% 
                     group_by(age_grp, sex, hcw, tch) %>% 
                     summarise(tested = sum(!duplicated(anon_id))) %>% 
                     ungroup())
cntrl_tst <- tibble(test_week = roll, data = cntrl_tst) %>% 
  unnest(data)

## Join the controls tested onto the controls dataset
## note that there are some controls in every stratum so allows re-weighting
## first join is to capture all weeks
cntrl_tst <- cntrls_all %>% 
  left_join(cntrl_tst %>% distinct(age_grp, sex, test_week)) %>%
  left_join(cntrl_tst) %>% 
  mutate_at(vars(tested), ~ if_else(is.na(.x), 0L, .x)) %>% 
  arrange(test_week)

## update summary with total number who tested positive in the POPULATION (ie the number of test positive cases)
popn_pos <- map(roll, ~ cc_data %>% 
                     filter(test_week == .x) %>% 
                     group_by(age_grp, sex, hcw, tch) %>% 
                     summarise(tested_pos = sum( as.integer(!duplicated(anon_id)) *
                                                   as.integer(ncov_result == "Positive"))) %>% 
                     ungroup())
popn_pos <- tibble(test_week = roll, data = popn_pos) %>% 
  unnest(data)
## Add those who tested positive in the population to those who were tested in the control
tstd_all <- cntrl_tst %>% 
  left_join(popn_pos) %>% 
  mutate_at(vars(tested_pos), ~ if_else(is.na(.x), 0L, .x)) %>% 
  arrange(test_week) %>% 
  rename(tested_pos_all = tested_pos,
         tested_cntrls = tested)

## add population denominator
denom <- readRDS("Scratch_data/denom_pop_tch_hcw.Rds")
denom <- denom %>% 
  mutate(hcw = if_else(hcw_teacher == "Healthcare Worker", 1L, 0L),
         tch = if_else(hcw_teacher == "Teacher", 1L, 0L)) %>% 
  select(-hcw_teacher) %>% 
  mutate(sex = if_else(male ==1, 1, 2)) %>% 
  select(-male)
tstd_all <- denom %>% 
  inner_join(tstd_all) 
denom_agg <- denom %>% 
  group_by(hcw, tch) %>% 
  summarise(pop_all = sum(pop)) %>% 
  ungroup()
## Impute the numbers tested in the population from the proportion tested in the controls
## then use this to estimate the proportion who tested positive in the population
## note expecting this to be NaN where there are no controls for tested_all
## and to be NaN where there are no tested_cntrs for prp_tstd
## Since this 
tstd_all <- tstd_all %>% 
  mutate(tested_all = pop * tested_cntrls/cntrls,
         prp_tstd = tested_pos_all/tested_all)
tstng_smry <- tstd_all

## Assign closed/open status based on testing date
lagdays <- 14
cc_data <- cc_data %>%
  mutate(
    closure_status = case_when(
      test_date <= as.Date("2020-08-12") + lagdays ~  "a_pre_reop",
      test_date >  as.Date("2020-08-12") + lagdays & test_date <= as.Date("2020-12-22") + lagdays ~  "b_reop2020",
      test_date >  as.Date("2020-12-22") + lagdays ~ "c_mixed"))

## Limit to strata where have at least one control and at least one cases
## Count number of cases and controls per stratum, some strata either 0 cases or 0 controls, and many strata only 1 control
cc_data_smry <- cc_data %>% 
  group_by(stratum_new, is_case) %>% 
  count() %>% 
  ungroup()
cc_data_smry <- cc_data_smry %>% 
  spread(is_case, n, fill = 0L) %>% 
  rename(case = `1`, control = `0`)
cc_data_smry_agg <- cc_data_smry %>% 
  count(case, control)
cc_data_smry_agg %>% 
  filter(case == 0 | control == 0) 
cc_data_smry_agg %>% 
  filter(!case == 0 & !control == 0) %>% 
  arrange(desc(n))

## Data with at least one case and one control in a stratum
cc_data_smry_agg %>% 
  filter(!case == 0 & !control == 0) %>% 
  summarise(tot = sum((case + control)*n),
            cases = sum(case),
            controls = sum(control))
## Data without at least one case and one control in a stratum, 12 cases and 31529 controls
cc_data_smry_agg %>% 
  filter(case == 0 | control == 0) %>% 
  summarise(tot = sum((case + control)*n),
            cases = sum(case*n),
            controls = sum(control*n))
cc_data %>% 
  semi_join(cc_data_smry %>% filter(control ==0 | case ==0)) %>% 
  count(is_case)
cc_data <- cc_data %>% 
  anti_join(cc_data_smry %>% filter(control ==0 | case ==0))

## Run for different time periods ----
if(!file.exists("output/test_pos_models_periods.Rds")){
cc_data_nst <- cc_data %>% 
  filter(test_week <= TESTWEEKCUT) %>% 
  group_by(closure_status) %>% 
  nest()

cc_data_nst$res <- map(cc_data_nst$data, function(cc_data_choose){
  ## Run regression models for all people for all periods, note that stratum new must be a character variable
  mod1 <- glm( is_case ~ as.character(stratum_new) + tch + hcw + hhd + hhd_tch + npf + hh18_cat + age_grp + I(age_year/10) + factor(sex),
               data = cc_data_choose %>% filter(anon_id %in% cc_data_limit), family = "binomial")
  summary(mod1)
  mod2 <- update(mod1, . ~ . + como_count + simd + white + hh18_cat)
  summary(mod2)

bind_rows(unad = broom::tidy(mod1),
          adj = broom::tidy(mod2),
          .id = "models")
})


cc_unnst <- cc_data_nst %>% 
  select(closure_status, res) %>% 
  unnest("res") %>% 
  mutate(est = exp(estimate),
         lci = exp(estimate - 1.96*std.error),
         uci = exp(estimate + 1.96*std.error)) %>% 
  mutate_at(vars(est, lci, uci), function(x)  formatC(round(x, 2), digits = 2, format = "f", flag = "0")) %>% 
  mutate(res = paste0(est, " (", lci, "-", uci, ")")) %>% 
  select(closure_status, models, term, res) %>% 
  filter(term %in% c("hcw", "tch")) %>% 
  spread(term, res) %>% 
  mutate(Neither = "1") %>% 
  select(closure_status, models, Neither, HCW = hcw, Teacher = tch) %>% 
  arrange(closure_status, desc(models))

cc_data_nst$cnts <- map(cc_data_nst$data, ~ .x %>% 
  mutate(status = case_when(
    tch == 1 ~ "Teacher",
    hcw ==1 ~ "HCW",
    TRUE ~ "Neither")) %>% 
  count(is_case, status) %>% 
  mutate(case = if_else(is_case ==1, "Cases", "Controls")) %>% 
  select(-is_case) %>% 
  spread(case, n, fill = 0L) %>% 
  mutate(cc = paste(Cases, Controls, sep = "/")) %>% 
  select(status, cc) %>% 
  spread(status, cc, fill = 0L) %>% 
  select(Neither, everything()))
cc_unnst2 <- cc_data_nst %>% 
  mutate(models = "") %>% 
  select(closure_status, models, cnts) %>% 
  unnest("cnts")

cc_unnst3 <- bind_rows(cc_unnst2,
                       cc_unnst) %>% 
  arrange(closure_status) %>% 
  ungroup()
cc_data_nst$data <- NULL
clos_period_res <- cc_unnst3
saveRDS(cc_data_nst, "output/test_pos_models_periods.Rds") 
saveRDS(cc_unnst3, "output/test_pos_models_periods_summary.RDS")
}
cc_data_nst <- readRDS("output/test_pos_models_periods.Rds")
clos_period_res <- readRDS("output/test_pos_models_periods_summary.RDS")


## stratify models by age and sex only for re-open period ----
cc_data_reopen <- cc_data %>% 
  filter(closure_status == "b_reop2020", test_week <= TESTWEEKCUT)

if(!file.exists("output/test_pos_models_reopen_agesex.Rds")){
cc_data_nst <- cc_data_reopen %>% 
  group_by(age_grp, sex) %>% 
  nest()
cc_data_nst$res <- map(cc_data_nst$data, ~ glm( is_case ~ as.character(stratum_new) + 
                                                  tch + hcw + I(age_year/10) + simd + white + 
                                                  hhd + hhd_tch + npf + hh18_cat + como_count, data = .x %>% filter(anon_id %in% cc_data_limit), family = "binomial") )
cc_data_nst$res2 <- map(cc_data_nst$res, broom::tidy)
cc_data_res <- cc_data_nst %>% 
  select(sex, age_grp, res2) %>% 
  unnest(res2) %>% 
  ungroup()
cc_data_res %>% 
  filter(term == "tch") %>% 
  arrange(sex, age_grp)
cc_data_res %>% 
  filter(term == "hcw") %>% 
  arrange(sex, age_grp)
saveRDS(cc_data_res, "output/test_pos_models_reopen_agesex.Rds")
}
cc_data_res <- readRDS("output/test_pos_models_reopen_agesex.Rds")

## Do rolling effect estimates for 3 week periods ----
rolldate <- cc_data %>% 
  filter(is_case == 1) %>% 
  distinct(test_week, stratum_new)

# limit test week
warning("Applied limit to end of data, need to remove this if update data")
rolls <- unique(rolldate$test_week) %>% sort()
minrol <- min(rolls) + 1
maxrol <- max(rolls) - 1
warning("Applied limit to end of data, in following line, need to remove this if update data")
maxrol <- TESTWEEKCUT
rolls <- rolls[rolls >= minrol & rolls <= maxrol]
rollall <- list(st = rolls - 1, end = rolls + 1)

if(!file.exists("output/model_rolling.Rds")){
nstd <- map2(rollall$st, rollall$end, ~ {
  # browser()
  rollchoose <- rolldate %>% 
         filter(test_week >= .x, test_week <= .y)
  cc_data %>% 
    filter(stratum_new %in% rollchoose$stratum_new)
  }
  )
nstd <- tibble(data = nstd, test_week = rolls)
nstd$hcw <- map_int(nstd$data, ~ sum(as.integer(.x$is_case)[.x$hcw ==1]))
nstd$tch <- map_int(nstd$data, ~ sum(as.integer(.x$is_case)[.x$tch ==1]))
nstd$nei <- map_int(nstd$data, ~ sum(as.integer(.x$is_case)[!.x$tch ==1 & !.x$hcw ==1]))

nstd$mdl <- map(nstd$data, ~ glm(is_case ~ as.character(stratum_new) + tch + hcw + age_grp + I(age_year/10) + simd + 
                                   white + hhd + hhd_tch + npf + hh18_cat + como_count + factor(sex),
                                data = .x %>% filter(anon_id %in% cc_data_limit), 
                                family = "binomial")
)
nstd$mdl2 <- map(nstd$data, ~ glm(is_case ~ as.character(stratum_new) + tch + hcw + npf + age_grp + I(age_year/10) + factor(sex),
                                 data = .x %>% filter(anon_id %in% cc_data_limit), 
                                 family = "binomial")
)

nstd$res <- map(nstd$mdl, broom::tidy)
nstd$res2 <- map(nstd$mdl2, broom::tidy)

nstd$mdl <- map(nstd$mdl, ~ stripGlmLR(.x))
nstd$mdl2 <- map(nstd$mdl2, ~ stripGlmLR(.x))

## Extract model results into a dataframe
unnstd <- nstd %>% 
  mutate(test_week = rolls) %>% 
  select(test_week, res) %>% 
  unnest(res) %>% 
  filter(term %in% c("tch", "hcw", "hhd", "hhd_tch")) %>% 
  select(test_week, term, estimate, std.error)

unnstd2 <- nstd %>% 
  mutate(test_week = rolls) %>% 
  select(test_week, res2) %>% 
  unnest(res2) %>% 
  filter(term %in% c("tch", "hcw", "hhd", "hhd_tch")) %>% 
  select(test_week, term, estimate, std.error)

unnstd_ref <- unnstd %>% 
  distinct(test_week) %>% 
  mutate(term = "nei",
         estimate = 0)
## Bind effect estimates from each model
unnstd <- bind_rows(unnstd,
                    unnstd_ref) %>% 
  arrange(test_week, term)
unnstd2 <- bind_rows(unnstd2,
                    unnstd_ref) %>% 
  arrange(test_week, term)
unnstd <- bind_rows(minim = unnstd2,
                    fully = unnstd, .id = "models")
saveRDS(unnstd, "output/model_rolling.Rds")
} else unnstd <- readRDS("output/model_rolling.Rds")


## Calculate proportion of positive tests at level of occupation (aggregate over age and sex) ----
tstng_smry2 <- tstng_smry %>% 
  group_by(hcw, tch, test_week) %>% 
  summarise_at(vars(tested_pos_all, tested_all), sum) %>% 
  ungroup() %>% 
  mutate(estimate = tested_pos_all/tested_all,
         term = case_when(hcw ==1 ~ "hcw",
                          tch ==1 ~ "tch",
                          TRUE ~ "nei")) %>% 
  select(test_week, term, estimate) 


## Use cowplot to join figures
tstng_smry3 <- tstng_smry %>% 
  mutate(occ = case_when(
    hcw == 1 ~ "Healthcare worker",
    tch == 1 ~ "Teacher",
    TRUE ~ "Neither"),
  wk_st = StartWeek(test_week)) %>% 
  group_by(occ, wk_st) %>% 
  summarise_at(vars(tested_pos_all, tested_all, pop), sum) %>% 
  ungroup() 
  
tstng_smry3 <- tstng_smry3  %>% 
  mutate(tested_pos_prp = tested_pos_all/tested_all)
    
denom_agg <- denom_agg %>% 
  mutate(occ = case_when(
    hcw == 1 ~ "Healthcare worker",
    tch == 1 ~ "Teacher",
    TRUE ~ "Neither")) %>% 
  mutate(pop_tch = min(pop_all))
tested_raw <- ggplot(tstng_smry3 %>% inner_join(denom_agg), 
                     aes(x = wk_st, y = 100*tested_all/pop_all, colour = occ)) +
  geom_point(size = 0.25) +
  geom_line() +
  scale_y_continuous("Individuals tested (%)") 
tested_raw

tested_rat <- ggplot(tstng_smry3, 
                 aes(x = wk_st, y = 100*tested_pos_all/pop, colour = occ)) +
  geom_point(size = 0.25) +
  geom_line() +
  scale_y_continuous("Indiv. tstd. positive / indiv. (%)")
tested_rat

tested_pos <- ggplot(tstng_smry3, 
                 aes(x = wk_st, y = 100*tested_pos_all/tested_all, colour = occ)) +
  geom_point(size = 0.25) +
  geom_line() +
  scale_y_continuous("Indiv. tstd. positive / indiv. tstd. (%)") +
  theme(axis.text.x = element_text(angle = 90)) 
tested_pos

unnstd2 <- unnstd %>% 
  mutate(wk_st = StartWeek(test_week),
         wave = case_when(
           wk_st < as.Date("2020-06-12") ~ "First wave",
           wk_st > as.Date("2020-08-01") ~ "Second wave",
           TRUE ~ "interwave"),
          occ = case_when(
            term == "hcw" ~ "Healthcare worker",
            term == "tch" ~ "Teacher",
            TRUE ~ "Neither")) #%>% 
  # filter(wave != "interwave")

for_roll_plot <- unnstd2 %>% 
  mutate(occ = case_when(
    term %in% c("hcw", "hhd") ~ "Healthcare worker",
    term %in% c("tch", "hhd_tch") ~ "Teacher",
    TRUE ~ "Neither"),
    hhd = factor(term %in% c("hhd", "hhd_tch"), 
                 levels = c(F,T), 
                 labels = c("Individual", "Household member"))) %>% 
  filter(wk_st <= StartWeek(TESTWEEKCUT),
         models == "fully",
         estimate > -5) %>% 
  mutate(uci= exp(estimate + 1.96*std.error),
         lci = exp(estimate - 1.96*std.error),
         estimate = exp(estimate))

test_or <- ggplot(for_roll_plot,  #%>%
                    # filter(hhd == "Individual"), 
                 aes(x = wk_st, y = estimate, ymin = lci, ymax = uci,
                     colour = occ,
                     linetype = hhd,
                     group = interaction(occ, hhd, wave))) +
  geom_point(size = 0.25) +
  geom_line() +
  scale_y_log10("Odds ratio", limits = c(1/10, 13)) +
  scale_color_discrete("")
test_or

## repeat for RR
test_rr_data <- readRDS("output/model_rolling_rate_ratio_cleaned.Rds")

nei <- test_rr_data %>% 
  distinct(models, test_week) %>% 
  mutate(params = "nei",
         est = 1)
test_rr_data <- bind_rows(test_rr_data,
                          nei)
test_rr_data <- test_rr_data %>%
  rename(term = params, estimate = est) %>%
  mutate(
    wk_st = StartWeek(test_week),
    wave = case_when(
      wk_st < as.Date("2020-06-12") ~ "First wave",
      wk_st > as.Date("2020-08-01") ~ "Second wave",
      TRUE ~ "interwave"
    ),
    occ = case_when(
      term %in% c("hcw", "hhd") ~ "Healthcare worker",
      term %in% c("tch", "hhd_tch") ~ "Teacher",
      TRUE ~ "Neither"
    ),
    hhd = factor(
      term %in% c("hhd", "hhd_tch"),
      levels = c(F, T),
      labels = c("Individual", "Household member")
    )
  ) %>%
  filter(
    wk_st <= StartWeek(TESTWEEKCUT),
    models == "fully",
    term %in% c("tch", "nei", "hhd", "hcw", "hhd_tch")
  ) 


test_rr <- ggplot(test_rr_data, 
                  aes(x = wk_st, y = estimate, 
                      colour = occ, 
                      group = interaction(occ, hhd),
                      linetype = hhd)) +
  geom_point(size = 0.25) +
  geom_line() +
  scale_y_log10("Rate ratio", limits = c(1/10, 13)) +
  scale_color_discrete("") +
  scale_linetype("") +
  theme(axis.text.x = element_text(angle = 90), legend.direction = "horizontal") 
test_rr

legend <- cowplot::get_legend(test_rr)
allplot <- list(
  tested_raw, tested_rat, tested_pos,
  NULL, test_rr, test_or)
allplot <- map(allplot, ~ .x + theme(legend.position='none'))
allplot <- map(allplot, ~ .x + 
                 scale_x_date("Date tested", date_breaks = "months", date_labels = "%b-%y", limits = c(as.Date("2020-03-01"), 
                                                                                                          as.Date("2021-05-01"))) +
                 theme(axis.text.x = element_text(angle = 90)))
allplot[[4]] <- legend
library(cowplot)
final_panel <- plot_grid(plotlist = allplot,labels = c("A", "B", "C", "" ,"D", "E"), label_x = 0.5,  nrow = 2)
final_panel
tiff("output/plot_test_outcomes.tiff", res = 600, compression = "lzw", height = 15, width = 30, unit = "cm")
final_panel
dev.off()

pdf("output/plot_test_outcomes.pdf", paper = "a4r")
final_panel 
dev.off()

saveRDS(final_panel, "output/plot_test_outcomes.Rds")


## Test positivity for MB
col_choose <- c("Not teacher or HCW" = "#d1ab75",
                "Teacher" = "#2d543d")

tested_pos_recent <- ggplot(tstng_smry3 %>% 
                              filter(wk_st >= as.Date("2020-07-01"), occ %in% c("Neither", "Teacher")) %>% 
                              mutate(occ = if_else(occ == "Neither", "Not teacher or HCW", "Teacher")), 
                     aes(x = wk_st, y = 100*tested_pos_all/tested_all, colour = occ)) +
  geom_point(size = 0.25) +
  geom_line() +
  scale_x_date("Week beginning", date_breaks = "week", date_labels = "%d-%b", limits = as.Date(c("2020-06-29", "2021-01-04"))) +
  scale_y_continuous("Indiv. tstd. positive / indiv. tstd. (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_colour_manual("", breaks = names(col_choose), 
                      values = col_choose)
tested_pos_recent

col_choose2 <- c("Teacher" = "#2d543d",
                "Teacher's household member" = "#8fd175")
odds_ratio_recent <- ggplot(for_roll_plot %>% 
                               filter(wk_st >= as.Date("2020-07-01"), term %in% c("tch", "hhd_tch")) %>% 
                               mutate(cmpr = case_when(
                                 # term == "nei" ~ "General population",
                                 term == "tch" ~ "Teacher",
                                 term == "hhd_tch" ~ "Teacher's household member")
                               ),
                               aes(x = wk_st, y = estimate, 
                                   colour = cmpr, fill = cmpr, ymin = lci, ymax = uci)) +
  geom_line() +
  geom_ribbon(alpha = 0.25, colour = NA) +
  scale_y_log10("Odds ratio (versus general population*)", breaks = c(0.25, 0.5, round(1/1.5,1), 1, 1.5, 2, 3, 4)) +
  scale_x_date("Week beginning", date_breaks = "2 week", date_labels = "%d-%b") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_manual("", breaks = names(col_choose2), values = col_choose2) +
  scale_fill_manual("", breaks = names(col_choose2), values = col_choose2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_cartesian(ylim = c(0.5, 5))
odds_ratio_recent


tiff("output/tch_mb1.tiff", res = 300, width = 20, height = 12, unit = "cm", compression = "lzw")
tested_pos_recent
dev.off()
tiff("output/tch_mb2.tiff", res = 300, width = 20, height = 12, unit = "cm", compression = "lzw")
odds_ratio_recent
dev.off()

