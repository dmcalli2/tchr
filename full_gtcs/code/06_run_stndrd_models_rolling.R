#6 rolling_effect_estimate_all_cases

TESTWEEKCUT <- 77

library(tidyverse)
library(survival)
library(lubridate)
library(janitor)


## Note unlike the analysis for "any case" shown in the main manuscript, this is limited to individuals testing positive 

## Depends
# code/full_gtcs/01_create_basefile.R

## Outputs
# "Outputs/model_rolling_rate_ratio.Rds"


## Clean and transform case control data ----
source("code/full_gtcs/01_create_basefile.R")
source("code/full_gtcs/02_process_ecoss.R")

## Limit data to relevant variables and collapse small healthboard into other
cc_data <- cc_data %>% 
  mutate(hb2019name = if_else(
    hb2019name %in% c("NHS Ayrshire and Arran", "NHS Borders", "NHS Dumfries and Galloway",
                      "NHS Fife", "NHS Forth Valley", "NHS Grampian", "NHS Greater Glasgow and Clyde",
                      "NHS Highland", "NHS Lanarkshire", "NHS Lothian", "NHS Tayside"),
    hb2019name,
    "Other")) %>%
  select(anon_id, is_case, hosp, stratum, 
         tch, hcw, hhd, hhd_tch, hh_tch_det,como_count, white, simd, 
         age_year, sex, age_grp,
         sector_grp, specimendate, severe, hb2019name, teaching_status, hh18_cat,
         npf) %>% 
  rename(test_date = specimendate)

## For speed exclude any strata with no teachers or healthcare workers, still a very large dataset
cc_data <- cc_data %>% 
  mutate(any_interest = tch==1 | hcw ==1 | hhd ==1 | hhd_tch ==1) %>% 
  group_by(stratum) %>% 
  mutate(strat_keep = any(any_interest)) %>% 
  ungroup() %>% 
  filter(strat_keep)

## Exclude stratum where case did not test positive
ecoss_vect <- ecoss %>% filter(ncov_result == "Positive") %>% 
  mutate(anon_test = paste(anon_id, test_date, sep = "_")) %>% 
  pull(anon_test) %>% 
  unique()
cc_data <- cc_data %>% 
  mutate(case_pos = paste(anon_id, test_date, sep = "_") %in% ecoss_vect) %>% 
  group_by(stratum) %>% 
  mutate(case_pos = any(case_pos)) %>% 
  ungroup()
cc_data %>% 
  count(is_case, case_pos)
cc_data <- cc_data %>% 
  filter(case_pos) %>% 
  select(-case_pos)

## Count number of teachers and hcws in case control study
cc_data %>% 
  count(tch, hcw)
cc_data %>% 
  count(hhd, hhd_tch)
cc_data %>% 
  count(hh18_cat, hhd)
cc_data %>% 
  count(hcw, hhd)
cc_data %>% 
  count(tch, hhd_tch)

## Count number of teachers and hcws in case control study
cc_data %>% 
  count(hh18_cat, hhd_tch)
## Count teaching status among teachers
cc_data %>% 
  count(tch, teaching_status)
## Run week by week analysis ----
cc_data <- cc_data %>% 
  mutate(test_week = Week2020(test_date))

rolldate <- cc_data %>% 
  filter(is_case ==1) %>% 
  distinct(test_week, stratum)


# 
warning("Applied limit to end of data, need to change this if update data")
rolls <- unique(rolldate$test_week) %>% sort()
minrol <- min(rolls) + 1
maxrol <- max(rolls) - 1
warning("Applied limit to end of data, in following line, need to change this if update data")
maxrol <- TESTWEEKCUT
rolls <- rolls[rolls >= minrol & rolls <= maxrol]
rollall <- list(st = rolls - 1, end = rolls + 1)

safeclogit <- safely(clogit)

  nstd <- map2(rollall$st, rollall$end, ~ {
    rollchoose <- rolldate %>% 
      filter(test_week >= .x, test_week <= .y)
    cc_data %>% 
      filter(stratum %in% rollchoose$stratum)
  }
  )
  nstd <- tibble(data = nstd, test_week = rolls)
  nstd$hcw <- map_dbl(nstd$data, ~ sum(.x$is_case[.x$hcw ==1]))
  nstd$tch <- map_dbl(nstd$data, ~ sum(.x$is_case[.x$tch ==1]))
  nstd$hhd <- map_dbl(nstd$data, ~ sum(.x$is_case[.x$hhd ==1]))
  nstd$nei <- map_dbl(nstd$data, ~ sum(.x$is_case[!.x$tch ==1 & !.x$hcw ==1 & !.x$hhd]))
  nstd$mdl <- map(nstd$data, ~ safeclogit(is_case ~ strata(stratum) + tch + hcw + hhd + hhd_tch + npf + 
                                            simd + white + como_count + hh18_cat,
                                   data = .x)$result %>% ChopSize()
  )
  nstd$mdl2 <- map(nstd$data, ~ safeclogit(is_case ~ strata(stratum) + tch + hcw + hhd + hhd_tch + npf,
                                    data = .x)$result %>% ChopSize()
  )
  nstd$data <- NULL
  saveRDS(nstd, "output/model_rolling_rate_ratio.Rds")
