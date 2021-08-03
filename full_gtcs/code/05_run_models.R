#05_run_models

## This allows you to run a sensitivity analysis with in-house HID for 2021 rather than uprn throughout
if(!exists("uprn_only")) uprn_only <- TRUE

# Runs clogit models for different combinations of data and formulae
## Runs mode quickly if exclude strata with no teachers, healthcare workers, household members of teachers
## or household members of healthcare workers
concisedata <- TRUE

# Depends
# code/full_gtcs/01_create_basefile.R

# Outputs
# cut-down model objects (data dropped) for subsequent tabulation
# "output/models_nodata/teach.Rds"
# output/models_nodata/severe.Rds
# "output/models_nodata/sector.Rds"

library(tidyverse)
library(survival)
library(lubridate)

for(sens_analysis in c("main_analysis",
                       "no1st",
                       "no2nd")[if(uprn_only == FALSE) 1 else 1:3]){

  if(sens_analysis == "main_analysis" & uprn_only == FALSE) sens_analysis = paste0(sens_analysis, "_uprn_hid")
  print(sens_analysis)
  
## Data ----
source("code/full_gtcs/01_create_basefile.R")

## Limit data to relevant variables 
cc_data <- cc_data %>% 
  select(anon_id, closure_status, is_case, hosp, stratum, 
         tch, hcw, hhd, hhd_tch, hh_tch_det, como_count, white, simd, 
         age_year, sex, age_grp,
         sector_grp_tv, specimendate, severe,
         hh04, hh1217, hh511, hh18,
         hid,
         hh18_cat,
         npf, vt, vt1, vt2,
         closure_status_clps) 
cc_data <- cc_data %>% 
  distinct()
## Limit data if running sens analysis
# if(sens_analysis == "main_analysis") {
# }
if(sens_analysis == "no1st") {
  cc_data <- cc_data %>% 
    filter(vt == "a_unvacc")
}
if(sens_analysis == "no2nd") {
  cc_data <- cc_data %>% 
    filter(vt %in% c("a_unvacc", "b_post1st"))
}

cc_data <- cc_data %>% 
  select(-anon_id)

## Model formulae ----
form_hosp <- 
  list(
    unad_hosp = hosp    ~ strata(stratum) + tch + hcw + hhd + hhd_tch + npf,
    adj_hosp  = hosp    ~ strata(stratum) + tch + hcw + hhd + hhd_tch + npf + hh18_cat + como_count + white + simd,
    unad_hosp_hh = hosp    ~ strata(stratum) + tch + hcw + hhd + hh_tch_det + npf,
    adj_hosp_hh  = hosp    ~ strata(stratum) + tch + hcw + hhd + hh_tch_det + npf + como_count + white + simd + hh18_cat)

form <- list(
  unad_case = is_case ~ strata(stratum) + tch + hcw + hhd + hhd_tch + npf ,
  adj_case  = is_case ~ strata(stratum) + tch + hcw + hhd + hhd_tch  + npf + hh18_cat + como_count + white + simd,
  unad_case_hh = is_case ~ strata(stratum) + tch + hcw + hhd + hh_tch_det + npf ,
  adj_case_hh  = is_case ~ strata(stratum) + tch + hcw + hhd + hh_tch_det + npf + como_count + white + simd  + hh18_cat)

form_sect_hosp <-
  list(
    unad_hosp = hosp    ~ strata(stratum) + sector_grp_tv + hcw + hhd + hhd_tch + npf ,
    adj_hosp  = hosp    ~ strata(stratum) + sector_grp_tv + hcw + hhd + hhd_tch + npf  + como_count + white + simd + hh18_cat,
    unad_hosp_hh = hosp    ~ strata(stratum) + sector_grp_tv + hcw + hhd + hh_tch_det + npf ,
    adj_hosp_hh  = hosp    ~ strata(stratum) + sector_grp_tv + hcw + hhd + hh_tch_det + npf  + como_count + white + simd + hh18_cat)

form_sect <- list(
  unad_case = is_case ~ strata(stratum) + sector_grp_tv + hcw + hhd + hhd_tch + npf ,
  adj_case  = is_case ~ strata(stratum) + sector_grp_tv + hcw + hhd + hhd_tch + npf + como_count + white + simd + hh18_cat,
  unad_case_hh = is_case ~ strata(stratum) + sector_grp_tv + hcw + hhd + hh_tch_det + npf,
  adj_case_hh  = is_case ~ strata(stratum) + sector_grp_tv + hcw + hhd + hh_tch_det + npf + como_count + white + simd + hh18_cat)

form_sev <-
  list(
    unad_sevr = severe  ~ strata(stratum) + tch + hcw + hhd + hhd_tch + npf,
    adj_sevr  = severe  ~ strata(stratum) + tch + hcw + hhd + hhd_tch + npf + como_count + white + simd + hh18_cat)

MakeCase <- function(mydf) {
  mydf %>% 
    mutate(case = if_else(is_case == TRUE, "Case", "Control")) %>% 
    count(tch, hcw, case) %>% 
    spread(case, n, fill = 0L)
}
MakeHosp <- function(mydf) {
  mydf %>% 
    mutate(case = if_else(hosp == 1, "Case", "Control")) %>% 
    count(tch, hcw, case) %>% 
    spread(case, n, fill = 0L)
}


## For speed in running exclude strata with no teachers or HCWs of either type ----
if(concisedata == TRUE){
  cc_data <- cc_data %>% 
    mutate(anyexp = (tch + hcw + hhd + hhd_tch + npf) >=1) %>% 
    group_by(stratum) %>% 
    mutate(any_tch_hcw = any(anyexp)) %>% 
    ungroup() %>% 
    filter(any_tch_hcw) %>% 
    select(-anyexp, -any_tch_hcw)
}

## Create data with different sub-groups ----
alldata <- cc_data %>% 
  mutate(allpt = "overall") %>% 
  group_by(allpt) %>% 
  nest()
closure <- cc_data %>% 
  group_by(closure_status) %>% 
  nest()
age_sex <- cc_data %>% 
  group_by(age_grp, sex) %>% 
  nest()
age_sex_closure <- cc_data %>% 
  group_by(age_grp, sex, closure_status) %>% 
  nest()

## Bind all data frames
cmbn <- bind_rows(
  univ = alldata,
  closure = closure,
  age_sex = age_sex,
  age_sex_closure = age_sex_closure,
  .id = "mdls"
)

# replace missing values to indicate includes all of these (eg age and sex)
cmbn <- cmbn %>% 
  mutate(allpt = if_else(is.na(allpt), "allof", as.character(allpt)),
         closure_status = if_else(is.na(closure_status), "allof", closure_status),
         sex = if_else(is.na(sex), "allof", as.character(sex)),
         age_grp = if_else(is.na(age_grp), "allof", age_grp)) %>% 
  ungroup()

cmbn_rv <- cmbn %>% 
  select(-data) %>% 
  group_by(mdls, allpt, closure_status, sex, age_grp) %>% 
  mutate(n = length(age_grp)) %>% 
  ungroup() 

rm(alldata, closure, age_sex, age_sex_closure)

cmbn$hosp_data <- map(cmbn$data, ~ .x %>% 
                        group_by(stratum) %>% 
                        mutate(hosp_strat = any(hosp ==1L)) %>% 
                        ungroup() %>% 
                        filter(hosp_strat) %>% 
                        select(-hosp_strat))

## Run models for all sub-divisions of data ----
safeclogit <- safely(clogit)

## Apply model for teacher status 
cmbn[names(form)] <- map(form, function(myform){
  print(myform)
  map(cmbn$data, ~ safeclogit(myform, data = .x)$result %>% ChopSize)
})
cmbn[names(form_hosp)] <- map(form_hosp, function(myform){
  print(myform)
  map(cmbn$hosp_data, ~ safeclogit(myform, data = .x)$result %>% ChopSize)
})
cmbn_teach <- cmbn %>%
  select(-data, -hosp_data)
saveRDS(cmbn_teach, paste0("output/", sens_analysis, "/teach.Rds"))

## Apply model for teacher types
cmbn[names(form_sect)] <- map(form_sect, function(myform){
  print(myform)
  map(cmbn$data, ~ safeclogit(myform, data = .x)$result %>% ChopSize)
})
cmbn[names(form_sect_hosp)] <- map(form_sect_hosp, function(myform){
  print(myform)
  map(cmbn$hosp_data, ~ safeclogit(myform, data = .x)$result %>% ChopSize)
})
cmbn_sect <- cmbn %>%
  select(-data, -hosp_data) 
saveRDS(cmbn_sect, paste0("output/", sens_analysis, "/sector.Rds"))

## Apply model for severe disease stratify by all cause and closure status only
cmbn_sevr <- cmbn %>% 
  filter(mdls %in%  c("univ", "closure"))

cmbn_sevr$data_sevr <- map(cmbn_sevr$data, ~ .x %>% 
                             group_by(stratum) %>% 
                             mutate(sevr_strat = any(severe ==1L)) %>% 
                             ungroup() %>% 
                             filter(sevr_strat) %>% 
                             select(-sevr_strat))

cmbn_sevr[names(form_sev)] <- map(form_sev, function(myform){
  print(myform)
  map(cmbn_sevr$data_sevr, ~ safeclogit(myform, data = .x)$result %>% ChopSize)
})
cmbn_sevr <- cmbn_sevr %>% 
  select(-data, -data_sevr, -hosp_data)
saveRDS(cmbn_sevr, paste0("output/",sens_analysis, "/severe.Rds"))

## Event counts post vaccination among teachers
hosp %>% 
  filter(vt != "a_unvacc", tch == 1) %>% 
  count(hosp, severe)
}
