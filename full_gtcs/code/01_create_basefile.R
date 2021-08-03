## Pulls together and cleans case-control, teacher and healthcare workers data for analysis.

# Is set-up to only run once per session as all subsequent scripts call this script. 
# Re-run by deleting "Scratch_data/basefile.Rds" or changing runanyway to TRUE

# Outputs of this script are
# Scratch_data/basefile.Rds - which is the cleaned case control data with the additional datasets added and
# "Scratch_data/cleaned_teachers.Rds" - a completely de-identified dataset to create a denominator table for teachers

# Code allows other scripts us to pass runanyway as TRUE before a source statement
if(!exists("runanyway")) runanyway <- FALSE
## Code to switch between UPRN only and UPRN/HID; use this to runa  sensitivity analysis using in-house household algorithm rather than UPRN
if(!exists("uprn_only")) uprn_only <- TRUE

## Packages ----
library(here)
library(janitor)
library(tidyverse)

## Functions ----
## Calculate the start of each "week" - note is in reference to the 1st of Jan 2020
StartWeek <- function(x) {
  as.Date("2020-01-01") + 7 * (x - 1)
}

#ConvertWeek
Week2020 <- function(x) {
  a <- (x - as.Date("2020-01-01")) / 7
  a <- round(as.integer(a))
  1 + a
}

ChopSize <- function(x) {
  if(!any(class(x) == "coxph")) x <- "" else{
    x$residuals <- NULL
    x$y <- NULL
    x$linear.predictors <- NULL
  }
  x
}

stripGlmLR = function(cm) {
  cm$y = c()
  cm$model = c()
  
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c()  
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  
  
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()
  
  cm
}


DropLt5 <- function(mtrx) {
  mtrx[parse_integer(mtrx) <=5 & !is.na(parse_integer(mtrx))] <- "-"
  mtrx
}

FormatP <- function(x) {
  case_when(
    x <0.001 ~ "<0.001",
    x <0.01 ~ formatC(x, format = "f", digits = 3),
    TRUE ~ formatC(x, format = "f", digits = 2)
  )
}
ExtractEst <- function(model){
  # Take estimate and 95% CI
  a <- summary(model)
  a <- a$conf.int
  a [, c(1, 3, 4), drop = FALSE]
}
CleanEst <- function(x){
  params <- rownames(x)
  x <- as_tibble(x)
  names(x) <- c("est", "lci", "uci")
  x$params <- params
  x %>% 
    select(params, everything())
}
ConvertListModelSmriesDf <- function(modellist){
  ## Convert model summaries into a single dataframe
  modeldfs <- map(modellist, function(x) {
    params <- rownames(x)
    x <- as_tibble(x)
    names(x) <- c("est", "lci", "uci")
    x$params <- params
    x %>% 
      select(params, everything())
  })
  modeldfs <- bind_rows(modeldfs, .id = "model")
  modeldfs_neat <- modeldfs %>% 
    mutate_at(vars(est, lci, uci), function(x)  formatC(round(x, 2), digits = 2, format = "f", flag = "0")) %>% 
    mutate(res = paste0(est, " (", lci, "-", uci, ")"))
  modeldfs$res <- modeldfs_neat$res
  rm(modeldfs_neat)
  modeldfs
}


##  Read data ----

## If already run and want UPRN only analysis, read in instead of re-running
if (runanyway == FALSE &
    uprn_only == FALSE &
    file.exists("Data/basefile.Rds")) {
  cc_data <-
    readRDS("Data/basefile.Rds")
## If already run and want UPRN/in-house analysis with read in instead of re-running
} else if (runanyway == FALSE & 
           uprn_only == TRUE &
           file.exists("Data/basefile_uprn_only.Rds")) {
  cc_data <- readRDS("Data/basefile_uprn_only.Rds")
} else {
  
## read in comorbidity data
cc_all <- readRDS("Scratch_data/comorbidity_count_jul2021.Rds")

# Read in teachers data
teachers <- readRDS("//TOPLEVEL/Case_control/data/2021-07-28/tchrs_linked_ANON_2021-07-28.rds")
teachers <- teachers %>% 
  distinct()

# Read in case control
cc_data <- readRDS(
"//TOPLEVEL/Case_control/data/2021-07-28/CC_linked_ANON_2021-07-28.rds"
) %>%
  clean_names()

## limit ages early to reduce data size
cc_data <- cc_data %>% 
  filter(age_years >=21, age_years <=65)

## Limit to events on or after February 1st 2020
cc_data <- cc_data %>% 
  rename(specimendate = specimen_date) %>% 
  filter(specimendate >= as.Date("2020-02-01"))

## Make selection of household variables based on UPRN/in-house selection
if(uprn_only == TRUE) {
cc_data[ , c("hid", "hh04", "hh511", "hh1217", "hh18")] <- NULL
cc_data <- cc_data %>% 
  rename(hid = hid_uprn,
         hh04 = uprn_0_4,
         hh511 = uprn_5_11,
         hh1217 = uprn_12_17,
         hh18 = uprn_adult) } else {
cc_data <- cc_data %>% 
  rename(inh_hid = hid, 
         inh_hh04 = hh04, 
         inh_hh511 = hh511, 
         inh_hh1217 = hh1217, 
         inh_hh18 = hh18)
cc_data[cc_data$specimendate <= as.Date("2020-12-31")   , c("hid", "hh04", "hh511", "hh1217", "hh18")] <- 
  cc_data[cc_data$specimendate <= as.Date("2020-12-31") , c("hid_uprn", "uprn_0_4", "uprn_5_11", "uprn_12_17", "uprn_adult")]

cc_data[cc_data$specimendate > as.Date("2020-12-31") , c("hid", "hh04", "hh511", "hh1217", "hh18")] <- 
  cc_data[cc_data$specimendate > as.Date("2020-12-31") , c("inh_hid", "inh_hh04", "inh_hh511", "inh_hh1217", "inh_hh18")]

cc_data <- cc_data[, setdiff(names(cc_data), c("hid_uprn", "uprn_0_4", "uprn_5_11", "uprn_12_17", "uprn_adult",
                                               "inh_hid", "inh_hh04", "inh_hh511", "inh_hh1217", "inh_hh18"))]
}

## remove 3.2% of households in those aged 21 to 65 with more than 10 adults
cc_data <- cc_data %>% 
  filter(hh18 <= 10)
  
# Read in HCW data
hcw <- readRDS("//TOPLEVEL/HCW/hcw2021/processed_data/hcw_for_tch_analysis20210728.Rds")
npf <- hcw$npf

## Identify healthcare workers who are undetermined or in the excluded non-patient facing (ie those who could be working in hospitals)
hcw_drop <- hcw$drop_hcw

## identify hcw household members (PF only)
pf_hhd <- hcw$pf_hhld
## Identify patient facing healthcare workers
hcw <- hcw$pf

# read in sicsag data
sicsag <-
  readRDS(
    "//TOPLEVEL/Case_control/data/2021-07-28/CC_SICSAG_ANON_2021-07-28.rds"
  ) %>%
  clean_names() 

sicsag <- sicsag %>% 
  filter(covid_icu_or_hdu %in% c(1,3))

# identify unique anon_id/admission dates
sicsag <- sicsag %>% 
  select(anon_id, admit_unit) %>% 
  mutate(admit_unit = as.Date(admit_unit))
sicsag <- sicsag %>% 
  distinct()

# Approx 4000 related to cases
sicsag <- sicsag %>% 
  inner_join(cc_data %>% filter(is_case == TRUE) %>% distinct(anon_id, specimendate))

## Read in hosp data
rapid <-
  readRDS(
    "//TOPLEVEL/Case_control/data/2021-07-28/CC_RAPID_ANON_2021-07-28.rds"
  ) %>%
  clean_names() 
rapid <- rapid %>% 
  select(-patient_ethnic_group_code) %>% 
  distinct()

## Limit rapid to cases, this definition of cases does not include the 897 without a known positive test
rapid <- rapid %>% 
  semi_join(cc_data %>%
              filter(is_case == TRUE) %>%
              select(anon_id))
## Join to cc_data to get specimendates
rapid <- rapid %>% 
  inner_join(cc_data %>%
               filter(is_case == TRUE) %>% 
               select(anon_id, specimendate))

rapid <- rapid %>% 
  mutate(inhosp = if_else(specimendate >= admission_date & 
                            (is.na(discharge_date) | (specimendate <= discharge_date)), 
                          1L, 0L),
         adm28  = if_else(specimendate < admission_date & 
                            (as.integer(admission_date - specimendate) <= 28 ),
                          1L, 0L)
         )
## limit to one result per specimendate, take any in-hospitalisation or admission during this period
rapid <- rapid %>% 
  group_by(anon_id, specimendate) %>% 
  summarise_at(vars(inhosp, adm28), max) %>% 
  ungroup()
 
cc_data <- cc_data %>% 
  left_join(rapid) %>% 
  mutate_at(vars(inhosp, adm28), ~ if_else(is.na(.x), 0L, .x))

## 
sicsag <- sicsag %>% 
  filter(admit_unit >= specimendate, (admit_unit - specimendate) <= 21 ) %>% 
  distinct(anon_id, specimendate) %>% 
  mutate(icu_hdu_flag = 1L)

vax <- readRDS("Data/vaccination_cleaned.Rds")
vax_dt <- vax %>% 
  mutate(vax_seq = paste0("vt", vax_seq)) %>% 
  select(anon_id, vax_seq, vax_date) %>% 
  spread(vax_seq, vax_date)
rm(vax)
##  Comorbidity Wrangling ----
# Change diagnosis flags to integers and count the number of comorbidities per person
cc_all <- cc_all %>%
  mutate_if(is.factor, function(x)
    x %>% as.character() %>% as.integer()) %>%
 # count the number of comorbidities each person has
  mutate(como_count = rowSums(across(ihd_any:esoph_stomach_duod)))

# Match comorbidity data into to case control
cc_data <- cc_data %>%
  left_join(cc_all)
rm(cc_all)

##  SICSAG Wrangling ----
## Drop dates prior to Feb 2020
cc_data <- cc_data %>%
  left_join(sicsag) %>%
  # Replace missing data with dummy values
  replace_na(list(icu_hdu_flag = 0L))

##  Teachers Wrangling ----
teachers <- teachers %>% 
  filter(!anon_id == "ANONNANA")
teachers_xtra <- read_csv("//TOPLEVEL/Case_control/data/2021-07-28/CC_TEACHER_2021-07-28.csv") %>% 
  filter(!anon_id == "ANONNANA")
teachers_xtra <- teachers_xtra %>% 
  select(-serial) %>% 
  distinct()
xtra_dups <- teachers_xtra$anon_id[duplicated(teachers_xtra$anon_id)] %>% unique()
## One duplicate with one row as primary and ones as secondary; set to "Primary & Secondary"
teachers_xtra %>% 
  filter(anon_id == xtra_dups)
teachers_xtra <- teachers_xtra %>% 
  mutate(sector = if_else(anon_id == xtra_dups & sector %in% c("Primary", "Secondary"),
                          "Primary & Secondary", sector)) %>% 
  distinct()
teachers <- teachers %>% 
  inner_join(teachers_xtra)
rm(teachers_xtra, xtra_dups)

teachers <- teachers %>%
  mutate(
    sector_grp = case_when(
      sector %in% c("Independent Secondary",
                    "Secondary") ~ "Secondary",
      sector %in% c("Nursery & Primary", "Nursery") ~
        "Nursery/Primary or Nursery",
      sector %in% c("Primary") ~
        "Primary",
      TRUE ~ "Other"
    ),
    teacher_flag = 1
  ) %>%
  rename(teaching_status = status)

saveRDS(
  teachers %>%
    select(sector, teaching_status, sex, age = age_years, sector_grp, simd = simd2020_sc_quintile),
  "Scratch_data/cleaned_teachers.Rds"
)
teachers <- teachers %>%
  select(-sex,-age_years, -simd2020_sc_quintile)

## Take teacher role and status within in household; preferentially take teacher "throughout"
tchs_hhd_inh <- teachers %>% 
  select(hid, teacher_in_hh = teaching_status, tch_hhd_sectors = sector_grp) %>% 
  arrange(desc(teacher_in_hh)) %>% 
  distinct(hid, .keep_all = TRUE)

tchs_hhd_uprn <- teachers %>% 
  select(hid_uprn, teacher_in_hh = teaching_status, tch_hhd_sectors = sector_grp) %>% 
  arrange(desc(teacher_in_hh)) %>% 
  distinct(hid_uprn, .keep_all = TRUE)

# Drop unneeded teacher vars
teachers <- teachers %>%
  select(anon_id, sector, sector_grp, teaching_status, teacher_flag)

# Link teachers to case_control; 
cc_data <- cc_data  %>%
  select(-teacher_flag) %>% 
  left_join(teachers) %>%
  replace_na(list(
    teacher_flag = 0,
    sector_grp = "Not a teacher",
    teaching_status = "Not a teacher"
  ))
  
##  HCW Wrangling ----

# Join HCW data to case control
## Treat HCW as non-TV as is not main exposure and is more straightforward to do so, treat teacher as time-varying
cc_data <- cc_data %>% 
  mutate(hcw_hhd = case_when(
    anon_id %in% hcw$anon_id ~ "hcw",
    anon_id %in% pf_hhd$anon_id ~ "hhd",
    anon_id %in% npf$anon_id ~ "npf",
    TRUE ~ "not a healthcare worker"
  ))

##  Final case control pull together ----
cc_data <- cc_data %>% 
  mutate(hh18_cat = case_when(
    is.na(hh18) ~ "1",
    hh18 >=3 ~ "3+",
    TRUE ~ as.character(hh18)))

##  define cases, time etc ----
cc_data <- cc_data %>%
  rename(is_case_orig = is_case) %>% 
  mutate(is_case = if_else(is_case_orig == T | diag_case == 1 | cod_case == 1, 1L, 0L),
         is_case = if_else(is_case %>% is.na(), 0L, is_case))

## rename variables to those used with previous data pipeline
cc_data <- cc_data %>%
  rename(age_year = age_years,
         carehome = care_home,
         nursinghome = nursing_home) 

cc_data <- cc_data %>% 
  # Flag people resident in care/nursing homes
  mutate(res_care = if_else(carehome == 1 |
                              nursinghome == 1, 1, 0)) %>%
  # Create flags for events
  mutate(
    tested_pos = if_else(is_case == TRUE, 1, 0),
    hosp = if_else((adm28 == 1 |
                      inhosp == 1 
                      ) & is_case == 1,
                   1, 0),
    severe = if_else( (icu_hdu_flag == 1 | covid_cod == 1 |
                       covid_ucod == 1 | dead28 == 1) & is_case == 1, 1, 0),
    covid_death = if_else(covid_cod == 1 | covid_ucod == 1 |
                            (dead28 == 1 & is_case ==1), 1, 0)
  ) %>%
  
  # Assign ethnicity
  mutate(
    white = if_else(grepl("^1[A-Z]$", ethnic_smr_last) == T, 1, 0),
    s_asian = if_else(grepl("^3[ABCFGH]$", ethnic_smr_last) == T, 1, 0),
    chinese = if_else(grepl("^3[EJ]$", ethnic_smr_last) == T, 1, 0),
    black = if_else(
      grepl("^4[ABCDEY]$", ethnic_smr_last) == T |
        grepl("^5[ABCDY]$", ethnic_smr_last) == T,
      1,
      0
    ),
    other = if_else(
      grepl("^2[A-Z]$", ethnic_smr_last) == T |
        grepl("^3[DZ]$", ethnic_smr_last) == T |
        grepl("^6[AZ]$", ethnic_smr_last) == T,
      1,
      0
    )
  ) %>%
  mutate(unknown = if_else(
    grepl("^9", ethnic_smr_last) | is.na(ethnic_smr_last) |
      (white == 0 &
         black == 0 & s_asian == 0 & chinese == 0 &
         other == 0),
    1,
    0
  )) %>%
  
  # Assign age groups
  mutate(
    age_grp = case_when(
      age_year >= 21 & 
        age_year <= 30 ~ "21 - 30",
      age_year >= 31 &
        age_year <= 40 ~ "31 - 40",
      age_year >= 41 &
        age_year <= 50 ~ "41 - 50",
      age_year >= 51 &
        age_year <= 65 ~ "51 - 65",
      TRUE ~ "Unknown"
    )
  ) %>%
  
  # drop unnecessary variables
  select(-ethnic_smr_last,-carehome,-nursinghome) 

cc_data <- cc_data %>%
  mutate(
    hcw_teacher = case_when(
      hcw_hhd == "hcw" ~ "Healthcare Worker, patient facing",
      teacher_flag == 1 ~ "Teacher",
      TRUE ~ "Neither"
    )
  ) %>%
  mutate(week_of_death = if_else(
    !is.na(Week2020(date_of_death)) &
      covid_cod != 1,
    Week2020(date_of_death),
    99999
  )) %>%
  replace_na(list(
    week_of_death = 99999,
    ethnic_smr_last = "Unknown",
    simd = "Unknown"
  ))

cc_data <- cc_data %>%
  rename(shielding = shielding_flag)

## Add in teacher household identifiers ----
## remove duplicates in teachers household. Are exact duplicates
tchs_hhd_inh <- tchs_hhd_inh %>% 
  mutate(hhd_tch = 1L)
tchs_hhd_uprn <- tchs_hhd_uprn %>% 
  mutate(hhd_tch = 1L)

if(uprn_only == TRUE){
  cc_data <- cc_data %>% 
    left_join(tchs_hhd_uprn %>%  rename(hid = hid_uprn)) %>%
    mutate(teacher_in_hh = if_else(teacher_in_hh %>% is.na(), "Not a teacher", tch_hhd_sectors))
} else {
  ## Note need different for early and late period for UPRN versus in house
  cc_data_rly <- cc_data %>% 
    filter(specimendate <= as.Date("2020-12-31"))
  cc_data_lte <- cc_data %>% 
    filter(specimendate  > as.Date("2020-12-31"))
  print(nrow(cc_data))
  rm(cc_data) 
  
  cc_data_rly <- cc_data_rly %>%
    left_join(tchs_hhd_uprn %>% rename(hid = hid_uprn)) %>%
    mutate(teacher_in_hh = if_else(teacher_in_hh %>% is.na(), "Not a teacher", tch_hhd_sectors))
  cc_data_lte <- cc_data_lte %>%
    left_join(tchs_hhd_inh %>% mutate(hid = as.character(hid))) %>%
    mutate(teacher_in_hh = if_else(teacher_in_hh %>% is.na(), "Not a teacher", tch_hhd_sectors))
  cc_data <- bind_rows(cc_data_rly,
                       cc_data_lte)
  print(nrow(cc_data))
  rm(cc_data_rly,
     cc_data_lte)
}

cc_data <- cc_data %>%
  mutate(hhd_tch = if_else(is.na(hhd_tch), 0L, hhd_tch))

## Create Primary, Secondary, other variable for household members of teachers
cc_data <- cc_data %>% 
  mutate(teacher_in_hh = str_to_lower(teacher_in_hh),
         hh_tch_det = case_when(
           str_detect(teacher_in_hh, "primary") ~ "Primary",
           str_detect(teacher_in_hh, "secondary") ~ "Secondary",
           hhd_tch ==1 ~ "Other",
           TRUE ~ "Not a teacher"
         ))

## Set teaching status according to joiner/leaver, ----
## set teacher status based on leaving and joining, use tch for time varying and sector_grp_tv for tv
## use the hcw_teacher variable for baseline characteristics
## hhd_tch is also a time varying variable
## create a mutually exclusive variable for teacher, hcw, hhd members and neither for plots
cc_data <- cc_data %>%
  mutate(
    tch = case_when(
      specimendate >= as.Date("2020-08-01") &
        teaching_status == "leaver" ~ 0L,
      specimendate <  as.Date("2020-08-01") &
        teaching_status == "joiner" ~ 0L,
      TRUE ~ as.integer(teacher_flag)
    ),
    hhd_tch = case_when(
      specimendate >= as.Date("2020-08-01") &
        teaching_status == "leaver" ~ 0L,
      specimendate <  as.Date("2020-08-01") &
        teaching_status == "joiner" ~ 0L,
      TRUE ~ as.integer(hhd_tch)),
    hh_tch_det = if_else(hhd_tch ==0, "Not a teacher", hh_tch_det),
    sector_grp_tv = if_else(tch == 0, "Not a teacher", sector_grp),
    hhd = if_else(hcw_hhd == "hhd", 1L, 0L),
    hcw = if_else(hcw_hhd == "hcw", 1L, 0L),
    grps = case_when(
      hcw ==1 ~ "Healthcare worker",
      hhd ==1 ~ "Household member of HCW",
      tch ==1 ~ "Teacher",
      hhd_tch==1 ~ "Houshold member of Tch.",
      TRUE ~ "General population"))


## Calculate closed/open periods for schools ----
## Note already includes 5 day lag.
# A.	1st March 2020 to 24th August 2020 – Spring and Summer 2020 Closed period
# B.	25th August 2020 to 23rd Dec 2020– Winter 2020 Open period
# C.	24th Dec 2020 to 26th February 2021 – Winter 2021 Closed period
# D.	27th Feb 2021 to 23rd April 2021 – Winter and Spring 2021 Partially open period
# E.	24th April 2021 to 30th June 2021 – Spring 2021 Open period

cc_data <- cc_data %>%
  filter(specimendate <= as.Date("2021-06-30")) %>% 
  mutate(
    closure_status = case_when(
      specimendate <= as.Date("2020-08-24") ~  "a_pre_reop",
      specimendate >= as.Date("2020-08-25") & specimendate <= as.Date("2020-12-23") ~  "b_reop2020",
      specimendate >= as.Date("2020-12-24") & specimendate <= as.Date("2021-02-26") ~ "c_closed",
      specimendate >= as.Date("2021-02-27") & specimendate <= as.Date("2021-04-23") ~ "d_mixed",
      specimendate >= as.Date("2021-04-24") & specimendate <= as.Date("2021-06-30") ~ "e_reop2021"))
  
cc_data <- cc_data %>%
  mutate(
    closure_status_clps = case_when(
      closure_status %in% c("a_pre_reop", "c_closed")   ~ "a_closed",
      closure_status %in% c("d_mixed")                  ~ "b_ignore",
      closure_status %in% c("b_reop2020", "e_reop2021") ~ "c_open"))

## Calculate test week
cc_data <- cc_data %>% 
  mutate(test_week = Week2020(specimendate),
         week_start = StartWeek(test_week))

## Drop small number of individuals who are assigned as both a teacher and a hcw ----
cc_data <- cc_data %>% 
  filter(! (hcw_teacher == "Healthcare Worker, patient facing" & teacher_flag ==1))

## By definition, only allow individuals to be in household members of teachers/hcws where they are not themselves teachers or household members of teachers
cc_data <- cc_data %>% 
  mutate(hhd_tch = if_else(tch ==1 |teacher_flag ==1, 0L, hhd_tch),
         hh_tch_det = if_else(tch ==1 |teacher_flag ==1, "Not a teacher", hh_tch_det),
         hhd = if_else(hcw ==1, 0L, hhd))

## Add in non-patient facing healthcare workers as per response to editorial comments
cc_data <- cc_data %>% 
  mutate(npf = if_else(anon_id %in% npf$anon_id & !(hcw ==1L|hhd==1L) , 1L, 0L))

cc_data <- cc_data %>% 
  mutate(hcw_teacher = if_else(npf ==1, "Healthcare worker, non-patient facing", hcw_teacher))

## Add in vaccination status as time-varying covariate ----
cc_data <- cc_data %>% 
  left_join(vax_dt)
cc_data <- cc_data %>% 
  mutate(vt1_dt = vt1,
         vt2_dt = vt2,
         vt1 = if_else(!is.na(vt1), as.integer(specimendate-vt1), 0L),
         vt2 = if_else(!is.na(vt2), as.integer(specimendate-vt2), 0L),
         vt = case_when(
           vt2 >= 14 ~ "c_post2nd",
           vt1 >= 14   ~ "b_post1st",
           TRUE ~ "a_unvacc"
         ))

## Drop undetermined and uncertain npf healthcare workers
cc_data <- cc_data %>% 
  mutate(uncertain_hcw = anon_id %in% hcw_drop$anon_id) %>% 
  group_by(hid) %>% 
  mutate(uncertain_hcw = any(uncertain_hcw)) %>% 
  ungroup()

saveRDS(cc_data %>% filter(uncertain_hcw ==TRUE), "uncertain_hcw_exclude.Rds")

cc_data <- cc_data %>% 
  filter(!uncertain_hcw) %>% 
  select(-uncertain_hcw)

## Select relevant variables for subsequent analyses
cc_data <- cc_data %>%
  select(
    anon_id,
    stratum,
    is_case,
    is_case_orig,
    sex,
    age_year,
    simd = simd2020_sc_quintile,
    hb2019name,
    dm_type,
    diab_reg,
    ihd_any,
    heart_other_any,
    circulatory_other,
    ckd_any,
    oad_any,
    neuro_any,
    liver_any,
    immune_any,
    neoplasm_any,
    neoplasm_last_year,
    esoph_stomach_duod,
    como_count,
    sector,
    sector_grp,
    sector_grp_tv,
    teaching_status,
    teacher_flag,
    hcw_teacher,
    hcw_hhd,
    hosp,
    severe,
    white,
    s_asian,
    chinese,
    black,
    other,
    unknown,
    age_grp,
    shielding,
    hhd_tch,
    hh_tch_det,
    tch,
    hhd,
    hcw,
    closure_status,
    closure_status_clps,
    simd2020_sc_quintile,
    specimendate,
    test_week,
    grps,
    hh04, hh1217, hh511, hh18,
    hid,
    hh18_cat,
    npf,
    vt1, 
    vt2,
    vt,
    vt1_dt,
    vt2_dt
  )

## Tidy env ----
rm(hcw,
   hcw_drop,
   npf,
   pf_hhd,
   runanyway,
   sicsag,
   tchs_hhd_inh,
   tchs_hhd_uprn,
   teachers,
   vax_dt)

  if (uprn_only == TRUE) saveRDS(cc_data, "Data/basefile_uprn_only.Rds") else{
    saveRDS(cc_data, "Data/basefile.Rds")
  }
 }

