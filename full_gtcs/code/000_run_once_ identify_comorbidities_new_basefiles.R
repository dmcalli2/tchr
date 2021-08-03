### 1 - Housekeeping ----
## Load Packages
library(magrittr)
library(dplyr)
library(here)
library(readr)
library(janitor)
library(stringr)
library(lubridate)
library(tidyr)
library(data.table)

# Case Control
## limit data to working-age
cc_all <- as.data.table(readRDS(
  "TOPLEVEL/Case_control/data/2021-07-28/CC_linked_ANON_2021-07-28.rds")%>% 
    filter(age_years >=21, age_years<=65) %>% 
    mutate(anon_id_sp_dt = paste0(anon_id, "_",specimen_date)), 
  key="anon_id_sp_dt") %>% 
  clean_names() %>% 
  rename(specimendate = specimen_date) 

## Restrict datasets to time prior to first specimen date until time of last specimendate
## note there are 4 duplicates anon ids across 8 duplicate rows, review and drop; all controls, all with missing specimendate
dups <- cc_all$anon_id_sp_dt[duplicated(cc_all$anon_id_sp_dt)]
# cc_all %>% 
#   filter(anon_id_sp_dt %in% dups) %>% View()
cc_all <- cc_all %>% 
  filter(!anon_id_sp_dt %in% dups)

## Read in data
# SMR01 - For ICD10 codes
diagnoses <- readRDS("TOPLEVEL/Case_control/data/2021-07-28/CC_SMR01_ICD10_2021-07-28.rds") %>% 
    clean_names() 
diagnoses <- diagnoses %>% 
  ## Limit data to 
  ## FROM 9 months prior to specimendate and 
  ## TO 25 days before speciemendate
  ## Note there are multiple anon_ids, but anon_id and specimendate together is unique
  ## create new variable anon_id_sp_dt based on anon_id and specimendate 
  inner_join(cc_all %>% distinct(anon_id, specimendate)) %>% 
  filter(discharge_date <   (specimendate - days(25)),
        discharge_date >=  (specimendate - years(5))) %>% 
  mutate(anon_id_sp_dt = paste0(anon_id, "_",specimendate))
diagnoses <- diagnoses %>% 
  as.data.table(key = "anon_id_sp_dt")

# SMR01 - For OPCS4 codes
procedures <- readRDS(
  "TOPLEVEL/Case_control/data/2021-07-28/CC_SMR01_OPCS4_MAIN2021-07-28.rds") %>%
  clean_names()%>% 
  semi_join(cc_all %>% select(anon_id))

procedures <- procedures %>% 
  ## Limit data to 
  ## FROM 9 months prior to specimendate and 
  ## TO 25 days before speciemendate
  ## Note there are multiple anon_ids, but anon_id and specimendate together is unique
  inner_join(cc_all %>% distinct(anon_id, specimendate)) %>% 
  filter(discharge_date <   (specimendate - days(25)),
         discharge_date >=  (specimendate - years(5)))%>% 
  mutate(anon_id_sp_dt = paste0(anon_id, "_",specimendate))
procedures <- procedures %>% 
  as.data.table(key = "anon_id_sp_dt")

# PIS - For Drugs
pis <- readRDS(
  "TOPLEVEL/Case_control/data/2021-07-28/CC_PIS_ANON_2021-07-28.rds") %>% 
  semi_join(cc_all %>% select(anon_id))
pis <- pis %>% 
  select(anon_id, dispensed_date, bnf_paragraph_code) 

## Limit data to 
## FROM 9 months prior to specimendate and 
## TO 15 days before speciemendate
## Note there are multiple anon_ids, but anon_id and specimendate together is unique
pis <- pis %>% 
  inner_join(cc_all %>% distinct(anon_id, specimendate)) %>% 
  filter(dispensed_date <   (specimendate - days(15)),
         dispensed_date >=  (specimendate - months(9)))%>% 
  mutate(anon_id_sp_dt = paste0(anon_id, "_",specimendate))
pis$specimendate <- NULL
pis <- pis %>% 
  as.data.table(key = "anon_id_sp_dt")
pis <- pis[ , c("anon_id_sp_dt", "bnf_paragraph_code")]
pis$bnf_section_code <- substr(pis$bnf_paragraph_code, 1, 4)


### Ischaemic Heart Disease ----
## Identify those prescribed drugs for this condition
# nitrates are BNF code 020601
ids.bnf.IHD <- unique(pis$anon_id_sp_dt[substr(as.character(pis$bnf_paragraph_code), 1, 6) == "020601" |
                                       substr(as.character(pis$bnf_paragraph_code), 1, 6) == "020603"])

## Identify those who have an ICD10 code for this condition
ids.icd.IHD <- unique(diagnoses$anon_id_sp_dt[grep("^I2[0-5]", diagnoses$icd10)])

## Identify those who have an OPCS4 (pocedure) code for CABG and PTCA
ids.procedures.IHD <- unique(procedures$anon_id_sp_dt[grep("^K4[012349]|^K50", procedures$main_operation)])

## Join the anon ids from above together and remove duplicates
ids.IHD <- unique(c(ids.icd.IHD, ids.bnf.IHD, ids.procedures.IHD))

## Create a flag for ischameic heart disease in case control data
cc_all$ihd_any <- as.factor(as.integer(cc_all$anon_id_sp_dt %in% ids.IHD))


### Other Heart Disease ----
## Identify those prescribed drugs for this condition
# anti-arrhythmics are BNF code 0203
ids.bnf.heart.other <-
  unique(pis$anon_id_sp_dt[substr(as.character(pis$bnf_paragraph_code), 1, 4) == "0203"])  

## Identify those who have an ICD10 code for this condition
# heart disease is I05 to I52
ids.icd.heart.other <- unique(diagnoses$anon_id_sp_dt[grep("^I0[01256789]|^I1[0-5]|^I2[6-8]|^I3[0-9]|^I4[0-9]|^I5[0-2]",
                                                     diagnoses$icd10)])

## Identify those who have had a relevant procedure
ids.procedures.heart.other <- unique(procedures$anon_id_sp_dt[grep("^K57", procedures$main_operation)])

## Join the anon ids from above together and remove duplicates
ids.heart.other <- unique(c(ids.icd.heart.other, ids.bnf.heart.other, ids.procedures.heart.other))

## Create a flag for other heart disease in case control data
cc_all$heart_other_any <- as.factor(as.integer(cc_all$anon_id_sp_dt %in% ids.heart.other))


### Other Circulatory Disease ----
## Identify those who have an ICD 10 code for this condition
# I60 to I99 
ids.icd.circulatory.other <-
  unique(diagnoses$anon_id_sp_dt[grep("^I[6-9]|^Z95", diagnoses$icd10)])

## Create a flag for other circulatory disease in case control data
cc_all$circulatory_other <-
  as.factor(as.integer(cc_all$anon_id_sp_dt %in% ids.icd.circulatory.other))


### Chronic Kidney Disease ----
## Identify those who have an ICD 10 code for this condition
# this includes CKD stage 4
ids.icd.ckd <- unique(diagnoses$anon_id_sp_dt[grep("^N18[45]|^Z49[0-2]|^Z94[02]",
                                             diagnoses$icd10)])

## Identify those who have had a relevant procedure
ids.kidneytransplant <- unique(procedures$anon_id_sp_dt[grep("^M01[1234589]",
                                                       procedures$main_operation)])

## Join the anon ids from above together and remove duplicates
ids.ckd.any <- unique(c(ids.icd.ckd, ids.kidneytransplant))

## Create a flag for chronic kidney disease in case control data
cc_all$ckd_any <-  as.factor(as.integer(cc_all$anon_id_sp_dt %in% ids.ckd.any))


### Asthma and chronic lower respiratory disease ----
## Identify those who have an ICD 10 code for this condition
ids.icd.asthma <- unique(diagnoses$anon_id_sp_dt[grep("^J4[56]", diagnoses$icd10)])
ids.icd.chronresp <- unique(diagnoses$anon_id_sp_dt[grep("^J4[012347]|^J6[0-9]|^J70|^J8[0-6]|^J9[0-9]|^G47\\.?3",
                                                   diagnoses$icd10)])

## Identify those prescribed drugs for this condition
ids.bnf.broncho <- unique(pis$anon_id_sp_dt[as.integer(pis$bnf_section_code) >= 301 &
                                           as.integer(pis$bnf_section_code) <= 303])

## Join the anon ids from above together and remove duplicates
ids.oad.any <- unique(c(ids.icd.asthma, ids.icd.chronresp, ids.bnf.broncho))

## Create a flag for asthma/other lower respiratory disease in case control data
cc_all$oad_any <- as.factor(as.integer(cc_all$anon_id_sp_dt %in% ids.oad.any))


### Neurological disorders ---- 
## Identify those who have an ICD 10 code for this condition
# Include all Nervous chapter except G40 "Episodic and Paroxysmal Disorders"
# Exclude G0 meningitis and encephalitis, and G5 local neuropathies
# also include F03 dementia NOS
ids.icd.neuro <- unique(diagnoses$anon_id_sp_dt[grep("^F03|^G[1236789]", diagnoses$icd10)])

## Identify those prescribed drugs for this condition
ids.bnf.neuro <- unique(pis$anon_id_sp_dt[as.integer(pis$bnf_section_code) == 409 |
                                         as.integer(pis$bnf_section_code) == 411 |
                                           as.integer(substr(pis$bnf_paragraph_code, 2, 5)) == 8020])

## Join the anon ids from above together and remove duplicates
ids.neuro.any <- unique(c(ids.icd.neuro, ids.bnf.neuro))

## Create a flag for asthma/other lower respiratory disease in case control data
cc_all$neuro_any <- as.factor(as.integer(cc_all$anon_id_sp_dt %in% ids.neuro.any))

### Liver disease ----
## Identify those who have an ICD 10 code for this condition
liver.grep.string <- "^C22\\.?0|^I85\\.?0|^I98\\.?3|^K70\\.?[234|^K71\\.?7|^K72\\.?[019]|^K72\\.?[019|^K73|^K74\\.?[023456]|^K76\\.?7|^R18"
ids.icd.liver <- unique(diagnoses$anon_id_sp_dt[grep(liver.grep.string, diagnoses$icd10)])

## Create a flag for liver disease in case control data
cc_all$liver_any <- as.factor(as.integer(cc_all$anon_id_sp_dt %in% ids.icd.liver))


### Immunodeficiency and immunosuppression ----
## Identify those who have an ICD 10 code for this condition
ids.icd.immune <- unique(diagnoses$anon_id_sp_dt[grep("^B2[0-3|^D8[0-9]", diagnoses$icd10)])

## Identify those prescribed drugs for this condition
# 802 other immunomodulating drugs
# Methotrexate and chloroquine appear in musculoskeletal chapter 
ids.bnf.immune <- unique(pis$anon_id_sp_dt[as.integer(pis$bnf_section_code) == 802 |
                                          as.integer(pis$bnf_paragraph_code) == 50301])

## Join the anon ids from above together and remove duplicates
# immune.any includes primary immunodeficiency and secondary immunosuppression
ids.immune.any <- unique(c(ids.icd.immune, ids.bnf.immune))

## Create a flag for immunodeficeincy/immunosuppression in case control data
cc_all$immune_any <- as.factor(as.integer(cc_all$anon_id_sp_dt %in% ids.immune.any))

### Neoplasms ---- 
## Identify those who have an ICD 10 code for this condition
ids.icd.neoplasm <- unique(diagnoses[grepl("^C[0-9]|^D[0-4]", icd10), anon_id_sp_dt])

ids.icd.neoplasm.lastyear <-
  unique(diagnoses[as.integer(specimendate - as.Date(discharge_date)) + 25 > 365 &
                     grepl("^C[0-9]|^D[0-4]", icd10), anon_id_sp_dt])

## Identify those prescribed drugs for this condition
ids.bnf.neoplasm <- unique(pis$anon_id_sp_dt[as.integer(pis$bnf_section_code) == 801])

## Join the anon ids from above together and remove duplicates
ids.neoplasm.any <- unique(c(ids.icd.neoplasm, ids.bnf.neoplasm))

ids.neoplasm.lastyear <- unique(c(ids.icd.neoplasm.lastyear, ids.bnf.neoplasm))

## Create a flag for neoplasms in case control data
cc_all[, neoplasm_any := as.factor(as.integer(anon_id_sp_dt %in% ids.neoplasm.any | 
                                                can_reg==1))]

cc_all[, neoplasm_last_year := as.factor(as.integer(anon_id_sp_dt %in% ids.neoplasm.lastyear))]


### Disorders of esophagus, stomach and duodenum ----
## Identify those who have an ICD 10 code for this condition
ids.esoph.stomach.duod <-  unique(diagnoses$anon_id_sp_dt[grep("^K2[0-9]|^K3[01]", diagnoses$icd10)])

## Create a flag for immunodeficiency/immunosuppression in case control data
cc_all$esoph_stomach_duod <-
  as.factor(as.integer(cc_all$anon_id_sp_dt %in% ids.esoph.stomach.duod))

### Limit dataset for matching ----
cc_all <- cc_all %>% 
  as_tibble() %>% 
  select(anon_id, specimendate, ihd_any:esoph_stomach_duod)
saveRDS(cc_all, "Scratch_data/comorbidity_count_jul2021.Rds")



