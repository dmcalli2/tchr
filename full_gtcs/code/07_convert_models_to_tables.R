#07_convert_models_to_tables

library(survival)
library(tidyverse)
library(hrbrthemes)
## Script converts models to tables and produces these for appendix and main manuscript. Also calculates counts of cases and controls for these tables


# Flip between output/main_analysis and other choices both in reading and saving

## read files and convert to model summaries ----
map(list("main_analysis",
         "main_analysis_uprn_hid",
         "no1st",
         "no2nd"), function(folder_choose){
       print(folder_choose)

if(folder_choose == "main_analysis_uprn_hid") uprn_only <- FALSE
       source("code/full_gtcs/01_create_basefile.R")
       
       
if(folder_choose == "no1st") cc_data <- cc_data %>% filter(vt == "a_unvacc")
if(folder_choose == "no2nd") cc_data <- cc_data %>% filter(!vt == "b_post1st")

teach <- readRDS(paste0("output/", folder_choose, "/", "teach.Rds"))
sectr <- readRDS(paste0("output/", folder_choose, "/", "sector.Rds"))
sever <- readRDS(paste0("output/", folder_choose, "/", "severe.Rds"))


teach <- bind_rows(teach = teach, 
                   sectr = sectr, 
                   .id = "forms")
teach <- bind_rows(teach, sever %>% mutate(forms = "sever"))
teach <- teach %>% 
  gather("model", "res", unad_hosp, adj_hosp, unad_case, adj_case, unad_case_hh, adj_case_hh, unad_hosp_hh, adj_hosp_hh, unad_sevr, adj_sevr)
teach <- teach %>% 
  mutate(model_n = 1:nrow(teach)) 

teach$res_valid <- map_lgl(teach$res, ~ ("clogit" == class(.x)[1]))
fail <- teach %>% 
  filter(!res_valid)
teach <- teach %>% 
  filter(res_valid)
ExtractEst2 <- function(model){
  # Take estimate and se
  a <- summary(model)$coefficients
  a <- a[, c(1, 3)] 
  colnames(a) <- c("b", "b_se")
  a %>% as_tibble()
}
teach$est <- map(teach$res, ExtractEst)
teach$est <- map(teach$est, CleanEst)
teach$raw <- map(teach$res, ExtractEst2)
teach_res <- teach %>% 
  select(-res, -res_valid) %>% 
  unnest(c(est, raw))

teach_res <- teach_res %>% 
  mutate_at(vars(est, lci, uci), function(x)  formatC(round(x, 2), digits = 2, format = "f", flag = "0")) %>% 
  mutate(res = paste0(est, " (", lci, "-", uci, ")"),
         res = if_else(str_detect(res, "Inf|NA"), "-", res)) %>% 
  select(-est, -lci, -uci)

## Identify model exposures 
teach_res <- teach_res %>% 
  group_by(model_n) %>% 
  mutate(Individual = if_else(any("tch" %in% params), "Teachers overall", "By sector"),
         Household  = if_else(any("hhd_tch" %in% params), "Teachers overall", "By sector")) %>% 
  ungroup()

## Create table for supplementary appendix
params_lkp <- c("tch" = "Teacher", 
                "hcw" = "Healthcare worker, patient facing",
                "npf" = "Healthcare worker, non-patient facing",
                "hhd" = "Household member of healthcare worker",
                "hhd_tch" = "Household member of teacher",
  "sector_grp_tvNursery/Primary or Nursery" = "Nursery/Primary or Nursery", 
  "sector_grp_tvOther" = "Teacher in other sector",
  "sector_grp_tvPrimary" = "Primary",
  "sector_grp_tvSecondary" = "Secondary", 
  "hh_tch_detOther" = "Household member other teacher", 
  "hh_tch_detPrimary" = "Household member of primary teacher", 
  "hh_tch_detSecondary" = "Household member of secondary teacher"
)

teach_res <- teach_res %>%
  filter(params %in% names(params_lkp))

s1 <- teach_res %>% 
  mutate(
    # exposure = case_when(
    # exposure == "sectr" ~ "By sector",
    # str_detect(model, "_hh$") ~ "By sector household",
    # TRUE ~ "Single category"),
         mdls = case_when(
           mdls == "univ" ~ "Full pandemic",
           mdls == "closure" ~ "By period",
           mdls == "age_sex" ~ "By age and sex",
           mdls == "age_sex_closure" ~ "By age, sex and period"),
         closure_status = if_else(closure_status =="allof", "Full period", closure_status),
         sex = case_when(sex ==1 ~ "Men",
                         sex ==2 ~ "Women",
                         TRUE ~ "Men and women"),
         age_grp = if_else(age_grp == "allof", " All ages", age_grp),
         Outcome = case_when(
           str_detect(model, "case") ~ "Any case", 
           str_detect(model, "hosp") ~ "Hospitalisation",
           str_detect(model, "sevr") ~ "Severe"),
         model = if_else(str_detect(model, "^unad"), "Unadjusted", "Adjusted") ,
         params = params_lkp[params]
         ) %>% 
  rename(Stratification = mdls,
         "Time period" = closure_status,
         Sex = sex,
         "Age group (years)" = age_grp,
         Adjustment = model,
         Parameter = params) %>% 
  select(Individual, Household, Outcome, Adjustment, everything()) %>% 
  # select(-allpt) %>% 
  distinct() %>% 
  mutate(res = if_else(res %>% is.na(), "-", res)) %>% 
  select(-model_n, -forms) %>% 
  distinct()

## Calcualte differences across time periods
cmpr_b <- s1 %>% 
  select(-res) %>% 
  filter(!`Time period` == "Full period") %>% 
  rename(est = b, se = b_se) %>% 
  pivot_wider(names_from = `Time period`, values_from = c(est, se))

cmpr_c <- cmpr_b %>% 
  mutate(est_b_reop2020 = est_b_reop2020 - est_a_pre_reop,
         est_e_reop2021 = est_e_reop2021 - est_a_pre_reop,
         se_b_reop2020  = (se_b_reop2020^2 + se_a_pre_reop^2)^0.5,
         se_e_reop2021  = (se_e_reop2021^2 + se_a_pre_reop^2)^0.5,
         est_reop2020 = exp(est_b_reop2020),
         est_reop2021 = exp(est_e_reop2021),
         lci_reop2020 = exp(est_b_reop2020 - 1.96*se_b_reop2020),
         uci_reop2020 = exp(est_b_reop2020 + 1.96*se_b_reop2020),
         lci_reop2021 = exp(est_e_reop2021 - 1.96*se_e_reop2021),
         uci_reop2021 = exp(est_e_reop2021 + 1.96*se_e_reop2021)
         ) %>% 
  mutate_at(vars(est_reop2020, est_reop2021, 
                 lci_reop2020, uci_reop2020,
                 lci_reop2021, uci_reop2021), ~ formatC(round(.x,2), digits = 2, flag = "0", format = "f")) %>% 
  mutate(reop2020 = paste0(est_reop2020, " (", lci_reop2020, "-", uci_reop2020, ")"),
         reop2021 = paste0(est_reop2021, " (", lci_reop2021, "-", uci_reop2021, ")")) %>% 
  select(Individual:Parameter, reop2020, reop2021)

cmpr_d <- cmpr_c %>% 
  filter(Stratification == "By period", Individual == "Teachers overall", Household == "Teachers overall") %>% 
  select(Outcome, Adjustment, Parameter, reop2020, reop2021) %>% 
  pivot_longer(c(reop2020, reop2021), names_to = "Period", values_to = "res") %>% 
  pivot_wider(names_from = Parameter, values_from = res) %>% 
  mutate(Period = if_else(Period == "reop2020", "Re-opened 2020", "Re-opened 2021")) %>% 
  arrange(Outcome, desc(Period), desc(Adjustment)) %>% 
  select(`Outcome`, `Adjustment`, `Period`, 
         `Teacher`, 
         `Household member of teacher`, 
         `Healthcare worker, patient facing`, 
         `Household member of healthcare worker`,
         `Healthcare worker, non-patient facing`)
saveRDS(cmpr_d, paste0("output/", folder_choose, "/", "compare_across_periods.Rds"))

s1_over <- s1 %>% 
  select(-b, -b_se) %>% 
  filter(allpt %in% c("overall")) %>% 
  mutate(`Time period` = "Full pandemic") %>% 
  spread(Outcome, res, fill = "-") 
s1_rest  <- s1 %>% 
  mutate(res = if_else(str_length(res)>50, "-", res)) %>% 
  select(-b, -b_se) %>% 
  filter(`Time period` %in% c("a_pre_reop", "b_reop2020", "c_closed", "d_mixed", "e_reop2021")) %>% 
  spread(Outcome, res, fill = "-") 

saveRDS(list(s1 = s1 %>% select(-b, -b_se),
             s1_over = s1_over,
             s1_rest = s1_rest), paste0("output/", folder_choose, "/", "hazard_ratio_supplementary_appendix.Rds"))


## results Tables 1 and 2
t1_2 <- bind_rows(s1_over %>% 
                    filter(`Time period` == "Full pandemic", Individual == "Teachers overall", Household == "Teachers overall"),
                  s1_rest %>% 
                    filter(Stratification == "By period", Individual == "Teachers overall", Household == "Teachers overall")
                  ) %>% 
  select(Adjustment, `Time period`, Parameter, `Any case`, Hospitalisation, Severe)

t1_2_plot <- t1_2 %>% 
  gather("Outcome", "value", `Any case`, Hospitalisation, Severe) %>% 
  separate(value, into = c("est", "lci", "uci"), sep = "\\(|\\-") %>% 
  mutate_at(vars(est, lci, uci), ~ .x %>% str_remove_all("\\)") %>% parse_double) %>% 
  mutate(x = factor(Parameter,
                    levels = c("Teacher", "Household member of teacher",  
                               "Healthcare worker, patient facing", "Household member of healthcare worker",
                               "Healthcare worker, non-patient facing"),
                    labels = c("Teacher", "Teacher\n(household)",  
                               "Hlthcare\nPF", "Hlthcare\nPF\n(household)",
                               "Hlthcare\nNPF")),
         tp = factor(`Time period`, 
                     levels = c("a_pre_reop", "b_reop2020", "c_closed", "d_mixed", "e_reop2021","Full pandemic"),
                     labels = c("Spring/Summer 2020 (closed)",  "Autumn term 2020 (open)", "Winter 2020/21 (closed)", "Spring term 2021 (phased)", "Summer term 2021 (open)", "All time periods") ),
         outclr = factor(Outcome,
                         levels = c("Severe", "Hospitalisation", "Any case")))

## All outcomes on plot
hosp_sev <- ggplot(t1_2_plot %>% 
                     filter(Adjustment == "Adjusted") %>% 
                     mutate(myshape = if_else(tp == "All time periods", 19L, 0L),
                            mysize = if_else(tp == "All time periods", 1.5, 1.5)),
                   aes(x = x, y = est, ymin = lci, ymax = uci, colour = Outcome, shape = myshape)) +
  geom_linerange(position = position_dodge2(width = 0.2)) +
  geom_point(position = position_dodge2(width = 0.2), mapping = aes(size = mysize)) +
  facet_wrap(~ tp, nrow = 2) +
  scale_y_log10("Rate ratio (95% CI)") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "darkgrey") +
  scale_colour_manual("",  
                      breaks = c("Any case",
                                 "Hospitalisation",
                                 "Severe"),
                      values = c("#8fd175",
                                 "#d18975",
                                 "#3f2d54")) + #values = (ipsum_pal()(3)[c(1,2,3)])) +
  scale_shape_identity("") +
  scale_size_identity() +
  scale_x_discrete("") +
  theme_ipsum() +
  coord_cartesian(ylim = c(1/15, 5))
hosp_sev

saveRDS(hosp_sev, paste0("output/", folder_choose, "/", "plot_hrs_outcomes_grps.Rds"))

## same plot for supplement with CIs and not excluding healthcare workers
hosp_sev_supp <- ggplot(t1_2_plot,
                   aes(x = x, y = est, ymin = lci, ymax = uci, colour = Outcome, linetype = Adjustment, shape = Adjustment)) +
  geom_point(position = position_dodge2(width = 0.2)) +
  geom_linerange(position = position_dodge2(width = 0.2)) +
  facet_wrap(~ tp, ncol = 2) +
  scale_y_log10("Rate ratio (95% CI)") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "darkgrey") +
  scale_colour_manual("",  
                      breaks = c("Any case",
                                 "Hospitalisation",
                                 "Severe"),
                      values = c("#8fd175",
                                 "#d18975",
                                 "#3f2d54")) + #values = (ipsum_pal()(3)[c(1,2,3)])) +
  scale_shape("") +
  scale_x_discrete("") +
  theme_ipsum() +
  coord_cartesian(ylim = c(1/15, 5))
hosp_sev_supp
saveRDS(hosp_sev_supp, paste0("output/", folder_choose, "/", "plot_hrs_outcomes_grps_supplement.Rds"))


## Tables
# Table 2 manuscript, teachers and HCWs and households
t2 <- t1_2 %>% 
  gather("Outcome", "value", `Any case`, Hospitalisation, Severe) %>% 
  spread(Parameter, value) %>% 
  filter(`Time period` %in% c("a_pre_reop", "b_reop2020", "c_closed", "d_mixed", "e_reop2021","Full pandemic")) %>% 
  mutate(`Time period` = factor(`Time period`, 
                                levels = c("a_pre_reop", "b_reop2020", "c_closed", "d_mixed", "e_reop2021","Full pandemic"),
                                labels = c("Spring/Summer 2020 (closed)",  "Autumn term 2020 (open)", "Winter 2020/21 (closed)", "Spring term 2021 (phased)", "Summer term 2021 (open)", "All time periods") ))
t2 <- t2 %>% 
  mutate(Neither = "1")
t2 <- t2[ , c("Time period", "Outcome", "Adjustment", "Neither", "Teacher", "Household member of teacher",
              "Healthcare worker, patient facing", 
              "Household member of healthcare worker",
              "Healthcare worker, non-patient facing"
              )] %>% 
  arrange(as.integer(`Time period`), Outcome, desc(Adjustment))

t3 <- bind_rows(s1_over %>% 
                       filter(Individual == "By sector", Household == "By sector"),
                s1_rest %>% 
                       filter(Stratification == "By period", Individual == "By sector", Household == "By sector")
) %>% 
  select(Adjustment, `Time period`, Parameter, `Any case`, Hospitalisation) %>% 
  mutate(`Time period` = factor(`Time period`, 
                                levels = c("a_pre_reop", "b_reop2020", "c_closed", "d_mixed", "e_reop2021","Full pandemic"),
                                labels = c("Spring/Summer 2020 (closed)",  "Autumn term 2020 (open)", "Winter 2020/21 (closed)", "Spring term 2021 (phased)", "Summer term 2021 (open)", "All time periods") ))
t3 <- t3 %>% 
  filter(Parameter %in% c(  hh_tch_detPrimary = "Household member of primary teacher", 
                            hh_tch_detSecondary = "Household member of secondary teacher", 
                            hh_tch_detOther = "Household member other teacher", 
                            `sector_grp_tvNursery/Primary or Nursery` = "Nursery/Primary or Nursery", 
                            sector_grp_tvPrimary = "Primary", 
                            sector_grp_tvSecondary = "Secondary", 
                            sector_grp_tvOther = "Teacher in other sector")) %>% 
  gather("Outcome", "value", `Any case`, Hospitalisation) %>% 
  group_by(Adjustment, `Time period`, Parameter, Outcome) %>% 
  summarise(value = value %>% 
              paste(collapse = "|") %>% 
              str_remove_all(fixed("-|")) %>% 
              str_remove_all(fixed("|-"))) %>% 
  ungroup() %>% 
  spread(Parameter, value)
t3 <- t3[ ,c("Time period", "Outcome", "Adjustment",
  "Nursery/Primary or Nursery", "Primary", "Secondary", "Teacher in other sector",
  "Household member of primary teacher",  "Household member of secondary teacher", "Household member other teacher")] %>% 
  arrange(`Time period`, Outcome, desc(Adjustment))


## simplify counts ----
## count individuals and household members overall
## restrict to data analysed, ie exclude strata excluded in regression models
cc_data <- cc_data %>% 
  mutate(anyexp = (tch + hcw + hhd + hhd_tch + npf) >=1) %>% 
  group_by(stratum) %>% 
  mutate(any_tch_hcw = any(anyexp)) %>% 
  ungroup() %>% 
  filter(any_tch_hcw) %>% 
  select(-anyexp, -any_tch_hcw)

cc_data <- cc_data %>% 
  mutate(hcw_teacher2 = case_when(
    hcw_teacher == "Teacher" & tch ==0 ~ "Neither",
    hcw_teacher == "Teacher" ~ "Teacher",
    hcw_teacher == "Healthcare Worker, patient facing"     ~ "Healthcare worker, patient facing",
    hcw_teacher == "Healthcare worker, non-patient facing" ~ "Healthcare worker, non-patient facing",
    hhd_tch ==1 ~ "Household member of teacher",
    hhd ==1 ~ "Household member of healthcare worker",
    TRUE ~ "Neither"),
    tp = factor(closure_status,
                levels = c("a_pre_reop", "b_reop2020", "c_closed", "d_mixed", "e_reop2021"),
                labels = c("Spring/Summer 2020 (closed)",  "Autumn term 2020 (open)", "Winter 2020/21 (closed)", "Spring term 2021 (phased)", "Summer term 2021 (open)") ))
library(data.table)
cc_data <- data.table(cc_data)
setkey(cc_data, anon_id)

CaseMakerDT <- function(x) {
  x <- setorder(x, -chs, anon_id)
  x <- x[, .SD[1], by = .(anon_id, hcw_teacher2, tp)] 
  x <- x[, .(cases = sum(chs),
             controls = sum(1-chs)),
         by = .(hcw_teacher2, tp)]
  x %>% 
    as_tibble() %>% 
    arrange(hcw_teacher2, tp) %>% 
    mutate(res = paste(cases, controls, sep = "/")) %>% 
    select(-cases, -controls)
}

CaseMakerDTOverall <- function(x) {
  x <- setorder(x, -chs, anon_id)
  x <- x[, .SD[1], by = .(anon_id, hcw_teacher2)] 
  x <- x[, .(cases = sum(chs),
             controls = sum(1-chs)),
         by = .(hcw_teacher2)]
  x %>% 
    as_tibble() %>% 
    arrange(hcw_teacher2) %>% 
    mutate(tp = "Full pandemic",
           res = paste(cases, controls, sep = "/")) %>% 
    select(-cases, -controls)
}

CaseMaker <- function(x) {
  x %>% 
    arrange(desc(chs), anon_id) %>% 
    distinct(anon_id, .keep_all = TRUE) %>% 
    summarise(cases = sum(chs ==1),
              controls = sum(chs ==0)) %>% 
    ungroup() %>% 
    mutate(res = paste(cases, controls, sep = "/")) %>% 
    select(-cases, -controls)
}

cnts_vera <- cc_data %>% 
  rename(chs = is_case) %>% 
  CaseMakerDT() %>% 
  rename(Anycase = res)
cnts_verb <- cc_data %>% 
  rename(chs = hosp) %>% 
  CaseMakerDT() %>% 
  rename(Hospitalisation = res)
cnts_verc <- cc_data  %>% 
  rename(chs = severe) %>% 
  CaseMakerDT() %>% 
  rename(Severe = res)

cnts_ver <- cnts_vera %>% 
                 inner_join(cnts_verb) %>% 
                 inner_join(cnts_verc) %>% 
  rename(hcw_teacher = hcw_teacher2) %>% 
  gather("Outcome", "res", -hcw_teacher, -tp) %>% 
  mutate(Model = "Counts")

cnts_ver_fpa <- cc_data %>% 
  rename(chs = is_case) %>% 
  CaseMakerDTOverall() %>% 
  rename(Anycase = res)
cnts_ver_fpb <- cc_data %>% 
  rename(chs = hosp) %>% 
  CaseMakerDTOverall() %>% 
  rename(Hospitalisation = res)
cnts_ver_fpc <- cc_data  %>% 
  rename(chs = severe) %>% 
  CaseMakerDTOverall() %>% 
  rename(Severe = res)

cnts_verfp <- cnts_ver_fpa %>% 
  inner_join(cnts_ver_fpb) %>% 
  inner_join(cnts_ver_fpc) %>% 
  rename(hcw_teacher = hcw_teacher2) %>% 
  gather("Outcome", "res", -hcw_teacher, -tp) %>% 
  mutate(Model = "Counts")

cnts_ver <- bind_rows(cnts_ver, cnts_verfp)

t2_final <- bind_rows(t2 %>% slice(0),
                      cnts_ver %>% 
                        rename("Time period" = tp, Adjustment = Model) %>% 
                        spread(hcw_teacher, res) %>% 
                        mutate(Adjustment = "Cases/Controls",
                               Outcome = if_else(Outcome == "Anycase", "Any case", Outcome)),
                      t2) %>% 
  arrange(`Time period`, Outcome)
t2_final <- t2_final %>% 
  mutate(Outcome = if_else(Outcome == "Anycase", "Any case", Outcome))

t2_final <- t2_final %>% 
  mutate(`Time period` = if_else(`Time period` == "Full pandemic", "All time periods", `Time period`)) %>% 
  arrange(factor(`Time period`, levels = c("All time periods",
                                           "Spring/Summer 2020 (closed)",
                                           "Autumn term 2020 (open)",
                                           "Winter 2020/21 (closed)",
                                           "Spring term 2021 (phased)",
                                           "Summer term 2021 (open)"
                                           )),
          factor(Outcome, levels = c("Any case", "Hospitalisation", "Severe")),
          factor(Adjustment, levels = c("Cases/Controls", "Unadjusted", "Adjusted")))
saveRDS(t2_final, paste0("output/", folder_choose, "/", "table2_main_manuscript.Rds"))


# count individuals by sector, teachers only
CaseMakerDTSector <- function(x) {
  x <- setorder(x, -chs, anon_id)
  x <- x[, .SD[1], by = .(anon_id, sector_grp_tv, tp)] 
  x <- x[, .(cases = sum(chs),
             controls = sum(1-chs)),
         by = .(sector_grp_tv, tp)]
  x %>% 
    as_tibble() %>% 
    arrange(sector_grp_tv, tp) %>% 
    mutate(res = paste(cases, controls, sep = "/")) %>% 
    select(-cases, -controls)
}

CaseMakerDTSectorOverall <- function(x) {
  x <- setorder(x, -chs, anon_id)
  x <- x[, .SD[1], by = .(anon_id, sector_grp_tv)] 
  x <- x[, .(cases = sum(chs),
             controls = sum(1-chs)),
         by = .(sector_grp_tv)]
  x %>% 
    as_tibble() %>% 
    arrange(sector_grp_tv) %>% 
    mutate(tp = "Full pandemic",
           res = paste(cases, controls, sep = "/")) %>% 
    select(-cases, -controls)
}


cc_data2 <-cc_data %>% 
  filter(hcw_teacher2  == "Teacher")
indiv_cnt2a <- cc_data2 %>% 
  rename(chs = is_case) %>% 
  CaseMakerDTSector() %>% 
  rename(Anycase = res)
indiv_cnt2b <- cc_data2 %>% 
  rename(chs = hosp) %>% 
  CaseMakerDTSector() %>% 
  rename(Hospitalisation = res)
indiv_cnt2c <- cc_data2 %>% 
  rename(chs = severe) %>% 
  CaseMakerDTSector() %>% 
  rename(Severe = res)
indiv_cnt2 <-  indiv_cnt2a %>% 
  inner_join(indiv_cnt2b) %>% 
  inner_join(indiv_cnt2c) %>% 
  gather("Outcome", "res", -sector_grp_tv, -tp) %>% 
  mutate(Model = "Counts")
indiv_cnt2 <- indiv_cnt2 %>% 
  mutate(Outcome = if_else(Outcome == "Anycase", "Any case", Outcome))
cnts_t3 <- indiv_cnt2 %>% 
  rename("Time period" = tp, Adjustment = Model) %>% 
  spread(sector_grp_tv, res)

indiv_cnt_fp2a <- cc_data2 %>% 
  rename(chs = is_case) %>% 
  CaseMakerDTSectorOverall() %>% 
  rename(Anycase = res)
indiv_cnt_fp2b <- cc_data2 %>% 
  rename(chs = hosp) %>% 
  CaseMakerDTSectorOverall() %>% 
  rename(Hospitalisation = res)
indiv_cnt_fp2c <- cc_data2 %>% 
  rename(chs = severe) %>% 
  CaseMakerDTSectorOverall() %>% 
  rename(Severe = res)
indiv_cnt_fp2 <-  indiv_cnt_fp2a %>% 
  inner_join(indiv_cnt_fp2b) %>% 
  inner_join(indiv_cnt_fp2c) %>% 
  gather("Outcome", "res", -sector_grp_tv, -tp) %>% 
  mutate(Model = "Counts")
indiv_cnt_fp2 <- indiv_cnt_fp2 %>% 
  mutate(Outcome = if_else(Outcome == "Anycase", "Any case", Outcome))
cnts_t3_fp <- indiv_cnt_fp2 %>% 
  rename("Time period" = tp, Adjustment = Model) %>% 
  spread(sector_grp_tv, res)

cnts_t3 <- bind_rows(cnts_t3,
                     cnts_t3_fp)

t3final <- t3 %>% 
  select(`Time period`, `Outcome`, Adjustment, `Nursery/Primary or Nursery`, Primary,  Secondary, Other = `Teacher in other sector`  )

cnts_t3 <- cnts_t3 %>% 
  select(`Time period`, `Outcome`, Adjustment, `Nursery/Primary or Nursery`, Primary,  Secondary, Other) %>% 
  mutate(`Time period` = if_else(`Time period` == "Full pandemic", "All time periods", `Time period`))

t3final2 <- bind_rows(cnts_t3,
                      t3final) %>% 
  filter(!Outcome == "Severe") %>% 
  mutate(Adjustment = if_else(Adjustment == "Counts", "Cases/Controls", Adjustment))

t3final2 <- t3final2 %>% 
  arrange(factor(`Time period`, levels = c("All time periods",
                                           "Spring/Summer 2020 (closed)",
                                           "Autumn term 2020 (open)",
                                           "Winter 2020/21 (closed)",
                                           "Spring term 2021 (phased)",
                                           "Summer term 2021 (open)"
  )),
  factor(Outcome, levels = c("Any case", "Hospitalisation")),
  factor(Adjustment, levels = c("Cases/Controls", "Unadjusted", "Adjusted"))) 

saveRDS(t3final2, paste0("output/", folder_choose, "/", "table3_main_manuscript.Rds"))
})


## Extract rolling effect estimates for HR for cases ----
rm(cc_data, cc_data2)

rr_roll <- readRDS("output/model_rolling_rate_ratio.Rds")
rr_roll$mdl_do <- (map_lgl(rr_roll$mdl, ~ "character" %in% class(.x)))
rr_roll$mdl2_do <- (map_lgl(rr_roll$mdl2, ~ "character" %in% class(.x)))
## one only which did not run, weeks 30 and 31
rr_roll %>% filter(mdl_do | mdl2_do)
rr_roll <- rr_roll %>%
  filter(!mdl_do, ! mdl2_do)  %>%
  select(-mdl_do, -mdl2_do)

rr_roll$est <- map(rr_roll$mdl, ExtractEst)
rr_roll$est <- map(rr_roll$est, CleanEst)

rr_roll$est2 <- map(rr_roll$mdl2, ExtractEst)
rr_roll$est2 <- map(rr_roll$est2, CleanEst)

rr_roll_res1 <- rr_roll %>%
  select(test_week, est) %>%
  unnest(est)
rr_roll_res2 <- rr_roll %>%
  select(test_week, est2) %>%
  unnest(est2)

mdls <- bind_rows(minim = rr_roll_res2,
                  fully = rr_roll_res1,
                  .id = "models")
saveRDS(mdls, "output/model_rolling_rate_ratio_cleaned.Rds")
