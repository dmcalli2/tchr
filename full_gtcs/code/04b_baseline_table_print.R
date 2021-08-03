# 04b_baseline_table_print
library(tidyverse)
library(tableone)

CombineNP <- function(baseline_chars, baseline_chars_p){
  baseline_chars_p2 <- baseline_chars_p 
  baseline_chars_p2[-(1:3),] <- paste0(str_trim(baseline_chars_p[-(1:3),]), "%", " (", baseline_chars[-(1:3),], ")")
  baseline_chars_p2[] <- baseline_chars_p2 %>% 
    str_replace_all(fixed("(%)"), "")
  baseline_chars_p2[1,] <- baseline_chars[1,]
  baseline_chars_p2
}

DropLt5 <- function(mtrx) {
  mtrx[parse_integer(mtrx) <=5 & !is.na(parse_integer(mtrx))] <- "-"
  mtrx
}

bc <- readRDS("output/table1_components.Rds")

baseline_chars <- bc$baseline_chars %>% 
  print(contDigits = 0, noSpaces = TRUE, dropEqual = TRUE, catDigits = 0, format = "f")
baseline_chars <- DropLt5(baseline_chars)
baseline_chars_hbge10 <- bc$baseline_chars_hb %>% 
  print(contDigits = 0, noSpaces = TRUE, dropEqual = TRUE, catDigits = 0, format = "p")
baseline_chars_hblt10 <- bc$baseline_chars_hb %>% 
  print(contDigits = 0, noSpaces = TRUE, dropEqual = TRUE, catDigits = 1, format = "p")
baseline_chars_hb <- baseline_chars_hbge10
baseline_chars_hb[!is.na( as.integer(baseline_chars_hb)) & as.integer(baseline_chars_hb) < 10] <- baseline_chars_hblt10[!is.na( as.integer(baseline_chars_hb)) & as.integer(baseline_chars_hb) < 10]

baseline_chars_hb <- DropLt5(baseline_chars_hb)
baseline_chars_hb <- CombineNP(baseline_chars, baseline_chars_hb)
baseline_chars_hb[baseline_chars_hb == "% ()"] <- ""
rm(baseline_chars)

baseline_sect <- bc$baseline_sect %>% 
  print(contDigits = 0, noSpaces = TRUE, dropEqual = TRUE, catDigits = 0, format = "f")
baseline_sect <- DropLt5(baseline_sect)
baseline_sect_hbge10 <- bc$baseline_sect_hb %>% 
  print(contDigits = 0, noSpaces = TRUE, dropEqual = TRUE, catDigits = 0, format = "p")
baseline_sect_hblt10 <- bc$baseline_sect_hb %>% 
  print(contDigits = 0, noSpaces = TRUE, dropEqual = TRUE, catDigits = 1, format = "p")
baseline_sect_hb <- baseline_sect_hbge10
baseline_sect_hb[!is.na( as.integer(baseline_sect_hb)) & as.integer(baseline_sect_hb) < 10] <- baseline_sect_hblt10[!is.na( as.integer(baseline_sect_hb)) & as.integer(baseline_sect_hb) < 10]

baseline_sect_hb <- DropLt5(baseline_sect_hb)
baseline_sect_hb <- CombineNP(baseline_sect, baseline_sect_hb)
baseline_sect_hb[baseline_sect_hb == "% ()"] <- ""
rm(baseline_sect)


Reformat <- function(baseline_chars, baseline_sect){
  final_table1 <- cbind(baseline_chars[,c("Neither", "Healthcare Worker, patient facing", 
                                          "Healthcare worker, non-patient facing", "Teacher")], baseline_sect[,c(2,4,5,3)])
  # final_table1[nrow(final_table1),2] <- "-"
  lbls <- c("n" = "Number of individuals",
            "age_year (mean (SD))" = "Age, mean(sd*) *standard deviation",
            "female (%)" = "Female, %",
            "simd (%)" = "SIMD % (n)",
            "   1 - most deprived" = "   1 - most deprived", 
            "   2" = "   2", 
            "   3" = "   3",
            "   4" = "   4",
            "   5 - least deprived" = "   5 - least deprived",
            "   NA" = "   Unknown SIMD", 
            "ethnicity (%)" = "Ethnicity, % (n)",
            "   Asian" = "   Asian", 
            "   Black" = "   Black",
            "   Chinese" = "   Chinese",
            "   Other" = "   Other", 
            "   Unknown" = "   Unknown",
            "   White" = "   White",
            "shielding (%)" = "Shielding % (n)",
            "como_count (%)" = "Comorbidity count % (n)", 
            "   0" = "None", 
            "   1" = "One", 
            "   2+" = "Two or more", 
            "hhd (%)" = "Household member of patient facing healthcare worker,% (n)", 
            "hhd_tch (%)" = "Household member of teacher, % (n)",
            "dose_ever (%) " = "Vaccinated, % (n)",
            "a_unvacc" = "Unvaccinated",
            "b_post1st" = "First dose",
            "c_post2nd" = "Second dose")
  
  rownames(final_table1) <- lbls
  
  colnames(final_table1) <- c("Neither", "Healthcare Worker, patient facing", 
                              "Healthcare worker, non-patient facing", "Teacher",
                              "Nursery or Comb", "Primary", "Secondary", "Other")
  final_table1
}
final_table_hb <- Reformat(baseline_chars_hb, baseline_sect_hb)
final_table_hb["Household member of patient facing healthcare worker,% (n)", c("Healthcare Worker, patient facing", "Healthcare worker, non-patient facing")] <- "-"
final_table_hb["Household member of teacher, % (n)" , c("Teacher",         "Nursery or Comb", "Primary",         "Secondary",       "Other")] <- "-"
final_table_hb["Female, %", ] <- paste0(final_table_hb["Female, %", ], "%")

saveRDS(final_table_hb, "output/final_table1_svy_hb_simd.Rds")
