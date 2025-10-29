library(tidyverse)

detect_pam_distal_double_mismatch <- library_alignment %>% 
  dplyr::rename(pam_start = pam_site) %>%
  mutate(
    pam_start = ifelse(strand == '+', pam_start, pam_start-2),
    pam_end = ifelse(strand == '+', pam_start+2, pam_start+2),
    start_pos = ifelse(strand == '+', pam_start-21, pam_start+3),
    end_pos = ifelse(strand == '+', pam_start-1, pam_start+23),
    cut_pos = ifelse(strand == '+', pam_start-4, pam_start+5)) %>%
  relocate(pam_start, .after = cut_pos) %>%
  relocate(start_pos, .before = cut_pos) %>%
  relocate(end_pos, .after = start_pos) %>%
  relocate(pam_end, .after = pam_start) %>%
  rowwise() %>%
  mutate(
    mismatch_loc = list(which(strsplit(spacer, '')[[1]] != strsplit(protospacer, '')[[1]]))) %>%
  mutate(
    type =
      case_when(is.na(mismatch_loc[1]) | mismatch_loc[1] == 0 ~ "perfect",
                strand == "+" & mismatch_loc[1] < 3 & is.na(mismatch_loc[2]) ~ "pam-distal single mismatch",
                strand == "-" & mismatch_loc[1] > 18 & is.na(mismatch_loc[2]) ~ "pam-distal single mismatch",
                strand == '+' & (mismatch_loc[1] < 3 & mismatch_loc[2] < 3) ~ "pam-distal double mismatch",
                strand == '-' & (mismatch_loc[1] > 18 & mismatch_loc[2] > 18) ~ "pam-distal double mismatch",
                length(mismatch_loc) < 2 ~ "single mismatch",
                length(mismatch_loc) > 1 ~ "double mismatch",
                TRUE ~ "Uncategorized"
                ))



