library(tidyverse)

# function_df_count
function_df_count <- function(df_count, count_sgrna, num_genes, num_alignments, notes){
  df_count <- rbind(df_count, c(count_sgrna, num_genes, num_alignments, notes))
  kable(df_count)
  return(df_count)
}

# function to keep alignment only for normal chromosome
alignment_normal_chr <- function(alignment_data) {
    list_chromosome <- paste0('chr', c(1:22, 'X', 'Y'))
    alignment_data <- alignment_data %>%
        filter(chr %in% list_chromosome | is.na(chr))

    return(alignment_data)    
}


# function to remove multi-target sgRNA
multi_target_sgrna <- function(alignment_data){
    
    multi_target_sgRNAs_count_dropped_guides_per_gene <- alignment_data %>%
        filter(n_mismatches == 0) %>%
        select(sgRNA, gene) %>%
        count(sgRNA, gene) %>%
        dplyr::rename('Frequency' = 'n') %>%
        filter(Frequency > 1)
    
    dropped_due_to_multi_target <- multi_target_sgRNAs_count_dropped_guides_per_gene %>%
        count(gene) %>%
        dplyr::rename('multi-target sgRNA' = 'n')

    multi_target_sgRNA_list_guides <- c(multi_target_sgRNAs_count_dropped_guides_per_gene$sgRNA)

    return(list(
        dropped_guides_num = multi_target_sgRNAs_count_dropped_guides_per_gene,
        dropped_genes_num = dropped_due_to_multi_target,
        list_guides = multi_target_sgRNA_list_guides
    ))
}

# function to remove single mismatch sgRNA
single_mismatch_sgrna <- function(alignment_data) {
    single_mismatch_sgRNA_count_dropped_guides_per_gene <- alignment_data %>%
        filter(n_mismatches == 1) %>%
        count(sgRNA, gene) %>%
        dplyr::rename('Frequency' = 'n')

    dropped_due_to_single_mismatch <- single_mismatch_sgRNA_count_dropped_guides_per_gene %>%
        count(gene) %>%
        dplyr::rename('single mismatch sgRNA' = 'n')

    single_mismatch_sgRNA_list_guides <- c(unique(single_mismatch_sgRNA_count_dropped_guides_per_gene$sgRNA))

    return(list(
        dropped_guides_num = single_mismatch_sgRNA_count_dropped_guides_per_gene,
        dropped_genes_num = dropped_due_to_single_mismatch,
        list_guides = single_mismatch_sgRNA_list_guides
    ))

}

# function to remove sgRNA aligned with double mismatches at the most PAM distal site

pam_distal_double_mismatch <- function(alignment_data){
    detect_pam_distal_double_mismatch <- alignment_data %>% 
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
            if_else(is.na(mismatch_loc[1]) | mismatch_loc[1] == 0,'perfect',
            if_else(strand == '+' & (mismatch_loc[1] < 3 & mismatch_loc[2] < 3),'pam-distal double mismatch',
            if_else(strand == '-' & (mismatch_loc[1] > 18 & mismatch_loc[2] > 18),'pam-distal double mismatch',
                'double mismatch'
                ))
            )
        )

    pam_distal_double_mismatch <- detect_pam_distal_double_mismatch %>%
        filter(type == 'pam-distal double mismatch')

    pam_distal_count_dropped_guides_per_gene <- pam_distal_double_mismatch %>%
        select(sgRNA, gene) %>%
        count(sgRNA, gene) %>%
        dplyr::rename('Frequency' = 'n')

    dropped_due_to_pam_distal_double_mismatch <- pam_distal_count_dropped_guides_per_gene %>%
        count(gene) %>%
        dplyr::rename('PAM distal double mismatch' = 'n')
    
    pam_distal_double_mismatch_list_guides <- c(unique(pam_distal_count_dropped_guides_per_gene$sgRNA))

    return(list(
        pam_distal_df = detect_pam_distal_double_mismatch,
        dropped_guides_num = pam_distal_count_dropped_guides_per_gene,
        dropped_genes_num = dropped_due_to_pam_distal_double_mismatch,
        list_guides = pam_distal_double_mismatch_list_guides
    ))

}