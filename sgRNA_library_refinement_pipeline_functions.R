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

read_ccds_data <- function(ccds_filename){
    ccds <- read_delim(ccds_filename, delim = '\t')

    ccds$cds_from <- as.integer(ccds$cds_from)
    ccds$cds_to <- as.integer(ccds$cds_to)

    list_chromosome <- paste0('chr', c(1:22, 'X', 'Y'))

    ccds <- ccds %>%
        dplyr::rename('chromosome' = '#chromosome') %>%
        mutate(chromosome = str_c("chr", chromosome)) %>%
        filter(ccds_status %in% c("Public", "Reviewed, update pending", "Under review, update")) %>%
        filter(chromosome %in% list_chromosome, !is.na(cds_from), !is.na(cds_to))

    return(ccds)

}

transform_gene_annotation_ccds <- function(imported_ccds_data){
    list_chromosome <- paste0('chr', c(1:22, 'X', 'Y'))
    genome_info <- GenomeInfoDb::Seqinfo(genome = "hg38")[list_chromosome]

    ccds_exon <- imported_ccds_data %>%
        mutate(cds_interval = str_replace_all(cds_locations, "[\\[\\]]", "") %>% str_split("\\s*,\\s*")) %>%
        tidyr::unnest(cds_interval) %>%
        group_by(gene, gene_id, cds_locations) %>% 
        mutate(exon_code = ifelse(cds_strand=="+", 1:dplyr::n(), dplyr::n():1)) %>% 
        ungroup() %>%
        dplyr::mutate(cds_start = str_extract(cds_interval, "^[0-9]+") %>% as.integer,
                        cds_end = str_extract(cds_interval, "[0-9]+$") %>% as.integer) %>%
        dplyr::select(gene, gene_id, chromosome, start=cds_start, end=cds_end, strand=cds_strand,
                        gene_start=cds_from, gene_end=cds_to, exon_code)

    gene_annot_granges <- ccds_exon %>% 
        GenomicRanges::makeGRangesFromDataFrame(seqinfo = genome_info, 
                                                keep.extra.columns = T, 
                                                na.rm = TRUE)
    return(list(
        ccds_exon_df = ccds_exon,
        gene_annot_granges_df = gene_annot_granges
    ))

}

# temporary function to meet the requirement of cut_pos for granges alignments
add_cut_pos_pam_pos <- function(alignment_data){
    alignment_data <- alignment_data %>%
        dplyr::rename(pam_start = pam_site) %>%
        mutate(
            pam_start = ifelse(strand == '+', pam_start, pam_start-2),
            pam_end = ifelse(strand == '+', pam_start+2, pam_start+2),
            start_pos = ifelse(strand == '+', pam_start-21, pam_start+3),
            end_pos = ifelse(strand == '+', pam_start-1, pam_start+23),
            cut_pos = ifelse(strand == '+', pam_start-4, pam_start+5),
            unique_aln_id = paste0('sgrna_aln_', row_number())) %>%
        relocate(pam_start, .after = cut_pos) %>%
        relocate(start_pos, .before = cut_pos) %>%
        relocate(end_pos, .after = start_pos) %>%
        relocate(pam_end, .after = pam_start) %>%
        relocate(unique_aln_id, .before = sgRNA)

    return(alignment_data)
}

make_granges_from_alignment_data <- function(alignment_data){
    list_chromosome <- paste0('chr', c(1:22, 'X', 'Y'))
    genome_info <- GenomeInfoDb::Seqinfo(genome = "hg38")[list_chromosome]

    guide_aln_granges <- alignment_data %>%
        select(unique_id = unique_aln_id,
                id = sgRNA, 
                Spacer = spacer, 
                Protospacer = protospacer, 
                Chr = chr, 
                Start = cut_pos, 
                strand, 
                mismatch = n_mismatches) %>%
        mutate(End = Start) %>%
        GenomicRanges::makeGRangesFromDataFrame(seqinfo = genome_info, 
                                                keep.extra.columns = T,
                                                na.rm = TRUE)
    
    return(guide_aln_granges)
}

find_overlaps_gene_annotation_and_alignment <- function(guide_aln_granges, gene_annot_granges) {
    hits <- GenomicRanges::findOverlaps(guide_aln_granges, gene_annot_granges, ignore.strand=T) %>% as_tibble

    gene_df <- hits %>%
    transmute(unique_aln_id = guide_aln_granges$unique_id[queryHits] %>% as.character(),
              sgrna = guide_aln_granges$id[queryHits],
              spacer = guide_aln_granges$Spacer[queryHits],
              protospacer = guide_aln_granges$Protospacer[queryHits],
              mismatches = guide_aln_granges$mismatch[queryHits],
              chr = GenomicRanges::seqnames(guide_aln_granges)[queryHits] %>% as.character() ,
              cut_pos = GenomicRanges::start(guide_aln_granges)[queryHits] %>% as.integer(),
              strand = GenomicRanges::strand(guide_aln_granges)[queryHits] %>% as.character(),
              gene = gene_annot_granges$gene[subjectHits],
              CDS_strand = GenomicRanges::strand(gene_annot_granges)[subjectHits] %>% as.character(),
              CDS_start = GenomicRanges::start(gene_annot_granges)[subjectHits] %>% as.integer(),
              CDS_end = GenomicRanges::end(gene_annot_granges)[subjectHits] %>% as.integer()) %>%
    distinct()

    return(gene_df)

}
