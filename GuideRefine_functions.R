library(tidyverse)
library(dplyr)
library(stringr)

#### MAIN FUNCTIONS

# function_df_count
function_df_count <- function(df_count, count_sgrna, num_genes, num_alignments, notes){
  df_count <- rbind(df_count, c(count_sgrna, num_genes, num_alignments, notes))
  kable(df_count)
  return(df_count)
}

# FUNCTION TO REMOVE MULTI-TARGET-/OFF-TARGET SGRNAS

# function to keep alignment only for normal chromosome
alignment_normal_chr <- function(alignment_data) {
    list_chromosome <- paste0('chr', c(1:22, 'X', 'Y'))
    alignment_data <- alignment_data %>%
        filter(chr %in% list_chromosome | is.na(chr))

    return(alignment_data)    
}


# function to remove multi-target sgRNA
multi_target_sgrna <- function(alignment_data){
    
    # Step 1: Identify sgRNAs that aligns perfectly with 0 mismatches to other location > 1 times
    multi_target_sgRNAs_count_dropped_guides_per_gene <- alignment_data %>%
            filter(n_mismatches == 0) %>%
            select(sgRNA, gene) %>%
            count(sgRNA, gene) %>%
            dplyr::rename('Frequency' = 'n') %>%
            filter(Frequency > 1)
      
    # Step 2: Identify sgRNAs with any non-zero mismatches
    # mismatched_sgRNAs <- alignment_data %>%
    #   filter(n_mismatches > 0) %>%
    #   pull(sgRNA) %>%
    #   unique()
    # 
    # Remove any sgRNA from Step 1 that includes n_mismatches > 0
    # multi_target_sgRNAs_count_dropped_guides_per_gene <- perfect_hits %>%
    #   filter(!sgRNA %in% mismatched_sgRNAs)
    
    # Number of dropped sgRNA due to multi-target
    dropped_due_to_multi_target <- multi_target_sgRNAs_count_dropped_guides_per_gene %>%
        count(gene) %>%
        dplyr::rename('multi-target sgRNA' = 'n')
    
    # List of multi-target sgRNAs
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

# GENE ANNOTATION FOR EACH SGRNA SEQUENCES (1)

# function to read NCBI Coding DNA sequences
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

# function to transform gene annotation CCDS
transform_gene_annotation_ccds <- function(imported_ccds_data, chosen_genome){
    list_chromosome <- paste0('chr', c(1:22, 'X', 'Y'))
    genome_info <- GenomeInfoDb::Seqinfo(genome = chosen_genome)[list_chromosome]
    
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

# function to make granges from alignment data
make_granges_from_alignment_data <- function(alignment_data, ref_bsgenome_selected){
  
    if(seqlevelsStyle(ref_bsgenome_selected) == "NCBI") {
      list_chromosome <- c(1:22, "X", "Y", "MT")
    } else if (seqlevelsStyle(ref_bsgenome_selected) == "UCSC") {
      list_chromosome <- paste0("chr", c(1:22, "X", "Y"))
    }
  
    genome_info <- GenomeInfoDb::Seqinfo(genome = unique(genome(ref_bsgenome_selected)))[list_chromosome]

    guide_aln_granges <- alignment_data %>%
        dplyr::select(unique_id = unique_aln_id,
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

# function to find the overlaps of gene annotation and the alignment
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


# function to visualize the gene symbols 

visualize_gene_symbol <- function(gene_df) {
  
  gene_df %>%
    filter(mismatches == 0) %>%
    mutate(notes = case_when(str_extract(sgrna, "(?<=sg)([A-Za-z0-9-]+)") == gene ~ "symbol match",
                             str_extract(sgrna, "(?<=sg)([A-Za-z0-9-]+)") != gene ~ "symbol mismatch")) %>%
    select(gene, notes) %>%
    distinct() %>%
    count(notes) %>%
    dplyr::rename("how_many" = "n") %>%
    ggplot(aes(notes, how_many, fill = notes)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = paste("Validation of sgRNA alignment with exon and gene symbol annotation"),
         x = "",
         y = "Number of genes") +
    theme(axis.text.x = element_text(size = 12, colour = "black"),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title.y = element_text(size = 14, colour = "black")) +
    geom_text(aes(label = how_many), vjust=-0.5) +
    scale_x_discrete(labels=c("symbol match" = "gene symbol match \nwith exon annotation",
                              "symbol mismatch" = "gene symbol not match \nwith exon annotation")) +
    scale_fill_manual(
      values = c("symbol match" = "#1f77b4", 
                 "symbol mismatch" = "#ff7f0e")) +
    ylim(0, length(unique(gene_df$gene))+1000)
  
}


# function to detect any previous/alias symbol based on the alignment and correct it
update_gene_symbol <- function(gene_df) {
  # categorize sgRNA with gene symbol match and sgRNA with gene symbol mismatch
  gene_df_check <- gene_df %>%
    filter(mismatches == 0) %>%
    mutate(notes = case_when(str_extract(sgrna, "(?<=sg)([A-Za-z0-9-]+)") == gene ~ "symbol match",
                             str_extract(sgrna, "(?<=sg)([A-Za-z0-9-]+)") != gene ~ "symbol mismatch"))
  
  # check for duplicate spacer in gene_df_match and discard it (sgRNA targeting the same gene but different exon location)
  gene_df_match <- gene_df_check %>%
    filter(notes == "symbol match") %>%
    mutate(duplicate_spacer_gene = duplicated(spacer, gene)) %>%
    arrange(sgrna) %>%
    filter(!duplicate_spacer_gene == TRUE)
  
  # separate sgRNA with gene symbol mismatch dataframe, discard duplicated spacer
  gene_df_mismatch <- gene_df_check %>%
    filter(notes == "symbol mismatch") %>%
    mutate(duplicate_spacer_gene = duplicated(spacer, gene)) %>%
    arrange(gene) %>%
    filter(!duplicate_spacer_gene == TRUE)
  
  # arrange sequential numbers here
  run_length <- rle(gene_df_mismatch$gene)$lengths
  sequential_numbers <- sequence(run_length)
  
  # correct sgRNA - gene name according to the alignment exon location and add new sgRNA number id alongside _SC (corrected gene symbol)
  gene_df_mismatch <- gene_df_mismatch %>%
    mutate(new_sgrna_id = paste("sg", gene, "_", sequential_numbers, "_SC" , sep ="")) %>%
    select(-sgrna) %>%
    dplyr::rename("sgrna" = "new_sgrna_id") %>%
    relocate(sgrna, .before = spacer) %>%
    mutate(notes = "symbol corrected")
  
  # find corrected genes that already exist in gene_df_match put it into a new dataframe called as corrected_exist_gene
  list_corrected_exist_gene <- intersect(gene_df_match$gene, gene_df_mismatch$gene)
  corrected_exist_gene <- gene_df_match %>%
    filter(gene %in% list_corrected_exist_gene)
  
  # filter out problematic gene (gene that originally targeted by sgRNA in the library, but also targeted by another sgRNA with corrected symbol)
  gene_df_match <- gene_df_match %>%
    filter(!gene %in% list_corrected_exist_gene)
  
  # join gene_df_mismatch with gene_df_match to find the duplicated sgRNA
  gene_df_mismatch <- gene_df_mismatch %>%
    full_join(corrected_exist_gene) %>%
    arrange(gene, desc(notes)) %>%
    mutate(duplicate_spacer_gene = duplicated(spacer, gene))
  
  # full_join the gene_df_mismatch with gene_df_match to become a gene_df_corrected
  gene_df_corrected <- gene_df_match %>%
    full_join(gene_df_mismatch)
  
  return(list("gene_symbol_check" = gene_df_check, 
              "gene_symbol_corrected" = gene_df_corrected))
}

#### TEMPORARY FUNCTION

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
