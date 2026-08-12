library(tidyverse)

STANDARD_CHROMOSOMES <- paste0("chr", c(1:22, "X", "Y"))

# Mirrors the inline filtering in GuideRefine.Rmd for use in downstream analysis.
alignment_normal_chr <- function(alignment) {
  alignment %>% filter(chr %in% STANDARD_CHROMOSOMES | is.na(chr))
}

function_df_count <- function(df_count, count_sgrna, num_genes, num_alignments, notes){
  df_count <- rbind(df_count, c(count_sgrna, num_genes, num_alignments, notes))
  return(df_count)
}

multi_target_sgrna <- function(alignment_data){
    
    # Identify sgRNAs with more than 1 perfect alignment (0 mismatches) anywhere in the genome
    multi_target_sgRNAs <- alignment_data %>%
        filter(n_mismatches == 0) %>%
        count(sgRNA) %>%
        dplyr::rename('Frequency' = 'n') %>%
        filter(Frequency > 1)

    # Join back to alignment data to get per-gene breakdown for reporting
    multi_target_sgRNAs_count_dropped_guides_per_gene <- alignment_data %>%
        filter(sgRNA %in% multi_target_sgRNAs$sgRNA) %>%
        select(sgRNA, gene) %>%
        distinct()

    dropped_due_to_multi_target <- multi_target_sgRNAs_count_dropped_guides_per_gene %>%
        count(gene) %>%
        dplyr::rename('multi-target sgRNA' = 'n')

    multi_target_sgRNA_list_guides <- unique(multi_target_sgRNAs$sgRNA)

    return(list(
        dropped_guides_num = multi_target_sgRNAs_count_dropped_guides_per_gene,
        dropped_genes_num = dropped_due_to_multi_target,
        list_guides = multi_target_sgRNA_list_guides
    ))
}

single_mismatch_sgrna <- function(alignment_data) {
    single_mismatch_sgRNA_count_dropped_guides_per_gene <- alignment_data %>%
        filter(n_mismatches == 1) %>%
        count(sgRNA, gene) %>%
        dplyr::rename('Frequency' = 'n')

    dropped_due_to_single_mismatch <- single_mismatch_sgRNA_count_dropped_guides_per_gene %>%
        count(gene) %>%
        dplyr::rename('single mismatch sgRNA' = 'n')

    single_mismatch_sgRNA_list_guides <- unique(single_mismatch_sgRNA_count_dropped_guides_per_gene$sgRNA)

    return(list(
        dropped_guides_num = single_mismatch_sgRNA_count_dropped_guides_per_gene,
        dropped_genes_num = dropped_due_to_single_mismatch,
        list_guides = single_mismatch_sgRNA_list_guides
    ))

}

# function to remove sgRNA aligned with mismatches at the most PAM distal site
pam_distal_mismatch <- function(alignment_data){
  
    annotate_cut_pos <- annotate_pam_positions(alignment_data)

    # There is a possibility that a single guide share the pam-distal single mismatch and pam-distal double mismatch
    # that is OK as long as they are detected in either one condition
    
    detect_pam_distal_mismatches <- annotate_cut_pos %>%
        mutate(
            mismatch_loc = mapply(
                function(s, p) which(strsplit(s, "")[[1]] != strsplit(p, "")[[1]]),
                spacer, protospacer,
                SIMPLIFY = FALSE
            ),
            guide_len = nchar(spacer),
            n_mm = lengths(mismatch_loc),
            pos1 = vapply(mismatch_loc, function(x) if (length(x) >= 1) x[1L] else NA_integer_, integer(1L)),
            pos2 = vapply(mismatch_loc, function(x) if (length(x) >= 2) x[2L] else NA_integer_, integer(1L)),
            type = case_when(
                n_mm == 0                                                      ~ "perfect",
                strand == "+" & pos1 < 3 & n_mm == 1                          ~ "pam-distal single mismatch",
                strand == "-" & pos1 > (guide_len - 2) & n_mm == 1            ~ "pam-distal single mismatch",
                strand == "+" & pos1 < 3 & pos2 < 3                           ~ "pam-distal double mismatch",
                strand == "-" & pos1 > (guide_len - 2) & pos2 > (guide_len - 2) ~ "pam-distal double mismatch",
                n_mm == 1                                                      ~ "single mismatch",
                n_mm > 1                                                       ~ "double mismatch",
                TRUE                                                           ~ "Uncategorized"
            )
        ) %>%
        select(-guide_len, -n_mm, -pos1, -pos2)

    # PAM-distal single mismatches
    pam_distal_single_mismatch <- detect_pam_distal_mismatches %>%
      filter(type == 'pam-distal single mismatch') %>%
      pull(sgRNA) %>%
      unique()
    
    pam_distal_single_mismatch_count_dropped_guides_per_gene <- detect_pam_distal_mismatches %>%
      filter(sgRNA %in% pam_distal_single_mismatch) %>%
      select(sgRNA, gene) %>%
      count(sgRNA, gene) %>%
      dplyr::rename('Frequency' = 'n')
    
    dropped_due_to_pam_distal_single_mismatch <- pam_distal_single_mismatch_count_dropped_guides_per_gene %>%
      count(gene) %>%
      dplyr::rename('PAM distal single mismatch' = 'n')
    
    
    # PAM-distal double mismatches
    pam_distal_double_mismatch <- detect_pam_distal_mismatches %>%
      filter(type == 'pam-distal double mismatch') %>%
      pull(sgRNA) %>%
      unique()
    
    pam_distal_double_mismatch_count_dropped_guides_per_gene <- detect_pam_distal_mismatches %>%
      filter(sgRNA %in% pam_distal_double_mismatch) %>%
      select(sgRNA, gene) %>%
      count(sgRNA, gene) %>%
      dplyr::rename('Frequency' = 'n')
    
    dropped_due_to_pam_distal_double_mismatch <- pam_distal_double_mismatch_count_dropped_guides_per_gene %>%
      count(gene) %>%
      dplyr::rename('PAM distal double mismatch' = 'n')
    
    return(list(
        pam_distal_df = detect_pam_distal_mismatches,
        dropped_guides_num_pam_single = pam_distal_single_mismatch_count_dropped_guides_per_gene,
        dropped_genes_num_pam_single = dropped_due_to_pam_distal_single_mismatch,
        dropped_guides_num_pam_double = pam_distal_double_mismatch_count_dropped_guides_per_gene,
        dropped_genes_num_pam_double = dropped_due_to_pam_distal_double_mismatch,
        list_guides_pam_single = pam_distal_single_mismatch,
        list_guides_pam_double = pam_distal_double_mismatch
    ))

}

# Parses an NCBI CCDS annotation file into one row per transcript,
# keeping only actively-maintained records on standard chromosomes.
read_ccds_data <- function(ccds_filename) {
    active_status <- c("Public", "Reviewed, update pending", "Under review, update")

    readr::read_tsv(ccds_filename, col_types = readr::cols(.default = "c")) %>%
        dplyr::rename(chromosome = `#chromosome`) %>%
        dplyr::mutate(
            chromosome = paste0("chr", chromosome),
            cds_from = as.integer(cds_from),
            cds_to = as.integer(cds_to)
        ) %>%
        dplyr::filter(
            chromosome %in% STANDARD_CHROMOSOMES,
            ccds_status %in% active_status,
            !is.na(cds_from),
            !is.na(cds_to)
        )
}

# Expands CCDS transcripts into one row per exon (parsed from the
# "start-end, start-end, ..." cds_locations column) and builds the
# matching GRanges gene annotation object.
transform_gene_annotation_ccds <- function(imported_ccds_data, chosen_genome) {
    exon_table <- imported_ccds_data %>%
        dplyr::mutate(exon_range = stringr::str_extract_all(cds_locations, "[0-9]+-[0-9]+")) %>%
        tidyr::unnest(exon_range) %>%
        tidyr::separate(exon_range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
        dplyr::group_by(gene, gene_id, cds_locations) %>%
        dplyr::mutate(
            exon_code = dplyr::row_number(),
            exon_code = dplyr::if_else(cds_strand == "-", dplyr::n() - exon_code + 1L, exon_code)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::transmute(
            gene, gene_id, chromosome, start, end,
            strand = cds_strand,
            gene_start = cds_from,
            gene_end = cds_to,
            exon_code
        )

    seqinfo <- GenomeInfoDb::Seqinfo(genome = chosen_genome)[STANDARD_CHROMOSOMES]
    gene_annot_granges <- GenomicRanges::makeGRangesFromDataFrame(
        exon_table, seqinfo = seqinfo, keep.extra.columns = TRUE, na.rm = TRUE
    )

    list(ccds_exon_df = exon_table, gene_annot_granges_df = gene_annot_granges)
}

make_granges_from_alignment_data <- function(alignment_data, ref_bsgenome_selected) {
    style <- GenomeInfoDb::seqlevelsStyle(ref_bsgenome_selected)
    chrom_names <- switch(style,
        NCBI = c(as.character(1:22), "X", "Y", "MT"),
        UCSC = paste0("chr", c(as.character(1:22), "X", "Y")),
        stop("Unsupported seqlevels style: ", style)
    )
    seqinfo <- GenomeInfoDb::Seqinfo(genome = unique(GenomeInfoDb::genome(ref_bsgenome_selected)))[chrom_names]

    alignment_data %>%
        dplyr::transmute(
            unique_id = unique_aln_id,
            id = sgRNA,
            Spacer = spacer,
            Protospacer = protospacer,
            Chr = chr,
            Start = cut_pos,
            End = cut_pos,
            strand,
            mismatch = n_mismatches
        ) %>%
        GenomicRanges::makeGRangesFromDataFrame(seqinfo = seqinfo, keep.extra.columns = TRUE, na.rm = TRUE)
}

# Overlaps guide cut-site GRanges with gene-exon GRanges and returns one
# row per (guide, overlapping exon) pair with both sets of annotations.
find_overlaps_gene_annotation_and_alignment <- function(guide_aln_granges, gene_annot_granges) {
    overlap_hits <- GenomicRanges::findOverlaps(guide_aln_granges, gene_annot_granges, ignore.strand = TRUE)
    guide_idx <- S4Vectors::queryHits(overlap_hits)
    gene_idx <- S4Vectors::subjectHits(overlap_hits)

    tibble::tibble(
        unique_aln_id = as.character(guide_aln_granges$unique_id[guide_idx]),
        sgrna         = guide_aln_granges$id[guide_idx],
        spacer        = guide_aln_granges$Spacer[guide_idx],
        protospacer   = guide_aln_granges$Protospacer[guide_idx],
        mismatches    = guide_aln_granges$mismatch[guide_idx],
        chr           = as.character(GenomicRanges::seqnames(guide_aln_granges))[guide_idx],
        cut_pos       = as.integer(GenomicRanges::start(guide_aln_granges))[guide_idx],
        strand        = as.character(GenomicRanges::strand(guide_aln_granges))[guide_idx],
        gene          = gene_annot_granges$gene[gene_idx],
        CDS_strand    = as.character(GenomicRanges::strand(gene_annot_granges))[gene_idx],
        CDS_start     = as.integer(GenomicRanges::start(gene_annot_granges))[gene_idx],
        CDS_end       = as.integer(GenomicRanges::end(gene_annot_granges))[gene_idx]
    ) %>%
        dplyr::distinct()
}


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
  gene_df_check <- gene_df %>%
    filter(mismatches == 0) %>%
    mutate(notes = case_when(str_extract(sgrna, "(?<=sg)([A-Za-z0-9-]+)") == gene ~ "symbol match",
                             str_extract(sgrna, "(?<=sg)([A-Za-z0-9-]+)") != gene ~ "symbol mismatch"))
  
  # check for duplicate spacer in gene_df_match and discard it (sgRNA targeting the same gene but different exon location)
  gene_df_match <- gene_df_check %>%
    filter(notes == "symbol match") %>%
    mutate(duplicate_spacer_gene = duplicated(data.frame(spacer, gene))) %>%
    arrange(sgrna) %>%
    filter(!duplicate_spacer_gene == TRUE)
  
  # separate sgRNA with gene symbol mismatch dataframe, discard duplicated spacer
  gene_df_mismatch <- gene_df_check %>%
    filter(notes == "symbol mismatch") %>%
    mutate(duplicate_spacer_gene = duplicated(data.frame(spacer, gene))) %>%
    arrange(gene) %>%
    filter(!duplicate_spacer_gene == TRUE)
  
  run_length <- rle(gene_df_mismatch$gene)$lengths
  sequential_numbers <- sequence(run_length)
  
  # correct sgRNA - gene name according to the alignment exon location and add new sgRNA number id alongside _SC (corrected gene symbol)
  gene_df_mismatch <- gene_df_mismatch %>%
    mutate(new_sgrna_id = paste("sg", gene, "_", sequential_numbers, "_SC" , sep ="")) %>%
    select(-sgrna) %>%
    dplyr::rename("sgrna" = "new_sgrna_id") %>%
    relocate(sgrna, .before = spacer) %>%
    mutate(notes = "symbol corrected")
  
  list_corrected_exist_gene <- intersect(gene_df_match$gene, gene_df_mismatch$gene)
  corrected_exist_gene <- gene_df_match %>%
    filter(gene %in% list_corrected_exist_gene)
  
  # filter out problematic gene (gene that originally targeted by sgRNA in the library, but also targeted by another sgRNA with corrected symbol)
  gene_df_match <- gene_df_match %>%
    filter(!gene %in% list_corrected_exist_gene)
  
  gene_df_mismatch <- gene_df_mismatch %>%
    full_join(corrected_exist_gene, by = names(corrected_exist_gene)) %>%
    arrange(gene, desc(notes)) %>%
    mutate(duplicate_spacer_gene = duplicated(data.frame(spacer, gene)))
  
  gene_df_corrected <- gene_df_match %>%
    full_join(gene_df_mismatch, by = names(gene_df_mismatch))
  
  return(list("gene_symbol_check" = gene_df_check, 
              "gene_symbol_corrected" = gene_df_corrected))
}

annotate_pam_positions <- function(alignment_data) {
    alignment_data %>%
        dplyr::rename(pam_start = pam_site) %>%
        mutate(
            pam_start = ifelse(strand == '+', pam_start, pam_start-2),
            pam_end   = ifelse(strand == '+', pam_start+2, pam_start+2),
            start_pos = ifelse(strand == '+', pam_start-21, pam_start+3),
            end_pos   = ifelse(strand == '+', pam_start-1, pam_start+23),
            cut_pos   = ifelse(strand == '+', pam_start-4, pam_start+5)) %>%
        relocate(pam_start, .after = cut_pos) %>%
        relocate(start_pos, .before = cut_pos) %>%
        relocate(end_pos, .after = start_pos) %>%
        relocate(pam_end, .after = pam_start)
}
