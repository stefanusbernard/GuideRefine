# sgRNA Scoring Pipeline
# by Stefanus Bernard

# SET CURRENT WORKING DIRECTORY
setwd(getwd())

# INSTALL PACKAGES WHEREAS REQUIRED
list_packages <- c('crisprVerse', 
                      'crisprBase',
                      'crisprBowtie',
                      'crisprDesign',
                      'crisprDesignData',
                      'crisprScore',
                      'BSgenome.Hsapiens.UCSC.hg38',
                      'tidyverse',
                      'dplyr',
                      'ggfortify',
                      'reshape',
                      'hgnc')
new_packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)


# LOAD REQUIRED PACKAGES
library(crisprVerse)
library(crisprBase)
library(crisprBowtie)
library(crisprDesign)
library(crisprDesignData)
library(crisprScore)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(dplyr)
library(ggfortify)
library(reshape)
library(hgnc)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

# LOAD SGRNA LIBRARY (SHOULD HAVE 3 COLUMNS CONSISTING OF 'sgRNA', 'gene', and 'spacer')

sgRNA_library_selected <- 'Cellecta'
sgRNA_library <- read_csv(paste('raw_data/', sgRNA_library_selected, '.csv', sep = ''))

print(paste('Unique sgRNA in the library: ', length(unique(sgRNA_library$sgRNA))))
print(paste('Unique gene in the library: ', length(unique(sgRNA_library$gene))))

# BUILD BOWTIE INDEX

# fastaFile <- 'bowtie_index/hg38.fa.gz'

list.files(path = 'bowtie_index/', pattern = '*.ebwt$')

check_index_file = 'hg38.1.ebwt'
file_test('-f', paste('bowtie_index/', check_index_file, sep='')) # returns TRUE


if (file_test('-f', paste('bowtie_index/', check_index_file, sep='')) == FALSE) {
  bowtie_build(fastaFile,
               outdir="./bowtie2_index",
               force=TRUE,
               prefix='hg38')
} else {
  print(paste(check_index_file, 'is in directory, indicating indexing by bowtie have been done previously'))
}


# ALIGNMENT OF SGRNA TO HG38 UCSC USING BOWTIE

data(SpCas9, package="crisprBase")
crisprNuclease <- SpCas9

spacers_sequence <- c(unique(sgRNA_library$spacer))

bowtie_index <- file.path(getwd(), '/bowtie_index', 'hg38')

# change which spacer would like to be used cellecta or avana
scoring_alignment <- runCrisprBowtie(spacers_sequence,
                                     crisprNuclease=crisprNuclease,
                                     n_mismatches=0,
                                     canonical=TRUE,
                                     bowtie_index=bowtie_index)

scoring_data <- sgRNA_library %>%
  left_join(scoring_alignment, by = join_by(spacer))

scoring_data

rm(scoring_alignment)
rm(SpCas9)
gc()

# DATA CLEANING BEFORE CONSTRUCTING GUIDESET

# DROP ROW WHERE NA DETECTED 
# ONLY SELECT STANDARD CHROMOSOME
# DISCARD SGRNA TARGETING INTRON CONTROL GENE

list_chromosome <- paste0('chr', c(1:22, 'X', 'Y'))

scoring_data <- scoring_data %>%
  drop_na() %>%
  filter(chr %in% list_chromosome) %>%
  filter(gene != 'INTRON_CTRL')

# IDENTIFY AND DROP MULTI-TARGET GUIDES

function_guides_filtering <- function(alignment_data) {
  
  target_guides_count <- alignment_data %>%
    count(sgRNA) %>%
    dplyr::rename('Frequency of alignments' = 'n')
  
  multi_target_guides <- target_guides_count %>%
    filter(`Frequency of alignments` > 1)
  
  top_10_multi_target_guides <- target_guides_count %>%
    filter(`Frequency of alignments` > 1) %>%
    arrange(desc(`Frequency of alignments`)) %>%
    top_n(n = 10)
  
  single_target_guides <- target_guides_count %>%
    filter(`Frequency of alignments` == 1)
  
  no_target_guides <- target_guides_count %>%
    filter(`Frequency of alignments` < 1)
  
  cleaned_alignment <- alignment_data %>%
    filter(sgRNA %in% single_target_guides$sgRNA)
  
  return(list(single_target_guides, multi_target_guides, top_10_multi_target_guides, cleaned_alignment))
  
}

data <- function_guides_filtering(scoring_data)
single_target_guides <- data[[1]]
multi_target_guides <- data[[2]]
top_10_multi_target_guides <- data[[3]]
cleaned_alignment <- data[[4]]

rm(scoring_data)
rm(list_chromosome)
gc()

# PLOT MULTI-TARGET GUIDES

hgnc_tbl <- import_hgnc_dataset()
sgRNA_read_counts <- read_csv('raw_data/All_sgRNA_Count.csv')
sgRNA_read_counts <- sgRNA_read_counts %>%
  select(-DAY24_XL413_GL)

ggplot_data <- sgRNA_read_counts %>%
  filter(sgRNA %in% top_10_multi_target_guides$sgRNA) %>%
  left_join(top_10_multi_target_guides, by = 'sgRNA') %>%
  mutate(sgRNA = paste(sgRNA, ' - Frequency:', Frequency)) %>%
  select(-Frequency) %>%
  pivot_longer(cols=c('DAY0_CNTR_A0', 'DAY20_DMSO_FD'),
               names_to='experiment',
               values_to='sgRNA_count')

final_plot <- ggplot_data %>%
  ggplot(aes(experiment, sgRNA_count, color=sgRNA, group=sgRNA)) +
  geom_point(size=1.5, shape=16) +
  geom_line(linewidth=0.5) +
  labs(x = 'Experiment',
       y = 'sgRNA Count',
       title = paste('Top 10 Multi-Target sgRNAs in ', sgRNA_library_selected, ' Library', sep=''),
       color = 'sgRNA') +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.1, face = 'bold', size = 15),
        plot.subtitle = element_text(hjust = 0.5, face = 'bold', size = 10),
        axis.title.x = element_text(vjust=0.1),
        axis.title.y = element_text(hjust=0.5,vjust = 3),
        panel.grid.minor.x = element_blank())

final_plot

ggsave(
  path = 'plot/',
  plot = final_plot,
  filename = paste(sgRNA_library_selected, '-multi-target-guides.png', sep=''))

# COUNT THE PROCESS

sgrna_origin <- length(sgRNA_library$sgRNA)
sgrna_multi_target <- length(multi_target_guides$Frequency)
sgrna_single_target <- length(single_target_guides$Frequency)

stage <- c('Raw sgRNA Data', 'Single-Target sgRNA', 'Multi-Target sgRNA')
number_of_sgrna <- c(sgrna_origin,
                     sgrna_single_target, 
                     sgrna_multi_target)

sgrna_count_table <- data.frame(stage, number_of_sgrna)
write.csv(sgrna_count_table, file = paste(sgRNA_library_selected, '-count-table.csv', sep=''), row.names = FALSE)

rm(single_target_guides)
rm(multi_target_guides)
rm(target_guides_count)

# CONSTRUCTING GUIDESET

cleaned_alignment <- cleaned_alignment[!is.na(cleaned_alignment$pam_site),,drop=FALSE]

ids <- paste0('sgRNA_', seq_len(nrow(cleaned_alignment)))
sgRNA <- c(cleaned_alignment$sgRNA)
gene <- c(cleaned_alignment$gene)
df_sgRNA <- data.frame(ids, sgRNA, gene)
write_csv(df_sgRNA, file = 'scoring_guideset_sgrna.csv')

scoring_guideset <- GuideSet(ids = ids,
                             sgRNA = cleaned_alignment$sgRNA,
                             protospacers = cleaned_alignment$spacer,
                             pams = cleaned_alignment$pam,
                             pam_site = cleaned_alignment$pam_site,
                             seqnames = cleaned_alignment$chr,
                             strand = cleaned_alignment$strand,
                             CrisprNuclease = crisprNuclease,
                             bsgenome = bsgenome)

scoring_guideset$gene <- cleaned_alignment$gene

rm(ids)
gc()

# CHARACTERIZATION OF ON AND OFF TARGET EFFECT

data(txdb_human, package="crisprDesignData")
txObject <- txdb_human

scoring_guideset <- addSpacerAlignments(scoring_guideset,
                                        txObject=txObject,
                                        aligner='bowtie',
                                        aligner_index=bowtie_index,
                                        bsgenome=bsgenome,
                                        n_mismatches=2)


scoring_guideset <- addOffTargetScores(scoring_guideset)
scoring_guideset <- addOnTargetScores(scoring_guideset,
                                      methods= c('ruleset3'),
                                      tracrRNA = c('Hsu2013', 'Chen2013'))

rm(txdb_human)
rm(txObject)
rm(bsgenome)
rm(crisprNuclease)

scoring_guideset_table <- data.frame(scoring_guideset)
scoring_alignment <- data.frame(scoring_guideset$alignments)

saveRDS(scoring_guideset, file = "scoring_guideset.rds")
saveRDS(scoring_guideset_table, file = 'scoring_guideset_table.rds')
saveRDS(scoring_alignment, file = 'scoring_alignment.rds')

# Import RDS file to avoid run out of memory

library(tidyverse)
library(dplyr)
library(reshape)

scoring_guideset <- readRDS('scoring_guideset.rds')
scoring_guideset_table <- readRDS('scoring_guideset_table.rds')
scoring_alignment <- readRDS('scoring_alignment.rds')
scoring_guideset_sgrna <- read_csv('scoring_guideset_sgrna.csv')

scoring_alignment <- scoring_alignment %>%
  dplyr::rename('ids' = 'group_name')

dim(scoring_guideset_table)
# [1] 145932     20

dim(scoring_alignment)
# [1] 290454     23

# Join scoring alignment table with scoring_guideset_sgrna to indicate which sgRNA have off-target effect

scoring_alignment_sgrna <- left_join(scoring_alignment, scoring_guideset_sgrna, by = 'ids')
scoring_alignment_sgrna <- scoring_alignment_sgrna %>%
  select(ids, sgRNA, gene, spacer, score_cfd, score_mit, 
         n_mismatches, cds, fiveUTRs, threeUTRs, exons,
         introns, intergenic)


# analysis of score proportion based on sgrna_guideset_table

scoring_guideset_count <- scoring_guideset_table %>%
  select(score_cfd, score_mit) %>%
  melt() %>%
  dplyr::rename('score_category' = 'variable',
                'score' = 'value')

ggplot(scoring_guideset_count, aes(x = score_category, y = score, group = score_category)) +
  geom_boxplot() + 
  labs(x = 'Score Category', y = 'Score')

function_data_filter <- function(guideset_table_df, score, mismatch) {
  guideset_table_df <- guideset_table_df %>%
    filter()
  
  
  
  
}

scoring_guideset_test <- scoring_guideset_table %>%
  # select(sgRNA, n1, n2, score_cfd, score_mit, score_ruleset3) %>%
  select(sgRNA, n1, n2, score_cfd, score_mit) %>%
  # filter(score_cfd == 1 & score_ruleset3 > 0)
  # filter(n1 == 0 & n2 == 0)
  pivot_longer(c(n1, n2), names_to = 'mismatch_indicator', values_to = 'possible_alignment') %>%
  # pivot_longer(c(score_cfd, score_mit, score_ruleset3), names_to = 'scoring_method', values_to = 'score')
  pivot_longer(c(score_cfd, score_mit), names_to = 'scoring_method', values_to = 'score')

ggplot(scoring_guideset_test, aes(x = scoring_method, y = score, group = scoring_method)) +
  geom_boxplot() + 
  labs(x = 'Scoring Method', y = 'Score') + 
  scale_x_discrete(
    # breaks = c('score_cfd', 'score_mit', 'score_ruleset3'), labels = c('Score CFD', 'Score MIT', 'Score Ruleset 3'))
    breaks = c('score_cfd', 'score_mit'), labels = c('Score CFD', 'Score MIT')) +
  theme(
    axis.title.x = element_text(size = 0),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25))


scoring_guideset_cfd <- scoring_guideset_test %>%
  filter(scoring_method == 'score_cfd' & score == 1)
scoring_guideset_mit <- scoring_guideset_test %>%
  filter(scoring_method == 'score_mit' & score == 1)
scoring_guideset_ruleset3 <- scoring_guideset_test %>%
  filter(scoring_method == 'score_ruleset3')

ggplot(scoring_guideset_cfd, aes(x = mismatch_indicator, y = score, group = mismatch_indicator)) +
  geom_boxplot() + 
  labs(x = 'Mismatch Indicator', y = 'CFD Score') + 
  scale_x_discrete(
    breaks = c('n1', 'n2'), labels = c('1 Mismatch', '2 Mismatch'))
  

ggplot(scoring_guideset_mit, aes(x = mismatch_indicator, y = score, group = mismatch_indicator)) +
  geom_boxplot() + 
  scale_x_discrete(
    breaks = c('n1', 'n2'), labels = c('1 Mismatch', '2 Mismatch'))













off_target_sgrna <- scoring_alignment_sgrna %>%
  filter(n_mismatches > 0 & !is.na(exons)) %>%
  arrange(gene) %>%
  mutate(sl_indicator = paste(gene, exons, sep = '_'))

# finding synthetic lethality interaction

sl_data <- read_csv('./raw_data/Table_S8.csv')
sl_data_interaction <- c(sl_data$sorted_gene_pair)
sl_interaction <- off_target_sgrna %>%
  filter(off_target_sgrna$sl_indicator %in% sl_data_interaction)

print(paste('number of unique sgRNAs: ', length(unique(sl_interaction$sgRNA))))
print(paste('number of unique genes with synthetic lethality pair: ', length(unique(sl_interaction$gene))))

write.csv(sl_interaction, file = 'Cellecta_off-target_SL_interaction.csv')

# checking guide abundance in predicted off-target and synthetic lethal genes

sl_interation_genes <- c(unique(sl_interaction$sgRNA))
sgRNA_read_counts <- read_csv('raw_data/All_sgRNA_Count.csv')
sgRNA_read_counts_test <- sgRNA_read_counts %>%
  filter(Gene %in% c('CST1','CST2')) %>%
  mutate(
    `DMSO/CONTROL` = round(DAY20_DMSO_FD/DAY0_CNTR_A0, 3),
  )


write_csv(sgRNA_read_counts, 'sgRNA_Cellecta_off-target_SL_interaction.csv')




print(length(unique(scoring_alignment_sgrna$gene)))

FBWX7_sl_genes <- read_csv('./raw_data/FBWX7_sl_genes.txt')

FBWX7_sl_interaction <- sl_interaction %>%
  filter(gene %in% FBWX7_sl_genes$FBWX7_sl_genes)

print(paste('number of unique sgRNAs synthetic lethal to FBWX7: ', length(unique(FBWX7_sl_interaction$sgRNA))))
print(paste('number of unique genes with synthetic lethality pair to FBWX7: ', length(unique(FBWX7_sl_interaction$gene))))





scoring_guideset_table$flank5 <- 'ACCG'
scoring_guideset_table$flank3 <- 'GTT'

# Manually input rule set 1 score

function_score_ruleset1 <- function(flank5, spacer, pam, flank3) {
  input  <- paste0(flank5, spacer, pam, flank3) 
  results <- getRuleSet1Scores(input)
  return(as.character(results$score))
}

score_ruleset1 <- c()

for (i in seq_along(scoring_guideset_table$sgRNA)) {
  value <- function_score_ruleset1(scoring_guideset_table$flank5[i], scoring_guideset_table$protospacer[i], scoring_guideset_table$pam[i], scoring_guideset_table$flank3[i])
  score_ruleset1 <- c(score_ruleset1, value)
}

scoring_guideset_table$score_ruleset1 <- score_ruleset1

# cross compare the score of cfd, mit, and rule set 1

scoring_guideset_test <- scoring_guideset_table %>%
  select(sgRNA, gene, score_cfd, score_mit, score_ruleset1) %>%
  filter(score_cfd == 1)

ggplot(scoring_guideset_test, aes(x = variable, y = value)) + 
  geom_boxplot(notch=TRUE)

write.csv(scoring_guideset_table, file = paste(sgRNA_library_selected, '_scoring_guideset.csv', sep=''), row.names = FALSE)



scoring_guideset_score <- scoring_guideset_table %>%
  select(sgRNA, gene, alignments, score_cfd, score_mit, seqnames)

write.csv(scoring_guideset_score, file = paste(sgRNA_library_selected, '_FULL_score.csv', sep=''), row.names = FALSE)














