library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(dplyr)
library(scales)
library(ggpubr)
library(tidyr)
library(readr)
library(rstatix)
library(gghalves)
library(ggbeeswarm)
source("/fs/ess/PCON0160/ben/bin/mighty.R")

install_and_load <- function(pkgs, repos = getOption("repos")) {
  # Ensure repos is set (avoid empty CRAN mirror)
  if (length(repos) == 0 || repos == "@CRAN@") {
    repos <- "https://cloud.r-project.org"
  }
  
  for (pkg in pkgs) {
    # Check if package is already installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing package '%s'...", pkg))
      install.packages(pkg, repos = repos)
    }
    
    # Load it, suppressing startup messages
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }
  
  invisible(NULL)
}

packages = c("ggplot2", "RColorBrewer", "ggrepel", "dplyr", "tidyr", "scales", "readr", "rstatix", "gghalves", "ggbeeswarm", "ggpubr")
install_and_load(packages)


prg1_turbo_id = read_delim("/fs/ess/PCON0160/ben/projects/2024_disl2/2025_DISL2_degrades_piRNAs/prg1_TURBOID_analysis.tsv", delim = "\t")
prg1_turbo_id %>% 
  head()
enriched = prg1_turbo_id %>% filter(prg1_pv < 0.05 & prg1_lfc >= 1)

dat = read_delim("./af3_output/prg1_screen/20250518/merged_pae_scores.tsv", delim = "\t")
dat = read_delim("/fs/ess/PCON0160/ben/pipelines/alphafold3_pipeline/af3_output/prg1_screen/merged_pae_scores.tsv", col_names = c("pae_min", "sample"), skip = 1)

pathways = read_delim("/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/pathways.tsv", delim = "\t", col_names = c("gene", "biotype", "class"))
pathways = pathways %>% 
  separate(gene, into = c("geneid", "seqid", "locusid"), sep = ",") %>% 
  mutate(tag = ifelse(is.na(locusid), seqid, locusid))

pathways

dat = dat %>% 
  mutate(sample = gsub("_data", "", sample)) %>% 
  rowwise() %>% 
  mutate(bait = strsplit(sample, "_")[[1]][1]) %>% 
  mutate(tag = strsplit(sample, "PRG1_")[[1]][2]) %>% 
  mutate(tag = gsub("CELE_", "", tag)) %>% 
  ungroup()

head(dat)


dat_grouped = dat %>% 
  group_by(bait, tag) %>% 
  summarise(MD = mean(pae_min), CV = ( sd(pae_min)/mean(pae_min) )*100) %>% 
  ungroup() %>% 
  mutate(turboID = ifelse(tag %in% enriched$Alternate.ID & MD < 5 & CV < 10, 3,
                          ifelse(tag %in% enriched$Alternate.ID & MD >= 5, 1, 
                                 ifelse(MD < 5 & CV < 10, 2, 0)))) %>% 
  mutate(is_in_turboID = ifelse(tag %in% enriched$Alternate.ID, T, F)) %>% 
  left_join(pathways, by = c("tag"))


write.table(dat_grouped, file = "prg1_screen_annotated_germline.tsv", row.names = F, col.names = T, quote = F, sep = "\t")

dat_grouped['turboID_f'] = factor(dat_grouped$turboID, levels = c(3, 2, 1, 0))

cols = c("1" = "magenta",
         "2" = "green", 
         "3" = "blue",
         "0" = "grey")

labs = c("Enriched in TurboID::PRG-1 & pae < 5 & CV < 10", 
         "pae < 5 & CV < 10",
           "Enriched in TurboID::PRG-1", 
           "Other")

dat_grouped %>% filter(grepl("Y54G11A.3", tag))

ggplot(data = dat_grouped %>% arrange(desc(turboID_f)), aes(x = MD, y = CV)) +
  geom_point(size = 1, shape = 19, color = "grey50", alpha = 0.5) + 
  geom_point(data = dat_grouped %>% filter(grepl("Y54G11A.3", tag)), aes(x = MD, y = CV), color = "magenta", size = 2) + 
  geom_label_repel(data = dat_grouped %>% filter(grepl("Y54G11A.3", tag)), aes(x = MD, y = CV, label = tag), color = "magenta", box.padding = 1.4) + 
  my_theme() + 
  labs(x = "median minimum pae score\n(average of 5 models)", y = "coefficient of variantion")

dat_grouped %>% 
  filter(MD < 5 & CV < 10)
dat_grouped %>% filter(MD < 5 & CV < 10) %>% arrange(desc(MD))

dat_grouped %>% filter(tag == "Y57A10A.31")

gene_ids = read_delim("/users/PAS1473/benpasto1/PCON/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.gene_ids.txt")

known_e3 <- data.frame(gene_id = c(
  "T01C3.3", "C12C8.3", "Y105E8A.14", "Y55F3AM.6", "F56D2.2", "C06A5.8", "Y71F9AL.10",
  "T24D1.2", "K12B6.8", "Y47G6A.14", "C53A5.6", "K01G5.1", "H05L14.2", "C01G6.4",
  "ZK287.5", "EEED8.16", "T05H10.5", "C06A5.9", "Y51F10.2", "Y52E8A.2", "C15F1.5",
  "T24D1.5", "K04C2.4", "C55A6.1", "F26F4.7", "B0393.6", "C32D5.11", "F35G12.9",
  "Y45G12B.2", "C17E4.3", "R10A10.2", "T20F5.7", "C36A4.8", "D2085.4", "C49H3.5",
  "R06F6.2", "B0432.13", "C32D5.10", "K02B12.8", "C34E10.4", "C56A3.4", "C09E7.9",
  "T10F2.4", "Y65B4BR.4"
))
length(known_e3$gene_id)

gene_ids_known_e3 = gene_ids %>% filter(seq_id %in% known_e3$gene_id)

gene_ids_known_e3 %>% filter(grepl("rike", locus_id))

dat_grouped = dat_grouped %>% 
  mutate(germline_e3 = ifelse(tag %in% gene_ids_known_e3$seq_id, 1, 
                              ifelse(tag %in% gene_ids_known_e3$locus_id, 1, 0)))

gene_ids_known_e3 %>% filter(!seq_id %in% dat_grouped$tag) %>% filter(!locus_id %in% dat_grouped$tag)


dat_grouped['germline_expressing_e3'] = factor(dat_grouped$germline_e3, levels = c(1, 0))

dat_grouped %>% filter(grepl("rike", tag))

cols = c("1" = "magenta", "0" = "grey70")

labs = c("germline expressing E3", "Other")

ggplot(data = dat_grouped %>% arrange(desc(germline_expressing_e3)), aes(x = MD, y = CV, color = germline_expressing_e3)) +
  geom_point(size = 1, shape = 21) + 
  geom_label_repel(data = dat_grouped %>% filter(germline_expressing_e3 == 1 & CV < 10 & MD < 5), aes(label = tag), box.padding = 1, size = 3) + 
  my_theme() + 
  labs(x = "median minimum pae score\n(average of 5 models)", y = "coefficient of variantion") + 
  scale_color_manual(values = cols, labels = labs)

dat_grouped %>% filter(germline_expressing_e3 == 1 & MD < 5)

write.table(dat_grouped, "prg1_screen_annotated.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

dat_grouped %>% filter(tag == "marc-3")
# Aringher 
# marc-3 = C17E4.3 = F05
# zhp-3 = B01

###################################
# CSR1


gene_ids = read_delim("/users/PAS1473/benpasto1/PCON/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.gene_ids.txt")

known_e3 <- data.frame(gene_id = c(
  "T01C3.3", "C12C8.3", "Y105E8A.14", "Y55F3AM.6", "F56D2.2", "C06A5.8", "Y71F9AL.10",
  "T24D1.2", "K12B6.8", "Y47G6A.14", "C53A5.6", "K01G5.1", "H05L14.2", "C01G6.4",
  "ZK287.5", "EEED8.16", "T05H10.5", "C06A5.9", "Y51F10.2", "Y52E8A.2", "C15F1.5",
  "T24D1.5", "K04C2.4", "C55A6.1", "F26F4.7", "B0393.6", "C32D5.11", "F35G12.9",
  "Y45G12B.2", "C17E4.3", "R10A10.2", "T20F5.7", "C36A4.8", "D2085.4", "C49H3.5",
  "R06F6.2", "B0432.13", "C32D5.10", "K02B12.8", "C34E10.4", "C56A3.4", "C09E7.9",
  "T10F2.4", "Y65B4BR.4"
))
length(known_e3$gene_id)

gene_ids_known_e3 = gene_ids %>% filter(seq_id %in% known_e3$gene_id)




dat = read_delim("/fs/ess/PCON0160/ben/pipelines/alphafold3_pipeline/af3_output/csr1_screen/merged_pae_scores.tsv", col_names = c("pae_min", "sample"), skip = 1)

dat = dat %>% 
  mutate(sample = gsub("_data", "", sample)) %>% 
  rowwise() %>% 
  mutate(bait = strsplit(sample, "_")[[1]][1]) %>% 
  mutate(tag = strsplit(sample, "CSR1_")[[1]][2]) %>% 
  mutate(tag = gsub("CELE_", "", tag)) %>% 
  ungroup()

head(dat)

dat_grouped = dat %>% 
  group_by(bait, tag) %>% 
  summarise(MD = mean(pae_min), CV = ( sd(pae_min)/mean(pae_min) )*100) %>% 
  ungroup() %>% 
  mutate(turboID = ifelse(MD < 5 & CV < 10,1,0))

write.table(dat_grouped, "csr1_screen_annotated.tsv", sep = "\t", col.names = T, row.names = F, quote = F)


dat_grouped['turboID_f'] = factor(dat_grouped$turboID, levels = c(1,0))

cols = c("1" = "magenta",
         "0" = "grey")

labs = c("pae < 5 & CV < 10",
         "Other")

ggplot(data = dat_grouped %>% arrange(desc(turboID_f)), aes(x = MD, y = CV, color = turboID_f)) +
  geom_point(size = 1, shape = 21) + 
  my_theme() + 
  labs(x = "median minimum pae score\n(average of 5 models)", y = "coefficient of variantion") + 
  scale_color_manual(values = cols, labels = labs)

dat_grouped %>% 
  filter(MD < 5 & CV < 10)

dat_grouped %>% filter(MD < 5 & CV < 10) %>% arrange((MD))


#
dat_grouped = dat_grouped %>% 
  mutate(germline_e3 = ifelse(tag %in% gene_ids_known_e3$seq_id, 1, 
                              ifelse(tag %in% gene_ids_known_e3$locus_id, 1, 0)))

gene_ids_known_e3 %>% filter(!seq_id %in% dat_grouped$tag) %>% filter(!locus_id %in% dat_grouped$tag)


dat_grouped['germline_expressing_e3'] = factor(dat_grouped$germline_e3, levels = c(1, 0))

dat_grouped %>% filter(grepl("rike", tag))

cols = c("1" = "magenta", "0" = "grey70")

labs = c("germline expressing E3", "Other")

ggplot(data = dat_grouped %>% arrange(desc(germline_expressing_e3)), aes(x = MD, y = CV, color = germline_expressing_e3)) +
  geom_point(size = 1, shape = 21) + 
  geom_label_repel(data = dat_grouped %>% filter(germline_expressing_e3 == 1 & CV < 10 & MD < 6), aes(label = tag), box.padding = 1, size = 3) + 
  my_theme() + 
  labs(x = "median minimum pae score\n(average of 5 models)", y = "coefficient of variantion") + 
  scale_color_manual(values = cols, labels = labs)


dat_grouped %>% filter(germline_expressing_e3 == 1 & CV < 10 & MD < 8)


###################################
# Y57A10A.13
dat = read_delim("/fs/ess/PCON0160/ben/pipelines/alphafold3_pipeline/af3_output/Y57A10A.13_screen/merged_pae_scores.tsv", col_names = c("pae_min", "sample"), skip = 1)

dat = dat %>% 
  mutate(sample = gsub("_data", "", sample)) %>% 
  rowwise() %>% 
  mutate(bait = strsplit(sample, "_")[[1]][1]) %>% 
  mutate(tag = strsplit(sample, "Y57A10A.13_")[[1]][2]) %>% 
  mutate(tag = gsub("CELE_", "", tag)) %>% 
  ungroup()

head(dat)

dat_grouped = dat %>% 
  group_by(bait, tag) %>% 
  summarise(MD = mean(pae_min), CV = ( sd(pae_min)/mean(pae_min) )*100) %>% 
  ungroup() %>% 
  mutate(turboID = ifelse(MD < 5 & CV < 10,1,0))

dat_grouped['turboID_f'] = factor(dat_grouped$turboID, levels = c(1,0))

cols = c("1" = "magenta",
         "0" = "grey")

labs = c("pae < 5 & CV < 10",
         "Other")

ggplot(data = dat_grouped %>% arrange(desc(turboID_f)), aes(x = MD, y = CV, color = turboID_f)) +
  geom_point(size = 1, shape = 21) + 
  my_theme() + 
  labs(x = "median minimum pae score\n(average of 5 models)", y = "coefficient of variantion") + 
  scale_color_manual(values = cols, labels = labs)

dat_grouped %>% 
  filter(MD < 5 & CV < 10)

dat_grouped %>% filter(MD < 5 & CV < 10) %>% arrange((MD)) %>% view()

write.table(dat_grouped, file = "Y57A10A13_screen.tsv", sep = "\t", col.names = T, row.names = F, quote = F)















