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
prg1_turbo_id = prg1_turbo_id %>%  mutate(tag = gsub("CELE_", "", Alternate.ID)) 
prg1_turbo_id = prg1_turbo_id %>% select(tag, prg1_lfc, prg1_pv)

dat = read_delim("./af3_output/prg1_hits_from_AF3_screen/prg1_hits_from_screen_against_proteome.tsv", delim = "\t")

pathways = read_delim("/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/pathways.tsv", delim = "\t", col_names = c("gene", "biotype", "class"))
pathways = pathways %>% 
  separate(gene, into = c("geneid", "seqid", "locusid"), sep = ",") %>% 
  mutate(tag = ifelse(is.na(locusid), seqid, locusid))

dat = dat %>% 
  mutate(sample = gsub("_data", "", sample)) %>% 
  rowwise() %>% 
  mutate(bait = strsplit(sample, "_")[[1]][1]) %>% 
  mutate(tag = strsplit(sample, "_")[[1]][2]) %>% 
  mutate(tag = gsub("CELE_", "", tag)) %>% 
  ungroup()

head(dat)

dat_grouped = dat %>% 
  group_by(bait, tag) %>% 
  summarise(MD = mean(pae_min), CV = ( sd(pae_min)/mean(pae_min) )*100) %>% 
  ungroup() 


dat_grouped = dat_grouped %>% 
  left_join(prg1_turbo_id, by = c("tag")) %>% 
  arrange(tag)

write.table(dat_grouped, file = "prg1_screen_reciporical_hit_screen.tsv", quote = F, sep = "\t", row.names = F, col.names = T)




