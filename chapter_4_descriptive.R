## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
##
## Script name: NEKSUS Descriptive molecular epidemiology 
##
## Purpose of script: Descriptive molecular epidemiology and visualisations for NEKSUS dataset
##                    including phylogenetic trees, metadata assocaitions, and plasmids.
##
## Author: Dorottya Nagy
##
## Date Created: 17/10/2026
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

# Set working directory
setwd("~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/")

# load/install packages
#install.packages("effectsize")

# Load required packages
library(ape)
library(colorspace)
library(data.table)
library(dplyr)
library(effectsize)
library(ggnewscale)
library(ggplot2)
library(ggtree)
library(gridExtra)
library(igraph)
library(lubridate)
library(mclust)
library(patchwork)
library(phytools)
library(purrr)
library(RColorBrewer)
library(readr)
library(readxl)
library(scales)
library(stringr)
library(tidyr)
library(zoo)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD METADATA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load metadata
samples_assembled_with_species_hcai_deduplicated <- read.csv("samples_assembled_with_species_hcai_deduplicated.csv")
#colnames(samples_assembled_with_species_hcai_deduplicated) # 2695
#nrow(samples_assembled_with_species_hcai_deduplicated) # 2695
#View(samples_assembled_with_species_hcai_deduplicated) # 2695

# load data
all_metadata_deduplicated <- read.csv("all_metadata_deduplicated.csv")
#View(all_metadata_deduplicated)

# load amrfinder data
amrfinder_metadata <- read.csv("amrfinder_metadata_with_NAs_updated.csv")
#nrow(amrfinder_metadata)
#colnames(amrfinder_metadata)
amrfinder_metadata_deduplicated <- amrfinder_metadata |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  mutate(seqid_corrected = case_when(run == "run0" ~ paste0("pilot_", sample),
                                     run != "run0" & sample != "AB90" & sample != "AF210" ~ gsub("AH", "AH-", sample),
                                     run != "run0" & sample == "AB90" ~ "AB90-A",
                                     run != "run0" & sample == "AF210" ~ "AF210-A",
                                     TRUE ~ sample),
         tip_label = paste0(run, "_", seqid_corrected)) |>
  filter(!is.na(Element.symbol)) |>
  filter(Type == "AMR")
amrfinder_metadata_deduplicated$region <- as.factor(amrfinder_metadata_deduplicated$region)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. MLST composition over time ####
# pre- and post- COVID-19 pandemic
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * E. coli BSIs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
denom <- nrow(all_metadata_deduplicated[all_metadata_deduplicated$species_complex == "E. coli" & all_metadata_deduplicated$sampletype_cat == "blood",])
top_ecoli_mlst <- all_metadata_deduplicated |>
  filter(species_complex == "E. coli" & sampletype_cat == "blood") |>
  group_by(kleborate_mlst_ST) |>
  summarise(count = n(),
            freq = count/denom * 100) |>
  arrange(desc(freq))
#View(top_ecoli_mlst)

top_mlst_names_ecoli <- c("ST131", "ST73", "ST69", "ST95", "ST12", "ST127", "ST1193", "ST10", "other")
ecoli_cols <- c(
  "ST131" = "#4477AA",
  "ST73"  = "#EE6677",
  "ST69"  = "#228833",
  "ST95"  = "#66CCEE",
  "ST12"  = "#AA3377",
  "ST127" = "#EE7733",
  "ST1193"= "#CCBB44",
  "ST10"  = "#44AA99",
  "other" = "#E6E6E6"
)

# make combined df of pre-and mpost-pandemic
bsac_counts_ecoli <- c(218, 222, 82, 145, 57, 37, 1, 33, 532)
oxf_counts_ecoli <- c(455, 567, 312, 317, 144, 113, 33, 53, 1467)
neksus_counts_ecoli <- c(251, 218, 114, 94, 76, 36, 36, 27, 619)
# add ST1193, ST127, ST10


names(neksus_counts_ecoli) <- top_mlst_names_ecoli
names(bsac_counts_ecoli) <- top_mlst_names_ecoli
names(oxf_counts_ecoli) <- top_mlst_names_ecoli
neksus_vs_bsac_ecoli <- rbind(neksus_counts_ecoli, bsac_counts_ecoli)
neksus_vs_oxf_ecoli <- rbind(neksus_counts_ecoli, oxf_counts_ecoli)
# check percentages
bsac_pct <- bsac_counts_ecoli / 1327 *100
oxf_pct <- oxf_counts_ecoli / 3461 *100
neksus_pct <- neksus_counts_ecoli / 1471 *100

# bind df
ecoli_df <- tibble(
  Dataset = c("BSAC\n(2003-2012)","Oxfordshire\n(2008-2018)","NEKSUS\n(2023-2024)"),
  ST131 = c(bsac_counts_ecoli["ST131"], oxf_counts_ecoli["ST131"], neksus_counts_ecoli["ST131"]),
  ST73  = c(bsac_counts_ecoli["ST73"],  oxf_counts_ecoli["ST73"],  neksus_counts_ecoli["ST73"]),
  ST69  = c(bsac_counts_ecoli["ST69"],  oxf_counts_ecoli["ST69"],  neksus_counts_ecoli["ST69"]),
  ST95  = c(bsac_counts_ecoli["ST95"],  oxf_counts_ecoli["ST95"],  neksus_counts_ecoli["ST95"]),
  ST12  = c(bsac_counts_ecoli["ST12"],  oxf_counts_ecoli["ST12"],  neksus_counts_ecoli["ST12"]),
  ST127  = c(bsac_counts_ecoli["ST127"],  oxf_counts_ecoli["ST127"],  neksus_counts_ecoli["ST127"]),
  ST1193  = c(bsac_counts_ecoli["ST1193"],  oxf_counts_ecoli["ST1193"],  neksus_counts_ecoli["ST1193"]),
  ST10  = c(bsac_counts_ecoli["ST10"],  oxf_counts_ecoli["ST10"],  neksus_counts_ecoli["ST10"]),
  other = c(bsac_counts_ecoli["other"], oxf_counts_ecoli["other"], neksus_counts_ecoli["other"])
) |>
  pivot_longer(
    -Dataset,
    names_to="ST",
    values_to="Count"
  ) |>
  group_by(Dataset) |>
  mutate(
    Prop = Count/sum(Count),
    ST = factor(ST, levels = rev(top_mlst_names_ecoli)),
    Dataset = factor(Dataset, levels = c("BSAC\n(2003-2012)","Oxfordshire\n(2008-2018)","NEKSUS\n(2023-2024)"))
  )

## Plot
p_ecoli <- ggplot(ecoli_df,  aes(Dataset,  Prop, fill = ST)) +
  geom_col(width=.75) +
  geom_vline(  xintercept=2.5, linetype="dashed" ) +
  scale_fill_manual(values = ecoli_cols) +
  scale_y_continuous( labels=scales::percent) +
  labs(title="Escherichia coli",
       x=NULL,
       y="Proportion"
  ) +
  guides(fill = guide_legend(ncol = 6, byrow = TRUE)) +
  theme_minimal() +
  theme(
    panel.grid=element_blank(),
    plot.title=element_text(face="bold"),
    legend.position="bottom"
  )
p_ecoli

#~~~~~~~~#
# chi-squared test
# helper function
# Calculate enrichment, adjusted residuals, and Bonferroni-corrected p-values
get_chi_enrichment <- function(chi_result) {
  
  observed <- chi_result$observed
  expected <- chi_result$expected
  
  # Observed / expected ratio
  oe_ratio <- observed / expected
  
  # Adjusted Pearson residuals
  adjusted_residuals <- chi_result$stdres
  
  # Two-sided p-values from standard normal distribution
  p_values <- 2 * pnorm(-abs(adjusted_residuals))
  
  # Bonferroni correction across ST categories
  # Only test NEKSUS row because the two rows are complementary
  p_adjusted <- p.adjust(
    p_values["neksus_counts_ecoli", ],
    method = "bonferroni"
  )
  
  # Return results for NEKSUS
  results <- data.frame(
    ST = colnames(observed),
    observed = observed["neksus_counts_ecoli", ],
    expected = expected["neksus_counts_ecoli", ],
    OE_ratio = oe_ratio["neksus_counts_ecoli", ],
    adjusted_residual = adjusted_residuals["neksus_counts_ecoli", ],
    p_value = p_values["neksus_counts_ecoli", ],
    p_bonferroni = p_adjusted,
    significant = p_adjusted < 0.05
  )
  
  # Identify enrichment or depletion
  results$direction <- ifelse(
    results$OE_ratio > 1,
    "Enriched",
    "Depleted"
  )
  
  return(results)
}

# run
# NEKSUS vs BSAC
neksus_bsac_chi <- chisq.test(neksus_vs_bsac_ecoli, simulate.p.value = TRUE,B = 100000)
bsac_enrichment <- get_chi_enrichment(neksus_bsac_chi)
print(bsac_enrichment)


# NEKSUS vs Oxford
neksus_oxf_chi <- chisq.test(neksus_vs_oxf_ecoli, simulate.p.value = TRUE,B = 100000)
oxf_enrichment <- get_chi_enrichment(neksus_oxf_chi)
print(oxf_enrichment)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Klb BSIs ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
denom <- nrow(all_metadata_deduplicated[grepl("K. ", all_metadata_deduplicated$species_complex) & all_metadata_deduplicated$sampletype_cat == "blood",])
top_kleb_mlst <- all_metadata_deduplicated |>
  filter(grepl("K. ", species_complex) & sampletype_cat == "blood") |>
  group_by(kleborate_mlst_ST) |>
  summarise(count = n(),
            freq = count/denom * 100) |>
  arrange(desc(freq))
#View(top_kleb_mlst)


# load Lipworth et al Oxfordshire Kleb data
kleborate_oxford <- read.delim("~/DPhil_Clin_Medicine/DPhil/NEKSUS/main_pipeline_v2/kleborate_oxford.txt")
#View(kleborate_oxford)
oxfordshire_denom <- nrow(kleborate_oxford) # 886
kleborate_oxford_mlsts <- kleborate_oxford |>
  mutate(mlst = paste0(species, "|", ST)) |>
  group_by(mlst) |>
  summarise(count = n()) |>
  mutate(freq = count / oxfordshire_denom) |>
  arrange(desc(freq))
#View(kleborate_oxford_mlsts)

top_mlst_names_kleb <- c("ST45", "ST14", "ST29", "ST17", "ST307", "ST54", "ST23", "ST20", "ST35", "ST37", "ST490", "other")
kleb_cols <- c(
  "ST45"  = "#4477AA",
  "ST14"  = "#EE6677",
  "ST29"  = "#228833",
  "ST23"  = "#CCBB44",
  "ST17"  = "#AA3377",
  "ST307" = "#EE7733",
  "ST54"  = "#44AA99",
  "ST20"  = "#999933",
  "ST35"  = "#CC6677",
  "ST37"  = "#882255",
  "ST490" = "#332288",
  "other" = "#E6E6E6"
)

neksus_counts_kleb <- c(15,12,10,8,8,8,8,6,5,4,0,384)
oxf_counts_kleb <- c(19,17,13,11,6,4,1,13,11,19,45,727)
names(neksus_counts_kleb) <- top_mlst_names_kleb
names(oxf_counts_kleb) <- top_mlst_names_kleb
# bind rows
neksus_vs_oxf_kleb <- rbind(neksus_counts_kleb, oxf_counts_kleb)

# make df
kleb_df <- tibble(
    Dataset=c("BSAC\n(2003-2012)","Oxfordshire\n(2008-2018)","NEKSUS\n(2023-2024)"),
    ST45  = c(NA,19,15),
    ST14  = c(NA,17,12),
    ST29  = c(NA,13,10),
    ST17  = c(NA,11,8),
    ST307 = c(NA,6,8),
    ST54  = c(NA,4,8),
    ST23  = c(NA,1,9),
    ST20  = c(NA,13,7),
    ST35  = c(NA,11,6),
    ST37  = c(NA,19,4),
    ST490 = c(NA,45,0),
    other = c(NA,727,384)
  ) |>
  pivot_longer(
    -Dataset,
    names_to="ST",
    values_to="Count"
  ) |>
  group_by(Dataset) |>
  mutate(
    Prop = Count/sum(Count),
    ST = factor(ST, levels = rev(top_mlst_names_kleb)),
    Dataset=factor(Dataset, levels = c("BSAC\n(2003-2012)","Oxfordshire\n(2008-2018)","NEKSUS\n(2023-2024)"))
  )

# plot kleb
p_kleb <-  ggplot(kleb_df, aes(Dataset, Prop, fill = ST)) +
  geom_col(width=.75, na.rm=TRUE) +
  geom_vline(xintercept=2.5, linetype="dashed") +
  scale_fill_manual(values = kleb_cols) +
  scale_y_continuous( labels=scales::percent) +
  labs(
    title="Klebsiella spp.",
    x=NULL,
    y="Proportion"
  ) +
  guides(fill = guide_legend(ncol = 6, byrow = TRUE)) +
  theme_minimal() +
  theme(
    panel.grid=element_blank(),
    plot.title=element_text(face="bold"),
    legend.position="bottom"
  )
p_kleb
# combined plot
combined_pre_post_pandemic_plot <- (p_ecoli / p_kleb) +
  plot_layout(guides = "keep")
combined_pre_post_pandemic_plot
ggsave("combined_pre_post_pandemic_plot.png",  combined_pre_post_pandemic_plot, units = "in", width = 8, height = 7, dpi = 300)

#~~~~~~~~~~~~~#
# chi-squared test
set.seed(123)
# function
get_chi_enrichment <- function(chi_result) {
  
  observed <- chi_result$observed
  expected <- chi_result$expected
  
  # Observed / expected ratio
  oe_ratio <- observed / expected
  
  # Adjusted Pearson residuals
  adjusted_residuals <- chi_result$stdres
  
  # Two-sided p-values from standard normal distribution
  p_values <- 2 * pnorm(-abs(adjusted_residuals))
  
  # Bonferroni correction across ST categories
  # Only test NEKSUS row because the two rows are complementary
  p_adjusted <- p.adjust(
    p_values["neksus_counts_kleb", ],
    method = "bonferroni"
  )
  
  # Return results for NEKSUS
  results <- data.frame(
    ST = colnames(observed),
    observed = observed["neksus_counts_kleb", ],
    expected = expected["neksus_counts_kleb", ],
    OE_ratio = oe_ratio["neksus_counts_kleb", ],
    adjusted_residual = adjusted_residuals["neksus_counts_kleb", ],
    p_value = p_values["neksus_counts_kleb", ],
    p_bonferroni = p_adjusted,
    significant = p_adjusted < 0.05
  )
  
  # Identify enrichment or depletion
  results$direction <- ifelse(
    results$OE_ratio > 1,
    "Enriched",
    "Depleted"
  )
  
  return(results)
}

# NEKSUS vs Oxf for kleb
neksus_oxf_chi_kleb <- chisq.test(neksus_vs_oxf_kleb, simulate.p.value = TRUE,B = 100000)
oxf_enrichment_keb <- get_chi_enrichment(neksus_oxf_chi_kleb)
print(oxf_enrichment_keb)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * ST23 cluster investigate ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
write.csv(kleborate_kpsc_cleaned, "kleborate_kpsc_cleaned.csv", row.names = FALSE)

kpsc_st23 <- kleborate_kpsc_cleaned |>
  filter(klebsiella_pneumo_complex__mlst__ST == "ST23")
nrow(kpsc_st23)
kpsc_st23_seqids <- c("AJ251" ,  "AB121" ,  "AJ157",   "302145L", "AB10",    "AB8"  ,   "AF65" ,   "AF66"  ,  "AF110"  )

# load amrfinder metadata
amrfinder_metadata_autocycler_contigs <- read.csv("amrfinder_metadata_autocycler_contigs.csv")

kpsc_st23_amr <- amrfinder_metadata_autocycler_contigs |>
  filter(sample %in% kpsc_st23$sample)
View(kpsc_st23_amr)

# summarise
kpsc_st23_by_contig <- kpsc_st23_amr |>
  group_by(sample, region, sampletype_cat, samplecollectiondate,  autocycler_contig_topology, autocycler_contig_id, autocycler_unitig_length, Type) |>
  arrange(Element.symbol) |>
  summarise(genes = paste(Element.symbol, collapse = ",")) |>
  filter(autocycler_unitig_length > 600) |>
  # make wide
  pivot_wider(id_cols = c(sample, region, sampletype_cat, samplecollectiondate, autocycler_contig_topology, autocycler_contig_id, autocycler_unitig_length),
              names_from = Type, 
              values_from = genes )

#View(kpsc_st23_by_contig)

# load autocycler contig data 
#autocycler_contigs_mobsuite_genomad <- read.csv("autocycler_contigs_mobsuite_genomad_with_twice_seq.csv")
#colnames(autocycler_contigs_mobsuite_genomad)
plasmid_subcommunities_meta <- read.csv("plasmid_subcommunities_meta.csv")
colnames(plasmid_subcommunities_meta)

plasmid_subcommunities_meta_seected <- plasmid_subcommunities_meta |>
  mutate(primary_secondary_cluster = paste0(combo_primary_cluster_id, ":", combo_secondary_cluster_id)) |>
  select(plasmid, pling_subcommunity, num_plasmids_in_subcomm, molecule_type_genomad_majority, combo_rep_types, primary_secondary_cluster )

# left join
kpsc_st23_by_contig_plasmids <- kpsc_st23_by_contig |>
  mutate(plasmid = paste0(sample, "_", autocycler_contig_id)) |>
  left_join(plasmid_subcommunities_meta_seected, by = c("plasmid" = "plasmid")) |>
  ungroup() |>
  mutate(molecule_type_genomad_majority = case_when(is.na(molecule_type_genomad_majority) ~ "chromosome",
                                                    TRUE ~ molecule_type_genomad_majority)) |>
  select(sample, region, sampletype_cat, samplecollectiondate, autocycler_unitig_length, autocycler_contig_topology, 
         molecule_type_genomad_majority, combo_rep_types, primary_secondary_cluster, pling_subcommunity, AMR, STRESS, VIRULENCE)|>
  arrange(samplecollectiondate, sample, desc(autocycler_unitig_length))
#View(kpsc_st23_by_contig_plasmids)
colnames(kpsc_st23_by_contig_plasmids)

#save as df
write.csv(kpsc_st23_by_contig_plasmids, "kpsc_st23_by_contig_plasmids.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Phylogenetic Trees ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data
all_ecoli_metadata_deduplicated <- read.csv("all_ecoli_metadata_deduplicated.csv")
all_kleb_metadata_deduplicated <- read.csv("all_kleb_metadata_deduplicated.csv")
all_citro_metadata_deduplicated <- read.csv("all_citro_metadata_deduplicated.csv")
all_entero_metadata_deduplicated <- read.csv("all_entero_metadata_deduplicated.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * E. coli phylo tree ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data <- all_ecoli_metadata_deduplicated
#nrow(data) #1627
#colnames(data)
#View(data)
# remove ecoli| prefix from kleborate mlst
data <- data |>
  mutate(kleborate_mlst_ST = gsub("ecoli\\|", "", kleborate_mlst_ST))

# get total number of isolates
n_isolates <- length(unique(data$tip_label)) # 1559

# make MLST categories and assign colours
big_mlst <- data.frame(table(data$kleborate_mlst_ST)) 
single_mlst <- big_mlst |> filter(Freq == 1) 
big_mlst <- filter(big_mlst,!Var1=='-') |> 
  filter(Freq > (n_isolates * 0.02) ) |>
  arrange(desc(Freq)) # more than 5 isolates of same ST

data <- data |>
  mutate(mlst_cat = case_when(is.na(kleborate_mlst_ST) ~ "NA",
                              kleborate_mlst_ST %in% big_mlst$Var1 ~ kleborate_mlst_ST,
                              kleborate_mlst_ST %in% single_mlst$Var1 ~ "Single",
                              TRUE ~ "Minor"))


# Fixed colours for Single and Minor
fixed_colours <- c(
  "Single" = "#BDBDBD",
  "Minor" = "#E6E6E6"
)

# Palette for the non-Minor/non-Single STs, ordered from most common to least common
st_palette <- c(
  "#4477AA", # blue
  "#EE6677", # coral red
  "#228833", # green
  "#66CCEE", # cyan
  "#AA3377", # purple
  "#EE7733", # orange
  "#CCBB44", # mustard yellow
  "#44AA99", # teal
  "#882255", # burgundy
  "#999933"  # olive
)

# Count STs, excluding Minor and Single, and rank from most to least common
st_counts <- data |>
  filter(!mlst_cat %in% c("Minor", "Single"), !is.na(mlst_cat)) |>
  dplyr::count(mlst_cat, sort = TRUE)

# Assign colours according to abundance
st_colours <- setNames(
  st_palette[seq_len(min(nrow(st_counts), length(st_palette)))],
  st_counts$mlst_cat[seq_len(min(nrow(st_counts), length(st_palette)))]
)

# Combine fixed and automatically assigned colours
mlst_colours <- c(
  st_colours,
  fixed_colours
)

# make sure colours are character
mlst_colours <- as.character(mlst_colours)

# Order factor levels from most common to least common
# while keeping Single and Minor at the beginning
mlst_order <- st_counts$mlst_cat
mlst_levels <- c(
  mlst_order,
  "Minor",
  "Single"
)

data$mlst_cat <- factor(
  data$mlst_cat,
  levels = mlst_levels
)

# Check
#mlst_colours
#Single   Minor   ST131   ST73   ST69   ST95   ST12   ST127  ST1193   ST10 
#"#BDBDBD" "#E6E6E6" "#4477AA" "#EE6677" "#228833" "#66CCEE" "#AA3377" "#EE7733" "#CCBB44" "#44AA99" 

#~~~~~~~~~~~~~~#
# Categorise poppunk clusters
big_cluster <- data.frame(table(data$poppunk_cluster)) |>
  filter(Var1 != "-") |>
  filter(Freq > (n_isolates * 0.02)) |>
  arrange(desc(Freq))

single_cluster <- data.frame(table(data$poppunk_cluster)) |>
  filter(Var1 != "-", Freq == 1)

data <- data |>
  mutate(
    poppunk_cluster_cat = case_when(
      is.na(poppunk_cluster) ~ "NA",
      poppunk_cluster %in% big_cluster$Var1 ~ paste0(poppunk_cluster),
      poppunk_cluster %in% single_cluster$Var1 ~ "Single",
      TRUE ~ "Minor"))

#table(data$poppunk_cluster_cat, data$mlst_cat, useNA = "ifany")


# assign popPUNK cluster colours matching main MLST in that group
# find dominnat MLST cat in each poppunk cluster
cluster_mlst <- data |>
  filter(
    grepl("^ecoli_", poppunk_cluster_cat),
    !is.na(mlst_cat)) |>
  dplyr::count(poppunk_cluster_cat, mlst_cat, name = "n") |>
  group_by(poppunk_cluster_cat) |>
  arrange(desc(n), mlst_cat) |>
  slice(1) |>
  ungroup()
# add poppunk cluster sizes
cluster_sizes <- data |>
  filter(grepl("^ecoli_", poppunk_cluster_cat)) |>
  dplyr::count(poppunk_cluster_cat, name = "cluster_size")
cluster_mlst <- cluster_mlst |>
  left_join(cluster_sizes, by = "poppunk_cluster_cat")


# assign base colour palette for fastBAPS cluster that have Minor/Single MLSTs predominating
minor_single_palette <- c(
  "#CCB974", # mustard
  "#C44E52", # red
  "#8172B2", # purple
  "#64B5CD", # blue-green
  "#55A868", # green
  "#E17C05", # orange
  "#937860", # brown
  "#DA8BC3"  # pink
)

minor_single_clusters <- cluster_mlst |>
  filter(mlst_cat %in% c("Minor", "Single")) |>
  arrange(desc(cluster_size)) |>
  mutate(
    new_colour = minor_single_palette[
      ((row_number() - 1) %% length(minor_single_palette)) + 1
    ]
  )

# Add the appropriate base colour
cluster_mlst <- cluster_mlst |>
  left_join(
    minor_single_clusters |>
      select(poppunk_cluster_cat, new_colour),
    by = "poppunk_cluster_cat"
  ) |>
  mutate(
    base_colour = case_when(
      mlst_cat %in% c("Minor", "Single") ~ new_colour,
      TRUE ~ unname(mlst_colours[mlst_cat])
    )
  )

# create paler colours and different colours if multiple clusters for same mlst cat
poppunk_colour_table <- cluster_mlst |>
  group_by(mlst_cat) |>
  arrange(desc(cluster_size), .by_group = TRUE) |>
  mutate(n_same_mlst = n(),
         # Largest cluster = darkest
         # Smaller clusters = progressively paler
         lighten_amount = if_else(
           n_same_mlst == 1,
           0.55,
           seq(
             from = 0.35,
             to = 0.75,
             length.out = n())),
         poppunk_cluster_colour = mapply(
           colorspace::lighten,
           col = base_colour,
           amount = lighten_amount)) |>
  ungroup()

poppunk_order <- poppunk_colour_table |>
  mutate(mlst_rank = match(as.character(mlst_cat), mlst_order)) |>
  arrange(is.na(mlst_rank),
          mlst_rank,
          desc(cluster_size)) |>
  pull(poppunk_cluster_cat)

poppunk_order <- c(
  poppunk_order[!poppunk_order %in% c("Minor", "Single")],
  "Minor",
  "Single"
)


# create cluster colour vector:
poppunk_cluster_colours <- c(
  setNames(
    poppunk_colour_table$poppunk_cluster_colour,
    poppunk_colour_table$poppunk_cluster_cat
  ),
  "Single" = "#BDBDBD",
  "Minor" = "#E6E6E6"
)

#ecoli_4  ecoli_25   ecoli_8  ecoli_10   ecoli_9  ecoli_17  ecoli_22   ecoli_5  ecoli_15    Single     Minor 
#"#7AA5D9" "#C4DCFF" "#FFBCC1" "#7CD285" "#B8E9FF" "#FF93C9" "#FFC2B0" "#86E1CF" "#F2DF9C" "#BDBDBD" "#E6E6E6" 
#~~~~~~~~~~~~~~~~~~~~~~#

# Create grouping and colour palette for fastBAPS clusters
big_cluster <- data.frame(table(data$fastbaps_L3)) |>
  filter(Var1 != "-") |>
  filter(Freq > (n_isolates * 0.02)) |>
  arrange(desc(Freq))

single_cluster <- data.frame(table(data$fastbaps_L3)) |>
  filter(Var1 != "-", Freq == 1)

data <- data |>
  mutate(
    fastbaps_L3_cat = case_when(
      is.na(fastbaps_L3) ~ "NA",
      fastbaps_L3 %in% big_cluster$Var1 ~ paste0(fastbaps_L3),
      fastbaps_L3 %in% c("ecoli_30", "ecoli_32", "ecoli_36", "ecoli_24") ~ paste0(fastbaps_L3), # add in 3 fastbaps clusters maing up ST95, and 1 for ST127
      fastbaps_L3 %in% single_cluster$Var1 ~ "Single",
      TRUE ~ "Minor"))

#table(data$fastbaps_L3_cat, data$mlst_cat, useNA = "ifany")

# assign fastbaps cluster colours matching main MLST in that group
# find dominnat MLST cat in each poppunk cluster
cluster_mlst <- data |>
  filter(
    grepl("^ecoli_", fastbaps_L3_cat),
    !is.na(mlst_cat)) |>
  dplyr::count(fastbaps_L3_cat, mlst_cat, name = "n") |>
  group_by(fastbaps_L3_cat) |>
  arrange(desc(n), mlst_cat) |>
  slice(1) |>
  ungroup()
# add cluster sizes
cluster_sizes <- data |>
  filter(grepl("^ecoli_", fastbaps_L3_cat)) |>
  dplyr::count(fastbaps_L3_cat, name = "cluster_size")

cluster_mlst <- cluster_mlst |>
  left_join(cluster_sizes, by = "fastbaps_L3_cat")

# assign base colour palette for fastBAPS cluster that have Minor/Single MLSTs predominating
minor_single_palette <- c(
  "#C44E52", # red
  "#8172B2", # purple
  "#64B5CD", # blue-green
  "#CCB974", # mustard
  "#55A868", # green
  "#E17C05", # orange
  "#937860", # brown
  "#DA8BC3"  # pink
)

minor_single_clusters <- cluster_mlst |>
  filter(mlst_cat %in% c("Minor", "Single")) |>
  arrange(desc(cluster_size)) |>
  mutate(
    new_colour = minor_single_palette[
      ((row_number() - 1) %% length(minor_single_palette)) + 1
    ]
  )

# Add the appropriate base colour
cluster_mlst <- cluster_mlst |>
  left_join(
    minor_single_clusters |>
      select(fastbaps_L3_cat, new_colour),
    by = "fastbaps_L3_cat"
  ) |>
  mutate(
    base_colour = case_when(
      mlst_cat %in% c("Minor", "Single") ~ new_colour,
      TRUE ~ unname(mlst_colours[mlst_cat])
    )
  )

# create paler colours and different colours if multiple clusters for same mlst cat
cluster_colour_table <- cluster_mlst |>
  group_by(mlst_cat) |>
  arrange(desc(cluster_size), .by_group = TRUE) |>
  mutate(n_same_mlst = n(),
         # Largest cluster = darkest
         # Smaller clusters = progressively paler
         lighten_amount = if_else(
           n_same_mlst == 1,
           0.55,
           seq(
             from = 0.35,
             to = 0.75,
             length.out = n())),
         fastbaps_cluster_colour = mapply(
           colorspace::lighten,
           col = base_colour,
           amount = lighten_amount)) |>
  ungroup()

fastbaps_order <- cluster_colour_table |>
  mutate(mlst_rank = match(as.character(mlst_cat), mlst_order)) |>
  arrange(is.na(mlst_rank),
          mlst_rank,
          desc(cluster_size)) |>
  pull(fastbaps_L3_cat)

fastbaps_order <- c(
  fastbaps_order[!fastbaps_order %in% c("Minor", "Single")],
  "Minor",
  "Single"
)


# create cluster colour vector:
fastbaps_cluster_colours <- c(
  setNames(
    cluster_colour_table$fastbaps_cluster_colour,
    cluster_colour_table$fastbaps_L3_cat
  ),
  "Single" = "#BDBDBD",
  "Minor" = "#E6E6E6"
)
#fastbaps_cluster_colours
#ecoli_114 ecoli_19 ecoli_70 ecoli_68 ecoli_71 ecoli_61 ecoli_32 ecoli_30 ecoli_36 ecoli_46 ecoli_24 ecoli_99  ecoli_2 ecoli_140 ecoli_163 ecoli_121  Single   Minor 
#"#7AA5D9" "#C4DCFF" "#FF9BA5" "#FFBCC1" "#FFDADD" "#7CD285" "#90E0FF" "#B8E9FF" "#D9F3FF" "#FF93C9" "#FFC2B0" "#F2E174" "#86E1CF" "#98E3FC" "#F87E82" "#E0D8FF" "#BDBDBD" "#E6E6E6" 

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# create metadata panels data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Region, healthcare exposure, sampletype
#table(data$sampletype_cat, useNA = "ifany")
#table(data$region, useNA = "ifany")
#table(data$healthcare_exposure, useNA = "ifany")
#table(data$healthcare_exposure_original, useNA = "ifany") # 3 less unknown. use this original column.

type_place_hcai <- data |>
  select(tip_label, sampletype_cat, region, healthcare_exposure_original) |>
  rename(sampletype = sampletype_cat, 
         healthcare_exposure = healthcare_exposure_original) |>
  pivot_longer(cols = c(sampletype, region, healthcare_exposure) , names_to = "var")
#View(type_place_hcai)
type_place_hcai$var <- factor(type_place_hcai$var, levels = c("sampletype", "region", "healthcare_exposure"))
type_place_hcai$value <- as.factor(type_place_hcai$value)

# sample type 
colours_type <- c("blood" = "darkred","screening" = "goldenrod1", "unknown" = "grey")

sampletype_df <- type_place_hcai |>
  filter(var == "sampletype") |>
  mutate(
    value = factor(value, levels = names(colours_type)),
    xlab = factor("")
  )

# region
colours_region <- c(
  "London"       = "#E41A1C", # strong red (warm, central)
  "South East A"    = "#FB8072", # light coral
  "South East B"    = "#FCAE91", # salmon
  "South East C"    = "#FEE0D2", # pale rose
  "South West"     = "#FD8D3C", # orange
  "Midlands"      = "#B3DE69", # warm yellow-green
  "East"        = "#66C2A5", # turquoise-green
  "North West"     = "#377EB8", # deep blue
  "North East A" = "#807DBA", # purple-blue
  "North East B" = "#B3B3E6", # light lavender
  "unknown" = "grey"
)


region_df <- type_place_hcai |>
  filter(var == "region") |>
  mutate(
    value = factor(value, levels = names(colours_region)),
    xlab = factor("")
  )

# HCAI
colours_hcai <- c(
  "community" = "#FFF3B0",
  "quasi-community" = "#FFD37F",
  "quasi-nosocomial" = "#FCA17D",
  "nosocomial" = "#C75197",
  "unknown" = "grey"
)

hcai_df <- type_place_hcai |>
  filter(var == "healthcare_exposure") |>
  mutate(
    value = factor(value, levels = names(colours_hcai)),
    xlab = factor("")
  )


#~~~~~~~~~~~~~~~#
# PopPUNK and fastBAPS cluster metadata
poppunk_fastbaps_clusters <- data |>
  select(tip_label, poppunk_cluster_cat, fastbaps_L3_cat)|>
  pivot_longer(cols = c(poppunk_cluster_cat, fastbaps_L3_cat), names_to = "var")

# set order
poppunk_fastbaps_clusters$var <- factor(poppunk_fastbaps_clusters$var, levels = c("poppunk_cluster_cat", "fastbaps_L3_cat"))


# poppunk
poppunk_df <- poppunk_fastbaps_clusters|>
  filter(var == "poppunk_cluster_cat") |>
  mutate(value = factor(value, levels = poppunk_order),
         xlab = factor(""))
#table(poppunk_df$value)

#fastbaps
fastbaps_df <- poppunk_fastbaps_clusters|>
  filter(var == "fastbaps_L3_cat") |>
  mutate(value = factor(value, levels = fastbaps_order),
         xlab = factor(""))  
#table(fastbaps_df$value, useNA = "ifany")


#~~~~~~~~~~~~~#
# AMR metadata
amrfinder_data <- amrfinder_metadata_deduplicated |>
  filter(tip_label %in% data$tip_label) |>
  filter(X..Coverage.of.reference > 95 & X..Identity.to.reference > 95 & Method != "INTERNAL_STOP")
# check
#setdiff(amrfinder_data$tip_label, tree$tip.label)
#setdiff(tree$tip.label, amrfinder_data$tip_label)
length(unique(amrfinder_data$sample)) # 1,555
#colnames(amrfinder_data)

# define presence/absence colours for metadata
presence_colours <- c("absent" = "#ffffff",
                      "plasmid" = "#0d0c3b", 
                      "chromosome" = "grey", 
                      both = "#377EB8")


# herlper function to prepare gene content data for plotting:
make_gene_plot_df <- function(amrfinder_data, genes, data) {
  # all isolates x all genes
  all_combinations <- expand_grid(
    tip_label = data$tip_label,
    Element.symbol = genes
  )
  # determine where each gene occurs
  gene_locations <- amrfinder_data |>
    filter(Element.symbol %in% genes) |>
    select(tip_label, Element.symbol, molecule_type) |>
    distinct() |>
    group_by(tip_label, Element.symbol) |>
    summarise(
      has_plasmid = any(molecule_type == "plasmid", na.rm = TRUE),
      has_chromosome = any(molecule_type == "chromosome", na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      presence = case_when(
        has_plasmid & has_chromosome ~ "both",
        has_plasmid ~ "plasmid",
        has_chromosome ~ "chromosome",
        TRUE ~ "absent"
      )
    ) |>
    select(tip_label, Element.symbol, presence)
  # add absent combinations
  plot_df <- all_combinations |>
    left_join(
      gene_locations,
      by = c("tip_label", "Element.symbol")
    ) |>
    mutate(
      presence = replace_na(presence, "absent"),
      presence = factor(
        presence,
        levels = c("absent", "plasmid", "chromosome", "both")
      )
    )
  return(plot_df)
}

# select genes occuring in > 5 isolates for plotting
gene_counts_over5 <- data.frame(table(amrfinder_data$Element.symbol)) |> 
  filter(Freq >= 5) #e.g. filter out if occur <5 times in population
gene_counts_over10 <- data.frame(table(amrfinder_data$Element.symbol)) |> 
  filter(Freq >= 10) #e.g. filter out if occur <5 times in population

#filter out gent res genes
carb <- filter(amrfinder_data, grepl('CARBAPENEM', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over5$Var1) 
carb_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = carb$Element.symbol,
                                   data = data)
#length(unique(carb_selected$Element.symbol))
unique(carb_selected$Element.symbol)

ceph <- filter(amrfinder_data,grepl('CEPHALOSPORIN',Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1) |> filter(!Element.symbol=='blaEC-5')
ceph_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = ceph$Element.symbol,
                                   data = data)
unique(ceph_selected$Element.symbol)
#ceph_selected$present <- as.factor(ceph_selected$present)

beta_lactam <- filter(amrfinder_data, grepl('BETA-LACTAM', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1) |> 
  filter(!Element.symbol == 'blaEC')
beta_lactam_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                          genes = beta_lactam$Element.symbol,
                                          data = data)
unique(beta_lactam_selected$Element.symbol)
#beta_lactam_selected$present <- as.factor(beta_lactam_selected$present)

gent <- filter(amrfinder_data, grepl('GENTAMICIN', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1)
gent_selected <- make_gene_plot_df( amrfinder_data = amrfinder_data,
                                    genes = gent$Element.symbol,
                                    data = data)
unique(gent_selected$Element.symbol)
#gent_selected$present <- as.factor(gent_selected$present)

quin <- filter(amrfinder_data, grepl('QUINOLONE', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1)
quin_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = quin$Element.symbol,
                                   data = data)
unique(quin_selected$Element.symbol)
#quin_selected$present <- as.factor(quin_selected$present)

#~~~~~~~~~~~~#
# build tree
tree <- read.tree("all_ecoli_verticall.newick")
# trip tree tips to remove duplicate isoaltes that are not in the metadata df
length(tree$tip.label) # 1627
removed_duplicates <- setdiff(tree$tip.label, data$tip_label)
tree <- drop.tip(tree, removed_duplicates)
length(tree$tip.label) # 1559
# remove negative edges, rescale and reroot for beter visualisation
tree$edge.length[tree$edge.length < 0] <- 0
tree <- midpoint_root(tree) #so does this
tree$edge.length <- (tree$edge.length)^(1/3) #this just makes the tree look better
tree_plot <- ggtree(tree, options(ignore.negative.edge=TRUE))

p <- tree_plot %<+% data + 
  geom_tippoint(aes(color = mlst_cat), size = 1.5) + 
  scale_color_manual(values = mlst_colours) + 
  guides(color = guide_legend(
    title = "MLST",
    order = 1,
    ncol = 3,
    override.aes = list(size = 4))) +
  geom_treescale(x = 0, y = length(tree$tip.label) * 0.90, offset = length(tree$tip.label) * 0.01, width = 0.1, color = "black") + 
  theme_tree() 
p

p2 <- p + 
  # add clermont phylogroup
  new_scale_fill() +
  # add poppunk clustering
  new_scale_fill() +
  geom_facet(panel = "Pop\nPUNK", data = poppunk_df, geom = geom_tile,
             aes( x = xlab, fill = value), width = 1) + 
  scale_fill_manual(values = poppunk_cluster_colours,
                    guide = guide_legend(title = "PopPUNK", order = 2, ncol = 2)) +
    # Type metadata
  new_scale_fill() +
  geom_facet(panel = "Sample\ntype", data = sampletype_df, geom = geom_tile,
             mapping = aes(x = xlab, fill = value), width = 1) +
  scale_fill_manual(values = colours_type,
                  name = "Sample type", guide = guide_legend(order = 3, ncol = 2)) +
  # Region
  new_scale_fill() +
  geom_facet(panel = "Region", data = region_df, geom = geom_tile,
             mapping = aes(x = xlab, fill = value), width = 1) +
  scale_fill_manual(values = colours_region,
                    name = "Region", guide = guide_legend(order = 4, ncol = 2)) +
  # HCAI
  new_scale_fill() +
  geom_facet(panel = "HCAI", data = hcai_df, geom = geom_tile,
             mapping = aes(x = xlab, fill = value), width = 1) +
  scale_fill_manual(values = colours_hcai,
                    name = "HCAI", guide = guide_legend( order = 5, ncol = 2)) +
  # AMR genes 
  new_scale_fill() + # reset
  geom_facet(panel = "CARB", data = carb_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) + 
  geom_facet(panel = "CEPH", data = ceph_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) +
  geom_facet(panel = "B-L", data = beta_lactam_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) +
  geom_facet(panel = "GENT", data = gent_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) + 
  geom_facet(panel = "QUINOLONE", data = quin_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) +
  scale_fill_manual(values = presence_colours, name = "Gene location",  
                    guide = guide_legend(order = 6,  ncol = 2)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "white", color = "black"),
        legend.position = "right",
        legend.box = "vertical") +
  scale_x_discrete() +
  labs(color="MLST")
p2


# adjust panel size
gt <- ggplotGrob(p2)
panel_widths <- c(45, 3, 3, 3, 3,
                  5, 9, 3, 3, 15) # same order as geom_facet() calls
panel_idx <- which(grepl("^panel", gt$layout$name))
gt$widths[gt$layout$l[panel_idx]] <- unit(panel_widths, "null")

# preview
grid.newpage()
grid.draw(gt)

# save
ggsave("all_ecoli_verticall_tree_minimal.png", plot = gt, width = 20, height = 12, dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Kleb Phylo tree ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
data <- all_kleb_metadata_deduplicated
#nrow(data) #518
#colnames(data)
table(data$kleborate_mlst_ST, useNA = "ifany")
#View(data)
# remove kleb| prefix from kleborate mlst
data <- data |>
  mutate(kleborate_mlst_ST = gsub("kleb\\|", "", kleborate_mlst_ST))

# get total number of isolates
n_isolates <- length(unique(data$tip_label)) # 518

# make MLST categories and assign colours
big_mlst <- data.frame(table(data$kleborate_mlst_ST)) 
single_mlst <- big_mlst |> filter(Freq == 1) 
big_mlst <- filter(big_mlst,!Var1=='-') |> 
  filter(Freq > (n_isolates * 0.01) ) |>
  arrange(desc(Freq)) # more than 5 isolates of same ST

data <- data |>
  mutate(mlst_cat = case_when(is.na(kleborate_mlst_ST) ~ "NA",
                              kleborate_mlst_ST %in% big_mlst$Var1 ~ kleborate_mlst_ST,
                              kleborate_mlst_ST %in% single_mlst$Var1 ~ "Single",
                              TRUE ~ "Minor"))


# Fixed colours for Single and Minor
fixed_colours <- c(
  "Single" = "#BDBDBD",
  "Minor" = "#E6E6E6"
)

# Palette for the non-Minor/non-Single STs, ordered from most common to least common
# longer palette than for E. coli
st_palette <- c(
  "#4477AA", # blue
  "#EE6677", # coral
  "#228833", # green
  "#CCBB44", # mustard
  "#66CCEE", # cyan
  "#AA3377", # purple
  "#EE7733", # orange
  "#44AA99", # teal
  "#882255", # burgundy
  "#999933", # olive
  
  "#332288", # indigo
  "#117733", # dark green
  "#DDCC77", # sand
  "#CC6677", # muted rose
  "#88CCEE"  # light blue
)

# Count STs, excluding Minor and Single, and rank from most to least common
st_counts <- data |>
  filter(!mlst_cat %in% c("Minor", "Single"), !is.na(mlst_cat)) |>
  dplyr::count(mlst_cat, sort = TRUE)

# Assign colours according to abundance
st_colours <- setNames(
  st_palette[seq_len(min(nrow(st_counts), length(st_palette)))],
  st_counts$mlst_cat[seq_len(min(nrow(st_counts), length(st_palette)))]
)

# Combine fixed and automatically assigned colours
mlst_colours <- c(
  st_colours,
  fixed_colours
)

# make sure colours are character
mlst_colours <- as.character(mlst_colours)

# Order factor levels from most common to least common
# while keeping Single and Minor at the beginning
mlst_order <- st_counts$mlst_cat
mlst_levels <- c(
  mlst_order,
  "Minor",
  "Single"
)

data$mlst_cat <- factor(
  data$mlst_cat,
  levels = mlst_levels
)

# Check
table(data$mlst_cat, useNA = "ifany")
#mlst_colours
#"#4477AA" "#EE6677" "#228833" "#CCBB44" "#66CCEE" "#AA3377" "#EE7733" "#44AA99" "#882255" "#999933" "#332288" "#117733" "#DDCC77" "#CC6677" "#BDBDBD" "#E6E6E6"

#~~~~~~~~~~~~~~#
# Categorise poppunk clusters
big_cluster <- data.frame(table(data$poppunk_cluster)) |>
  filter(Var1 != "-") |>
  filter(Freq > (n_isolates * 0.01)) |>
  arrange(desc(Freq))

single_cluster <- data.frame(table(data$poppunk_cluster)) |>
  filter(Var1 != "-", Freq == 1)

data <- data |>
  mutate(
    poppunk_cluster_cat = case_when(
      is.na(poppunk_cluster) ~ "NA",
      poppunk_cluster %in% big_cluster$Var1 ~ paste0(poppunk_cluster),
      poppunk_cluster %in% single_cluster$Var1 ~ "Single",
      TRUE ~ "Minor"))

table(data$poppunk_cluster_cat, data$mlst_cat, useNA = "ifany")


# assign popPUNK cluster colours matching main MLST in that group
# find dominnat MLST cat in each poppunk cluster
cluster_mlst <- data |>
  filter(
    grepl("^k", poppunk_cluster_cat),
    !is.na(mlst_cat)) |>
  dplyr::count(poppunk_cluster_cat, mlst_cat, name = "n") |>
  group_by(poppunk_cluster_cat) |>
  arrange(desc(n), mlst_cat) |>
  slice(1) |>
  ungroup()
# add poppunk cluster sizes
cluster_sizes <- data |>
  filter(grepl("^k", poppunk_cluster_cat)) |>
  dplyr::count(poppunk_cluster_cat, name = "cluster_size")
cluster_mlst <- cluster_mlst |>
  left_join(cluster_sizes, by = "poppunk_cluster_cat")


# assign base colour palette for fastBAPS cluster that have Minor/Single MLSTs predominating
minor_single_palette <- c(
  "#CCB974", # mustard
  "#C44E52", # red
  "#8172B2", # purple
  "#64B5CD", # blue-green
  "#55A868", # green
  "#E17C05", # orange
  "#937860", # brown
  "#DA8BC3"  # pink
)

minor_single_clusters <- cluster_mlst |>
  filter(mlst_cat %in% c("Minor", "Single")) |>
  arrange(desc(cluster_size)) |>
  mutate(
    new_colour = minor_single_palette[
      ((row_number() - 1) %% length(minor_single_palette)) + 1
    ]
  )

# Add the appropriate base colour
cluster_mlst <- cluster_mlst |>
  left_join(
    minor_single_clusters |>
      select(poppunk_cluster_cat, new_colour),
    by = "poppunk_cluster_cat"
  ) |>
  mutate(
    base_colour = case_when(
      mlst_cat %in% c("Minor", "Single") ~ new_colour,
      TRUE ~ unname(mlst_colours[mlst_cat])
    )
  )

# create paler colours and different colours if multiple clusters for same mlst cat
poppunk_colour_table <- cluster_mlst |>
  group_by(mlst_cat) |>
  arrange(desc(cluster_size), .by_group = TRUE) |>
  mutate(n_same_mlst = n(),
         # Largest cluster = darkest
         # Smaller clusters = progressively paler
         lighten_amount = if_else(
           n_same_mlst == 1,
           0.55,
           seq(
             from = 0.35,
             to = 0.75,
             length.out = n())),
         poppunk_cluster_colour = mapply(
           colorspace::lighten,
           col = base_colour,
           amount = lighten_amount)) |>
  ungroup()

poppunk_order <- poppunk_colour_table |>
  mutate(mlst_rank = match(as.character(mlst_cat), mlst_order)) |>
  arrange(is.na(mlst_rank),
          mlst_rank,
          desc(cluster_size)) |>
  pull(poppunk_cluster_cat)

poppunk_order <- c(
  poppunk_order[!poppunk_order %in% c("Minor", "Single")],
  "Minor",
  "Single"
)


# create cluster colour vector:
poppunk_cluster_colours <- c(
  setNames(
    poppunk_colour_table$poppunk_cluster_colour,
    poppunk_colour_table$poppunk_cluster_cat
  ),
  "Single" = "#BDBDBD",
  "Minor" = "#E6E6E6"
)
poppunk_cluster_colours
#kpsc_3         kpsc_12     kpsc_48_128         kpsc_37        kosc_148 kpsc_11_108_277          kpsc_8        kpsc_297         kpsc_29 
#"#97C1F5"       "#FFBCC1"       "#7CD285"       "#F2E174"       "#B8E9FF"       "#FF93C9"       "#FFC2B0"       "#86E1CF"       "#EC8CB3" 
#kpsc_16         kpsc_36          kpsc_4         kpsc_10         kpsc_13         kosc_12         kosc_30          Single           Minor 
#"#D3D37A"       "#9D95E3"       "#7AC98A"       "#FFB2BD"       "#E4D18D"       "#FFA7A9"       "#E0D8FF"       "#BDBDBD"       "#E6E6E6" 

#~~~~~~~~~~~~~~~~~~~~~~#

# Create grouping and colour palette for fastBAPS clusters
big_cluster <- data.frame(table(data$fastbaps_L3)) |>
  filter(Var1 != "-") |>
  filter(Freq > (n_isolates * 0.01)) |>
  arrange(desc(Freq))

single_cluster <- data.frame(table(data$fastbaps_L3)) |>
  filter(Var1 != "-", Freq == 1)

data <- data |>
  mutate(
    fastbaps_L3_cat = case_when(
      is.na(fastbaps_L3) ~ "NA",
      fastbaps_L3 %in% big_cluster$Var1 ~ paste0(fastbaps_L3),
      fastbaps_L3 %in% c("kleb_30", "kleb_32", "kleb_36", "kleb_24") ~ paste0(fastbaps_L3), # add in 3 fastbaps clusters maing up ST95, and 1 for ST127
      fastbaps_L3 %in% single_cluster$Var1 ~ "Single",
      TRUE ~ "Minor"))

#table(data$fastbaps_L3_cat, data$mlst_cat, useNA = "ifany")

# assign fastbaps cluster colours matching main MLST in that group
# find dominnat MLST cat in each poppunk cluster
cluster_mlst <- data |>
  filter(
    # grepl("^kleb_", fastbaps_L3_cat),
    !is.na(mlst_cat)) |>
  dplyr::count(fastbaps_L3_cat, mlst_cat, name = "n") |>
  group_by(fastbaps_L3_cat) |>
  arrange(desc(n), mlst_cat) |>
  slice(1) |>
  ungroup()
# add cluster sizes
cluster_sizes <- data |>
  # filter(grepl("^kleb_", fastbaps_L3_cat)) |>
  dplyr::count(fastbaps_L3_cat, name = "cluster_size")
cluster_mlst <- cluster_mlst |>
  left_join(cluster_sizes, by = "fastbaps_L3_cat")

# assign base colour palette for fastBAPS cluster that have Minor/Single MLSTs predominating
minor_single_palette <- c(
  "#C44E52", # red
  "#8172B2", # purple
  "#64B5CD", # blue-green
  "#CCB974", # mustard
  "#55A868", # green
  "#E17C05", # orange
  "#937860", # brown
  "#DA8BC3"  # pink
)

minor_single_clusters <- cluster_mlst |>
  filter(mlst_cat %in% c("Minor", "Single")) |>
  arrange(desc(cluster_size)) |>
  mutate(
    new_colour = minor_single_palette[
      ((row_number() - 1) %% length(minor_single_palette)) + 1
    ]
  )

# Add the appropriate base colour
cluster_mlst <- cluster_mlst |>
  left_join(
    minor_single_clusters |>
      select(fastbaps_L3_cat, new_colour),
    by = "fastbaps_L3_cat"
  ) |>
  mutate(
    base_colour = case_when(
      mlst_cat %in% c("Minor", "Single") ~ new_colour,
      TRUE ~ unname(mlst_colours[mlst_cat])
    )
  )

# create paler colours and different colours if multiple clusters for same mlst cat
cluster_colour_table <- cluster_mlst |>
  group_by(mlst_cat) |>
  arrange(desc(cluster_size), .by_group = TRUE) |>
  mutate(n_same_mlst = n(),
         # Largest cluster = darkest
         # Smaller clusters = progressively paler
         lighten_amount = if_else(
           n_same_mlst == 1,
           0.55,
           seq(
             from = 0.35,
             to = 0.75,
             length.out = n())),
         fastbaps_cluster_colour = mapply(
           colorspace::lighten,
           col = base_colour,
           amount = lighten_amount)) |>
  ungroup()

fastbaps_order <- cluster_colour_table |>
  mutate(mlst_rank = match(as.character(mlst_cat), mlst_order)) |>
  arrange(is.na(mlst_rank),
          mlst_rank,
          desc(cluster_size)) |>
  pull(fastbaps_L3_cat)

fastbaps_order <- c(
  fastbaps_order[!fastbaps_order %in% c("Minor", "Single")],
  "Minor",
  "Single"
)


# create cluster colour vector:
fastbaps_cluster_colours <- c(
  setNames(
    cluster_colour_table$fastbaps_cluster_colour,
    cluster_colour_table$fastbaps_L3_cat
  ),
  "Single" = "#BDBDBD",
  "Minor" = "#E6E6E6"
)
#fastbaps_cluster_colours

#kpsc_7    kpsc_6    kosc_4   kpsc_12   kpsc_13   kpsc_15    Minor   kosc_14   kpsc_14    kpsc_4    kosc_7 kaerogenes_2   kosc_15    kosc_3    Single   kosc_13    kpsc_2    Single 
#"#7AA5D9"  "#C4DCFF"  "#B8E9FF"  "#FF93C9"  "#D3D37A"  "#85D2EA"  "#91E1A1"  "#FFEDAF"  "#F87E82"  "#B3A5E5"  "#EBD895"  "#FFBA96"  "#DABDA6"  "#FFCBEE"  "#FFBCBD"  "#DAD0FF"  "#C8F0FF"  "#BDBDBD" 
#Minor 
#"#E6E6E6" 

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# create metadata panels data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Region, healthcare exposure, sampletype
#table(data$sampletype_cat, useNA = "ifany")
#table(data$region, useNA = "ifany")
#table(data$healthcare_exposure, useNA = "ifany")
#table(data$healthcare_exposure_original, useNA = "ifany") # 3 less unknown. use this original column.

type_place_hcai <- data |>
  select(tip_label, sampletype_cat, region, healthcare_exposure_original) |>
  rename(sampletype = sampletype_cat, 
         healthcare_exposure = healthcare_exposure_original) |>
  pivot_longer(cols = c(sampletype, region, healthcare_exposure) , names_to = "var")
#View(type_place_hcai)
type_place_hcai$var <- factor(type_place_hcai$var, levels = c("sampletype", "region", "healthcare_exposure"))
type_place_hcai$value <- as.factor(type_place_hcai$value)

# sample type 
colours_type <- c("blood" = "darkred","screening" = "goldenrod1", "unknown" = "grey")

sampletype_df <- type_place_hcai |>
  filter(var == "sampletype") |>
  mutate(
    value = factor(value, levels = names(colours_type)),
    xlab = factor("")
  )

# region
colours_region <- c(
  "London"       = "#E41A1C", # strong red (warm, central)
  "South East A"    = "#FB8072", # light coral
  "South East B"    = "#FCAE91", # salmon
  "South East C"    = "#FEE0D2", # pale rose
  "South West"     = "#FD8D3C", # orange
  "Midlands"      = "#B3DE69", # warm yellow-green
  "East"        = "#66C2A5", # turquoise-green
  "North West"     = "#377EB8", # deep blue
  "North East A" = "#807DBA", # purple-blue
  "North East B" = "#B3B3E6", # light lavender
  "unknown" = "grey"
)


region_df <- type_place_hcai |>
  filter(var == "region") |>
  mutate(
    value = factor(value, levels = names(colours_region)),
    xlab = factor("")
  )

#~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~#
#  define phylo colours
kleb_colours <- c(
  # pneumoniae / variicola group — warm family
  "Klebsiella pneumoniae subsp. pneumoniae" =
    "#F3B61F",  # golden orange
  "Klebsiella quasipneumoniae subsp. quasipneumoniae" =
    "#E69F00",  # orange
  "Klebsiella quasipneumoniae subsp. similipneumoniae" =
    "#D55E00",  # strong vermillion
  "Klebsiella variicola subsp. variicola" =
    "#CC79A7",  # reddish purple
  "Klebsiella quasivariicola" =
    "#9C5A8E",  # deeper purple
  # aerogenes — separate teal/green family
  "Klebsiella aerogenes" =
    "#009E73",  # bluish green
  # oxytoca complex — blue family
  "Klebsiella oxytoca" =
    "#0072B2",  # strong blue
  "Klebsiella michiganensis" =
    "#56B4E9",  # sky blue
  "Klebsiella grimontii" =
    "#3B8BC2",  # medium blue
  "Klebsiella pasteurii" =
    "#8EC9E6"   # pale blue
)


# PopPUNK and fastBAPS cluster metadata
subspp_poppunk_fastbaps_clusters <- data |>
  select(tip_label, kraken2_subspecies, poppunk_cluster_cat, fastbaps_L3_cat)|>
  pivot_longer(cols = c(kraken2_subspecies, poppunk_cluster_cat, fastbaps_L3_cat), names_to = "var")

# set order
subspp_poppunk_fastbaps_clusters$var <- factor(subspp_poppunk_fastbaps_clusters$var, levels = c("kraken2_subspecies", "poppunk_cluster_cat", "fastbaps_L3_cat"))
#separate dfs
subspp_df <- subspp_poppunk_fastbaps_clusters |>
  filter(var == "kraken2_subspecies") |>
  mutate(value = factor(value),
         xlab = factor(""))

# poppunk
poppunk_df <- subspp_poppunk_fastbaps_clusters|>
  filter(var == "poppunk_cluster_cat") |>
  mutate(value = factor(value, levels = poppunk_order),
         xlab = factor(""))
#table(poppunk_df$value)


#~~~~~~~~~~~~~~~~~~~~~~#
# AMR metadata
amrfinder_data <- amrfinder_metadata_deduplicated |>
  filter(tip_label %in% data$tip_label) |>
  filter(X..Coverage.of.reference > 95 & X..Identity.to.reference > 95 & Method != "INTERNAL_STOP")

# check
#setdiff(amrfinder_data$tip_label, tree$tip.label)
#setdiff(tree$tip.label, amrfinder_data$tip_label)
length(unique(amrfinder_data$sample)) # 518
#colnames(amrfinder_data)

# define presence/absence colours for metadata
presence_colours <- c("absent" = "#ffffff",
                      "plasmid" = "#0d0c3b", 
                      "chromosome" = "grey", 
                      both = "#377EB8")



# select genes occuring in > 5 isolates for plotting
gene_counts_over5 <- data.frame(table(amrfinder_data$Element.symbol)) |> 
  filter(Freq >= 5) #e.g. filter out if occur <5 times in population
gene_counts_over10 <- data.frame(table(amrfinder_data$Element.symbol)) |> 
  filter(Freq >= 10) #e.g. filter out if occur <5 times in population

#filter out gent res genes
carb <- filter(amrfinder_data, grepl('CARBAPENEM', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over5$Var1)
carb_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = carb$Element.symbol,
                                   data = data)
unique(carb_selected$Element.symbol)

ceph <- filter(amrfinder_data,grepl('CEPHALOSPORIN',Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1) |> filter(!Element.symbol=='blaEC-5')
ceph_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = ceph$Element.symbol,
                                   data = data)
unique(ceph_selected$Element.symbol)
#ceph_selected$present <- as.factor(ceph_selected$present)

beta_lactam <- filter(amrfinder_data, grepl('BETA-LACTAM', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1) |> 
  filter(!Element.symbol == 'blaEC')
beta_lactam_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data, 
                                          genes = beta_lactam$Element.symbol,
                                          data = data)
unique(beta_lactam_selected$Element.symbol)
#beta_lactam_selected$present <- as.factor(beta_lactam_selected$present)

gent <- filter(amrfinder_data, grepl('GENTAMICIN', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1)
gent_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = gent$Element.symbol,
                                   data = data)
unique(gent_selected$Element.symbol)
#gent_selected$present <- as.factor(gent_selected$present)

quin <- filter(amrfinder_data, grepl('QUINOLONE', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1)
quin_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = quin$Element.symbol,
                                   data = data)
unique(quin_selected$Element.symbol)
#quin_selected$present <- as.factor(quin_selected$present)

#~~~~~~~~~~~~#
# build tree
tree <- read.tree("all_kleb_verticall.newick")
# trip tree tips to remove duplicate isoaltes that are not in the metadata df
length(tree$tip.label) # 1627
removed_duplicates <- setdiff(tree$tip.label, data$tip_label)
tree <- drop.tip(tree, removed_duplicates)
length(tree$tip.label) # 1559
# remove negative edges, rescale and reroot for beter visualisation
tree$edge.length[tree$edge.length < 0] <- 0
tree <- midpoint_root(tree) #so does this
tree$edge.length <- (tree$edge.length)^(1/3) #this just makes the tree look better
tree_plot <- ggtree(tree, options(ignore.negative.edge=TRUE))

p <- tree_plot %<+% data + 
  geom_tippoint(aes(color = mlst_cat), size = 2) + 
  scale_color_manual(values = mlst_colours) + 
  guides(color = guide_legend(title = "MLST", order = 1, ncol = 3,
                              override.aes = list(size = 4))) +
  geom_treescale(x = 0, y = length(tree$tip.label) * 0.90, offset = length(tree$tip.label) * 0.01, width = 0.1, color = "black") + 
  theme_tree() 
p

p2 <- p + 
  new_scale_fill() +
  # add subspecies
  new_scale_fill() +
  geom_facet(panel = "Sub\nspp.", data = subspp_df, geom = geom_tile,
             aes( x = xlab, fill = value), width = 1) + 
  scale_fill_manual(values = kleb_colours,
                    guide = guide_legend(title = "Species-subspecies", order = 2, ncol = 1)) +
  
  # add poppunk clustering 
  new_scale_fill() +
  geom_facet(panel = "Pop\nPUNK", data = poppunk_df, geom = geom_tile,
             aes( x = xlab, fill = value), width = 1) + 
  scale_fill_manual(values = poppunk_cluster_colours,
                    guide = guide_legend(title = "PopPUNK", order = 3, ncol = 2)) +
  
  # Type metadata
  new_scale_fill() +
  geom_facet(panel = "Sample\ntype", data = sampletype_df, geom = geom_tile,
             mapping = aes(x = xlab, fill = value), width = 1) +
  scale_fill_manual(values = colours_type,
                    name = "Sample type", guide = guide_legend(order = 4, ncol = 2)) +
  # Region
  new_scale_fill() +
  geom_facet(panel = "Region", data = region_df, geom = geom_tile,
             mapping = aes(x = xlab, fill = value), width = 1) +
  scale_fill_manual(values = colours_region,
                    name = "Region", guide = guide_legend(order = 5, ncol = 2)) +
  # AMR genes 
  new_scale_fill() + # reset
  geom_facet(panel = "CARB", data = carb_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) + 
  geom_facet(panel = "CEPH", data = ceph_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) +
  geom_facet(panel = "B-L", data = beta_lactam_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) +
  geom_facet(panel = "GENT", data = gent_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) + 
  geom_facet(panel = "QUINOLONE", data = quin_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) +
  scale_fill_manual(values = presence_colours, name = "Gene location",  
                    guide = guide_legend(order = 6,  ncol = 2)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "white", color = "black"),
        legend.position = "right",
        legend.box = "vertical") +
  scale_x_discrete() +
  labs(color="MLST")
p2


# adjust panel size
gt <- ggplotGrob(p2)
panel_widths <- c(45, 2.5, 2.5, 2.5, 2.5,
                  8, 6, 6, 3, 15) # same order as geom_facet() calls
panel_idx <- which(grepl("^panel", gt$layout$name))
gt$widths[gt$layout$l[panel_idx]] <- unit(panel_widths, "null")

# preview
grid.newpage()
grid.draw(gt)

# save
ggsave("all_kleb_verticall_tree_minimal.png", plot = gt, width = 21.5, height = 12, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * All CPE phylo tree ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# use original metadata table to include the serratia marcescens isolate
data <- samples_assembled_with_species_hcai_deduplicated |>
  filter(duplicate_assembly_qc_pass != "duplicate") |>
  filter(sampletype_cat == "screening" | grepl("Citrobacter", kraken2_species) | grepl("Enterobacter", kraken2_species))
nrow(data) #179
#View(data) #179
#colnames(data)

# fix tip labels
data <- data |>
  mutate(tip_label = sequencing_id) |>
  select(tip_label, everything())

data <- data |> 
  mutate(mlst = case_when(mlst != "-" & typing_scheme == "ecoli_achtman_4" ~ paste0("E. coli ST", mlst),
                          mlst != "-" & typing_scheme == "klebsiella" ~ paste0("K. pneumoniae ST", mlst),
                          mlst != "-" & typing_scheme == "koxytoca" ~ paste0("K. oxytoca ST", mlst),
                          mlst != "-" & typing_scheme == "kaerogenes" ~ paste0("K. aerogenes ST", mlst),
                          mlst != "-" & typing_scheme == "cfreundii" ~ paste0("C. freundii ST", mlst),
                          mlst != "-" & typing_scheme == "ecloacae" ~ paste0("E. cloacae ST", mlst),
                          mlst == "-" & kraken2_species == "Serratia marcescens" ~ "S. marcescens",
                          mlst == "-" & kraken2_species != "Serratia marcescens" ~ kleborate_mlst_ST,
                           TRUE ~ mlst))
table(data$mlst, useNA = "ifany")
#View(data)

# get total number of isolates
n_isolates <- length(unique(data$tip_label)) # 180

# make MLST categories and assign colours
big_mlst <- data.frame(table(data$mlst)) 
single_mlst <- big_mlst |> filter(Freq == 1) 
big_mlst <- filter(big_mlst,!Var1=='-') |> 
  filter(Freq > (2) ) |>
  arrange(desc(Freq)) # more than 5 isolates of same ST

data <- data |>
  mutate(mlst_cat = case_when(is.na(mlst) ~ "NA",
                              mlst %in% big_mlst$Var1 ~ mlst,
                              mlst %in% single_mlst$Var1 ~ "Single",
                              TRUE ~ "Two isolates"))
table(data$mlst_cat, useNA = "ifany")

# Fixed colours for Single and Minor
cpe_colours <- c(
  "K. oxytoca ST31"     = "#CC79A7",  # pink
  
  # E. coli - distinct greens, blues, teals and lime
  "E. coli ST38"        = "#009E73",  # emerald green
  "E. coli ST10"        = "#0072B2",  # strong blue
  "E. coli ST167"       = "#66A61E",  # olive green
  "E. coli ST399"       = "#56B4E9",  # sky blue
  "E. coli ST409"       = "#006D2C",  # forest green
  "E. coli ST69"        = "#17BECF",  # cyan
  
  # Other species complexes - unchanged
  "K. pneumoniae ST147" = "#D55E00",  # orange
  "C. freundii ST22"    = "#E69F00",  # golden yellow
  
  # E. coli - additional distinct shades
  "E. coli ST405"       = "#1F4E79",  # navy blue
  "E. coli ST410"       = "#7AD151",  # bright lime green
  "E. coli ST448"       = "#2A9D8F",  # muted turquoise
  "E. coli ST648"       = "#3B528B",  # indigo blue
  
  # Other species complexes - unchanged
  "K. pneumoniae ST16"  = "#E76F51",  # coral
  "K. pneumoniae ST231" = "#A63A00",  # burnt orange
  
  # Neutral categories
  "Two isolates"              = "#BDBDBD",
  "Single"              = "#E6E6E6"
)

# Match MLST factor order
cpe_colours <- cpe_colours[mlst_levels]

# make sure colours are character
cpe_colours <- as.character(cpe_colours)

st_counts <- data |>
  group_by(mlst_cat) |>
  summarise(count = n()) |>
  arrange(desc(count))

# Order factor levels from most common to least common
# while keeping Single and Minor at the beginning
mlst_order <- setdiff(st_counts$mlst_cat, c("Two isolates", "Single"))
mlst_levels <- c(mlst_order, "Two isolates", "Single")

data$mlst_cat <- factor(data$mlst_cat, levels = mlst_levels)

# Check
table(data$mlst_cat, useNA = "ifany")
cpe_colours
#~~~~~~~~~~~~~~#
data <- data |>
  mutate(species_complex = case_when(kraken2_species_complex == "Escherichia coli" ~ "E. coli",
                                     kraken2_species_complex == "Klebsiella pneumoniae species complex" ~ "K. pneumoniae complex",
                                     kraken2_species_complex == "Klebsiella oxytoca species complex" ~ "K. oxytoca complex",
                                     kraken2_species_complex == "Klebsiella aerogenes" ~ "K. aerogenes",
                                     kraken2_species_complex == "Citrobacter freundii" ~ "C. freundii complex",
                                     kraken2_species_complex == "Enterobacter cloacae complex" ~ "E. cloacae complex",
                                     kraken2_species_complex == "Serratia marcescens" ~ "S. marcescens",
                                     TRUE ~ NA))

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# create metadata panels data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Region, healthcare exposure, sampletype
#table(data$sampletype_cat, useNA = "ifany")
#table(data$region, useNA = "ifany")
#table(data$healthcare_exposure, useNA = "ifany")
#table(data$healthcare_exposure_original, useNA = "ifany") # 3 less unknown. use this original column.

type_place_hcai <- data |>
  select(tip_label, sampletype_cat, region, healthcare_exposure_original) |>
  rename(sampletype = sampletype_cat, 
         healthcare_exposure = healthcare_exposure_original) |>
  pivot_longer(cols = c(sampletype, region, healthcare_exposure) , names_to = "var")
#View(type_place_hcai)
type_place_hcai$var <- factor(type_place_hcai$var, levels = c("sampletype", "region", "healthcare_exposure"))
type_place_hcai$value <- as.factor(type_place_hcai$value)

# sample type 
colours_type <- c("blood" = "darkred","screening" = "goldenrod1", "unknown" = "grey")

sampletype_df <- type_place_hcai |>
  filter(var == "sampletype") |>
  mutate(value = factor(value, levels = names(colours_type)),
         xlab = factor("")
  )

# region
colours_region <- c(
  "London"       = "#E41A1C", # strong red (warm, central)
  "South East A"    = "#FB8072", # light coral
  "South East B"    = "#FCAE91", # salmon
  "South East C"    = "#FEE0D2", # pale rose
  "South West"     = "#FD8D3C", # orange
  "Midlands"      = "#B3DE69", # warm yellow-green
  "East"        = "#66C2A5", # turquoise-green
  "North West"     = "#377EB8", # deep blue
  "North East A" = "#807DBA", # purple-blue
  "North East B" = "#B3B3E6", # light lavender
  "unknown" = "grey"
)


region_df <- type_place_hcai |>
  filter(var == "region") |>
  mutate(
    value = factor(value, levels = names(colours_region)),
    xlab = factor("")
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
phylogroup_poppunk_fastbaps_clusters <- data |>
  select(tip_label, species_complex, poppunk_cluster)|>
  pivot_longer(cols = c(species_complex, poppunk_cluster), names_to = "var")

# set order
phylogroup_poppunk_fastbaps_clusters$var <- factor(phylogroup_poppunk_fastbaps_clusters$var, levels = c("species_complex", "poppunk_cluster"))
#separate dfs
species_complex_df <- phylogroup_poppunk_fastbaps_clusters |>
  filter(var == "species_complex") |>
  mutate(value = factor(value),
         xlab = factor(""))

species_colours <- c(
  "E. coli" = "#009E73",
  "K. pneumoniae complex" = "#D55E00",
  "K. oxytoca complex" = "#CC79A7",
  "K. aerogenes" = "#7F3C8D",
  "C. freundii complex" = "#E69F00",
  "E. cloacae complex" = "#0072B2",
  "S. marcescens" = "#8B1A1A"
)


#~~~~~~~~~~~~~#
# AMR metadata
amrfinder_data <- amrfinder_metadata_deduplicated |>
  mutate(tip_label = sample) |>
  filter(tip_label %in% data$tip_label) |>
  filter(X..Coverage.of.reference >95 & X..Identity.to.reference >95 & Method != "INTERNAL_STOP") 

# check
#setdiff(amrfinder_data$tip_label, tree$tip.label)
#setdiff(tree$tip.label, amrfinder_data$tip_label)
length(unique(amrfinder_data$sample)) #180
#colnames(amrfinder_data)

# define presence/absence colours for metadata
presence_colours <- c("absent" = "#ffffff",
                      "plasmid" = "#0d0c3b", 
                      "chromosome" = "grey", 
                      both = "#377EB8")


# select genes occuring in > 5 isolates for plotting
gene_counts_over5 <- data.frame(table(amrfinder_data$Element.symbol)) |> 
  filter(Freq >= 5) #e.g. filter out if occur <5 times in population
gene_counts_over10 <- data.frame(table(amrfinder_data$Element.symbol)) |> 
  filter(Freq >= 10) #e.g. filter out if occur <5 times in population


#filter out gent res genes
carb <- filter(amrfinder_data, grepl('CARBAPENEM', Subclass)) |> 
  distinct(Element.symbol) #|> filter(Element.symbol %in% gene_counts_over5$Var1) 
carb_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = carb$Element.symbol,
                                   data = data)
unique(carb_selected$Element.symbol)

ceph <- filter(amrfinder_data,grepl('CEPHALOSPORIN',Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1) |> filter(!Element.symbol=='blaEC-5')
ceph_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = ceph$Element.symbol,
                                   data = data)
unique(ceph_selected$Element.symbol)
#ceph_selected$present <- as.factor(ceph_selected$present)

beta_lactam <- filter(amrfinder_data, grepl('BETA-LACTAM', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1) |> 
  filter(!Element.symbol == 'blaEC')
beta_lactam_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                          genes = beta_lactam$Element.symbol,
                                          data = data)
unique(beta_lactam_selected$Element.symbol)
#beta_lactam_selected$present <- as.factor(beta_lactam_selected$present)

gent <- filter(amrfinder_data, grepl('GENTAMICIN', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1)
gent_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = gent$Element.symbol,
                                   data = data)
unique(gent_selected$Element.symbol)
#gent_selected$present <- as.factor(gent_selected$present)

quin <- filter(amrfinder_data, grepl('QUINOLONE', Subclass)) |> 
  distinct(Element.symbol) |> filter(Element.symbol %in% gene_counts_over10$Var1)
quin_selected <- make_gene_plot_df(amrfinder_data = amrfinder_data,
                                   genes = quin$Element.symbol,
                                   data = data)
unique(quin_selected$Element.symbol)
#quin_selected$present <- as.factor(quin_selected$present)

#~~~~~~~~~~~~~~#
# build tree
tree <- read.tree("cpe_mashtree.dnd")
# trip tree tips to remove duplicate isoaltes that are not in the metadata df
length(tree$tip.label) # 210


#intersect(tree$tip.label, data$tip_label)
#setdiff(data$tip_label, tree$tip.label)
removed_duplicates <- setdiff(tree$tip.label, data$tip_label)
tree <- drop.tip(tree, removed_duplicates)
length(tree$tip.label) # 179
# remove negative edges, rescale and reroot for beter visualisation
tree$edge.length[tree$edge.length < 0] <- 0
tree <- midpoint_root(tree) #so does this
tree$edge.length <- (tree$edge.length)^(1/3) #this just makes the tree look better
tree_plot <- ggtree(tree, options(ignore.negative.edge=TRUE))

p <- tree_plot %<+% data + 
  geom_tippoint(aes(color = mlst_cat), size = 3) + 
  scale_color_manual(values = cpe_colours,
                     limits = mlst_levels,
                     breaks = mlst_levels,
                     drop = FALSE) + 
  guides(color = guide_legend(title = "MLST",
    order = 1, ncol = 2,
    override.aes = list(size = 4))) +
  geom_treescale(x = 0, y = length(tree$tip.label) * 0.90, offset = length(tree$tip.label) * 0.01, width = 0.1, color = "black") + 
  theme_tree() 
p

p2 <- p + 
  # add clermont phylogroup
  new_scale_fill() +
  geom_facet(panel = "Spp.", data = species_complex_df, geom = geom_tile,
             aes( x = xlab, fill = value), width = 1) + 
  scale_fill_manual(values = species_colours,
                    guide = guide_legend(title = "Species complex", order = 2, ncol = 1)) +
 # Type metadata
  new_scale_fill() +
  geom_facet(panel = "Sample\ntype", data = sampletype_df, geom = geom_tile,
             mapping = aes(x = xlab, fill = value), width = 0.8) +
  scale_fill_manual(values = colours_type,
                    name = "Sample type", guide = guide_legend(order = 5, ncol = 3)) +
  # Region
  new_scale_fill() +
  geom_facet(panel = "Region", data = region_df, geom = geom_tile,
             mapping = aes(x = xlab, fill = value), width = 0.8) +
  scale_fill_manual(values = colours_region,
                    name = "Region", guide = guide_legend(order = 6, ncol = 2)) +
   # AMR genes 
  new_scale_fill() + # reset
  geom_facet(panel = "CARB", data = carb_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) + 
  geom_facet(panel = "CEPH", data = ceph_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) +
  geom_facet(panel = "B-L", data = beta_lactam_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) +
  geom_facet(panel = "GENT", data = gent_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) + 
  geom_facet(panel = "QUINOLONE", data = quin_selected, geom = geom_tile, 
             aes(x = Element.symbol, fill = presence), width = 1) +
  scale_fill_manual(values = presence_colours, name = "Gene location",  
                    guide = guide_legend(order = 6,  ncol = 2)) +
  
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        strip.text = element_text(size = 7),
        strip.background = element_rect(fill = "white", color = "black"),
        legend.position = "right",
        legend.box = "vertical") +
  scale_x_discrete() +
  labs(color="MLST")
p2


# adjust panel size
gt <- ggplotGrob(p2)
panel_widths <- c(45, 3, 3, 3,
                  30, 7, 3, 5, 13) # same order as geom_facet() calls
panel_idx <- which(grepl("^panel", gt$layout$name))
gt$widths[gt$layout$l[panel_idx]] <- unit(panel_widths, "null")

# preview
grid.newpage()
grid.draw(gt)

# save
ggsave("all_cpe_verticall_tree.png", plot = gt, width = 20, height = 12, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# MLSt composition across epidemiological compartment ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * E. coli associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# metadata associations:
# create categorical variable for number of ARG classes and for BL carriage
# load data if not already loaded
#amrfinder_metadata <- read.csv("amrfinder_metadata_with_NAs_updated.csv")

# heler function to separate ot classes
get_amr_classes <- function(x) {
  # Split combined classes
  classes <- str_split(x, "/") |>
    unlist() |>
    str_trim()
   # Remove missing/empty values
  classes <- unique(
    classes[
      !is.na(classes) &
        classes != ""
    ]
  )
  
  # No classes
  if (length(classes) == 0) {
    return(character(0))
  }
  
  # Remove broad EFFLUX/MULTIDRUG categories
  specific_classes <- setdiff(
    classes,
    c("EFFLUX", "MULTIDRUG")
  )
  
  # If specific classes exist, use only those
  if (length(specific_classes) > 0) {
    return(sort(specific_classes))
  }
  
  # Otherwise EFFLUX/MULTIDRUG were the only classes present,
  # so treat them collectively as one class
  return("EFFLUX/MULTIDRUG")
}

# call function
amr_categories <- amrfinder_metadata |>
  filter(sample %in% all_metadata_deduplicated$sequencing_id,
         Type == "AMR" | is.na(Type),
         X..Coverage.of.reference >= 95 & Method != "INTERNAL_STOP" & X..Identity.to.reference >= 95,
         Element.symbol != "blaEC") |>
  group_by(sample) |>
  summarise(
    amr_count = sum(!is.na(Element.symbol)),
    
    amr_count_unique = n_distinct(
      Element.symbol,
      na.rm = TRUE
    ),
    
    amr_classes = {
      classes <- get_amr_classes(Class)
      
      if (length(classes) == 0) {
        NA_character_
      } else {
        paste(classes, collapse = "; ")
      }
    },
    
    amr_num_classes = {
      classes <- get_amr_classes(Class)
      length(classes)
    },
    
    amr_beta_lactam = {
      
      subclasses <- str_to_upper(
        str_trim(Subclass)
      )
      
      subclasses <- subclasses[
        !is.na(subclasses) &
          subclasses != ""
      ]
      
      if (length(subclasses) == 0) {
        "no_amr"
        
      } else if (any(
        str_detect(subclasses, "CARBAPENEM")
      )) {
        "carbapenem"
        
      } else if (any(
        str_detect(subclasses, "CEPHALOSPORIN")
      )) {
        "cephalosporin"
        
      } else if (any(
        str_detect(subclasses, "BETA-LACTAM")
      )) {
        "other_beta_lactam"
        
      } else {
        "non_BL_amr"
      }
    },
    
    .groups = "drop"
  ) |>
  mutate(
    amr_presence_binary = case_when(
      amr_count == 0 ~ "absent",
      amr_count > 0 ~ "present",
      TRUE ~ NA_character_
    ),
    
    amr_count_category = case_when(
      amr_count == 0 ~ "absent",
      amr_count <= 2 ~ "1-2",
      amr_count == 3 ~ "3",
      amr_count == 4 ~ "4",
      amr_count == 5 ~ "5",
      amr_count <= 7 ~ "6-7",
      amr_count > 7 ~ "8+",
      TRUE ~ NA_character_
    ),
    
    amr_count_unique_category = case_when(
      amr_count_unique == 0 ~ "absent",
      amr_count_unique <= 2 ~ "1-2",
      amr_count_unique == 3 ~ "3",
      amr_count_unique == 4 ~ "4",
      amr_count_unique == 5 ~ "5",
      amr_count_unique <= 7 ~ "6-7",
      amr_count_unique > 7 ~ "8+",
      TRUE ~ NA_character_
    ),
    
    amr_num_classes_category = case_when(
      amr_num_classes == 0 ~ "absent",
      amr_num_classes == 1 ~ "1",
      amr_num_classes == 2 ~ "2",
      amr_num_classes == 3 ~ "3",
      amr_num_classes >= 4 ~ "4+",
      TRUE ~ NA_character_
    ),
    
    amr_beta_lactam = factor(
      amr_beta_lactam,
      levels = c(
        "no_amr",
        "non_BL_amr",
        "other_beta_lactam",
        "cephalosporin",
        "carbapenem"
      )
    )
  )
#View(amr_categories)
table(amr_categories$amr_count)
table(amr_categories$amr_count_unique)
table(amr_categories$amr_num_classes )

table(amr_categories$amr_beta_lactam )
table(amr_categories$amr_count_category )
table(amr_categories$amr_count_unique_category )
table(amr_categories$amr_num_classes_category, useNA = "ifany")

# save
write.csv(amr_categories, "amr_categories.csv", row.names = FALSE)
# reload if needed
#amr_categories <- read.csv("amr_categories.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# repeat by autocycler contig
contigs_info_merged <- read.csv("autocycler_amrfinder_unitig_contig_ids.csv")
nrow(contigs_info_merged)
colnames(contigs_info_merged)
# check
#setdiff(amrfinder_metadata$sample, contigs_info_merged$sample) # 0 in amrfidner data that is not in contigs info
#setdiff(contigs_info_merged$sample, amrfinder_metadata$sample) # twice sequenced samples are extra

#left join onto amrfinder metadata
amrfinder_metadata_autocycler_contigs <- amrfinder_metadata |>
  left_join(contigs_info_merged, by = c("sample" = "sample", "Contig.id" = "amrfinder_Contig.id"))
sum(is.na(amrfinder_metadata_autocycler_contigs$autocycler_contig_id))

# save if needed
write.csv(amrfinder_metadata_autocycler_contigs, "amrfinder_metadata_autocycler_contigs.csv", row.names = FALSE)
# reload if needed
#amrfinder_metadata_autocycler_contigs <- read.csv("amrfinder_metadata_autocycler_contigs.csv")

# summarise AMRfinder hits per autocycler contig
amr_categories_contigs <- amrfinder_metadata_autocycler_contigs |>
  filter(sample %in% all_metadata_deduplicated$sequencing_id,
         Type == "AMR" | is.na(Type),
         X..Coverage.of.reference >= 95 & Method != "INTERNAL_STOP" & X..Identity.to.reference >= 95,
         Element.symbol != "blaEC") |>
  group_by(sample, autocycler_contig_id) |>
  summarise(amr_count = sum(!is.na(Element.symbol)),
            amr_count_unique = n_distinct(Element.symbol, na.rm = TRUE),
            amr_classes = {classes <- get_amr_classes(Class)
            
            if (length(classes) == 0) {
              NA_character_
            } else {
              paste(classes, collapse = "; ")
            }
            },
            
            amr_num_classes = {
              classes <- get_amr_classes(Class)
              length(classes)
            },
            
            amr_beta_lactam = {
              subclasses <- str_to_upper(  str_trim(Subclass) )
              subclasses <- subclasses[
                !is.na(subclasses) &
                  subclasses != ""
              ]
              
              if (length(subclasses) == 0) {
                "no_amr"
                
              } else if (any(
                str_detect(subclasses, "CARBAPENEM")
              )) {
                "carbapenem"
                
              } else if (any(
                str_detect(subclasses, "CEPHALOSPORIN")
              )) {
                "cephalosporin"
                
              } else if (any(
                str_detect(subclasses, "BETA-LACTAM")
              )) {
                "other_beta_lactam"
                
              } else {
                "non_BL_amr"
              }
            },
            
            .groups = "drop"
  ) |>
  mutate(
    amr_presence_binary = case_when(
      amr_count == 0 ~ "absent",
      amr_count > 0 ~ "present",
      TRUE ~ NA_character_
    ),
    
    amr_count_category = case_when(
      amr_count == 0 ~ "absent",
      amr_count <= 2 ~ "1-2",
      amr_count == 3 ~ "3",
      amr_count == 4 ~ "4",
      amr_count == 5 ~ "5",
      amr_count <= 7 ~ "6-7",
      amr_count > 7 ~ "8+",
      TRUE ~ NA_character_
    ),
    
    amr_count_unique_category = case_when(
      amr_count_unique == 0 ~ "absent",
      amr_count_unique <= 2 ~ "1-2",
      amr_count_unique == 3 ~ "3",
      amr_count_unique == 4 ~ "4",
      amr_count_unique == 5 ~ "5",
      amr_count_unique <= 7 ~ "6-7",
      amr_count_unique > 7 ~ "8+",
      TRUE ~ NA_character_
    ),
    
    amr_num_classes_category = case_when(
      amr_num_classes == 0 ~ "absent",
      amr_num_classes == 1 ~ "1",
      amr_num_classes == 2 ~ "2",
      amr_num_classes == 3 ~ "3",
      amr_num_classes >= 4 ~ "4+",
      TRUE ~ NA_character_
    ),
    
    amr_beta_lactam = factor(
      amr_beta_lactam,
      levels = c(
        "no_amr",
        "non_BL_amr",
        "other_beta_lactam",
        "cephalosporin",
        "carbapenem"
      )
    )
  )
#View(amr_categories_contigs)
table(amr_categories_contigs$amr_count)
table(amr_categories_contigs$amr_count_unique)
table(amr_categories_contigs$amr_num_classes )

table(amr_categories_contigs$amr_beta_lactam )
table(amr_categories_contigs$amr_count_category )
table(amr_categories_contigs$amr_count_unique_category )
table(amr_categories_contigs$amr_num_classes_category, useNA = "ifany")

# save
write.csv(amr_categories_contigs, "amr_categories_autocycler_contigs.csv", row.names = FALSE)
# reload if needed
#amr_categories_contigs <- read.csv("amr_categories_autocycler_contigs.csv")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# make categories of MLST, fastbaps, poppunk, and verticall clusters so colours are not overwhelming;
# make categories
all_metadata_deduplicated <- all_metadata_deduplicated |>
  group_by(species_complex) |>
  mutate(num_in_complex = n()) |>
  ungroup() |>
  group_by(kleborate_mlst_ST, species_complex) |>
  mutate(num_in_mlst = n()) |>
  ungroup() |>
  group_by(poppunk_cluster, species_complex) |> 
  mutate(num_in_poppunk = n()) |>
  ungroup() |>
  mutate(prop_in_mlst = num_in_mlst / num_in_complex,
         prop_in_poppunk = num_in_poppunk / num_in_complex
         ) |>
  mutate(freq_threshold = case_when(species_complex == "E. coli" ~ 0.02,
                                    species_complex == "K. pneumoniae complex" ~ 0.01,
                                    species_complex == "K. oxytoca complex" ~ 0.02,
                                    species_complex == "K. aerogenes" ~ 0.05,
                                    species_complex == "C. freundii complex" ~ 0.05,
                                    species_complex == "E. cloacae complex" ~ 0.05,
                                    TRUE ~ NA)) |>
  mutate(mlst_cat2 = case_when(prop_in_mlst > freq_threshold ~ kleborate_mlst_ST,
                               num_in_mlst == 1 ~ paste0(species_complex, "|Single"),
                               prop_in_mlst <= freq_threshold & num_in_mlst > 1 ~ paste0(species_complex, "|Minor")),
         poppunk_cat2 = case_when(prop_in_poppunk > freq_threshold ~ poppunk_cluster,
                                  num_in_poppunk == 1 ~ paste0(species_complex, "|Single"),
                                  prop_in_poppunk <= freq_threshold & num_in_poppunk > 1 ~ paste0(species_complex, "|Minor"))) |>
  mutate(mlst_cat2 = case_when(mlst_cat2 == "kaerogenes|4|dnaA(3)|fusA(3)|gyrB(1)|leuS(4)|pryG(3)|rplB(1)|rpoB(2)" ~ "kaerogenes|ST4",
                               mlst_cat2 == "kaerogenes|93|dnaA(3)|fusA(3)|gyrB(2)|leuS(18)|pryG(3)|rplB(1)|rpoB(10)" ~ "kaerogenes|ST93",
                               TRUE ~ mlst_cat2))

table(all_metadata_deduplicated$mlst_cat2, useNA = "ifany")
table(all_metadata_deduplicated$poppunk_cat2, useNA = "ifany")
sum(is.na(all_metadata_deduplicated$kleborate_mlst_ST))
#View(all_metadata_deduplicated_cat)
# investigate the NA MLSTs

# left join amrfinder data to all_metadata
all_metadata_deduplicated_amr <- all_metadata_deduplicated |>
  left_join(amr_categories, by = c("sequencing_id" = "sample")) 
nrow(all_metadata_deduplicated_amr) #2218
sum(is.na(all_metadata_deduplicated_amr$amr_count)) #2218
table(amr_categories$amr_num_classes_category, useNA = "ifany")
colnames(all_metadata_deduplicated_amr)

# re-categorise healthcare exposure
all_metadata_deduplicated_amr <- all_metadata_deduplicated_amr |>
  mutate(healthcare_exposure_cat = case_when(sampletype_cat == "screening" ~ "screening",
                                             TRUE ~ healthcare_exposure))

table(all_metadata_deduplicated_amr$healthcare_exposure_cat, useNA = "ifany")
table(all_metadata_deduplicated_amr$healthcare_exposure_cat, all_metadata_deduplicated_amr$sampletype_cat, useNA = "ifany")
table(all_metadata_deduplicated_amr$healthcare_exposure_cat, all_metadata_deduplicated_amr$species_complex, useNA = "ifany")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# statistical test uncategorised using: sampletype, region, healthcare exposure, and AMR gene cariage
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# define metadata and cluster variables
metadata_vars <- c(
  "sampletype_cat",
  "region",
  "healthcare_exposure_cat",
  "amr_num_classes_category",
  "amr_beta_lactam"
)

cluster_vars_cat <- c(
  "mlst_cat2",
  "poppunk_cat2")

species_complexes <- all_metadata_deduplicated_amr |>
  filter(!is.na(species_complex)) |>
  distinct(species_complex) |>
  pull(species_complex)

# Funcion to run one monte-carlo chi-squared test:
run_chisq_test <- function(
    data,
    species_complex_value,
    cluster_var,
    metadata_var,
    B = 100000
) {
  
  d <- data |>
    filter(
      species_complex == .env$species_complex_value,
      !is.na(.data[[cluster_var]]),
      !is.na(.data[[metadata_var]])
    )
  
  tab <- table(
    d[[cluster_var]],
    d[[metadata_var]]
  )
  
  # Remove empty rows and columns
  tab <- tab[
    rowSums(tab) > 0,
    colSums(tab) > 0,
    drop = FALSE
  ]
  
  # Cannot perform chi-squared if only one level remains
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    
    return(
      tibble(
        species_complex = species_complex_value,
        cluster_var = cluster_var,
        metadata_var = metadata_var,
        n = sum(tab),
        n_cluster_levels = nrow(tab),
        n_metadata_levels = ncol(tab),
        chi_squared = NA_real_,
        p_value = NA_real_,
        cramers_v = NA_real_
      )
    )
  }
  
  test <- chisq.test(
    tab,
    simulate.p.value = TRUE,
    B = B
  )
  
  n <- sum(tab)
  
  k <- min(
    nrow(tab),
    ncol(tab)
  )
  
  cramers_v <- sqrt(
    as.numeric(test$statistic) /
      (n * (k - 1))
  )
  
  tibble(
    species_complex = species_complex_value,
    cluster_var = cluster_var,
    metadata_var = metadata_var,
    n = n,
    n_cluster_levels = nrow(tab),
    n_metadata_levels = ncol(tab),
    chi_squared = unname(test$statistic),
    p_value = test$p.value,
    cramers_v = cramers_v
  )
}



#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# run all categorised cluster test
set.seed(123)

chisq_cat <- tidyr::crossing(
  species_complex_value = species_complexes,
  cluster_var = cluster_vars_cat,
  metadata_var = metadata_vars
) |>
  purrr::pmap_dfr(
    \(species_complex_value, cluster_var, metadata_var)
    run_chisq_test(
      data = all_metadata_deduplicated_amr,
      species_complex_value = species_complex_value,
      cluster_var = cluster_var,
      metadata_var = metadata_var,
      B = 100000
    )
  )

# bonferroni correction separately:
chisq_cat <- chisq_cat |>
  # group_by(species_complex) |> # remove if don;t want separate bonferroni correction per sc
  mutate(
    p_bonferroni = p.adjust(p_value, method = "bonferroni"),
    significant_bonferroni = p_bonferroni < 0.05) |>
  ungroup() |>
  arrange(species_complex, p_bonferroni)
#View(chisq_cat)

# save
write.csv(chisq_cat, "neksus_cluster_metadata_chisq_categorised.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~#
# observed, expected and standardised residuals:
# function to return category-level info - which features over-represented (species-complex aware)
get_chisq_cells_species <- function(
    data,
    species_complex_value,
    cluster_var,
    metadata_var
) {
  
  d <- data |>
    filter(
      species_complex == .env$species_complex_value,
      !is.na(.data[[cluster_var]]),
      !is.na(.data[[metadata_var]])
    )
  
  tab <- table(
    d[[cluster_var]],
    d[[metadata_var]]
  )
  
  # Remove empty rows / columns
  tab <- tab[
    rowSums(tab) > 0,
    colSums(tab) > 0,
    drop = FALSE
  ]
  
  # Need at least 2 x 2 table
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(tibble())
  }
  
  # Standard chi-square object used here to extract
  # expected counts and residuals
  test <- suppressWarnings(
    chisq.test(tab)
  )
  
  observed_df <- as.data.frame(as.table(tab)) |>
    rename(
      cluster = Var1,
      metadata_value = Var2,
      observed = Freq
    )
  
  expected_df <- as.data.frame(as.table(test$expected)) |>
    rename(
      cluster = Var1,
      metadata_value = Var2,
      expected = Freq
    )
  
  residual_df <- as.data.frame(as.table(test$residuals)) |>
    rename(
      cluster = Var1,
      metadata_value = Var2,
      pearson_residual = Freq
    )
  
  stdres_df <- as.data.frame(as.table(test$stdres)) |>
    rename(
      cluster = Var1,
      metadata_value = Var2,
      std_residual = Freq
    )
  
  observed_df |>
    left_join(
      expected_df,
      by = c("cluster", "metadata_value")
    ) |>
    left_join(
      residual_df,
      by = c("cluster", "metadata_value")
    ) |>
    left_join(
      stdres_df,
      by = c("cluster", "metadata_value")
    ) |>
    mutate(
      species_complex = species_complex_value,
      cluster_var = cluster_var,
      metadata_var = metadata_var,
      
      observed_expected_ratio =
        observed / expected,
      
      percent_difference_expected =
        100 * (observed - expected) / expected,
      
      residual_p =
        2 * pnorm(-abs(std_residual))
    ) |>
    select(
      species_complex,
      cluster_var,
      metadata_var,
      cluster,
      metadata_value,
      observed,
      expected,
      pearson_residual,
      std_residual,
      observed_expected_ratio,
      percent_difference_expected,
      residual_p
    )
}


#~~~~~~~~~~~~~~~~#
# call function
# run for categorical cluster variables
chisq_cells_cat <- tidyr::crossing(
  species_complex_value = species_complexes,
  cluster_var = cluster_vars_cat,
  metadata_var = metadata_vars
) |>
  purrr::pmap_dfr(
    \(species_complex_value, cluster_var, metadata_var)
    get_chisq_cells_species(
      data = all_metadata_deduplicated_amr,
      species_complex_value = species_complex_value,
      cluster_var = cluster_var,
      metadata_var = metadata_var
    )
  )
# separete bonferroni correction for cell-level residual p-vales
chisq_cells_cat <- chisq_cells_cat |>
  group_by(species_complex) |>
  mutate(residual_p_bonferroni = p.adjust(residual_p, method = "bonferroni"),
         residual_significant_bonferroni = residual_p_bonferroni < 0.05) |>
  ungroup()
#View(chisq_cells_cat)  # 103/2,198 significant

# save
write.csv(chisq_cells_cat, "neksus_cluster_metadata_chisq_categorised_residuals.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plotting df
# reshape to long
plot_meta_long <- all_metadata_deduplicated_amr |>
  pivot_longer(
    cols = all_of(metadata_vars),
    names_to = "metadata_variable",
    values_to = "metadata_value"
  ) |>
  filter(
    !is.na(species_complex),
    !is.na(metadata_value)
  )

# separately reshape cluster type
plot_df_long <- plot_meta_long |>
  pivot_longer(
    cols = all_of(cluster_vars_cat),
    names_to = "cluster_method",
    values_to = "cluster"
  ) |>
  filter(
    !is.na(cluster)
  ) |>
  mutate(
    cluster_method = recode(
      cluster_method,
      mlst_cat2 = "MLST",
      poppunk_cat2 = "PopPUNK"),
    
    metadata_variable = recode(
      metadata_variable,
      sampletype_cat = "Sample type",
      region = "Region",
      healthcare_exposure_cat = "Healthcare exposure",
      amr_num_classes_category = "AMR classes",
      amr_beta_lactam = "β-lactam AMR"
    )
  )

# calculate proportions
plot_df_summary <- plot_df_long |>
  dplyr::count(
    species_complex,
    cluster_method,
    metadata_variable,
    metadata_value,
    cluster,
    name = "n"
  ) |>
  dplyr::group_by(
    species_complex,
    cluster_method,
    metadata_variable,
    metadata_value
  ) |>
  dplyr::mutate(
    proportion = n / sum(n)
  ) |>
  dplyr::ungroup() |>
  # remove species prefix from cluster value
  mutate(cluster = str_remove(cluster, "^.*?\\|")) |>
  mutate(cluster = gsub("STST", "ST", cluster))
#View(plot_df_summary)

# set order of metadata variables
plot_df_summary$metadata_variable <- factor(plot_df_summary$metadata_variable, levels = c("Sample type", "Region", "Healthcare exposure", "AMR classes", "β-lactam AMR"))
#colnames(plot_df_summary)
unique(plot_df_summary$metadata_value)

# fix BL carriage category
plot_df_summary <- plot_df_summary |>
  mutate(metadata_value = case_when(metadata_value == "non_BL_amr" ~ "non-BL-AMR",
                                    metadata_value == "other_beta_lactam" ~ "other-BL",
                                    metadata_value == "" ~ "non-BL-AMR",
                                    TRUE ~ metadata_value))

# set metadata values order
hcai_order <- c("community", "quasi-community", "quasi-nosocomial", "nosocomial", "unknown", "screening")
bl_carriage_order <- c("non-BL-AMR", "other-BL", "cephalosporin", "carbapenem")
amr_class_num_order <- c("1", "2", "3", "4+")
hcai_bl_amr <- c(hcai_order, bl_carriage_order, amr_class_num_order)

metadata_values_no_hcai_bl <- setdiff(unique(plot_df_summary$metadata_value), hcai_bl_amr)
metadata_values_order <- c(hcai_bl_amr, metadata_values_no_hcai_bl)
#set order
plot_df_summary$metadata_value <- factor(plot_df_summary$metadata_value, levels = metadata_values_order)

# set order of summary df to have minor and single before all else
cluster_values_no_single <- setdiff(unique(plot_df_summary$cluster), c("Minor", "Single"))
cluster_values_ordered <- c("Minor", "Single", cluster_values_no_single )
plot_df_summary$cluster <- factor(plot_df_summary$cluster, levels = cluster_values_ordered)


# function to make one plot:
make_cluster_plot <- function(
    data,
    species,
    method,
    colours,
    legend_ncol = 8
) {
  
  d <- data |>
    filter(
      species_complex == species,
      cluster_method == method
    )
  
  # Keep only colours actually used in this panel
  used_clusters <- unique(d$cluster)
  colours <- colours[names(colours) %in% used_clusters]
  
  # Add a fallback for any categories missing from the palette
  missing_clusters <- setdiff(
    used_clusters,
    names(colours)
  )
  
  if (length(missing_clusters) > 0) {
    fallback <- setNames(
      rep("#999999", length(missing_clusters)),
      missing_clusters
    )
    
    colours <- c(colours, fallback)
  }
  
  ggplot(d,  aes(x = metadata_value, y = proportion, fill = cluster)) +
    geom_col(
      width = 0.8
    ) +
    facet_wrap(
      ~ metadata_variable,
      scales = "free_x",
      nrow = 1
    ) +
    scale_fill_manual(
      values = colours,
      drop = FALSE
    ) +
    guides(fill = guide_legend(ncol = legend_ncol, byrow = TRUE)) +
    scale_y_continuous(
      labels = percent_format(),
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = paste(species, "-", method),
      x = NULL,
      y = "Proportion",
      fill = method
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(
        face = "bold"
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      plot.title = element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    )
}

#~~~~~~~~~~~~~~~~~~~~~~#
# combined plot 
# for MLST

# define mlst colours:
# species specific palette handling
ecoli_mlst_colours <- c(
  "ST131"  = "#4477AA",
  "ST73"   = "#EE6677",
  "ST69"   = "#228833",
  "ST95"   = "#66CCEE",
  "ST12"   = "#AA3377",
  "ST127"  = "#EE7733",
  "ST1193" = "#CCBB44",
  "ST10"   = "#44AA99",
  "Minor"  = "#BDBDBD",
  "Single" = "#E6E6E6"
)

#kleb
kpsc_mlst_colours <- c(
  "ST45" = "#4477AA",
  "ST14" = "#EE6677",
  "ST29" = "#228833",
  "ST23" = "#CCBB44",
  "ST17" = "#AA3377",
  "ST307"= "#EE7733",
  "ST54" = "#44AA99",
  "ST13" = "#882255",
  "ST20" = "#999933",
  "ST870" ="#332288",
  "ST147" ="#117733",
  "ST25" = "#DDCC77",
  "ST35" = "#CC6677",
  "ST133" ="#66CCEE",  
  "ST461" ="#937860", 
  "ST834" ="#88CCEE",
  
  "Minor" ="#BDBDBD",
  "Single"="#E6E6E6"
)

kosc_mlst_colours <- c(
  "ST31"     = "#66CCEE",
  "ST37"     = "#4477AA",
  "ST157"    = "#EE6677",
  "ST88"     = "#228833",
  "ST108"    = "#CCBB44",
  "ST20"     = "#AA3377",
  "ST214-2LV"= "#EE7733",
  "ST237"    = "#44AA99",
  "ST282"    = "#882255",
  "ST362"    = "#999933",
  "ST40"      ="#332288",
  "ST58"      ="#117733",
  "Single"    ="#E6E6E6"
)

# call plots
p_ecoli <- make_cluster_plot( plot_df_summary,
                              "E. coli", "MLST",  
                              ecoli_mlst_colours, 
                              legend_ncol = 15
)
p_ecoli
ggsave("ecoli_mlst_composition_plot.png", p_ecoli, units = "in", width = 14, height = 5, dpi = 300)


p_kpsc <- make_cluster_plot(plot_df_summary,
                            "K. pneumoniae complex",
                            "MLST",
                            kpsc_mlst_colours,
                            legend_ncol = 18
)
p_kpsc
ggsave("kpsc_mlst_composition_plot.png", p_kpsc, units = "in", width = 14, height = 5.5, dpi = 300)

p_kosc <- make_cluster_plot(plot_df_summary,
  "K. oxytoca complex",
  "MLST",
  kosc_mlst_colours,
  legend_ncol = 15
)
p_kosc
ggsave("kosc_mlst_composition_plot.png", p_kosc, units = "in", width = 14, height = 5, dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# for popunk
# define poppunk colours
ecoli_poppunk_colours <- c(
  "ecoli_4" = "#7AA5D9",
  "ecoli_25" = "#C4DCFF",
  "ecoli_8" = "#FFBCC1",
  "ecoli_10" = "#7CD285",
  "ecoli_9" = "#B8E9FF",
  "ecoli_17" = "#FF93C9",
  "ecoli_22" = "#FFC2B0",
  "ecoli_5" = "#86E1CF",
  "ecoli_15" = "#F2DF9C",
  "Minor" = "#BDBDBD",
  "Single" = "#E6E6E6"
)

kpsc_poppunk_colours <- c(
  "kpsc_3" = "#97C1F5",
  "kpsc_12" = "#FFBCC1",
  "kpsc_10" = "#7CD285",
  "kpsc_11_108_277" =   "#F2E174",
  "kpsc_37" = "#B8E9FF",
  "kpsc_297" = "#FF93C9",
  "kpsc_48_128" = "#FFC2B0",
  "kpsc_8" = "#86E1CF",
  "kpsc_13" = "#EC8CB3" ,
  "kpsc_16" =  "#D3D37A",
  "kpsc_298" =  "#9D95E3",
  "kpsc_29" =  "#7AC98A" ,
  "kpsc_36" =  "#F2DF9C",
  "kpsc_53" =  "#7AA5D9",
  "kpsc_4" = "#FFB2BD",
  "Minor" = "#BDBDBD",
  "Single" = "#E6E6E6"
)


kosc_poppunk_colours <- c(
  "kosc_148"= "#97C1F5",
  "kosc_12" = "#FFBCC1",
  "kosc_30" = "#7CD285",
  "kosc_134" =   "#F2E174",
  "kosc_10" = "#B8E9FF",
  "kosc_149" = "#FF93C9",
  "kosc_150" = "#FFC2B0",
  "kosc_81" = "#86E1CF",
  "kosc_51" = "#EC8CB3" ,
  "kosc_98" =  "#D3D37A",
  "kosc_40" =   "#9D95E3",
  "Minor" = "#BDBDBD",
  "Single" = "#E6E6E6"
)

#call plots
p_ecoli <- make_cluster_plot(
  plot_df_summary,
  "E. coli",
  "PopPUNK",
  ecoli_poppunk_colours,
  legend_ncol = 15
)
p_ecoli
ggsave("ecoli_poppunk_composition_plot.png", p_ecoli, units = "in", width = 15, height = 5, dpi = 300)


p_kpsc <- make_cluster_plot(
  plot_df_summary,
  "K. pneumoniae complex",
  "PopPUNK",
  kpsc_poppunk_colours,
  legend_ncol = 15
)
p_kpsc
ggsave("kpsc_poppunk_composition_plot.png", p_kpsc, units = "in", width = 15, height = 5.5, dpi = 300)

p_kosc <- make_cluster_plot(
  plot_df_summary,
  "K. oxytoca complex",
  "PopPUNK",
  kosc_poppunk_colours,
  legend_ncol = 15
)
p_kosc
ggsave("kosc_poppunk_composition_plot.png", p_kosc, units = "in", width = 15, height = 5, dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Indovidual cluster investigation ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ecoli ST1193
st1193 <- samples_assembled_with_species_hcai_deduplicated |>
  filter(kleborate_mlst_ST == "ecoli|ST1193")
#View(st1193)  #36, linked to SEA reigon, visually appear tobe nosocomiallly-associated


#  ST834 associations 
st834 <- samples_assembled_with_species_hcai_deduplicated |>
  filter(kleborate_mlst_ST == "kpsc|ST834")
#View(st834)

#  Kosc ST31 associations 
st31 <- samples_assembled_with_species_hcai_deduplicated |>
  filter(kleborate_mlst_ST == "kosc|ST31")
#View(st31)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Duplicate comparisons ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data
all_duplicates_snp_distances <- read.csv("all_duplicates_snp_distances.csv")
#View(all_duplicates_snp_distances)
nrow(all_duplicates_snp_distances) # 117
table(all_duplicates_snp_distances$duplicate_type) # 31 and 86

# summarise cgSNV distances for both re-cultures and re-sampled duplicates
table(all_duplicates_snp_distances$core_snp_dist)
table(all_duplicates_snp_distances[all_duplicates_snp_distances$duplicate_type == "re-cultured",]$core_snp_dist)
table(all_duplicates_snp_distances[all_duplicates_snp_distances$duplicate_type == "re-sampled",]$core_snp_dist)


# make empirical tail threhsolds
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# calculate X% tail probability thresholds (SNV value below which X% of duplicates lie) - SNV count and rate

close_clusters_duplicate_pairs <- read.csv("close_clusters_duplicate_pairs.csv")

duplicate_thresholds <- close_clusters_duplicate_pairs |>
  filter(!is.na(snvs_per_time_mean)) |>
  group_by(species_complex) |>
  summarise(
    # absolute counts
    snvs_02.5 = quantile(core_snp_dist, 0.025,  na.rm = TRUE, type = 1 ),
    snvs_05 = quantile(core_snp_dist,  0.05, na.rm = TRUE, type = 1 ),
    snvs_25 = quantile(core_snp_dist, 0.25, na.rm = TRUE, type = 1 ),
    snvs_50 = quantile(core_snp_dist, 0.50,na.rm = TRUE, type = 1 ),
    snvs_75 = quantile(core_snp_dist,  0.75,na.rm = TRUE, type = 1 ),
    snvs_95 = quantile(core_snp_dist,  0.95,na.rm = TRUE,type = 1 ),
    snvs_97.5 = quantile(core_snp_dist, 0.975,na.rm = TRUE, type = 1 ),
    n_duplicates = n(),
    .groups = "drop"
  )
#View(duplicate_thresholds)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# species order
species_order <- c("E. coli",
                   "K. pneumoniae complex",
                   "K. oxytoca complex",
                   "K. aerogenes" ,
                   "C. freundii complex" ,
                   "E. cloacae complex",
                   "S. marcescens")

species_colours <- c(
  "E. coli" = "#009E73",
  "K. pneumoniae complex" = "#D55E00",
  "K. oxytoca complex" = "#CC79A7",
  "K. aerogenes" = "#7F3C8D",
  "C. freundii complex" = "#E69F00",
  "E. cloacae complex" = "#0072B2",
  "S. marcescens" = "#8B1A1A"
)

all_duplicates_snp_distances$duplicate_type <- factor(all_duplicates_snp_distances$duplicate_type, 
                                                      levels = c("re-cultured", "re-sampled"))
all_duplicates_snp_distances$species_complex <- factor(all_duplicates_snp_distances$species_complex, 
                                                       levels = species_order)
duplicate_thresholds$species_complex <- factor(duplicate_thresholds$species_complex, 
                                               levels = species_order)


# plot
duplicate_snv_dist_hist <- ggplot(data = all_duplicates_snp_distances, aes(x = core_snp_dist, colour = species_complex, fill = species_complex, alpha = duplicate_type)) +
  geom_histogram(binwidth = 1, colour = NA) +
  geom_vline(data = duplicate_thresholds, aes(xintercept = (snvs_95 + 0.5)),
             linetype = "dashed",inherit.aes = FALSE, linewidth = 0.3) +
  facet_wrap(duplicate_type ~ species_complex, scale = "free_y", ncol = 7) +
  scale_fill_manual(values = species_colours, na.value = "grey70") +
  scale_alpha_manual(values = c(0.3, 0.8)) +
  labs(x = "Core genome SNVs",
       y = "Count",
       fill = "Speices complex", 
       alpha = "Duplicate type") +
  theme_minimal()
duplicate_snv_dist_hist
ggsave("duplicate_snv_count_dist_hist_with_technical_duplicates.png", duplicate_snv_dist_hist, units = "in", width = 15, height = 5, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Match conitgs and plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# read in matched duplicate contig data
technical_duplicate_contigs_merged <- read.csv("technical_duplicate_contigs_merged.csv")
#~~~~~~~~~~~~~~~~~~#

genomad_molecule_type_colours <- c("virus" = "#F0E442",
                                   "plasmid" = "#E69F00",
                                   "chromosome" = "#0072B2")

technical_duplicate_contigs_merged$molecule_type_genomad_majority <- factor(technical_duplicate_contigs_merged$molecule_type_genomad_majority,
                                                                            levels = c("virus", "plasmid", "chromosome"))
# plot
technical_duplicate_contigs_plot <- ggplot(technical_duplicate_contigs_merged, aes(x = autocycler_contig_length.y , y = autocycler_contig_length.x, colour = molecule_type_genomad_majority)) +
  geom_point(size = 3, alpha = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, alpha = 0.5, linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = genomad_molecule_type_colours) +
  labs(x = "Retained assembly contig length",
       y = "Duplicate assembly contig length",
       colour = "geNomad\nmolecule type"
  ) +
  theme_minimal()
technical_duplicate_contigs_plot
ggsave("technical_duplicate_contigs_plot.png", technical_duplicate_contigs_plot, units = "in", width = 7, height = 6, dpi = 300)


#~~~~~~~~~~~~~~#
# read in technical duplicates plasmids distances
technical_duplicate_plasmids_distances <- read.csv("technical_duplicate_plasmids_distances.csv")

# plot histogram of distances
technical_duplicate_plasmid_distances <- ggplot(technical_duplicate_plasmids_distances, aes(x = distance, fill = molecule_type_genomad_majority)) +
  geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 4.5, alpha = 0.5, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = genomad_molecule_type_colours) +
  labs(x = "PLING DCJ/indel distance between main and duplicate contig pair",
       y = "Count",
       fill = "geNomad\nmolecule type"
  ) +
  theme_minimal()
technical_duplicate_plasmid_distances
ggsave("technical_duplicate_plasmid_distances_hist.png", technical_duplicate_plasmid_distances, units = "in", width = 6, height = 3, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Biological duplicate plots
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data if needed
biological_duplicate_contigs_merged <-  read.csv("biological_duplicate_contigs_merged.csv")

#plot
genomad_molecule_type_colours <- c("virus" = "#F0E442",
                                   "plasmid" = "#E69F00",
                                   "chromosome" = "#0072B2")

# plot
biological_duplicate_contigs_plot <- ggplot(biological_duplicate_contigs_merged, aes(x = autocycler_contig_length.y , y = autocycler_contig_length.x, colour = molecule_type_genomad_majority)) +
  geom_point(size = 3, alpha = 0.6) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, alpha = 0.5, linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = genomad_molecule_type_colours) +
  labs(x = "Retained assembly contig length",
       y = "Duplicate assembly contig length",
       colour = "geNomad\nmolecule type"
  ) +
  theme_minimal()
biological_duplicate_contigs_plot
ggsave("biological_duplicate_contigs_plot.png", biological_duplicate_contigs_plot, units = "in", width = 7, height = 6, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot histogram of plasmid/ virus pling distances
biological_duplicate_plasmid_distances <- ggplot(biological_duplicate_contigs_merged, aes(x = distance, fill = molecule_type_genomad_majority)) +
  geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 4.5, alpha = 0.5, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = genomad_molecule_type_colours) +
  labs(x = "PLING DCJ/indel distance between main and duplicate contig pair",
       y = "Count",
       fill = "geNomad\nmolecule type"
  ) +
  theme_minimal()
biological_duplicate_plasmid_distances
ggsave("biological_duplicate_plasmid_distances_hist.png", biological_duplicate_plasmid_distances, units = "in", width = 6, height = 3, dpi = 300)

# reload sumary data
biological_duplicate_isolate_summary <- read.csv("biological_duplicate_isolate_summary.csv")

# plot core genome SNP vs dcj indel distance
cg_snv_vs_dcj_scatter <- ggplot(biological_duplicate_isolate_summary, aes(x = max_core_snv_dist, y = max_dcj_dist)) +
  geom_point(position = position_jitter(width = 0.1, height = 0.1), alpha = 0.7) +
  labs(x = "core genome SNV distance",
       y = "Largest DCJ/indel distance"
  ) +
  theme_minimal()
cg_snv_vs_dcj_scatter
ggsave("biological_duplicate_cgSNV_vs_DCJ_dist_scatter_plot.png", cg_snv_vs_dcj_scatter, units = "in", width = 5, height = 4, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Chi-squared test and tabulation to compare composition of duplicates to main NEKSUS data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~#
# tabulate duplicates to see if they are representative of 
#chi-squared test
table(all_duplicates_snp_distances$species_complex, all_duplicates_snp_distances$duplicate_type)

# tabulate main df species complexes
neksus_main_species_complex <- samples_assembled_with_species_hcai_deduplicated |>
  mutate(species_complex = case_when(grepl("ecoli_achtman_4", typing_scheme) ~ "E. coli",
                                     grepl("klebsiella", typing_scheme) ~ "K. pneumoniae complex",
                                     grepl("koxytoca", typing_scheme) ~ "K. oxytoca complex",
                                     grepl("kaerogenes", typing_scheme) ~ "K. aerogenes",
                                     grepl("cfreundii", typing_scheme) ~ "C. freundii complex",
                                     grepl("ecloacae", typing_scheme) ~ "E. cloacae complex",
                                     sequencing_id == "AF31" ~ "S. marcescens",
                                     TRUE ~ NA
  )) |>
  group_by(sequencing_id) |> 
  mutate(count = n()) |>
  filter(!(count >1 & !grepl("failed", run))) |>
  ungroup() |>
  filter(assembly_qc_pass == TRUE,
         duplicate_assembly_qc_pass != "duplicate") 
nrow(neksus_main_species_complex) # 2119
table(neksus_main_species_complex$species_complex, useNA = "ifany")


tech_dup_freq <- c(21,7,2,0,0,0,1)
biol_dup_freq <- c(59,20,2,1,1,3,0)
neksus_main_freq <- c(1559,406,91,21,20,21,1)
species_complex_names <- c("E. coli",
                           "K. pneumoniae complex",
                           "K. oxytoca complex",
                           "K. aerogenes",
                           "C. freundii complex",
                           "E. cloacae complex",
                           "S. marcescens"
)



tech_dup_freq <- c(21,9,1)
biol_dup_freq <- c(59,23,4)
neksus_main_freq <- c(1559,518,42)

species_complex_names <- c("E. coli", "Klebsiella spp.", "Other")
names(tech_dup_freq) <- species_complex_names
names(biol_dup_freq) <- species_complex_names
names(neksus_main_freq) <- species_complex_names

set.seed(123)
tech_vs_main <- rbind(tech_dup = tech_dup_freq,  neksus_main = neksus_main_freq)
chisq.test(tech_vs_main, simulate.p.value = TRUE, B = 100000)
biol_vs_main <- rbind(biol_dup = biol_dup_freq,  neksus_main = neksus_main_freq)
chisq.test(biol_vs_main, simulate.p.value = TRUE, B = 100000)

# test all 3
tech_biol_main <- rbind(tech_dup = tech_dup_freq, biol_dup = biol_dup_freq, neksus_main_freq)
chisq.test(tech_biol_main, simulate.p.value = TRUE, B = 100000)

# sampletype distribution
table(neksus_main_species_complex$sampletype_cat, useNA = "ifany")

twice_seq <- samples_assembled_with_species_hcai_deduplicated |>
  filter(sequencing_id %in% technical_duplicate_snp_pairs$genome_a)
twice_seq$sequencing_id
table(twice_seq$sampletype_cat) # 16 blood, 15 screening

biol_dup <- samples_assembled_with_species_hcai_deduplicated |>
  filter(sequencing_id %in% biological_duplicate_snp_pairs$genome_a,
         assembly_qc_pass == TRUE) |>
  select(sequencing_id, sampletype)
#View(biol_dup)
nrow(biol_dup)
table(biol_dup$sampletype, useNA = "ifany")
#View(duplicates_check)
biol_dup_merged <- biological_duplicate_snp_pairs |>
  left_join(biol_dup, by = c("genome_a" = "sequencing_id"))
table(biol_dup_merged$sampletype)
# 67 ; 19


tech_sampletype_freq <- c(16,15)
biol_sampletype_freq <- c(67,19)
neksus_main_sampletype_freq <- c(1944,175)

sampletype_names <- c("blood", "screening")
names(tech_sampletype_freq) <- sampletype_names
names(biol_sampletype_freq) <- sampletype_names
names(neksus_main_sampletype_freq) <- sampletype_names

tech_vs_main <- rbind(tech_dup = tech_sampletype_freq,  neksus_main = neksus_main_sampletype_freq)
chisq.test(tech_vs_main, simulate.p.value = TRUE, B = 100000)
biol_vs_main <- rbind(biol_dup = biol_sampletype_freq,  neksus_main = neksus_main_sampletype_freq)
chisq.test(biol_vs_main, simulate.p.value = TRUE, B = 100000)


# test all 3
tech_biol_main <- rbind(tech_dup = tech_sampletype_freq, biol_dup = biol_sampletype_freq, neksus_main_sampletype_freq)
chisq.test(tech_biol_main, simulate.p.value = TRUE, B = 100000)
#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Nearest neighbour cgSNV analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data
nearest_neighbour_snv_count <- read.csv("nearest_neighbour_snv_count.csv")
nrow(nearest_neighbour_snv_count) # 2133, as includes the ties in nn distance
colnames(nearest_neighbour_snv_count)
# filter to keep only main TS for plott
nn_snv_count_filt <- nearest_neighbour_snv_count |>
  filter(species_complex %in% c("E. coli", "K. pneumoniae complex", "K. oxytoca complex"))

table(nn_snv_count_filt$species_complex, useNA = "ifany")

# make plot dfs (use just the count df)
plot_df_count <- nn_snv_count_filt |>
  mutate(nn_snv_count_plot = case_when(nn_snv_count == 0 ~ 0.5,
                                       TRUE ~ nn_snv_count))
plot_df_rate <- nn_snv_count_filt |>
  filter(!is.na(nn_species_complex)) |>
  mutate(nn_snv_per_day_plot = case_when(nn_snv_per_day == 0 ~ 0.0025,
                                         TRUE ~ nn_snv_per_day))

# set factor orders:
plot_df_count$pair_same_region <- factor(plot_df_count$pair_same_region, levels = c(TRUE, FALSE))
plot_df_count$pair_sampletype_combo <- factor(plot_df_count$pair_sampletype_combo, levels = c("blood|blood","screening|screening", "blood|screening"))
plot_df_count$species_complex <- factor(plot_df_count$species_complex, levels = c("E. coli",
                                                                                  "K. pneumoniae complex",
                                                                                  "K. oxytoca complex"))
plot_df_rate$pair_same_region <- factor(plot_df_rate$pair_same_region, levels = c(TRUE, FALSE))
plot_df_rate$pair_sampletype_combo <- factor(plot_df_rate$pair_sampletype_combo, levels = c("blood|blood","screening|screening", "blood|screening"))
plot_df_rate$species_complex <- factor(plot_df_rate$species_complex, 
                                       levels = c("E. coli",
                                                  "K. pneumoniae complex",
                                                  "K. oxytoca complex"))
duplicate_thresholds_filt <- duplicate_thresholds |>
  filter(species_complex %in% c("E. coli", "K. pneumoniae complex", "K. oxytoca complex"))
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SAMPLE TYPE
sampletype_cols <- c(
  "blood|blood" = "darkred",
  "screening|screening" = "goldenrod1",
  "blood|screening" = "#7B3294",
  "Unknown" = "grey70"
)

# SNV count
close_clusters_hist <- ggplot(data = plot_df_count, aes(x = nn_snv_count_plot, fill = pair_sampletype_combo)) +
  geom_histogram(binwidth = 0.05, colour = NA, alpha = 0.8) +
  #geom_vline(xintercept = 20, linetype = "dashed", linewidth = 0.3) +
  geom_vline(data = duplicate_thresholds_filt, aes(xintercept = snvs_95), linetype = "dashed", inherit.aes = FALSE, linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = sampletype_cols) +
  scale_x_log10(limits = c(-0.05, 500)) +
  labs(x = "Core genome SNVs",
       y = "Count",
       fill = "Sample type") +
  theme_minimal()
print(close_clusters_hist)
ggsave("nearest_neighbour_cgSNV_dist_hist_by_sampletype_log_trunc500.png", close_clusters_hist, units = "in", width = 14, height = 3, dpi = 300)


# SNV/time
close_clusters_hist <- ggplot(data = plot_df_rate, aes(x = nn_snv_per_day_plot, fill = pair_sampletype_combo)) +
  geom_histogram(binwidth = 0.025, colour = NA, alpha = 0.8) +
  geom_vline(xintercept = 0.003, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0.4, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = sampletype_cols) +
  scale_x_log10(limits = c(-0.006, 5)) +
  labs(x = "Core genome SNVs/day",
       y = "Count",
       fill = "Sample type") +
  theme_minimal()
close_clusters_hist
ggsave("nearest_neighbour_cgSNV_per_day_hist_by_sampletype_log.png", close_clusters_hist, units = "in", width = 14, height = 3, dpi = 300)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# REGION 
region_pair_cols = c("TRUE" = "#DDCC77", 
                     "FALSE" = "#332288")
# SNV count
close_clusters_hist <- ggplot(data = plot_df_count, aes(x = nn_snv_count_plot, fill = pair_same_region)) +
  geom_histogram(binwidth = 0.05, colour = NA, alpha = 1) +
  geom_vline(data = duplicate_thresholds_filt, aes(xintercept = snvs_95), linetype = "dashed", inherit.aes = FALSE, linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = region_pair_cols) +
  scale_x_log10(limits = c(-0.05, 500)) +
  labs(x = "Core genome SNVs",
       y = "Count",
       fill = "Same region") +
  theme_minimal()
close_clusters_hist
ggsave("nearest_neighbour_cgSNV_dist_hist_by_same_region_log_trunc500.png", close_clusters_hist, units = "in", width = 14, height = 3,dpi = 300)

# SNV/time
close_clusters_hist <- ggplot(data = plot_df_rate, aes(x = nn_snv_per_day_plot, fill = pair_same_region)) +
  geom_histogram(binwidth = 0.025, colour = NA, alpha = 1) +
  geom_vline(xintercept = 0.003, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0.4, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = region_pair_cols) +
  scale_x_log10(limits = c(-0.006, 5)) +
  labs(x = "Core genome SNVs/day",
       y = "Count",
       fill = "Same region") +
  theme_minimal()
close_clusters_hist
ggsave("nearest_neighbour_cgSNV_per_day_hist_by_same_region_log.png", close_clusters_hist, units = "in", width = 14, height = 3, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# HCAI
colours_hcai <- c(
  "community" = "#FFF3B0",
  "quasi-community" = "#FFD37F",
  "quasi-nosocomial" = "#FCA17D",
  "nosocomial" = "#C75197",
  "unknown" = "grey"
)
# colour by HEALTHCARE EXPOSURE
# SNV count
close_clusters_hist <- ggplot(data = plot_df_count, aes(x = nn_snv_count_plot, fill = nn_hcai)) +
  geom_histogram(binwidth = 0.05, colour = NA, alpha = 0.95) +
  geom_vline(data = duplicate_thresholds_filt , aes(xintercept = snvs_95), linetype = "dashed", inherit.aes = FALSE, linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = colours_hcai) +
  scale_x_log10(limits = c(-0.05, 500)) +
  labs(x = "Core genome SNVs",
       y = "Count",
       fill = "Nearest neighbour\nhealthcare exposure") +
  theme_minimal()
print(close_clusters_hist)
ggsave("nearest_neighbour_cgSNV_dist_hist_by_hcai_log_trunc500.png", close_clusters_hist, units = "in", width = 14, height = 3, dpi = 300)


# SNV/time
close_clusters_hist <- ggplot(data = plot_df_rate, aes(x = nn_snv_per_day_plot, fill = nn_hcai)) +
  geom_histogram(binwidth = 0.025, colour = NA, alpha = 0.95) +
  geom_vline(xintercept = 0.003, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0.4, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = colours_hcai) +
  scale_x_log10(limits = c(-0.006, 5)) +
  labs(x = "Core genome SNVs/day",
       y = "Count",
       fill = "Nearest neighbour\nhealthcare exposure") +
  theme_minimal()
close_clusters_hist
ggsave("nearest_neighbour_cgSNV_per_day_hist_by_hcai_log.png", close_clusters_hist, units = "in", width = 14, height = 3, dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot coloured by number of ABx classes resistant to
#table(plot_df_count$amr_num_classes_category, useNA = "ifany")

colours_amr_classes <- c(
  "1"  = "#DCEAF7",  # very light blue
  "2"  = "#9ECAE1",  # medium light blue
  "3"  = "#6A51A3",  # purple
  "4+" = "#3F007D"   # deep purple
)

# SNV count
close_clusters_hist <- ggplot(data = plot_df_count, aes(x = nn_snv_count_plot, fill = amr_num_classes_category)) +
  geom_histogram(binwidth = 0.05, colour = NA, alpha = 1) +
  geom_vline(data = duplicate_thresholds_filt, aes(xintercept = snvs_95), linetype = "dashed", inherit.aes = FALSE, linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = colours_amr_classes) +
  scale_x_log10(limits = c(-0.05, 500)) +
  labs(x = "Core genome SNVs",
       y = "Count",
       fill = "Nearest neighbour\nnumber of antibiotic\nclasses resistant to") +
  theme_minimal()
close_clusters_hist
ggsave("nearest_neighbour_cgSNV_dist_hist_by_num_amr_classes_log_trunc500.png", close_clusters_hist, units = "in", width = 14, height = 3,dpi = 300)

# SNV/time
close_clusters_hist <- ggplot(data = plot_df_rate, aes(x = nn_snv_per_day_plot, fill = amr_num_classes_category)) +
  geom_histogram(binwidth = 0.025, colour = NA, alpha = 0.95) +
  geom_vline(xintercept = 0.003, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0.4, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = colours_amr_classes) +
  scale_x_log10(limits = c(-0.006, 5)) +
  labs(x = "Core genome SNVs/day",
       y = "Count",
       fill = "Nearest neighbour\nnumber of antibiotic\nclasses resistant to") +
  theme_minimal()
close_clusters_hist
ggsave("nearest_neighbour_cgSNV_per_day_hist_by_num_amr_classes_log.png", close_clusters_hist, units = "in", width = 14, height = 3, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot coloured by beta-lactam resistance category
#table(plot_df_count$amr_beta_lactam, useNA = "ifany")

colours_amr_beta_lactam <- c(
  "no_amr"            = "#F2F0F7",  # very pale lavender
  "non_BL_amr"        = "#CBC9E2",  # light lavender
  "other_beta_lactam" = "#9E9AC8",  # medium purple
  "cephalosporin"     = "#756BB1",  # stronger purple
  "carbapenem"        = "#54278F"   # dark alarm purple
)

# SNV count
close_clusters_hist <- ggplot(data = plot_df_count, aes(x = nn_snv_count_plot, fill = amr_beta_lactam)) +
  geom_histogram(binwidth = 0.05, colour = NA, alpha = 1) +
  geom_vline(data = duplicate_thresholds_filt, aes(xintercept = snvs_95), linetype = "dashed", inherit.aes = FALSE, linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = colours_amr_beta_lactam) +
  scale_x_log10(limits = c(-0.05, 500)) +
  labs(x = "Core genome SNVs",
       y = "Count",
       fill = "Nearest neighbour\nBeta-lactam resistance\ncategory") +
  theme_minimal()
close_clusters_hist
ggsave("nearest_neighbour_cgSNV_dist_hist_by_BL_res_cat_log_trunc500.png", close_clusters_hist, units = "in", width = 14, height = 3,dpi = 300)


# cgSNVs/day
close_clusters_hist <- ggplot(data = plot_df_rate, aes(x = nn_snv_per_day_plot, fill = amr_beta_lactam)) +
  geom_histogram(binwidth = 0.025, colour = NA, alpha = 0.95) +
  geom_vline(xintercept = 0.003, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = 0.4, linetype = "dashed", linewidth = 0.3) +
  facet_wrap(~ species_complex, scale = "free", ncol = 6) +
  scale_fill_manual(values = colours_amr_beta_lactam) +
  scale_x_log10(limits = c(-0.006, 5)) +
  labs(x = "Core genome SNVs/day",
       y = "Count",
       fill = "Nearest neighbour\nBeta-lactam resistance\ncategory") +
  theme_minimal()
close_clusters_hist
ggsave("nearest_neighbour_cgSNV_per_day_hist_by_BL_res_cat_log.png", close_clusters_hist, units = "in", width = 14, height = 3, dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Close clusters analysis ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load data if needed
close_clusters_filtered_meta <- read.csv("neksus_close_clusters_filtered_meta.csv")

# add species specific threshold to pairwise df
pairs <- close_clusters_filtered_meta |>
  left_join(duplicate_thresholds, by = "species_complex")
nrow(pairs) # 80,480

# define 4 edge definitions
pairs <- pairs |>
  mutate(criterion_empirical95_count = !is.na(core_snp_dist) & !is.na(snvs_95) & core_snp_dist <= snvs_95, 
         criterion_rate003 = !is.na(snvs_per_time_mean) & snvs_per_time_mean <= 0.003,
         criterion_rate04 =  !is.na(snvs_per_time_mean) & snvs_per_time_mean <= 0.4)

table(pairs$criterion_empirical95_count, useNA = "ifany") # 212 links
table(pairs$criterion_rate04, useNA = "ifany") # 437 links
table(pairs$criterion_rate003, useNA = "ifany") # 44 links

# ensure each pair is undirected edge:
pairs <- pairs |>
  mutate(vertex_1 = pmin(genome_a, genome_b),
         vertex_2 = pmax(genome_a, genome_b)
  ) |>
  filter(vertex_1 != vertex_2)


# deduplicate if same pair repeated
pairs <- pairs |>
  group_by(vertex_1, vertex_2) |>
  slice(1) |>
  ungroup()

nrow(pairs) # still 80,483

# function to extract connected components:
get_components <- function(pair_df, criterion, all_isolates) {
  edges <- pair_df |>
    filter(.data[[criterion]]) |>
    select(vertex_1, vertex_2)
  
  # All isolates are vertices, including isolates with no qualifying edge
  vertices <- tibble(
    sequencing_id = unique(all_isolates)
  )
  
  if (nrow(edges) == 0) {
    return(
      vertices |>
        mutate(
          component_id = NA_integer_,
          component_size = 1L
        )
    )
  }
  
  g <- graph_from_data_frame(
    d = edges,
    directed = FALSE,
    vertices = vertices |>
      rename(name = sequencing_id)
  )
  
  comp <- components(g)
  
  tibble(
    sequencing_id = names(comp$membership),
    component_id = as.integer(comp$membership),
    component_size = as.integer(comp$csize[comp$membership])
  ) |>
    mutate(
      component_id = if_else(
        component_size >= 2,
        component_id,
        NA_integer_
      )
    )
}

# complete isolate universe:
all_isolates <- all_metadata_deduplicated |>
  distinct(sequencing_id) |>
  pull(sequencing_id)


# get component IDs for all definitions
components_empirical95_count <- get_components(pairs,"criterion_empirical95_count", all_isolates) |>
  mutate(clustering_definition = "snvs_per_time_empirical95_count")
components_rate003 <- get_components(pairs,"criterion_rate003", all_isolates) |>
  mutate(clustering_definition = "snvs_per_time_le0.003")
components_rate04 <- get_components(pairs,"criterion_rate04",all_isolates) |>
  mutate( clustering_definition = "snvs_per_time_le0.4")


# combine:
component_membership <- bind_rows(
  components_empirical95_count,
  components_rate003,
  components_rate04
)
#View(component_membership)
length(unique(component_membership$sequencing_id))

# reasonable cluster IDs
component_membership <- component_membership |>
  group_by(clustering_definition) |>
  mutate(
    component_id_new = if_else(
      !is.na(component_id),
      dense_rank(component_id),
      NA_integer_
    ),
    very_close_cluster_id = if_else(
      !is.na(component_id_new),
      paste0(
        clustering_definition,
        "_",
        str_pad(component_id_new, 3, pad = "0")
      ),
      NA_character_
    )
  ) |>
  ungroup() |>
  select(
    sequencing_id,
    clustering_definition,
    component_id = component_id_new,
    component_size,
    very_close_cluster_id
  )
# long format df saying which isolate belong to which component under which definition

# pivot into wider df
component_membership_wide <- component_membership |>
  select(
    sequencing_id,
    clustering_definition,
    very_close_cluster_id,
    component_size
  ) |>
  pivot_wider(
    names_from = clustering_definition,
    values_from = c(
      very_close_cluster_id,
      component_size
    ),
    names_sep = "__"
  )

#View(component_membership_wide)
nrow(component_membership_wide) # 2118
colnames(component_membership_wide)
table(component_membership_wide$component_size__snvs_per_time_empirical95_count, useNA = "ifany")
table(component_membership_wide$component_size__snvs_per_time_le0.003 , useNA = "ifany")
table(component_membership_wide$component_size__snvs_per_time_le0.4 , useNA = "ifany")

# save
write.csv(component_membership_wide, "close_cluster_component_membership_wide_3_definitions.csv", row.names = FALSE)
# read back in
#component_membership_wide <- read.csv("close_cluster_component_membership_wide_3_definitions.csv")

# get only qualifying pariwise connections:
edges_long <- bind_rows(pairs |>
                          filter(criterion_empirical95_count) |>
                          mutate(clustering_definition = "snvs_per_time_empirical95_count"),
                        pairs |>
                          filter(criterion_rate003) |>
                          mutate(clustering_definition = "snvs_per_time_le0.003"),
                        pairs |>
                          filter(criterion_rate04) |>
                          mutate(clustering_definition = "snvs_per_time_le0.4")
)

# assign qualifying edge to connected components
edges_with_components <- edges_long |>
  select(
    clustering_definition,
    vertex_1,
    vertex_2,
    core_snp_dist,
    snvs_per_time_mean,
    species_complex,
    everything()
  ) |>
  left_join(
    component_membership |>
      select(
        sequencing_id,
        clustering_definition,
        very_close_cluster_id
      ),
    by = c(
      "clustering_definition",
      "vertex_1" = "sequencing_id"
    )
  ) |>
  rename(
    component_id_1 = very_close_cluster_id
  ) |>
  left_join(
    component_membership |>
      select(
        sequencing_id,
        clustering_definition,
        very_close_cluster_id
      ),
    by = c(
      "clustering_definition",
      "vertex_2" = "sequencing_id"
    )
  ) |>
  rename(
    component_id_2 = very_close_cluster_id
  ) |>
  filter(
    !is.na(component_id_1),
    component_id_1 == component_id_2
  ) |>
  mutate(
    very_close_cluster_id = component_id_1
  )

# get species complex for isolates
isolate_species <- all_metadata_deduplicated |>
  select(sequencing_id, species_complex) |>
  distinct()


# join
component_membership <- component_membership |>
  left_join(isolate_species, by = "sequencing_id")

# define helpers
# define helpers
safe_median <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_quantile <- function(x, p) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, p, na.rm = TRUE))
}

safe_min <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

safe_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  max(x)
}

# define species specific labels and denominators
species_denominators <- all_metadata_deduplicated_amr |>
  dplyr::filter(!is.na(species_complex)) |>
  dplyr::distinct(sequencing_id, species_complex) |>
  dplyr::count(
    species_complex,
    name = "n_all_isolates"
  )

overall_denominator <- all_metadata_deduplicated_amr |>
  dplyr::distinct(sequencing_id) |>
  nrow()

# cluster sizes
cluster_sizes_species <- component_membership |>
  dplyr::filter(
    !is.na(very_close_cluster_id),
    !is.na(species_complex)
  ) |>
  dplyr::distinct(
    clustering_definition,
    species_complex,
    very_close_cluster_id,
    sequencing_id
  ) |>
  dplyr::count(
    clustering_definition,
    species_complex,
    very_close_cluster_id,
    name = "cluster_size"
  )

# overall cluster sizes
cluster_sizes_overall <- component_membership |>
  dplyr::filter(
    !is.na(very_close_cluster_id)
  ) |>
  dplyr::distinct(
    clustering_definition,
    very_close_cluster_id,
    sequencing_id
  ) |>
  dplyr::count(
    clustering_definition,
    very_close_cluster_id,
    name = "cluster_size"
  ) |>
  dplyr::mutate(
    species_complex = "Overall"
  )

# combine
cluster_sizes_all <- dplyr::bind_rows(
  cluster_sizes_species,
  cluster_sizes_overall
)

# summarise cluster sizes
cluster_size_summary <- cluster_sizes_all |>
  dplyr::group_by(
    clustering_definition,
    species_complex
  ) |>
  dplyr::summarise(
    n_non_singleton_clusters = dplyr::n(),
    
    cluster_size_median = safe_median(cluster_size),
    cluster_size_q1 = safe_quantile(cluster_size, 0.25),
    cluster_size_q3 = safe_quantile(cluster_size, 0.75),
    cluster_size_min = safe_min(cluster_size),
    cluster_size_max = safe_max(cluster_size),
    
    .groups = "drop"
  )


# number of isolates in non-singleton clusters
# species specific
isolates_clustered_species <- component_membership |>
  dplyr::filter(
    !is.na(very_close_cluster_id),
    !is.na(species_complex)
  ) |>
  dplyr::distinct(
    clustering_definition,
    species_complex,
    sequencing_id
  ) |>
  dplyr::count(
    clustering_definition,
    species_complex,
    name = "n_isolates_in_non_singleton_clusters"
  ) |>
  dplyr::left_join(
    species_denominators,
    by = "species_complex"
  )

# overall
isolates_clustered_overall <- component_membership |>
  dplyr::filter(
    !is.na(very_close_cluster_id)
  ) |>
  dplyr::distinct(
    clustering_definition,
    sequencing_id
  ) |>
  dplyr::count(
    clustering_definition,
    name = "n_isolates_in_non_singleton_clusters"
  ) |>
  dplyr::mutate(
    species_complex = "Overall",
    n_all_isolates = overall_denominator
  )

# combine and calculate percentages
isolates_clustered_summary <- dplyr::bind_rows(
  isolates_clustered_species,
  isolates_clustered_overall
) |>
  dplyr::mutate(
    pct_isolates_in_non_singleton_clusters =
      100 * n_isolates_in_non_singleton_clusters / n_all_isolates
  )



# pairwise links and SNV distributions
link_summary_species <- edges_with_components |>
  dplyr::filter(!is.na(species_complex)) |>
  dplyr::group_by(
    clustering_definition,
    species_complex
  ) |>
  dplyr::summarise(
    n_pairwise_links = dplyr::n(),
    
    snv_median = safe_median(core_snp_dist),
    snv_q1 = safe_quantile(core_snp_dist, 0.25),
    snv_q3 = safe_quantile(core_snp_dist, 0.75),
    snv_min = safe_min(core_snp_dist),
    snv_max = safe_max(core_snp_dist),
    
    snvs_per_day_median = safe_median(snvs_per_time_mean),
    snvs_per_day_q1 = safe_quantile(snvs_per_time_mean, 0.25),
    snvs_per_day_q3 = safe_quantile(snvs_per_time_mean, 0.75),
    snvs_per_day_min = safe_min(snvs_per_time_mean),
    snvs_per_day_max = safe_max(snvs_per_time_mean),
    
    .groups = "drop"
  )

# overall
link_summary_overall <- edges_with_components |>
  dplyr::group_by(
    clustering_definition
  ) |>
  dplyr::summarise(
    n_pairwise_links = dplyr::n(),
    
    snv_median = safe_median(core_snp_dist),
    snv_q1 = safe_quantile(core_snp_dist, 0.25),
    snv_q3 = safe_quantile(core_snp_dist, 0.75),
    snv_min = safe_min(core_snp_dist),
    snv_max = safe_max(core_snp_dist),
    
    snvs_per_day_median = safe_median(snvs_per_time_mean),
    snvs_per_day_q1 = safe_quantile(snvs_per_time_mean, 0.25),
    snvs_per_day_q3 = safe_quantile(snvs_per_time_mean, 0.75),
    snvs_per_day_min = safe_min(snvs_per_time_mean),
    snvs_per_day_max = safe_max(snvs_per_time_mean),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    species_complex = "Overall"
  )

# combine
link_summary <- dplyr::bind_rows(
  link_summary_species,
  link_summary_overall
)


# build complete grid
cluster_definitions <- unique(
  component_membership$clustering_definition
)

species_complexes <- all_metadata_deduplicated_amr |>
  dplyr::filter(!is.na(species_complex)) |>
  dplyr::distinct(species_complex) |>
  dplyr::pull(species_complex)

summary_grid <- tidyr::crossing(
  clustering_definition = cluster_definitions,
  species_complex = c(species_complexes, "Overall")
)

# combine everything
cluster_summary_final <- summary_grid |>
  dplyr::left_join(
    cluster_size_summary,
    by = c("clustering_definition",
      "species_complex"  )) |>
  dplyr::left_join(  isolates_clustered_summary,
    by = c(
      "clustering_definition",
      "species_complex"
    )
  ) |>
  dplyr::left_join( link_summary,
    by = c(
      "clustering_definition",
      "species_complex"
    )
  ) |>
  dplyr::mutate(
    n_non_singleton_clusters =  tidyr::replace_na(n_non_singleton_clusters, 0L),
    
    n_isolates_in_non_singleton_clusters =
      tidyr::replace_na(
        n_isolates_in_non_singleton_clusters,
        0L
      ),
    
    n_pairwise_links =  tidyr::replace_na(n_pairwise_links, 0L),
    
    pct_isolates_in_non_singleton_clusters =
      tidyr::replace_na(
        pct_isolates_in_non_singleton_clusters,
        0
      ),
    
    cluster_size_med_iqr = if_else(
      is.na(cluster_size_median),
      NA_character_,
      paste0(
        round(cluster_size_median, 1),
        " (",
        round(cluster_size_q1, 1),
        "-",
        round(cluster_size_q3, 1),
        ")"
      )
    ),
    
    cluster_size_range = if_else(
      is.na(cluster_size_min),
      NA_character_,
      paste0(
        sprintf("%02.0f", cluster_size_min),
        "-",
        sprintf("%02.0f", cluster_size_max)
      )
    ),
    
    isolates_in_clusters_n_pct = paste0(
      n_isolates_in_non_singleton_clusters,
      " (",
      sprintf(
        "%.1f",
        pct_isolates_in_non_singleton_clusters
      ),
      "%)"
    ),
    
    snv_med_iqr = if_else(
      is.na(snv_median),
      NA_character_,
      paste0(
        round(snv_median, 1),
        " (",
        round(snv_q1, 1),
        "-",
        round(snv_q3, 1),
        ")"
      )
    ),
    
    snv_range = if_else(
      is.na(snv_min),
      NA_character_,
      paste0(
        sprintf("%02.0f", snv_min),
        "-",
        sprintf("%02.0f", snv_max)
      )
    ),
    
    snvs_per_day_med_iqr = if_else(
      is.na(snvs_per_day_median),
      NA_character_,
      paste0(
        signif(snvs_per_day_median, 3),
        " (",
        signif(snvs_per_day_q1, 3),
        "-",
        signif(snvs_per_day_q3, 3),
        ")"
      )
    ),
    
    snvs_per_day_range = if_else(
      is.na(snvs_per_day_min),
      NA_character_,
      paste0(
        signif(snvs_per_day_min, 3),
        "-",
        signif(snvs_per_day_max, 3)
      )
    )
  ) |>
  dplyr::select(
    clustering_definition,
    species_complex,
    
    n_non_singleton_clusters,
    cluster_size_med_iqr,
    cluster_size_range,
    
    isolates_in_clusters_n_pct,
    
    n_pairwise_links,
    
    snv_med_iqr,
    snv_range,
    
    snvs_per_day_med_iqr,
    snvs_per_day_range,
    
    dplyr::everything()
  )

species_order <- c(species_order, "Overall")
cluster_summary_final$species_complex <- factor(cluster_summary_final$species_complex, levels = species_order)
#table(cluster_summary_final$species_complex)
write.csv(cluster_summary_final, "very_close_cluster_summary_final_3_defs.csv", row.names = FALSE)
cluster_summary_final <- read.csv("very_close_cluster_summary_final_3_defs.csv")
#View(cluster_summary_final)
#~~~~~~~~~~~~~~~~~~~#
 # make per-cluster summary
table(samples_assembled_with_species_hcai_deduplicated$component_size__snvs_per_time_le0.003, useNA = "ifany")
table(samples_assembled_with_species_hcai_deduplicated$component_size__snvs_per_time_le0.4, useNA = "ifany")
table(samples_assembled_with_species_hcai_deduplicated$component_size__snvs_per_time_empirical95_count, useNA = "ifany")


# Most frequent non-missing category
get_main_category <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return(NA_character_)
  
  names(sort(table(x), decreasing = TRUE))[1]
}

# Proportion belonging to the most frequent category
get_main_prop <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return(NA_real_)
  
  max(table(x)) / length(x)
}


# cluster summary function
summarise_close_clusters <- function(data,
                                     cluster_col,
                                     component_size_col,
                                     membership_col) {
  
  data |>
    group_by(sequencing_id) |>
    mutate(duplicate = n()) |>
    filter(!(duplicate > 2 & !grepl("fail", run))) |>
    ungroup() |>
    filter(duplicate_assembly_qc_pass != "duplicate") |>
    
    # Threshold-specific filtering
    filter(
      !is.na(.data[[membership_col]]),
      .data[[component_size_col]] > 2,
      !is.na(.data[[cluster_col]])
    ) |>
    
    group_by(.data[[cluster_col]]) |>
    
    summarise(
      close_cluster_size = n(),
      
      # MLST
      n_mlst = n_distinct(kleborate_mlst_ST, na.rm = TRUE),
      main_mlst = get_main_category(kleborate_mlst_ST),
      prop_main_mlst = get_main_prop(kleborate_mlst_ST),
      
      # Sample type
      n_sampletype = n_distinct(sampletype_cat, na.rm = TRUE),
      main_sample_type = get_main_category(sampletype_cat),
      prop_main_sample_type = get_main_prop(sampletype_cat),
      
      # Region
      n_region = n_distinct(region, na.rm = TRUE),
      main_region = get_main_category(region),
      prop_max_region = get_main_prop(region),
      
      # Healthcare exposure
      n_healthcare_exposure = n_distinct(
        healthcare_exposure_cat,
        na.rm = TRUE
      ),
      main_healthcare_exposure = get_main_category(
        healthcare_exposure_cat
      ),
      prop_main_healthcare_exposure = get_main_prop(
        healthcare_exposure_cat
      ),
      
      .groups = "drop"
    ) |>
    mutate(main_mlst_prop = paste0(main_mlst, " (", round(prop_main_mlst*100) , "%)"),
           main_sampletype_prop = paste0(main_sample_type, " (", round(prop_main_sample_type*100) , "%)"),
           
           main_region_prop = paste0(main_region, " (", round(prop_max_region *100) , "%)"),
           main_healthcare_exposure_prop = paste0(main_healthcare_exposure, " (", round(prop_main_healthcare_exposure *100) , "%)")
           ) |>
        rename(very_close_cluster_id = all_of(cluster_col)) |>
    select(very_close_cluster_id, close_cluster_size, 
           main_mlst_prop, main_sampletype_prop, 
           n_region, main_region_prop,
           n_healthcare_exposure, main_healthcare_exposure_prop, 
           everything() ) |>
    arrange(desc(close_cluster_size))
}


# add healthcare exposrue cate varibles
samples_assembled_with_species_hcai_deduplicated <- samples_assembled_with_species_hcai_deduplicated |>
  mutate(healthcare_exposure_cat = case_when(sampletype_cat == "screening" ~ "screening",
                                             TRUE ~ healthcare_exposure))
# generate summaries for all 3 thresholds:
# 95th percentile empirical
close_cluster_level_summary_empirical95_count <-
  summarise_close_clusters(
    data = samples_assembled_with_species_hcai_deduplicated,
    cluster_col = "very_close_cluster_id__snvs_per_time_empirical95_count",
    component_size_col = "component_size__snvs_per_time_empirical95_count",
    membership_col = "close_cluster_member_empirical95_count"
  )
#View(close_cluster_level_summary_empirical95_count)

# 0.4 cgSNVs /day
close_cluster_level_summary_0.4 <-
  summarise_close_clusters(
    data = samples_assembled_with_species_hcai_deduplicated,
    cluster_col = "very_close_cluster_id__snvs_per_time_le0.4",
    component_size_col = "component_size__snvs_per_time_le0.4",
    membership_col = "close_cluster_member_0.4_SNVs_per_day"
  )

# 0.003 cgSNVS/day
close_cluster_level_summary_0.003 <-
  summarise_close_clusters(
    data = samples_assembled_with_species_hcai_deduplicated,
    cluster_col = "very_close_cluster_id__snvs_per_time_le0.003",
    component_size_col = "component_size__snvs_per_time_le0.003",
    membership_col = "close_cluster_member_0.003_SNVs_per_day"
  )

write.csv(close_cluster_level_summary_empirical95_count, "close_cluster_level_summary_empirical95_count.csv", row.names = FALSE)
write.csv(close_cluster_level_summary_0.4, "close_cluster_level_summary_04.csv", row.names = FALSE)
write.csv(close_cluster_level_summary_0.003, "close_cluster_level_summary_003.csv", row.names = FALSE)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * LOGISTIC REGRESSION ####
# of close clusters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# add patient age to metadata 
age_cat <- samples_assembled_with_species_hcai_deduplicated |>
  mutate(meansamplecollectiondate = as.Date(meansamplecollectiondate),
         patientdob = as.Date(patientdob),
         
         # Age at sample collection
         age_at_collection_days = as.numeric(
           meansamplecollectiondate - patientdob
         ),
         
         age_at_collection_years = age_at_collection_days / 365.25,
         
         # Completed years of age
         age_at_collection_years_integer = floor(
           age_at_collection_years
         ),
         
         # Flag implausible ages
         age_invalid = !is.na(age_at_collection_years) &
           (age_at_collection_years < 0 |
              age_at_collection_years > 110),
         
         # Remove implausible ages from analysis
         age_at_collection_years = if_else(
           age_invalid,
           NA_real_,
           age_at_collection_years
         )
  )


summary(age_cat$age_at_collection_days)
summary(age_cat$age_at_collection_years)
summary(age_cat$age_at_collection_years_integer)

# create categories
age_cat <- age_cat |>
  mutate(
    age_category = case_when(
      age_at_collection_years < 5 ~ "0–4",
      age_at_collection_years < 18 ~ "5–17",
      age_at_collection_years < 40 ~ "18–39",
      age_at_collection_years < 60 ~ "40–59",
      age_at_collection_years < 75 ~ "60–74",
      age_at_collection_years < 85 ~ "75–84",
      age_at_collection_years >= 85 ~ "85+",
      TRUE ~ NA_character_
    ),
    
    age_category = factor(
      age_category,
      levels = c(
        "0–4", "5–17", "18–39",
        "40–59", "60–74", "75–84", "85+"
      ),
      ordered = TRUE
    )
  )

# vis
age_cat |>
  filter(!is.na(age_at_collection_years)) |>
  ggplot(aes(x = age_at_collection_years)) +
  geom_histogram(
    binwidth = 1,
    boundary = 0,
    fill = "#0072B2",
    colour = "white",
    alpha = 0.85
  ) +
  scale_x_continuous(
    breaks = seq(0, 110, 10),
    limits = c(0, 110)
  ) +
  labs(
    x = "Age at sample collection (years)",
    y = "Number of samples",
    title = "Distribution of patient ages"
  ) +
  theme_bw(base_size = 12)

# vis ge cat
age_cat |>
  filter(!is.na(age_category)) |>
  dplyr::count(age_category) |>
  mutate(
    proportion = n / sum(n)
  ) |>
  ggplot(aes(x = age_category, y = proportion)) +
  geom_col(
    fill = "#0072B2",
    width = 0.75
  ) +
  geom_text(
    aes(label = scales::percent(proportion, accuracy = 0.1)),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_y_continuous(
    labels = scales::percent,
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(
    x = "Patient age category (years)",
    y = "Percentage of samples",
    title = "Age distribution by category"
  ) +
  theme_bw(base_size = 12)

# summarise ages and check missing
# Overall age summary
age_cat |>
  summarise(
    n_samples = n(),
    n_missing_age = sum(is.na(age_at_collection_years)),
    n_invalid_age = sum(age_invalid, na.rm = TRUE),
    median_age = median(age_at_collection_years, na.rm = TRUE),
    mean_age = mean(age_at_collection_years, na.rm = TRUE),
    sd_age = sd(age_at_collection_years, na.rm = TRUE),
    min_age = min(age_at_collection_years, na.rm = TRUE),
    max_age = max(age_at_collection_years, na.rm = TRUE),
    q1_age = quantile(age_at_collection_years, 0.25, na.rm = TRUE),
    q3_age = quantile(age_at_collection_years, 0.75, na.rm = TRUE)
  )

# Age category counts and proportions
age_cat |>
  dplyr::count(age_category, .drop = FALSE) |>
  mutate(
    proportion = n / sum(n)
  )

# left join onyl metadta
age_cat_selected <- age_cat |>
  select(sequencing_id, meansamplecollectiondate, age_at_collection_days, 
         age_at_collection_years, age_at_collection_years_integer, age_category)

all_metadata_deduplicated_amr_age <- all_metadata_deduplicated_amr |>
  left_join(age_cat_selected, by = c("sequencing_id" = "sequencing_id"))  |>
  mutate(age_category = case_when(is.na(age_category) ~ "DOB unknown", 
                                  TRUE ~ age_category))

sum(is.na(all_metadata_deduplicated_amr_age$age_category)) # 236 NA
sum(is.na(all_metadata_deduplicated_amr_age$patientdob)) # 236 NA
table(all_metadata_deduplicated_amr_age$age_category, all_metadata_deduplicated_amr_age$sampletype_cat, useNA = "ifany")

# build modelling df
predictors <- c(
  "age_category",
  "sampletype_cat",
  "region",
  "healthcare_exposure_cat",
  "amr_num_classes_category",
  "amr_beta_lactam"
)

model_df <- component_membership |>
  select(
    sequencing_id,
    clustering_definition,
    very_close_cluster_id,
    component_size
  ) |>
  left_join(
    all_metadata_deduplicated_amr_age |>
      select(
        sequencing_id,
        species_complex,
        all_of(predictors)
      ) |>
      distinct(),
    by = "sequencing_id"
  ) |>
  mutate(
    close_cluster = if_else(
      !is.na(very_close_cluster_id),
      1L,
      0L
    )
  )

# check
model_df_check <- model_df |>
  dplyr::count(
    clustering_definition,
    species_complex,
    close_cluster
  ) |>
  pivot_wider(
    names_from = close_cluster,
    values_from = n,
    values_fill = 0
  )
#View(model_df_check)

# 2. set biologically meaningful reference categories
model_df <- model_df |>
  mutate(
    age_category = relevel(
      factor(age_category),
      ref = "60–74"),
    
    sampletype_cat = relevel(
      factor(sampletype_cat),
      ref = "blood"),
    
    region = relevel(
      factor(region),
      ref = "South East A"),
    
    healthcare_exposure_cat = relevel(
      factor(healthcare_exposure_cat),
      ref = "community"),
    
    amr_num_classes_category = relevel(
      factor(amr_num_classes_category),
      ref = "1"),
    
    amr_beta_lactam = relevel(
      factor(amr_beta_lactam),
      ref = "non_BL_amr")
  )
#View(model_df) # 1 row per isolate per definition


# 3. Descriptive statistics
descriptive_long <- model_df |>
  pivot_longer(
    cols = all_of(predictors),
    names_to = "metadata_variable",
    values_to = "metadata_category"
  ) |>
  filter(
    !is.na(species_complex),
    !is.na(metadata_category)
  ) |>
  mutate(
    metadata_category = as.character(metadata_category)
  ) |>
  dplyr::count(
    clustering_definition,
    species_complex,
    metadata_variable,
    metadata_category,
    close_cluster,
    name = "n"
  ) |>
  group_by(
    clustering_definition,
    species_complex,
    metadata_variable,
    metadata_category
  ) |>
  mutate(
    category_total = sum(n),
    percentage = 100 * n / category_total,
    
    n_pct = paste0(
      n,
      " (",
      sprintf("%.1f", percentage),
      "%)"
    )
  ) |>
  ungroup()
#View(descriptive_long)

# then make wide
descriptive_wide <- descriptive_long |>
  mutate(
    close_cluster = if_else(
      close_cluster == 1,
      "included",
      "not_included"
    )
  ) |>
  select(
    clustering_definition,
    species_complex,
    metadata_variable,
    metadata_category,
    close_cluster,
    n_pct
  ) |>
  pivot_wider(
    names_from = close_cluster,
    values_from = n_pct,
    values_fill = "0 (0.0%)"
  )
#View(descriptive_wide)

#~~~~~~~~~~~~~~~~~~~~~#
# 4. fit one adjusted model per specues x clustering definition
# fail-safe wrapper thatn can tolerate some species complexes missing outcomes
# record whether glm had convergence problems

fit_species_logistic <- function(
    data,
    clustering_definition_value,
    species_complex_value
) {
  
  d <- data |>
    dplyr::filter(
      clustering_definition ==
        .env$clustering_definition_value,
      
      species_complex ==
        .env$species_complex_value
    ) |>
    dplyr::select(
      close_cluster,
      dplyr::all_of(predictors)
    ) |>
    tidyr::drop_na() |>
    droplevels()
  
  n_events <- sum(d$close_cluster == 1)
  n_non_events <- sum(d$close_cluster == 0)
  
  if (
    nrow(d) == 0 ||
    length(unique(d$close_cluster)) < 2
  ) {
    
    return(
      list(
        model = NULL,
        
        diagnostics = tibble(
          clustering_definition =
            clustering_definition_value,
          
          species_complex =
            species_complex_value,
          
          n_model = nrow(d),
          n_clustered = n_events,
          n_not_clustered = n_non_events,
          n_parameters = NA_integer_,
          converged = FALSE,
          
          model_status =
            "not_estimable_no_outcome_variation"
        )
      )
    )
  }
  
  usable_predictors <- predictors[
    purrr::map_lgl(
      d[predictors],
      ~ dplyr::n_distinct(.x) > 1
    )
  ]
  
  form <- reformulate(
    usable_predictors,
    response = "close_cluster"
  )
  
  fit <- glm(
    form,
    data = d,
    family = binomial()
  )
  
  # Identify extreme fitted probabilities
  extreme_probs <- any(
    fitted(fit) < 1e-6 |
      fitted(fit) > 1 - 1e-6
  )
  
  list(
    model = fit,
    
    diagnostics = tibble(
      clustering_definition =
        clustering_definition_value,
      
      species_complex =
        species_complex_value,
      
      n_model = nrow(d),
      n_clustered = n_events,
      n_not_clustered = n_non_events,
      
      n_parameters =
        length(coef(fit)),
      
      converged =
        fit$converged,
      
      extreme_fitted_probabilities =
        extreme_probs,
      
      model_status = case_when(
        !fit$converged ~ "did_not_converge",
        extreme_probs ~ "possible_separation",
        TRUE ~ "fitted"
      )
    )
  )
}


# create grid
model_grid <- tidyr::crossing(
  clustering_definition_value = unique(model_df$clustering_definition),
  species_complex_value = unique(na.omit(model_df$species_complex))
)

# fit everything 
species_models <- model_grid |>
  mutate(result = map2(
    clustering_definition_value,
    species_complex_value,
    ~ fit_species_logistic(
      model_df,
      clustering_definition_value = .x,
      species_complex_value = .y
    )
  )
  )
#print(species_models)
#View(species_models)

# diagnostics
model_diagnostics <- species_models |>
  dplyr::transmute(result = purrr::map(result, "diagnostics")) |>
  tidyr::unnest(result)

#View(model_diagnostics)

# 5. Adjusted ORs, 95% CIs, Wald p-values
# extract function - using wald confidence intervals
extract_or <- function(
    fit,
    clustering_definition_value,
    species_complex_value
) {
  
  if (is.null(fit)) {
    return(tibble())
  }
  
  x <- broom::tidy(
    fit,
    conf.int = FALSE,
    exponentiate = FALSE
  ) |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::mutate(
      # Wald CI on log-odds scale
      conf.low_log = estimate - 1.96 * std.error,
      conf.high_log = estimate + 1.96 * std.error,
      
      # Transform back to odds-ratio scale
      adjusted_or = exp(estimate),
      ci_low = exp(conf.low_log),
      ci_high = exp(conf.high_log)
    ) |>
    dplyr::transmute(
      clustering_definition =
        clustering_definition_value,
      
      species_complex =
        species_complex_value,
      
      term,
      
      adjusted_or,
      ci_low,
      ci_high,
      
      wald_p = p.value,
      
      log_odds = estimate,
      std_error = std.error,
      
      coefficient_finite =
        is.finite(estimate) &
        is.finite(std.error) &
        is.finite(adjusted_or) &
        is.finite(ci_low) &
        is.finite(ci_high)
    )
  
  x
}

# run
adjusted_or_results <- species_models |>
  dplyr::transmute(
    clustering_definition_value,
    species_complex_value,
    model = purrr::map(
      result,
      ~ .x$model
    )
  ) |>
  purrr::pmap_dfr(
    \(clustering_definition_value,
      species_complex_value,
      model)
    
    extract_or(
      model,
      clustering_definition_value,
      species_complex_value
    )
  )

#View(adjusted_or_results)

# check for extreme or non-finite coefficients:
adjusted_or_results |>
  dplyr::filter(
    !coefficient_finite |
      adjusted_or > 100 |
      ci_high > 1000
  ) |>
  dplyr::arrange(
    species_complex,
    clustering_definition
  )

# check warnings
warnings()


# overall likelihood ration and p-value for each variable
# lrt function helper
get_lrt_results <- function(
    fit,
    clustering_definition_value,
    species_complex_value
) {
  
  if (is.null(fit)) {
    return(tibble())
  }
  
  model_terms <- attr(
    terms(fit),
    "term.labels"
  )
  
  map_dfr(
    model_terms,
    function(variable) {
      
      reduced_formula <- update(
        formula(fit),
        paste(". ~ . -", variable)
      )
      
      reduced_fit <- glm(
        reduced_formula,
        data = model.frame(fit),
        family = binomial()
      )
      
      lr <- anova(
        reduced_fit,
        fit,
        test = "LRT"
      )
      
      tibble(
        clustering_definition =   clustering_definition_value,
        
        species_complex =
          species_complex_value,
        
        metadata_variable = variable,
        
        lrt_chisq =
          lr$Deviance[2],
        
        lrt_df =
          lr$Df[2],
        
        lrt_p =
          lr$`Pr(>Chi)`[2]
      )
    }
  )
}

# call function
lrt_results <- species_models |>
  dplyr::transmute(
    clustering_definition_value,
    species_complex_value,
    
    model = purrr::map(
      result,
      "model"
    )
  ) |>
  pmap_dfr(
    \(clustering_definition_value,
      species_complex_value,
      model)
    
    get_lrt_results(
      model,
      clustering_definition_value,
      species_complex_value
    )
  )

#View(lrt_results)

# 7. Map coefficient sback to variables/ categories
# helper to cretae a lookup
make_term_lookup <- function(fit) {
  
  if (is.null(fit)) {
    return(tibble())
  }
  
  mf <- model.frame(fit)
  
  map_dfr(
    intersect(
      predictors,
      names(mf)
    ),
    function(variable) {
      
      x <- mf[[variable]]
      
      if (!is.factor(x)) {
        return(tibble())
      }
      
      levs <- levels(x)
      
      tibble(
        metadata_variable = variable,
        metadata_category = levs,
        reference = seq_along(levs) == 1
      )
    }
  )
}

# inspect coefficient names: 
adjusted_or_results |>
  distinct(term)

# make final cross tabulation with ORs and p-values
# split 'term#
adjusted_or_results_clean <- adjusted_or_results |>
  mutate(
    metadata_variable = case_when(
      str_starts(term, "age_category") ~ "age_category",
      str_starts(term, "sampletype_cat") ~ "sampletype_cat",
      str_starts(term, "region") ~ "region",
      str_starts(term, "healthcare_exposure_cat") ~ "healthcare_exposure_cat",
      str_starts(term, "amr_num_classes_category") ~ "amr_num_classes_category",
      str_starts(term, "amr_beta_lactam") ~ "amr_beta_lactam",
      TRUE ~ NA_character_
    ),
    
    metadata_category = case_when(
      metadata_variable == "age_category" ~
        str_remove(term, "^age_category"),
      metadata_variable == "sampletype_cat" ~
        str_remove(term, "^sampletype_cat"),
      metadata_variable == "region" ~
        str_remove(term, "^region"),
      metadata_variable == "healthcare_exposure_cat" ~
        str_remove(term, "^healthcare_exposure_cat"),
      metadata_variable == "amr_num_classes_category" ~
        str_remove(term, "^amr_num_classes_category"),
      metadata_variable == "amr_beta_lactam" ~
        str_remove(term, "^amr_beta_lactam"),
      TRUE ~ NA_character_
    ),
    
    OR_95CI = paste0(
      sprintf("%.2f", adjusted_or),
      " (",
      sprintf("%.2f", ci_low),
      "-",
      sprintf("%.2f", ci_high),
      ")"
    )
  )

# join descriptive table
descriptive_model_results <- descriptive_wide |>
  left_join(
    adjusted_or_results_clean |>
      select(
        clustering_definition,
        species_complex,
        metadata_variable,
        metadata_category,
        OR_95CI,
        wald_p
      ),
    by = c(
      "clustering_definition",
      "species_complex",
      "metadata_variable",
      "metadata_category"
    )
  )

# format p-value for reporting
descriptive_model_results <- descriptive_model_results |>
  mutate(
    wald_p = case_when(
      is.na(wald_p) ~ NA_character_,
      wald_p < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", wald_p)
    )
  )

# reorder columns
descriptive_model_results <- descriptive_model_results |>
  select(
    clustering_definition,
    species_complex,
    metadata_variable,
    metadata_category,
    included,
    not_included,
    OR_95CI,
    wald_p
  )
#View(descriptive_model_results)
# set variable orders
species_order <- c("E. coli",
                   "K. pneumoniae complex",
                   "K. oxytoca complex",
                   "K. aerogenes" ,
                   "C. freundii complex" ,
                   "E. cloacae complex",
                   "S. marcescens",
                   "Overall")

metadata_variable_order <- c("age_category",
                             "sampletype_cat",
                             "region",
                             "healthcare_exposure_cat",
                             "amr_num_classes_category",
                             "amr_beta_lactam"
)

age_cat_order <- c("0–4",
                   "5–17",
                   "18–39",
                   "40–59",
                   "60–74",
                   "75–84",
                   "85+",
                   "DOB unknown"
)
unique(descriptive_model_results$metadata_category)

# set metadata category value order
hcai_order <- c("community", "quasi-community", "quasi-nosocomial", "nosocomial", "unknown", "screening")
bl_carriage_order <- c("non_BL_amr", "other_beta_lactam", "cephalosporin", "carbapenem")
amr_class_num_order <- c("1", "2", "3", "4+")
hcai_bl_amr_age <- c(hcai_order, bl_carriage_order, amr_class_num_order, age_cat_order)
metadata_values_no_hcai_bl_age <- setdiff(unique(descriptive_model_results$metadata_category), hcai_bl_amr_age)
metadata_category_order <- c(hcai_bl_amr_age, metadata_values_no_hcai_bl)

descriptive_model_results$species_complex <- factor(descriptive_model_results$species_complex, 
                                                    levels = species_order)
descriptive_model_results$metadata_variable <- factor(descriptive_model_results$metadata_variable,  
                                                      levels = metadata_variable_order)
descriptive_model_results$metadata_category <- factor(descriptive_model_results$metadata_category, 
                                                      levels = metadata_category_order)

write.csv(descriptive_model_results, "close_cluster_LR_descriptive_model_results_3_def.csv", row.names = FALSE)
#View(descriptive_model_results)

# add in reference categories in a complete category universe
category_universe <- model_df |>
  pivot_longer(
    cols = all_of(predictors),
    names_to = "metadata_variable",
    values_to = "metadata_category"
  ) |>
  filter(!is.na(metadata_category)) |>
  mutate(
    metadata_category =
      as.character(metadata_category)
  ) |>
  distinct(
    species_complex,
    metadata_variable,
    metadata_category
  )


#~~~~~~~~~~~~#
# forest plot - create df
forest_df <- adjusted_or_results_clean |>
  left_join(
    lrt_results,
    by = c(
      "clustering_definition",
      "species_complex",
      "metadata_variable"
    )
  ) |>
  mutate(
    category_label = metadata_category,
    
    # Larger point = smaller Wald p
    significance_size = -log10(
      pmax(wald_p, 1e-10)
    ),
    
    # Prevent a tiny number of extreme p-values
    # dominating the point-size scale
    significance_size =
      pmin(significance_size, 10)
  )


# fuction for one clustering definition
make_or_forest_plot <- function(
    data,
    clustering_definition_value
) {
  
  d <- data |>
    dplyr::filter(
      clustering_definition ==
        .env$clustering_definition_value,
      
      species_complex %in% c(
        "E. coli",
        "K. pneumoniae complex"
      ),
      
      is.finite(adjusted_or),
      is.finite(ci_low),
      is.finite(ci_high)
    ) |>
    dplyr::mutate(
      
      metadata_variable = factor(
        metadata_variable,
        levels = c(
          "sampletype_cat",
          "region",
          "healthcare_exposure",
          "amr_num_classes_category",
          "amr_beta_lactam"
        ),
        labels = c(
          "Sample type",
          "Region",
          "Healthcare\nexposure",
          "Number of\nantibiotic class\nresistances",
          "Beta-lactam\nresistance category"
        )
      ),
      
      # Reverse order so the first variable appears at the top
      metadata_variable =
        forcats::fct_rev(metadata_variable),
      
      # Discrete p-value categories for point size
      p_size = dplyr::case_when(
        wald_p < 0.001 ~ "<0.001",
        wald_p < 0.01  ~ "0.001-0.009",
        wald_p < 0.05  ~ "0.01-0.049",
        wald_p < 0.10  ~ "0.05-0.099",
        TRUE           ~ ">=0.10"
      )
    )
  
  ggplot(
    d,
    aes(
      x = adjusted_or,
      y = category_label
    )
  ) +
    
    geom_vline(
      xintercept = 1,
      linetype = "dotted",
      linewidth = 0.4
    ) +
    
    geom_errorbar(
      aes(
        xmin = ci_low,
        xmax = ci_high
      ),
      orientation = "y",
      height = 0.15,
      linewidth = 0.4
    ) +
    
    geom_point(
      aes(
        size = p_size
      )
    ) +
    
    facet_grid(
      metadata_variable ~ species_complex,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    
    scale_x_log10() +
    
    scale_size_manual(
      values = c(
        "<0.001"    = 6,
        "0.001-0.009" = 5,
        "0.01-0.049"  = 4,
        "0.05-0.099"  = 3,
        ">=0.10"      = 2
      ),
      name = "Wald P-value",
      breaks = c(
        "<0.001",
        "0.001-0.009",
        "0.01-0.049",
        "0.05-0.099",
        ">=0.10"
      )
    ) +
    
    labs(
      title = clustering_definition_value,
      x = "Adjusted odds ratio (95% CI)",
      y = NULL
    ) +
    
    theme_minimal() +
    
    theme(
      panel.grid.minor = element_blank(),
      
      strip.text = element_text(
        face = "bold"
      ),
      
      strip.placement = "outside",
      
      # Put facet labels on the left
      strip.text.y.left = element_text(
        angle = 0,
        hjust = 1
      ),
      
      axis.text.y = element_text(
        hjust = 1
      ),
      
      legend.position = "bottom"
    )
}

# plotting helpers

threshold_levels <- c(
  "snvs_per_time_le0.4",
  "snvs_per_time_empirical95_count"#,
 # "snvs_per_time_le0.003"
)

threshold_labels <- c(
  "≤0.4 cgSNVs/day",
  "≤95th percentile cgSNVs"#,
  #"≤0.003 cgSNVs/day"
)

threshold_colours <- c(
  "≤0.4 cgSNVs/day"          = "#66B2A8",
  "≤95th percentile cgSNVs"  = "#005F73"   
 # "≤0.003 cgSNVs/day"           = "#B7DCD3"
)


# make self-contained df
forest_df <- adjusted_or_results |>
  mutate(
    metadata_variable = case_when(
      str_starts(term, "age_category") ~ "age_category",
      str_starts(term, "sampletype_cat") ~ "sampletype_cat",
      str_starts(term, "region") ~ "region",
      str_starts(term, "healthcare_exposure_cat") ~ "healthcare_exposure",
      str_starts(term, "amr_num_classes_category") ~ "amr_num_classes_category",
      str_starts(term, "amr_beta_lactam") ~ "amr_beta_lactam",
      TRUE ~ NA_character_
    ),
    
    metadata_category = case_when(
      metadata_variable == "age_category" ~        str_remove(term, "^age_category"),
      metadata_variable == "sampletype_cat" ~        str_remove(term, "^sampletype_cat"),
      metadata_variable == "region" ~        str_remove(term, "^region"),
      metadata_variable == "healthcare_exposure" ~        str_remove(term, "^healthcare_exposure_cat"),
      metadata_variable == "amr_num_classes_category" ~        str_remove(term, "^amr_num_classes_category"),
      metadata_variable == "amr_beta_lactam" ~        str_remove(term, "^amr_beta_lactam"),
      TRUE ~ NA_character_
    ),
    
    threshold = factor(
      clustering_definition,
      levels = threshold_levels,
      labels = threshold_labels
    )
  ) |>
  filter(
    species_complex %in% c(
      "E. coli",
      "K. pneumoniae complex"
    ),
    is.finite(adjusted_or),
    is.finite(ci_low),
    is.finite(ci_high),
    !is.na(metadata_variable)
  ) |>
  mutate(
    metadata_variable = factor(
      metadata_variable,
      levels = c(
        "age_category",
        "sampletype_cat",
        "region",
        "healthcare_exposure",
        "amr_num_classes_category",
        "amr_beta_lactam"
      ),
      labels = c(
        "Age\ncategory",
        "Sample\ntype",
        "Region",
        "Healthcare\nExposure",
        "Number of\nantibiotic class\nresistances",
        "Beta-lactam\nresistance\ncategory"
      )
    ),
    
    # metadata_variable = forcats::fct_rev(
    #  metadata_variable
    #),
    
    p_size = case_when(
      wald_p < 0.001 ~ "<0.001",
      wald_p < 0.01  ~ "0.001-0.009",
      wald_p < 0.05  ~ "0.01-0.049",
      wald_p < 0.10  ~ "0.05-0.099",
      TRUE           ~ ">=0.10"
    )
  )

# set order of categories
forest_df$metadata_category <- factor(forest_df$metadata_category, levels = metadata_category_order)

forest_df_filt <- forest_df |>
  filter(clustering_definition %in% c("snvs_per_time_empirical95_count", "snvs_per_time_le0.4"),
         !(species_complex == "K. pneumoniae complex" & metadata_variable == "Healthcare\nExposure"))

# plot
p_all_thresholds <- ggplot(
  forest_df_filt, # filter out 95th percentile rate threshold as this produce huge implausible clusters
  aes(
    x = adjusted_or,
    y = metadata_category,
    colour = threshold,
    size = p_size )) +
  geom_vline(
    xintercept = 1,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_errorbar(aes(xmin = ci_low,  xmax = ci_high),
    orientation = "y",
    height = 0.15,
    linewidth = 0.45,
    position = position_dodge(
      width = 0.6
    )
  ) +
  geom_point(
    position = position_dodge(
      width = 0.6
    )
  ) +
  facet_grid(
    metadata_variable ~ species_complex,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_x_log10() +
  scale_colour_manual(
    values = threshold_colours,
    drop = FALSE
  ) +
  scale_size_manual(values = c(
      "<0.001" = 6,
      "0.001-0.009" = 5,
      "0.01-0.049" = 4,
      "0.05-0.099" = 3,
      ">=0.10" = 2
    ),
    name = "Wald P-value"
  ) +
  labs(
    x = "Adjusted odds ratio (95% CI)",
    y = NULL,
    colour = "Close-cluster definition"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    
    strip.placement = "outside",
    
    strip.text.y.left = element_text(
      angle = 0,
      hjust = 1,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      hjust = 1
    ),
    
    legend.position = "bottom"
  )

p_all_thresholds
# save
ggsave("all_thresholds_forest_plot_close_clusters_LR_2_definitions.png", p_all_thresholds, units = "in", width = 14, height = 10, dpi = 300)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plasmids ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
all_plasmids_distances_meta <- read.csv("all_plasmids_distances_meta_longest_contig_plasmid_definition.csv")
#View(all_plasmids_distances_meta)
plasmid_subcommunities_meta <- read.csv("plasmid_subcommunities_meta.csv")
#View(plasmid_subcommunities_meta)
#nrow(plasmid_subcommunities_meta)
#colnames(plasmid_subcommunities_meta)


# split plasmid subcommunities with bimodal length distributions
# Function to split one subcommunity by plasmid length
split_subcommunity_by_length <- function(
    df,
    min_n_mixture = 5,
    bic_improvement = 5,
    min_fold_difference = 1.5,
    singleton_fold_difference = 5
) {
  
  lengths <- as.numeric(df$autocycler_contig_length)
  valid <- which(is.finite(lengths) & lengths > 0)
  
  # Default: do not split
  df$pling_subcommunity_length_split <- as.character(df$pling_subcommunity)
  df$length_group <- NA_character_
  df$length_split_method <- "none"
  if (length(valid) < 2) return(df)
    log_lengths <- log10(lengths[valid])
  n <- length(log_lengths)
  
  # 1. Check for extreme singleton outliers
  sorted_lengths <- sort(lengths[valid])
    # Smallest plasmid is much shorter than all others
  small_outlier <- (
    sorted_lengths[2] / sorted_lengths[1]
  ) >= singleton_fold_difference
  
  # Largest plasmid is much longer than all others
  large_outlier <- (
    sorted_lengths[n] / sorted_lengths[n - 1]
  ) >= singleton_fold_difference
  
  singleton_groups <- NULL
  
  if (small_outlier || large_outlier) {
    
    # Choose the stronger separation if both are present
    small_ratio <- sorted_lengths[2] / sorted_lengths[1]
    large_ratio <- sorted_lengths[n] / sorted_lengths[n - 1]
    
    if (small_outlier && small_ratio >= large_ratio) {
      singleton_groups <- ifelse(
        lengths[valid] == sorted_lengths[1], 2L, 1L
      )
    } else {
      singleton_groups <- ifelse(
        lengths[valid] == sorted_lengths[n], 1L, 2L
      )
    }
  }
  
  # 2. Gaussian mixture model
  mixture_groups <- NULL
  
  if (n >= min_n_mixture &&
      length(unique(log_lengths)) >= 3) {
    
    fit <- tryCatch(
      Mclust(
        log_lengths,
        G = 1:2,
        modelNames = "E",
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    
    if (!is.null(fit) && fit$G == 2) {
      
      # Compare BIC against a one-component model
      fit1 <- tryCatch(
        Mclust(
          log_lengths,
          G = 1,
          modelNames = "E",
          verbose = FALSE
        ),
        error = function(e) NULL
      )
      
      if (!is.null(fit1)) {
        
        delta_bic <- fit$bic - fit1$bic
        
        group_medians <- tapply(
          lengths[valid],
          fit$classification,
          median
        )
        
        fold_difference <- max(group_medians) /
          min(group_medians)
        
        if (
          delta_bic >= bic_improvement &&
          fold_difference >= min_fold_difference
        ) {
          mixture_groups <- fit$classification
        }
      }
    }
  }
  
  # 3. Choose splitting method
  if (!is.null(mixture_groups)) {
    
    groups <- mixture_groups
    method <- "mixture"
    
  } else if (!is.null(singleton_groups)) {
    
    groups <- singleton_groups
    method <- "singleton_outlier"
    
  } else {
    return(df)
  }
  
  # 4. Name groups by median length
  group_medians <- tapply(
    lengths[valid],
    groups,
    median
  )
  
  group_order <- names(sort(group_medians, decreasing = TRUE))
  
  group_labels <- setNames(
    c("a", "b"),
    group_order
  )
  
  labels <- unname(group_labels[as.character(groups)])
  
  df$length_group[valid] <- labels
  
  df$pling_subcommunity_length_split[valid] <- paste0(
    as.character(df$pling_subcommunity[valid]),
    "_",
    labels
  )
  
  df$length_split_method[valid] <- method
  
  df
}

# apply to df
plasmid_subcommunities_meta_split <- plasmid_subcommunities_meta |>
  group_split(pling_subcommunity, .keep = TRUE) |>
  map_dfr(split_subcommunity_by_length)

# summarise splittogn results
length_split_summary <-  plasmid_subcommunities_meta_split |>
  group_by(
    pling_subcommunity,
    pling_subcommunity_length_split,
    length_split_method
  ) |>
  summarise(
    n_plasmids = n(),
    median_length = median(autocycler_contig_length, na.rm = TRUE),
    min_length = min(autocycler_contig_length, na.rm = TRUE),
    max_length = max(autocycler_contig_length, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(pling_subcommunity, desc(median_length))

# Visualise length distributions
pling_subcommunity_split_plot <- plasmid_subcommunities_meta_split |>
  filter(length_split_method != "none") |>
  ggplot(
    aes(
      x = autocycler_contig_length,
      fill = length_group
    )
  ) +
  geom_histogram(
    bins = 25,
    alpha = 0.8,
    colour = NA,
    linewidth = 0.2,
    position = "identity"
  ) +
  scale_fill_manual(
    values = c(
      "a" = "#0072B2",  # Blue - larger plasmids
      "b" = "#E69F00"   # Orange - smaller plasmids
    ),
    labels = c(
      "a" = "a (larger)",
      "b" = "b (smaller)"
    )
  ) +
  scale_x_log10(
    labels = scales::label_number(big.mark = ",")
  ) +
  facet_wrap(~ pling_subcommunity, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "Plasmid length (bp, log10 scale)",
    y = "Number of plasmids",
    fill = "Length group"
  )

pling_subcommunity_split_plot
ggsave("pling_subcommunity_split_plot_bic_improvement_2.png", pling_subcommunity_split_plot, width = 10, height = 8, dpi = 300)


# save split df
colnames(plasmid_subcommunities_meta_split)
table(plasmid_subcommunities_meta_split$length_split_method)
table(plasmid_subcommunities_meta_split$pling_subcommunity_length_split)

write.csv(plasmid_subcommunities_meta_split, "plasmid_subcommunities_meta_split_bic_improvement_2.csv", row.names = FALSE)


# summarise plasmid split cluster size
length(unique(plasmid_subcommunities_meta_split$pling_subcommunity_length_split)) # 1567
length(unique(plasmid_subcommunities_meta_split$pling_subcommunity_length_split)) # 1567
split_cluster_summary <- plasmid_subcommunities_meta_split |>
  group_by(pling_subcommunity_length_split) |>
  summarise(count = n())
View(split_cluster_summary)
table(split_cluster_summary$count)
sum(split_cluster_summary$count) # 5596
sum(split_cluster_summary[split_cluster_summary$count >1, ]$count) # 4404
summary(split_cluster_summary[split_cluster_summary$count >1, ]$count) # 4404
length(split_cluster_summary[split_cluster_summary$count >1, ]$count) # 375


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# summarise plasmids per species complex and # plasmids per isolate
table(plasmid_subcommunities_meta$molecule_type_genomad_majority, useNA = "ifany")
table(plasmid_subcommunities_meta$sampletype_cat, useNA = "ifany")
table(plasmid_subcommunities_meta$species_complex, useNA = "ifany")
table(plasmid_subcommunities_meta$species_complex, plasmid_subcommunities_meta$molecule_type_genomad_majority, useNA = "ifany")
table(plasmid_subcommunities_meta$species_complex, plasmid_subcommunities_meta$sampletype_cat, useNA = "ifany")

# overall plasmid/virus per isolate
num_plasmids_per_isolate <- plasmid_subcommunities_meta_split |>
  filter(!grepl("run", unique_sample)) |>
  group_by(unique_sample, species_complex) |>
  summarise(num_non_chrom_in_isolate = n(),
            num_plasmids_in_isolate = sum(molecule_type_genomad_majority == "plasmid"),
            num_virus_in_isolate = sum(molecule_type_genomad_majority == "virus"))

num_plasmids_per_isolate_with_zero <- samples_assembled_with_species_hcai_deduplicated |>
  mutate(species_complex = case_when(grepl("ecoli_achtman_4", typing_scheme) ~ "E. coli",
                                     grepl("klebsiella", typing_scheme) ~ "K. pneumoniae complex",
                                     grepl("koxytoca", typing_scheme) ~ "K. oxytoca complex",
                                     grepl("kaerogenes", typing_scheme) ~ "K. aerogenes",
                                     grepl("cfreundii", typing_scheme) ~ "C. freundii complex",
                                     grepl("ecloacae", typing_scheme) ~ "E. cloacae complex",
                                     sequencing_id == "AF31" ~ "S. marcescens",
                                     TRUE ~ NA)) |>
  filter(assembly_qc_pass == TRUE,
         duplicate_assembly_qc_pass != "duplicate") |>
  select(sequencing_id, species_complex, sampletype_cat) |>
  left_join(num_plasmids_per_isolate |> select(-c(species_complex)), by = c("sequencing_id" = "unique_sample")) |>
  mutate(num_non_chrom_in_isolate = replace_na(num_non_chrom_in_isolate, 0),
         num_plasmids_in_isolate  = replace_na(num_plasmids_in_isolate, 0),
         num_virus_in_isolate = replace_na(num_virus_in_isolate, 0),
         genus = case_when(grepl("K.", species_complex) ~ "Klebsiella",
                           TRUE ~ species_complex))

colnames(num_plasmids_per_isolate_with_zero)
nrow(num_plasmids_per_isolate_with_zero)
sum(is.na(num_plasmids_per_isolate_with_zero$species_complex))

summary(num_plasmids_per_isolate_with_zero$num_non_chrom_in_isolate)
summary(num_plasmids_per_isolate_with_zero$num_plasmids_in_isolate)
summary(num_plasmids_per_isolate_with_zero$num_virus_in_isolate)

# per species complex summary
species_complex_plasmids_num_summary <- num_plasmids_per_isolate_with_zero |>
  group_by(species_complex) |>
  summarise(min_non_chrom = min(num_non_chrom_in_isolate),
            q1_non_chrom = quantile(num_non_chrom_in_isolate, 0.25),
            q2_non_chrom = quantile(num_non_chrom_in_isolate, 0.5),
            q3_non_chrom = quantile(num_non_chrom_in_isolate, 0.75),
            max_non_chrom = max(num_non_chrom_in_isolate),
            
            min_plasmid = min(num_plasmids_in_isolate),
            q1_plasmid = quantile(num_plasmids_in_isolate, 0.25),
            q2_plasmid = quantile(num_plasmids_in_isolate, 0.5),
            q3_plasmid = quantile(num_plasmids_in_isolate, 0.75),
            max_plasmid = max(num_plasmids_in_isolate),
            
            min_virus = min(num_virus_in_isolate),
            q1_virus = quantile(num_virus_in_isolate, 0.25),
            q2_virus = quantile(num_virus_in_isolate, 0.5),
            q3_virus = quantile(num_virus_in_isolate, 0.75),
            max_virus = max(num_virus_in_isolate)
  )
#View(species_complex_plasmids_num_summary)
# pairwise wilcox test
pairwise.wilcox.test(num_plasmids_per_isolate_with_zero$num_non_chrom_in_isolate, 
                     num_plasmids_per_isolate_with_zero$species_complex,
                     p.adjust.method =  "bonf",
                     paired = FALSE
)
pairwise.wilcox.test(num_plasmids_per_isolate_with_zero$num_non_chrom_in_isolate, 
                     num_plasmids_per_isolate_with_zero$genus,
                     p.adjust.method =  "bonf",
                     paired = FALSE
)

#sampletype summary
sampletype_plasmids_num_summary <- num_plasmids_per_isolate_with_zero |>
  group_by(sampletype_cat) |>
  summarise(min_non_chrom = min(num_non_chrom_in_isolate),
            q1_non_chrom = quantile(num_non_chrom_in_isolate, 0.25),
            q2_non_chrom = quantile(num_non_chrom_in_isolate, 0.5),
            q3_non_chrom = quantile(num_non_chrom_in_isolate, 0.75),
            max_non_chrom = max(num_non_chrom_in_isolate),
            
            min_plasmid = min(num_plasmids_in_isolate),
            q1_plasmid = quantile(num_plasmids_in_isolate, 0.25),
            q2_plasmid = quantile(num_plasmids_in_isolate, 0.5),
            q3_plasmid = quantile(num_plasmids_in_isolate, 0.75),
            max_plasmid = max(num_plasmids_in_isolate),
            
            min_virus = min(num_virus_in_isolate),
            q1_virus = quantile(num_virus_in_isolate, 0.25),
            q2_virus = quantile(num_virus_in_isolate, 0.5),
            q3_virus = quantile(num_virus_in_isolate, 0.75),
            max_virus = max(num_virus_in_isolate)
  )
#View(sampletype_plasmids_num_summary)

pairwise.wilcox.test(num_plasmids_per_isolate_with_zero$num_non_chrom_in_isolate, 
                     num_plasmids_per_isolate_with_zero$sampletype_cat,
                     p.adjust.method =  "bonf",
                     paired = FALSE
)

#per sampletype summary

#View(num_plasmids_per_isolate_with_zero)
# overall plasmid/ virus per isolate stratified by species complex
table(autocycler_contigs_mobsuite_genomad$molecule_type_genomad_majority)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * summarise split clusters ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~#
# subcommunity (split) level summary
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# create subcommunity-level summary

# Helper functions
# Modal value and its proportion
mode_prop <- function(x) {
  x <- x[!is.na(x) & x != ""]
  
  if (length(x) == 0) {
    return(
      list(
        main = NA_character_,
        prop = NA_real_
      )
    )
  }
  
  tab <- table(x)
  main <- names(tab)[which.max(tab)]
  
  list(
    main = main,
    prop = unname(max(tab) / length(x))
  )
}


# Number of unique comma-separated components across a cluster
# and modal component + proportion of plasmids containing it
comma_feature_summary <- function(x, n_plasmids) {
  
  x <- as.character(x)
  x <- x[!is.na(x) & x != "" & x != "-"]
  
  if (length(x) == 0) {
    return(list(
      n_unique = 0L,
      main = NA_character_,
      prop = NA_real_
    ))
  }
  
  # Split comma-separated values
  split_values <- stringr::str_split(x, ",")
  
  # Trim whitespace and remove empty entries
  split_values <- purrr::map(
    split_values,
    ~ {
      vals <- stringr::str_trim(.x)
      vals[vals != "" & vals != "-"]
    }
  )
  
  all_values <- unlist(split_values)
  
  if (length(all_values) == 0) {
    return(list(
      n_unique = 0L,
      main = NA_character_,
      prop = NA_real_
    ))
  }
  
  # Count plasmids containing each replicon type
  # (each replicon counted once per plasmid)
  all_values_unique <- unlist(
    purrr::map(split_values, unique)
  )
  
  tab <- table(all_values_unique)
  
  main <- names(tab)[which.max(tab)]
  
  n_with_main <- sum(
    purrr::map_lgl(
      split_values,
      ~ main %in% .x
    )
  )
  
  list(
    n_unique = dplyr::n_distinct(all_values),
    main = main,
    prop = n_with_main / n_plasmids
  )
}

# Number of comma-separated replicons per plasmid
n_replicons <- function(x) {
  
  ifelse(
    is.na(x) | x == "" | x == "-",
    0L,
    str_count(x, ",") + 1L
  )
}


# Descriptive statistics
num_stats <- function(x, prefix) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) {
    return(
      tibble(
        !!paste0(prefix, "_median") := NA_real_,
        !!paste0(prefix, "_mean") := NA_real_,
        !!paste0(prefix, "_min") := NA_real_,
        !!paste0(prefix, "_max") := NA_real_,
        !!paste0(prefix, "_q1") := NA_real_,
        !!paste0(prefix, "_q3") := NA_real_
      )
    )
  }
  
  tibble(
    !!paste0(prefix, "_median") := median(x),
    !!paste0(prefix, "_mean") := mean(x),
    !!paste0(prefix, "_min") := min(x),
    !!paste0(prefix, "_max") := max(x),
    !!paste0(prefix, "_q1") := quantile(x, 0.25, names = FALSE),
    !!paste0(prefix, "_q3") := quantile(x, 0.75, names = FALSE)
  )
}


subcommunity_summary <- plasmid_subcommunities_meta_split |>
  group_by(pling_subcommunity_length_split) |>
  group_modify(~ {
    
    d <- .x
    n_plasmids <- nrow(d)
    
    # Simple categorical variables
    species_complex <- mode_prop(d$kraken2_species)
    
    mlst <- mode_prop(d$kleborate_mlst_ST)
    poppunk <- mode_prop(d$poppunk_cluster)
    region <- mode_prop(d$region)
    sampletype <- mode_prop(d$sampletype_cat)
    hcai <- mode_prop(d$healthcare_exposure)
    genomad <- mode_prop(d$molecule_type_genomad_majority)
    fastbaps <- mode_prop(d$fastbaps_L3)
    
    # COMMA-SEPARATED FEATURE SUMMARY
    # combo_rep_types
    rep_summary <- comma_feature_summary(
      d$combo_rep_types,
      n_plasmids
    )
    
    # Number of replicons per plasmid
    replicon_counts <- n_replicons(
      d$combo_rep_types
    )
    
    replicon_stats <- num_stats(
      replicon_counts,
      "n_replicons"
    )
    
    # Exact-value categorical summaries
    relaxase <- mode_prop(d$combo_relaxase_types)
    mpf <- mode_prop(d$combo_mpf_types)
    orit <- mode_prop(d$combo_orit_types)
    mobility <- mode_prop(d$combo_predicted_mobility)
    primary_cluster <- mode_prop(d$combo_primary_cluster_id)
    secondary_cluster <- mode_prop(d$combo_secondary_cluster_id)
    
    # Date summary
    dates <- as.Date(d$meansamplecollectiondate)
    
    date_min <- if (all(is.na(dates))) {
      as.Date(NA)
    } else {
      min(dates, na.rm = TRUE)
    }
    
    date_max <- if (all(is.na(dates))) {
      as.Date(NA)
    } else {
      max(dates, na.rm = TRUE)
    }
    
    date_range_days <- if (
      is.na(date_min) |
      is.na(date_max)
    ) {
      NA_real_
    } else {
      as.numeric(date_max - date_min)
    }
    
    # Contig length
    length_stats <- num_stats(
      d$autocycler_contig_length,
      "contig_length"
    )
    
    # Count variables
    count_vars <- c(
      "n_amrfidner_hits",
      "n_amr_genes",
      "n_virulence_genes",
      "n_stress_genes",
      "n_beta_lactam",
      "n_ceph",
      "n_carb",
      "n_aminoglycoside",
      "n_gent",
      "n_quin",
      "n_phenicol",
      "n_tet",
      "n_trimethoprim",
      "n_sulphonamide"
    )
    
    count_summary <- purrr::map(
      count_vars,
      function(v) {
        
        x <- d[[v]]
        
        stats <- num_stats(x, v)
        
        tibble(
          !!paste0(v, "_n_plasmids_any") :=
            sum(!is.na(x) & x > 0)
        ) |>
          dplyr::bind_cols(stats)
      }
    ) |>
      dplyr::bind_cols()
    
    # Build one-row summary
    tibble(
      pling_cluster_size = n_plasmids,
      # Contig length
      length_stats,
      # combo_rep_types - exact values
      n_unique_rep_type_values =
        n_distinct(
          d$combo_rep_types[
            !is.na(d$combo_rep_types) &
              d$combo_rep_types != ""
          ]
        ),
      main_rep_type_value = mode_prop(d$combo_rep_types)$main,
      main_rep_type_value_prop = mode_prop(d$combo_rep_types)$prop,
      
      # combo_rep_types - split comma-separated contents
      n_unique_replicon_types = rep_summary$n_unique,
      main_replicon_type = rep_summary$main,
      main_replicon_type_prop = rep_summary$prop,
      
      # Number replicons per plasmid
      replicon_stats,
      # Other categorical features
      n_unique_relaxase_types =
        n_distinct(
          d$combo_relaxase_types[
            !is.na(d$combo_relaxase_types) &
              d$combo_relaxase_types != ""
          ]
        ),
      
      main_relaxase_type = relaxase$main,
      main_relaxase_type_prop = relaxase$prop,
      
      n_unique_mpf_types = n_distinct(
        d$combo_mpf_types[
          !is.na(d$combo_mpf_types) &
            d$combo_mpf_types != ""
        ]
      ),
      
      main_mpf_type = mpf$main,
      main_mpf_type_prop = mpf$prop,
      n_unique_orit_types =  n_distinct(
        d$combo_orit_types[
          !is.na(d$combo_orit_types) &
            d$combo_orit_types != ""
        ]
      ),
      
      main_orit_type =  orit$main,
      main_orit_type_prop =   orit$prop,
      n_unique_predicted_mobility = n_distinct(
        d$combo_predicted_mobility[
          !is.na(d$combo_predicted_mobility) &
            d$combo_predicted_mobility != ""
        ]
      ),
      
      main_predicted_mobility = mobility$main,
      main_predicted_mobility_prop = mobility$prop,
      n_unique_primary_clusters = n_distinct(
        d$combo_primary_cluster_id[
          !is.na(d$combo_primary_cluster_id)
        ]
      ),
      
      main_primary_cluster = primary_cluster$main,
      main_primary_cluster_prop = primary_cluster$prop,
      
      n_unique_secondary_clusters =
        n_distinct(
          d$combo_secondary_cluster_id[
            !is.na(d$combo_secondary_cluster_id)
          ]
        ),
      
      main_secondary_cluster = secondary_cluster$main,
      main_secondary_cluster_prop = secondary_cluster$prop,
      # Species
      num_species_complexes = n_distinct(d$kraken2_species, na.rm = TRUE),
      main_species_complex = species_complex$main,
      main_species_complex_prop =species_complex$prop,
      
      # MLST / PopPUNK / region / sample type / HCAI
      num_mlst =  n_distinct(d$kleborate_mlst_ST, na.rm = TRUE),
      main_mlst = mlst$main,
      main_mlst_prop =  mlst$prop,
      num_poppunk_clusters = n_distinct(d$poppunk_cluster, na.rm = TRUE),
      main_poppunk_cluster = poppunk$main,
      main_poppunk_cluster_prop = poppunk$prop,
      num_regions = n_distinct(d$region, na.rm = TRUE),
      main_region = region$main,
      main_region_prop = region$prop,
      num_sample_types =  n_distinct(d$sampletype_cat, na.rm = TRUE),
      main_sample_type = sampletype$main,
      main_sample_type_prop = sampletype$prop,
      num_healthcare_exposure_types = n_distinct(d$healthcare_exposure, na.rm = TRUE),
      main_healthcare_exposure = hcai$main,
      main_healthcare_exposure_prop =  hcai$prop,
      # Genomad
      num_genomad_molecule_types = n_distinct(d$molecule_type_genomad_majority, na.rm = TRUE),
      main_genomad_molecule_type = genomad$main,
      main_genomad_molecule_type_prop = genomad$prop,
      # FastBAPS
      num_fastbaps_L3 = n_distinct(d$fastbaps_L3, na.rm = TRUE),
      main_fastbaps_L3 =fastbaps$main,
      main_fastbaps_L3_prop = fastbaps$prop,
      # Dates
      meansamplecollectiondate_min =  date_min,
      meansamplecollectiondate_max = date_max,
      meansamplecollectiondate_range_days = date_range_days
    ) |>
      bind_cols(count_summary)
  }) |>
  ungroup()

# add length range and length IQR
subcommunity_summary <- subcommunity_summary |>
  mutate(length_range = contig_length_max - contig_length_min,
         length_IQR = contig_length_q3 - contig_length_q1)
View(subcommunity_summary)

# save
write.csv(subcommunity_summary, "pling_length_split_subcommunity_summary_longest_contig.csv", row.names = FALSE)
#~~~~~~~~~~~~~#
# read in
#subcommunity_summary <- read.csv("pling_subcommunity_summary_longest_contig.csv")
#View(subcommunity_summary)
#~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Nearest neighbour plasmid plots ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#nearest_neighbours_ties <- read.csv("nearest_neighbours_with_ties_longest_contig_definition.csv")

#~~~~~~~~~~~~~~~~~~~~~~~#
# plots
#~~~~~~~~~~~~~~~~~~~~~~~#
# plot cluster size
pling_cluster_size_hist <- ggplot(subcommunity_summary, aes(x = pling_cluster_size)) +
  geom_histogram(alpha = 0.6, binwidth = 0.05) +
  scale_x_log10() +
  labs(
    x = "Pling plasmid subcommunity size",
    y = "Count") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_cluster_size_hist
ggsave("pling_cluster_size_hist.png", pling_cluster_size_hist, units = "in", width = 6, height = 5, dpi = 300)

# summarise clusters by various metadata
# plot network of non-singleton clusters

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot distance pairs by metadata
# define colours
region_pair_cols = c("TRUE"  = "#DDCC77", 
                     "FALSE" = "#332288")

genomad_molecule_type_colours <- c("virus"      = "#F0E442",
                                   "plasmid"    = "#E69F00",
                                   "chromosome" = "#0072B2")


sampletype_cols <- c(
  "blood|blood" = "darkred",
  "screening|screening" = "goldenrod1",
  "blood|screening" = "#7B3294",
  "Unknown" = "grey70"
)

sampletype_individual_cols <- c(
  "blood"     = "darkred",
  "screening" = "goldenrod1",
  "Unknown"   = "grey70"
)

colours_hcai <- c(
  "community"        = "#FFF3B0",
  "quasi-community"  = "#FFD37F",
  "quasi-nosocomial" = "#FCA17D",
  "nosocomial"       = "#C75197",
  "unknown"          = "grey"
)

colours_region <- c(
  "London"       = "#E41A1C",
  "South East A" = "#FB8072",
  "South East B" = "#FCAE91",
  "South East C" = "#FEE0D2",
  "South West"   = "#FD8D3C",
  "Midlands"     = "#B3DE69",
  "East"         = "#66C2A5",
  "North West"   = "#377EB8",
  "North East A" = "#807DBA",
  "North East B" = "#B3B3E6",
  "unknown"      = "grey"
)

species_colours <- c(
  "E. coli"                = "#009E73",
  "K. pneumoniae complex"  = "#D55E00",
  "K. oxytoca complex"     = "#CC79A7",
  "K. aerogenes"           = "#7F3C8D",
  "C. freundii complex"    = "#E69F00",
  "E. cloacae complex"     = "#0072B2",
  "S. marcescens"          = "#8B1A1A"
)

# load subcommunity summary

# define important subcommunities
# filter out important subcommunities
important_subcommunities <- subcommunity_summary |>
  filter(pling_cluster_size > 1,
         !(n_amrfidner_hits_n_plasmids_any == 0 &  contig_length_median < 20000)) 
nrow(subcommunity_summary)
nrow(important_subcommunities) # 282

# make plotting df
within_subcomm_pling_distances <- all_plasmids_distances_meta |>
  filter(pling_subcommunity_1 == pling_subcommunity_2) |>
  filter(pling_subcommunity_1 %in% important_subcommunities$pling_subcommunity)
nrow(within_subcomm_pling_distances)



# plot 0 distance pairs by region, bacterial spp, ST, HCAI, sampletype
pling_dist_hist <- ggplot(within_subcomm_pling_distances, aes(x = distance, fill = same_region)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_x_log10() +
  scale_fill_manual(values = region_pair_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Same region") +
  #scale_y_log10() +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_dist_by_region_within_subcomm.png", pling_dist_hist, units = "in", width = 6, height = 5, dpi = 300)


# sample type
pling_dist_hist <- ggplot(within_subcomm_pling_distances, aes(x = distance, fill = sampletype_pair)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_x_log10() +
  scale_fill_manual(values = sampletype_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Sample type") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_dist_by_sampletype_within_subcomm.png", pling_dist_hist, units = "in", width = 6, height = 5, dpi = 300)



# same species complex
pling_dist_hist <- ggplot(within_subcomm_pling_distances, aes(x = distance, fill = same_species_complex)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_x_log10() +
  scale_fill_manual(values = region_pair_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Same species complex") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_dist_by_same_species_complex_within_subcomm.png", pling_dist_hist, units = "in", width = 6, height = 5, dpi = 300)

# same species
pling_dist_hist <- ggplot(within_subcomm_pling_distances, aes(x = distance, fill = same_species)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_x_log10() +
  scale_fill_manual(values = region_pair_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Same species") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_dist_by_same_species_within_subcomm.png", pling_dist_hist, units = "in", width = 6, height = 5, dpi = 300)


# same MLST
pling_dist_hist <- ggplot(within_subcomm_pling_distances, aes(x = distance, fill = same_mlst)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_x_log10() +
  scale_fill_manual(values = region_pair_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Same MLST") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_dist_by_same_mlst_within_subcomm.png", pling_dist_hist, units = "in", width = 6, height = 5, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~#
# plot nearest neighbour by metadata
# start nn without ties

# colour by region
pling_dist_hist <- ggplot(nearest_neighbours, aes(x = distance, fill = region_2)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = colours_region) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Nearest neighbour\nregion") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_no_ties_by_region.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

# colour by region - with ties
pling_dist_hist <- ggplot(nearest_neighbours_ties, aes(x = distance, fill = region_2)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = colours_region) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Nearest neighbour\nregion") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_with_ties_by_region.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

#~~~~~~~~~~~~~~~~~~~#
# colour by same region
pling_dist_hist <- ggplot(nearest_neighbours, aes(x = distance, fill = same_region)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = region_pair_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Nearest neighbour\nsame region") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_no_ties_by_same_region.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

# colour by same region - with ties
pling_dist_hist <- ggplot(nearest_neighbours_ties, aes(x = distance, fill = same_region)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = region_pair_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Nearest neighbour\nsame region") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_with_ties_by_same_region.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~#
# sample type
pling_dist_hist <- ggplot(nearest_neighbours, aes(x = distance, fill = sampletype_pair)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = sampletype_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Sample type of plasmid pair") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_no_ties_by_sampletype.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

# colour by region - with ties
pling_dist_hist <- ggplot(nearest_neighbours_ties, aes(x = distance, fill = sampletype_pair)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = sampletype_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Sample type of plasmid pair") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_with_ties_by_sampletype.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~#
# species of pair
pling_dist_hist <- ggplot(nearest_neighbours, aes(x = distance, fill = same_species)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = region_pair_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Same species") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_no_ties_by_same_species.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

# colour by region - with ties
pling_dist_hist <- ggplot(nearest_neighbours_ties, aes(x = distance, fill = same_species)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = region_pair_cols) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Sample species") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_with_ties_by_same_species.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

#~~~~~~~~~~~~~~~~~#
# species- value of nn
pling_dist_hist <- ggplot(nearest_neighbours, aes(x = distance, fill = species_complex_2)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = species_colours) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Nearest neighbour\nspecies") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_no_ties_by_species_complex.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

# colour by region - with ties
pling_dist_hist <- ggplot(nearest_neighbours_ties, aes(x = distance, fill = species_complex_2)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = species_colours) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Nearest neighbour\nspecies") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_with_ties_by_species_complex.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

# HCAI
pling_dist_hist <- ggplot(nearest_neighbours, aes(x = distance, fill = healthcare_exposure_2)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  #scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = colours_hcai) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Nearest neighbour\nhealthcare exposure") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_no_ties_by_hcai.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

# colour by region - with ties
pling_dist_hist <- ggplot(nearest_neighbours_ties, aes(x = distance, fill = healthcare_exposure_2)) +
  geom_histogram(alpha = 0.6, binwidth = 1, colour = NA) +
  scale_y_log10() +
  geom_vline(xintercept = 4.5, linetype = "dashed", alpha = 0.6, linewidth = 0.8 ) +
  scale_fill_manual(values = colours_hcai) +
  labs(
    x = "Pling DCJ/indel distance",
    y = "Count",
    fill = "Nearest neighbour\nhealthcare exposure") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
pling_dist_hist
ggsave("pling_dcj_nearest_neighbour_with_ties_by_healthcare_exposure.png", pling_dist_hist, units = "in", width = 8, height = 5, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Metadata associations plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 1. colours defined above
colours_hcai <- c(
  "community"        = "#FFF3B0",
  "quasi-community"  = "#FFD37F",
  "quasi-nosocomial" = "#FCA17D",
  "nosocomial"       = "#C75197",
  "unknown"          = "grey",
  "screening" = "grey60"
)

# 2. Filter singleton subcommunities and set global ordering
subcommunity_order <- plasmid_subcommunities_meta_split |>
  group_by(pling_subcommunity_length_split) |>
  summarise(count = n()) |>
  filter(count > 10) |>
  arrange(desc(count)) |>
  pull(pling_subcommunity_length_split)
length(subcommunity_order) # 71
# consider removing hub plasmids?
subcomm_order_no_hub <- setdiff(subcommunity_order, "hub_plasmid")
subcommunity_order <- c(subcomm_order_no_hub, "hub_plasmid")

dat <- plasmid_subcommunities_meta_split |>
  filter(pling_subcommunity_length_split %in% subcommunity_order) |>
  mutate(pling_subcommunity_length_split = factor(pling_subcommunity_length_split,
                                                  levels = subcommunity_order),
         healthcare_exposure = case_when(sampletype_cat == "screening" ~ "screening",
                                         TRUE ~ healthcare_exposure)
  )


# add majority community molecule type
dat <- dat |>
  group_by(pling_subcommunity_length_split) |>
  arrange(molecule_type_genomad_majority) |>
  mutate(molecule_type_majority_cluster = min(molecule_type_genomad_majority)) |>
  ungroup()

# 3. Common theme
theme_subcomm <- theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 11),
    plot.margin = margin(3, 5, 3, 5)
  ) 


# 4. Top panel: absolute counts by sample type
molecule_cols <- genomad_molecule_type_colours

p_count <- ggplot(dat, aes(x = pling_subcommunity_length_split, fill = molecule_type_majority_cluster)) +
  geom_bar(width = 0.9) +
  scale_fill_manual(values = molecule_cols, na.value = "grey70") +
  scale_y_log10() +
  labs(
    y = "Number of plasmids",
    fill = "geNomad\nmolecule type",
    title = "Plasmid subcommunity size"
  ) +
  theme_subcomm
p_count

# 5. Function for proportional categorical panels
plot_prop <- function(data, variable, colours, title, legend_title = title) {
  ggplot(data, aes(x = pling_subcommunity, fill = .data[[variable]])) +
    geom_bar(
      position = "fill",
      width = 0.9
    ) +
    scale_y_continuous(
      labels = percent_format(),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(
      values = colours,
      na.value = "grey80"
    ) +
    labs(
      y = "Proportion",
      fill = legend_title,
      title = title
    ) +
    theme_subcomm
}


# 6. Categorical metadata panels
p_sampletype <- plot_prop(
  dat,
  "sampletype_cat",
  sampletype_individual_cols,
  "Sample type"
)

p_region <- plot_prop(
  dat,
  "region",
  colours_region,
  "Region"
)

p_species <- plot_prop(
  dat,
  "species_complex",
  species_colours,
  "Species complex"
)

p_hcai <- plot_prop(
  dat,
  "healthcare_exposure",
  colours_hcai,
  "Healthcare exposure"
)


# 8. Continuous-variable plotting function
# Boxplot + individual jittered plasmids
plot_continuous <- function(data, variable, title, ylab = title, log_y = FALSE) {
  
  p <- ggplot(data, aes(x = pling_subcommunity, y = .data[[variable]])) +
    geom_boxplot(
      width = 0.65,
      outlier.shape = NA,
      linewidth = 0.3
    ) +
    geom_jitter(
      width = 0.18,
      height = 0,
      alpha = 0.35,
      size = 0.6
    ) +
    labs(
      y = ylab,
      title = title
    ) +
    theme_subcomm
  
  if (log_y) {
    p <- p + scale_y_log10(labels = label_number())
  }
  
  p
}


# 9. Contig length
p_length <- plot_continuous(
  dat,
  "autocycler_contig_length",
  "Plasmid contig length",
  "Contig length (bp)",
  log_y = TRUE
)


# 10. AMR / virulence / stress genes
gene_variables <- c(
  "n_amr_genes",
  "n_virulence_genes",
  "n_stress_genes"
)

gene_labels <- c(
  n_amr_genes       = "AMR genes",
  n_virulence_genes = "Virulence genes",
  n_stress_genes    = "Stress-response genes"
)

gene_plots <- purrr::map(
  gene_variables,
  \(v) plot_continuous(
    dat,
    v,
    gene_labels[[v]],
    "Gene count"
  )
)


# 13. Combine everything into one vertical figure
all_plots <- c(
  list(
    p_count,
    p_sampletype,
    p_region,
    p_species,
    p_hcai,
    p_length
  ),
  gene_plots
)

final_plot <- wrap_plots(
  all_plots,
  ncol = 1,
  guides = "collect"
)

final_plot
ggsave("plasmid_cluster_summary_faceted.png", final_plot, units = "in", width = 8, height = 9, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * metadata chi-squared associations ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plasmid-level chi-squared test
metadata_vars <- c(
  "sampletype_cat",
  "region",
  "species_complex",
  "kleborate_mlst_ST",
  "healthcare_exposure"
)
# function to analyse 1 metadata variable
run_subcommunity_association <- function(data, metadata_var, B = 100000) {
  # Keep complete observations
  dat <- data |>
    select(
      pling_subcommunity_length_split,
      metadata = all_of(metadata_var)
    ) |>
    drop_na() |>
    mutate(
      pling_subcommunity_length_split = as.factor(pling_subcommunity_length_split),
      metadata = as.factor(metadata)
    ) |>
    droplevels()
  
  # Contingency table
  tab <- table(
    dat$pling_subcommunity_length_split,
    dat$metadata
  )
  
  # Remove empty rows/columns
  tab <- tab[
    rowSums(tab) > 0,
    colSums(tab) > 0,
    drop = FALSE
  ]
  
  # Chi-squared test
  #
  # simulated p-value is useful because some subcommunities
  # will have small expected counts
  chi <- chisq.test(
    tab,
    simulate.p.value = TRUE,
    B = B
  )
  # Overall effect size: Cramer's V
  n <- sum(tab)
  
  cramers_v <- sqrt(
    unname(chi$statistic) /
      (
        n *
          min(
            nrow(tab) - 1,
            ncol(tab) - 1
          )
      )
  )
  
  # Overall result
  overall <- tibble(
    metadata_variable = metadata_var,
    n = n,
    n_subcommunities = nrow(tab),
    n_metadata_categories = ncol(tab),
    chi_squared = unname(chi$statistic),
    df = (nrow(tab) - 1) * (ncol(tab) - 1),
    p_value = chi$p.value,
    cramers_v = cramers_v
  )
  # Cell-level results
  observed <- as.data.frame(as.table(chi$observed))
  expected <- as.data.frame(as.table(chi$expected))
  stdres   <- as.data.frame(as.table(chi$stdres))
  
  names(observed) <- c(
    "pling_subcommunity_length_split",
    "metadata_category",
    "observed"
  )
  
  names(expected) <- c(
    "pling_subcommunity_length_split",
    "metadata_category",
    "expected"
  )
  
  names(stdres) <- c(
    "pling_subcommunity_length_split",
    "metadata_category",
    "std_residual"
  )
  cells <- observed |>
    left_join(expected,
              by = c(
                "pling_subcommunity_length_split",
                "metadata_category"
              )
    ) |>
    left_join(stdres,
              by = c(
                "pling_subcommunity_length_split",
                "metadata_category"
              )
    ) |>
    mutate(
      metadata_variable = metadata_var,
      # Fold enrichment relative to independence
      observed_expected_ratio = observed / expected,
      # log2 enrichment is often easier to plot
      log2_observed_expected = log2(
        observed_expected_ratio
      ),
      # Approximate two-sided p-value associated
      # with standardized Pearson residual
      residual_p = 2 * pnorm(
        -abs(std_residual)
      )
    ) |>
    select(
      metadata_variable,
      pling_subcommunity_length_split,
      metadata_category,
      observed,
      expected,
      observed_expected_ratio,
      log2_observed_expected,
      std_residual,
      residual_p
    )
  
  
  return(list(overall = overall, cells = cells))
}



plasmid_subcommunities_meta_filtered_1_per_ST_with_meta <- plasmid_subcommunities_meta_split |>
  mutate(healthcare_exposure = case_when(sampletype_cat == "screening" ~ "screening",
                                         TRUE ~ healthcare_exposure)) |>
  group_by(pling_subcommunity_length_split) |>
  mutate(num_in_split_group = n()) |>
  filter(num_in_split_group >=10) |>
  group_by(pling_subcommunity_length_split, kleborate_mlst_ST, sampletype_cat, 
           region, species_complex, healthcare_exposure) |>
  slice_head() |>
  ungroup() 
nrow(plasmid_subcommunities_meta_filtered_1_per_ST_with_meta) # 2,647
nrow(plasmid_subcommunities_meta_filtered_1_per_ST) # 1,408

# run all
set.seed(123)
association_results <- purrr::map(
  metadata_vars,
  ~ run_subcommunity_association(
    plasmid_subcommunities_meta_filtered_1_per_ST_with_meta,
    metadata_var = .x,
    B = 100000
  )
)

# overall significance and effect size
overall_results <- map_dfr(association_results, "overall") |>
  mutate(p_bonferroni = p.adjust(p_value, method = "bonferroni"),# Conservative correction across the six tests
         p_fdr = p.adjust( p_value, method = "BH")) |>   # Also calculate FDR
  arrange(p_value)
overall_results
#write.csv(overall_results, "pling_subcommunity_dcj_4_metadata_assoc_chisq_overall_results.csv", row.names = FALSE)
write.csv(overall_results, "pling_split_subcommunity_dcj_4_metadata_assoc_chisq_overall_results_1_per_ST_with_meta.csv", row.names = FALSE)

# examine cell-level values
cell_results <- map_dfr(association_results,"cells")

# adjust within each metadata varibale
cell_results <- cell_results |>
  group_by(metadata_variable) |>
  mutate( residual_p_fdr = p.adjust(residual_p, method = "BH"),
          residual_p_bonferroni = p.adjust(residual_p, method = "bonferroni") ) |>
  ungroup() |>
  mutate(direction = case_when(
    std_residual > 0 ~ "Enriched",
    std_residual < 0 ~ "Depleted",
    TRUE ~ "Expected"
  )
  )
#View(cell_results)
write.csv(cell_results, "pling_split_subcommunity_dcj_4_metadata_assoc_chisq_cell_results_1_per_ST_with_meta.csv", row.names = FALSE)
# !! care statistical dependency (isolates carrying multiple plasmids are non-independent)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * Network plot ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# fix amr_beta-lactam clasification
plasmid_subcommunities_meta_split <- plasmid_subcommunities_meta_split |>
  mutate(amr_beta_lactam = case_when(amr_beta_lactam == "0" ~ "no AMR",
                                     amr_beta_lactam == "no_amr" ~ "no AMR",
                                     amr_beta_lactam == "na AMR" ~ "no AMR",
                                     amr_beta_lactam == "non_BL_amr" ~ "non beta-lactam AMR",
                                     amr_beta_lactam == "other_beta_lactam" ~ "other beta-lactam",
                                     TRUE ~ amr_beta_lactam)) 
table(plasmid_subcommunities_meta_split$amr_beta_lactam, useNA = "ifany")

# keep one plasmid per row
# One metadata row per plasmid
plasmid_metadata <- plasmid_subcommunities_meta_split |>
  distinct(plasmid, .keep_all = TRUE)


# Add length-split IDs for both plasmids in each edge
all_plasmids_distances_meta_split <-all_plasmids_distances_meta |>
  left_join(plasmid_metadata |>
              select(
                plasmid,
                pling_subcommunity_length_split
              ) |>
              rename(
                pling_subcommunity_length_split_1 =
                  pling_subcommunity_length_split
              ),
            by = c("plasmid_1" = "plasmid")
  ) |>
  left_join(plasmid_metadata |>
              select(
                plasmid,
                pling_subcommunity_length_split
              ) |>
              rename(
                pling_subcommunity_length_split_2 =
                  pling_subcommunity_length_split
              ),
            by = c("plasmid_2" = "plasmid")
  )

# count plasmids in each new split subcommunity
subcommunity_counts <- plasmid_metadata |>
  filter(
    !is.na(pling_subcommunity_length_split),
    pling_subcommunity_length_split != "hub_plasmid"
  ) |>
  dplyr::count(
    pling_subcommunity_length_split,
    name = "n_plasmids"
  ) |>
  filter(n_plasmids > 1)


# Network edges: only connect plasmids in the SAME
# length-based subcommunity
network_edges <- all_plasmids_distances_meta_split |>
  filter(
    distance <= 4,
    plasmid_1 != plasmid_2,
    !is.na(pling_subcommunity_length_split_1),
    !is.na(pling_subcommunity_length_split_2),
    pling_subcommunity_length_split_1 ==
      pling_subcommunity_length_split_2
  ) |>
  semi_join(
    subcommunity_counts,
    by = c(
      "pling_subcommunity_length_split_1" =
        "pling_subcommunity_length_split"
    )
  ) |>
  transmute(
    from = plasmid_1,
    to = plasmid_2,
    distance = distance,
    pling_subcommunity_length_split =
      pling_subcommunity_length_split_1
  ) |>
  distinct()

# Network nodes
network_nodes <- plasmid_metadata |>
  semi_join(
    subcommunity_counts,
    by = "pling_subcommunity_length_split"
  ) |>
  transmute(
    plasmid,
    pling_subcommunity,
    pling_subcommunity_length_split,
    plasmid_length = autocycler_contig_length,
    n_amr = n_amrfidner_hits,
    n_amr_classes = amr_num_classes_category,
    amr_bl_category = amr_beta_lactam
  ) |>
  mutate(
    AMR = case_when(
      n_amr == 0 ~ "0",
      n_amr == 1 ~ "1",
      n_amr >= 2 ~ "2+",
      TRUE ~ NA_character_
    ),
    AMR = factor(AMR, levels = c("0", "1", "2+"))
  )

# efine colours
colours_amr_beta_lactam <- c(
  "no AMR"            = "grey60",  # grey
  "non beta-lactam AMR"="#B8D8F8",  # pale blue
  "other beta-lactam" = "#9E9AC8",  # medium purple
  "cephalosporin"     = "#756BB1",  # stronger purple
  "carbapenem"        = "#54278F"   # dark alarm purple
)

colours_amr_classes <- c(
  "0"  = "grey50", # "#9EB8C8","#E8E8E8",  # very pale grey
  "1"  = "#B8D8F8",  # pale blue
  "2"  = "#6BAED6",  # stronger medium blue
  "3"  = "#6A51A3",  # purple
  "4+" = "#3F007D"   # deep purple
)

# LAYOUT

# Circular layout for length-split subcommunities
subcommunity_layout <- subcommunity_counts |>
  mutate(
    radius = pmax(1.5, n_plasmids^0.7 * 0.35)
  ) |>
  arrange(desc(n_plasmids))

# Number of columns
n_cols <- ceiling(sqrt(nrow(subcommunity_layout)))

# Assign circles to grid cells
subcommunity_layout <- subcommunity_layout |>
  mutate(
    grid_col = (row_number() - 1) %% n_cols + 1,
    grid_row = (row_number() - 1) %/% n_cols + 1
  )

# Column sizes
column_sizes <- subcommunity_layout |>
  group_by(grid_col) |>
  summarise(
    max_radius = max(radius),
    .groups = "drop"
  ) |>
  arrange(grid_col) |>
  mutate(
    diameter = 2 * max_radius
  )

# Row sizes
row_sizes <- subcommunity_layout |>
  group_by(grid_row) |>
  summarise(
    max_radius = max(radius),
    .groups = "drop"
  ) |>
  arrange(grid_row) |>
  mutate(
    diameter = 2 * max_radius
  )

# Gap between circles
circle_gap <- 3

# Circle centres
column_sizes <- column_sizes |>
  mutate(
    centre_x = cumsum(diameter + circle_gap) -
      diameter / 2
  )

row_sizes <- row_sizes |>
  mutate(
    centre_y = -(
      cumsum(diameter + circle_gap) -
        diameter / 2
    )
  )

# Join coordinates
subcommunity_layout <- subcommunity_layout |>
  left_join(
    column_sizes |>
      select(grid_col, centre_x),
    by = "grid_col"
  ) |>
  left_join(
    row_sizes |>
      select(grid_row, centre_y),
    by = "grid_row"
  )


# Position plasmids within their NEW subcommunities
network_nodes <- network_nodes |>
  left_join(
    subcommunity_layout |>
      select(
        pling_subcommunity_length_split,
        centre_x,
        centre_y,
        radius,
        n_plasmids
      ),
    by = "pling_subcommunity_length_split"
  ) |>
  group_by(pling_subcommunity_length_split) |>
  arrange(plasmid, .by_group = TRUE) |>
  mutate(
    angle = seq(
      0,
      2 * pi,
      length.out = n() + 1
    )[-(n() + 1)],
    x = centre_x + radius * cos(angle),
    y = centre_y + radius * sin(angle)
  ) |>
  ungroup()

# Edge coordinates
network_edges_plot <- network_edges |>
  inner_join(
    network_nodes |>
      select(plasmid, x, y),
    by = c("from" = "plasmid")
  ) |>
  rename(
    x_from = x,
    y_from = y
  ) |>
  inner_join(
    network_nodes |>
      select(plasmid, x, y),
    by = c("to" = "plasmid")
  ) |>
  rename(
    x_to = x,
    y_to = y
  )




# make plot
plasmid_network_circle <- ggplot() +
  geom_segment(data = network_edges_plot,
               aes(x = x_from,y = y_from,
                   xend = x_to, yend = y_to,
                   alpha = distance), 
               colour = "grey50", linewidth = 0.05) +
  geom_point(data = network_nodes,  aes(x = x,  y = y,size = plasmid_length, colour = amr_bl_category)) +
  scale_colour_manual(values = colours_amr_beta_lactam,   name = "Beta-lactam\nresistance\nclass") +
  scale_size_continuous(
    trans = "sqrt",
    range = c(0.5, 4),
    breaks = c(
      1000,
      10000,
      50000,
      100000,
      300000
    ), labels = scales::label_number(scale = 1e-3, suffix = " kb"),name = "Plasmid length") +
  
  # Distance 0 = strong edge; distance 4 = faint edge
  scale_alpha_continuous(
    range = c(0.25, 0.02),
    limits = c(0, 4),
    breaks = 0:4,
    name = "DCJ/indel\ndistance"
  ) +
  coord_equal() +
  theme_void(base_family = "sans") +
  theme(legend.position = "right")

plasmid_network_circle

# save
ggsave(filename = "plasmid_network_amr_classes_split_subcommunities.png", plot = plasmid_network_circle,  width = 24,  height = 12,  units = "in",  dpi = 300, bg = "white")
ggsave("plasmid_network_highres_amr_classes_split_subcommunities.png", plasmid_network_circle, width = 40,  height = 20,  units = "in",  dpi = 600,  bg = "white")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# * DCJ/indel distance 0 clusters ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# extract connected components at 0 DCJ-indels (within existing subcommunities)

# 1. Get plasmid -> existing subcommunity mapping
plasmid_communities <- plasmid_subcommunities_meta |>
  select(plasmid,pling_subcommunity) |>
  distinct() |>
  filter(
    !is.na(plasmid),
    !is.na(pling_subcommunity)
  )

# buils distnce 0 edges
dcj0_edges <- all_plasmids_distances_meta |>
  filter(
    distance == 0,
    plasmid_1 != plasmid_2
  ) |>
  
  select(
    plasmid_1,
    plasmid_2
  ) |>
  
  # Community of plasmid 1
  left_join(
    plasmid_communities,
    by = c("plasmid_1" = "plasmid")
  ) |>
  rename(
    subcommunity_1 = pling_subcommunity
  ) |>
  
  # Community of plasmid 2
  left_join(
    plasmid_communities,
    by = c("plasmid_2" = "plasmid")
  ) |>
  rename(
    subcommunity_2 = pling_subcommunity
  ) |>
  
  # Only relationships within the same existing subcommunity
  filter(
    subcommunity_1 == subcommunity_2
  ) |>
  
  transmute(
    pling_subcommunity = subcommunity_1,
    from = plasmid_1,
    to = plasmid_2
  ) |>
  distinct()

# calculate connected components
get_dcj0_components <- function(subcommunity_name) {
  # All plasmids belonging to this subcommunity
  vertices <- plasmid_communities |>
    filter( pling_subcommunity == subcommunity_name ) |>
    transmute(name = plasmid)
  
  # Distance-0 edges in this subcommunity
  edges <- dcj0_edges |>
    filter(pling_subcommunity == subcommunity_name) |>
    select( from,to )
  
  # Build graph
  g <- graph_from_data_frame(
    d = edges,
    vertices = vertices,
    directed = FALSE
  )
  
  # Connected component membership
  membership_df <- tibble(
    plasmid = names(components(g)$membership),
    component = unname(components(g)$membership)
  )
  
  # Determine component sizes and rank largest -> smallest
  component_sizes <- membership_df |>
    dplyr::count(component,  name = "component_size") |>
    arrange(  desc(component_size),  component) |>
    mutate( component_rank = row_number() )
  
  # Add names
  membership_df |>
    left_join(component_sizes,by = "component") |>
    mutate(pling_subcommunity = subcommunity_name,
           dcj_0_cluster = paste0( subcommunity_name, "_",component_rank )) |>
    select( plasmid,
            pling_subcommunity,
            dcj_0_cluster,
            dcj_0_cluster_size = component_size,
            dcj_0_cluster_rank = component_rank
    )
}

# call function for every subcommunity
dcj_0_clusters <- map_dfr(
  unique(plasmid_communities$pling_subcommunity),
  get_dcj0_components
)

nrow(dcj_0_clusters) # 5596
nrow(dcj0_edges) # 77,905
colnames(dcj0_edges) # 77,905
#View(dcj_0_clusters)
#View(dcj0_edges)

# 3. merge onto metadata
plasmid_subcommunities_meta <- plasmid_subcommunities_meta |>
  left_join(dcj_0_clusters |> select(
    plasmid,
    dcj_0_cluster,
    dcj_0_cluster_size,
    dcj_0_cluster_rank
  ), 
  by = "plasmid")


# check
table(is.na(plasmid_subcommunities_meta$dcj_0_cluster))
#View(plasmid_subcommunities_meta)
write.csv(plasmid_subcommunities_meta, "plasmid_subcommunities_meta.csv", row.names = FALSE)

## left join onto plasmid distances
all_plasmids_distances_meta_test <- all_plasmids_distances_meta |>
  left_join(dcj_0_clusters |> select(-c(dcj_0_cluster_rank)), 
            by = c("plasmid_1" = "plasmid", "pling_subcommunity_1" = "pling_subcommunity")) |>
  rename(dcj_0_cluster_1 = dcj_0_cluster,  dcj_0_cluster_size_1 = dcj_0_cluster_sizes) |>
  left_join(dcj_0_clusters |> select(-c(dcj_0_cluster_rank)), 
            by = c("plasmid_2" = "plasmid", "pling_subcommunity_2" = "pling_subcommunity")) |>
  rename(dcj_0_cluster_2 = dcj_0_cluster,
         dcj_0_cluster_size_2 = dcj_0_cluster_sizes)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# look at only non-COl like, >20kbp plasmid close clusters, with any AMR, stress, viruelnce genes. 

# plot as graph

# see if any are across different STs/ look at chromosomal relatedness


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# filter to exclude small col-like plasmids carrying no aMR genes
nrow(plasmid_subcommunities_meta) # 5,596
length(unique(plasmid_subcommunities_meta$pling_subcommunity)) # 1547 (includes hub plasmids)
table(plasmid_subcommunities_meta$num_plasmids_in_subcomm)

plasmids_non_singletons <- plasmid_subcommunities_meta |>
  group_by(pling_subcommunity) |>
  mutate(num_plasmids_in_subcomm = n()) |>
  ungroup() |>
  filter(num_plasmids_in_subcomm > 1)
nrow(plasmids_non_singletons) # 4,416 (includes 72 hub plasmids)
length(unique(plasmids_non_singletons$pling_subcommunity)) # 367 (includes hub plasmids)

# remove subcommunites with no amrfinder hits
plasmids_non_singletons_no_small_non_amr <- plasmids_non_singletons |>
  group_by(pling_subcommunity) |>
  mutate(max_plasmid_length_in_subcomm = max(autocycler_contig_length),
         sum_amrfinder_subcomm = sum(n_amrfidner_hits),
         sum_amr_subcomm = sum(n_amr_genes),
         sum_stress_subcomm = sum(n_stress_genes),
         sum_virulence_subcomm = sum(n_virulence_genes)) |>
  ungroup() |>
  filter(!(max_plasmid_length_in_subcomm < 20000 & sum_amrfinder_subcomm == 0 ))
#View(plasmids_non_singletons_no_small_non_amr) # 2853 (includes 72 hub plasmids)
nrow(plasmids_non_singletons_no_small_non_amr) # 3280 (includes 72 hub plasmids)
length(unique(plasmids_non_singletons_no_small_non_amr$pling_subcommunity)) # 282 subcommunities (includes hub plasmids)
table(plasmids_non_singletons_no_small_non_amr$num_plasmids_in_subcomm)

# now look at how may close dcj 0 clustrs
plasmids_non_singletons_no_small_non_amr_dcj_0 <- plasmids_non_singletons_no_small_non_amr |>
  filter(dcj_0_cluster_size > 1)
#View(plasmids_non_singletons_no_small_non_amr_dcj_0) 
nrow(plasmids_non_singletons_no_small_non_amr_dcj_0) # 1819 plasmids in close clusters
length(unique(plasmids_non_singletons_no_small_non_amr_dcj_0$pling_subcommunity)) # 155 subcommunities (includes hub plasmids)
length(unique(plasmids_non_singletons_no_small_non_amr_dcj_0$dcj_0_cluster)) # 324 subcommunities (includes hub plasmids)
table(plasmids_non_singletons_no_small_non_amr_dcj_0$num_plasmids_in_subcomm)


# compare to numbers if I filter the cdj clusters
dcj_clusters_filtered <- plasmid_subcommunities_meta |>
  filter(dcj_0_cluster_size > 1) |>
  group_by(dcj_0_cluster) |>
  mutate(max_plasmid_length_in_dcj_0_cluster = max(autocycler_contig_length),
         sum_amrfinder_dcj_0_cluster = sum(n_amrfidner_hits),
         sum_amr_dcj_0_cluster = sum(n_amr_genes),
         sum_stress__dcj_0_cluster = sum(n_stress_genes),
         sum_virulence__dcj_0_cluster = sum(n_virulence_genes)) |>
  ungroup() |>
  filter(!(max_plasmid_length_in_dcj_0_cluster < 20000 & sum_amrfinder_dcj_0_cluster == 0)) |>
  group_by(dcj_0_cluster) |>
  mutate(num_mlsts_in_dcj_0_cluster = n_distinct(kleborate_mlst_ST) ) |>
  ungroup()

nrow(dcj_clusters_filtered) # 737 (2879 with the clo-like removed)
#View(dcj_clusters_filtered) 
length(unique(dcj_clusters_filtered$pling_subcommunity)) # 149 subcommunities (includes hub plasmids)
length(unique(dcj_clusters_filtered$dcj_0_cluster)) # 227 close clusters (includes hub plasmids)
table(dcj_clusters_filtered$dcj_0_cluster_size)
#table(dcj_clusters_filtered$num_plasmids_in_subcomm) # the 44 rows with 1 in subcommunity are the hub plasmids

dcj_0_cluster_sizes <- dcj_clusters_filtered |>
  group_by(dcj_0_cluster) |>
  summarise(size = n())
summary(dcj_0_cluster_sizes$size)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# plot
dcj_clusters_filtered <- dcj_clusters_filtered |>
  mutate(amr_beta_lactam = case_when(amr_beta_lactam == "0" ~ "no_amr",
                                     amr_beta_lactam == "no_amr" ~ "no AMR",
                                     amr_beta_lactam == "na AMR" ~ "no AMR",
                                     amr_beta_lactam == "non_BL_amr" ~ "non beta-lactam AMR",
                                     amr_beta_lactam == "other_beta_lactam" ~ "other beta-lactam",
                                     TRUE ~ amr_beta_lactam),
         multi_ST = case_when(num_mlsts_in_dcj_0_cluster == 1 ~ "single-ST",
                              num_mlsts_in_dcj_0_cluster > 1 ~ "multi-ST")) 
table(dcj_clusters_filtered$amr_beta_lactam, useNA = "ifany")
table(dcj_clusters_filtered$multi_ST, useNA = "ifany")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot code
# 1. Identify DCJ-0 clusters to plot
dcj_cluster_counts <- dcj_clusters_filtered |>
  distinct(dcj_0_cluster, plasmid, multi_ST) |>
  dplyr::count( dcj_0_cluster, multi_ST, name = "n_plasmids") |>
  filter(n_plasmids > 1)

# Inspect cluster sizes
dcj_cluster_counts |>
  arrange(desc(n_plasmids))


# 2. create nodes
network_nodes <- dcj_clusters_filtered |>
  semi_join(
    dcj_cluster_counts,
    by = "dcj_0_cluster"
  ) |>
  transmute(
    plasmid,
    dcj_0_cluster,
    plasmid_length = autocycler_contig_length,
    n_amr = n_amrfidner_hits,
    n_amr_classes = amr_num_classes_category,
    amr_beta_lactam = amr_beta_lactam
  ) |>
  distinct(
    plasmid,
    .keep_all = TRUE
  )

# set AMR factor:
network_nodes <- network_nodes |>
  mutate(n_amr_classes = factor(n_amr_classes,
                                levels = c("0", "1", "2", "3", "4+")))

# create edges using only those plasmids
dcj_cluster_membership <- dcj_clusters_filtered |>
  distinct( plasmid, dcj_0_cluster)

# network edges
network_edges <- all_plasmids_distances_meta |>
  filter(
    distance == 0,
    plasmid_1 != plasmid_2
  ) |>
  
  # Only retain plasmid_1 if it occurs in dcj_clusters_filtered
  inner_join(
    dcj_cluster_membership,
    by = c("plasmid_1" = "plasmid")
  ) |>
  rename(
    dcj_0_cluster_1 = dcj_0_cluster
  ) |>
  
  # Only retain plasmid_2 if it occurs in dcj_clusters_filtered
  inner_join(
    dcj_cluster_membership,
    by = c("plasmid_2" = "plasmid")
  ) |>
  rename(
    dcj_0_cluster_2 = dcj_0_cluster
  ) |>
  
  # This should be true by definition, but enforce it
  filter(
    dcj_0_cluster_1 == dcj_0_cluster_2
  ) |>
  
  transmute(
    from = plasmid_1,
    to = plasmid_2,
    distance,
    dcj_0_cluster = dcj_0_cluster_1
  ) |>
  distinct()


# layout code
dcj_cluster_layout <- dcj_cluster_counts |>
  mutate(
    radius = pmax(
      1.5,
      n_plasmids^0.7 * 0.35
    )
  ) |>
  arrange(desc(n_plasmids))


# Number of columns in grid
n_cols <- ceiling(sqrt(nrow(dcj_cluster_layout)))


# Assign each cluster to a grid position
dcj_cluster_layout <- dcj_cluster_layout |>
  mutate(
    grid_col = (row_number() - 1) %% n_cols + 1,
    grid_row = (row_number() - 1) %/% n_cols + 1
  )


# Work out how much horizontal space each column needs
column_sizes <- dcj_cluster_layout |>
  group_by(grid_col) |>
  summarise(
    max_radius = max(radius),
    .groups = "drop"
  ) |>
  arrange(grid_col) |>
  mutate(
    diameter = 2 * max_radius
  )


# Work out how much vertical space each row needs
row_sizes <- dcj_cluster_layout |>
  group_by(grid_row) |>
  summarise(
    max_radius = max(radius),
    .groups = "drop"
  ) |>
  arrange(grid_row) |>
  mutate(
    diameter = 2 * max_radius
  )


# Calculate circle centres
# Space between neighbouring cluster circles
circle_gap <- 3


column_sizes <- column_sizes |>
  mutate(
    centre_x =
      cumsum(diameter + circle_gap) -
      diameter / 2
  )


row_sizes <- row_sizes |>
  mutate(
    centre_y = -(
      cumsum(diameter + circle_gap) -
        diameter / 2
    )
  )


# Add centre coordinates to cluster layout
dcj_cluster_layout <- dcj_cluster_layout |>
  left_join(
    column_sizes |>
      select(
        grid_col,
        centre_x
      ),
    by = "grid_col"
  ) |>
  left_join(
    row_sizes |>
      select(
        grid_row,
        centre_y
      ),
    by = "grid_row"
  )

# check
dcj_cluster_layout |>
  select( dcj_0_cluster, n_plasmids,radius,  centre_x, centre_y)

# make df to add cluster outlines
cluster_backgrounds <- dcj_cluster_layout |>
  transmute(
    dcj_0_cluster,
    centre_x,
    centre_y,
    radius = radius + 0.8,
    multi_ST
  )



# position nodes
network_nodes <- network_nodes |>
  left_join(
    dcj_cluster_layout |>
      select(
        dcj_0_cluster,
        centre_x,
        centre_y,
        radius,
        n_plasmids
      ),
    by = "dcj_0_cluster"
  ) |>
  group_by(dcj_0_cluster) |>
  arrange(plasmid) |>
  mutate(
    angle = seq(
      0,
      2 * pi,
      length.out = n() + 1
    )[-(n() + 1)],
    
    x = centre_x + radius * cos(angle),
    y = centre_y + radius * sin(angle)
  ) |>
  ungroup()


# check
c(
  filtered_plasmids = n_distinct(dcj_clusters_filtered$plasmid),
  network_nodes = n_distinct(network_nodes$plasmid),
  clusters = n_distinct(network_nodes$dcj_0_cluster),
  edges = nrow(network_edges)
)


# 11. Give edges node coordinates
network_edges_plot <- network_edges |>
  
  # Coordinates of first plasmid
  inner_join(
    network_nodes |>
      select(
        plasmid,
        x,
        y
      ),
    by = c(
      "from" = "plasmid"
    )
  ) |>
  
  rename(
    x_from = x,
    y_from = y
  ) |>
  
  # Coordinates of second plasmid
  inner_join(
    network_nodes |>
      select(
        plasmid,
        x,
        y
      ),
    by = c(
      "to" = "plasmid"
    )
  ) |>
  
  rename(
    x_to = x,
    y_to = y
  )

library(ggforce)

# use subtle colours
cluster_outline_colours <- c(
  "single-ST" = "grey75",
  "multi-ST"  = "#D55E00"
)
cluster_fill <- c(
  "single-ST" = "grey96",
  "multi-ST" = "#FEE8C8"
)

colours_amr_beta_lactam <- c(
  "no AMR"            =  "grey60",  # grey
  "non beta-lactam AMR"        = "#B8D8F8",  # pale blue
  "other beta-lactam" = "#9E9AC8",  # medium purple
  "cephalosporin"     = "#756BB1",  # stronger purple
  "carbapenem"        = "#54278F"   # dark alarm purple
)

# 13. Plot
plasmid_network_circle <- ggplot() +
  # add circles
  geom_circle(data = cluster_backgrounds,
              aes(x0 = centre_x,  y0 = centre_y, r = radius, fill = multi_ST, colour = multi_ST),
              linewidth = 0.2,   alpha = 0.35,  inherit.aes = FALSE) +
  scale_colour_manual(values = cluster_outline_colours,
                      name = "Cluster lineage") +
  scale_fill_manual(values = cluster_fill,
                    name = "Cluster lineage") +
  ggnewscale::new_scale_colour() +
  ggnewscale::new_scale_fill() +
  # Distance-0 connections
  geom_segment( data = network_edges_plot,
                aes(  x = x_from,   y = y_from, xend = x_to,yend = y_to),
                colour = "grey60",linewidth = 0.04,   alpha = 0.15) +
  
  # Plasmids
  geom_point(   data = network_nodes,
                aes(
                  x = x,
                  y = y,
                  size = plasmid_length,
                  colour = amr_beta_lactam
                )
  ) +
  scale_colour_manual(values = colours_amr_beta_lactam,
                      name = "Beta-lactam\nresistance\ncategory",  drop = FALSE ) +
  scale_size_continuous(
    trans = "sqrt",
    range = c(0.5, 4),
    breaks = c(
      1000,
      10000,
      50000,
      100000,
      300000
    ),
    labels = scales::label_number(
      scale = 1e-3,
      suffix = " kb"
    ),  name = "Plasmid length") +
  coord_equal() +
  theme_void( base_family = "sans") +
  theme(  legend.position = "right" )
plasmid_network_circle


# 14. Save
ggsave( filename =  "plasmid_dcj0_cluster_network_amr_classes.png",
        plot = plasmid_network_circle, width = 24,height = 12, units = "in",dpi = 300,  bg = "white")


ggsave( filename =  "plasmid_dcj0_cluster_network_highres_amr_classes.png",
        plot = plasmid_network_circle,  width = 40, height = 20,units = "in",  dpi = 600,  bg = "white")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# summarise dcj-0 clustres
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# Most common value
mode_value <- function(x){
  
  x <- x[!is.na(x)]
  
  if(length(x) == 0)
    return(NA_character_)
  
  names(sort(table(x), decreasing = TRUE))[1]
}


# Proportion represented by the most common value
mode_prop <- function(x){
  
  x <- x[!is.na(x)]
  
  if(length(x) == 0)
    return(NA_real_)
  
  max(table(x))/length(x)
}

# 1 summary row per dcj
dcj_cluster_summary <- dcj_clusters_filtered |>
  
  group_by(dcj_0_cluster) |>
  
  summarise(
    
    
    ## Size
    
    
    cluster_size = n(),
    
    n_unique_samples =
      n_distinct(unique_sample),
    
    
    ## Lineages
    
    
    n_ST =
      n_distinct(kleborate_mlst_ST),
    
    dominant_ST =
      mode_value(kleborate_mlst_ST),
    
    dominant_ST_prop =
      mode_prop(kleborate_mlst_ST),
    
    n_species =
      n_distinct(species_complex),
    
    dominant_species =
      mode_value(species_complex),
    
    dominant_species_prop =
      mode_prop(species_complex),
    
    
    ## Sample metadata
    
    
    n_regions =
      n_distinct(region),
    
    dominant_region =
      mode_value(region),
    
    dominant_region_prop =
      mode_prop(region),
    
    n_sampletypes =
      n_distinct(sampletype_cat),
    
    dominant_sampletype =
      mode_value(sampletype_cat),
    
    dominant_sampletype_prop =
      mode_prop(sampletype_cat),
    
    n_healthcare =
      n_distinct(healthcare_exposure),
    
    dominant_healthcare =
      mode_value(healthcare_exposure),
    
    dominant_healthcare_prop =
      mode_prop(healthcare_exposure),
    
    
    ## Plasmid biology
    
    
    dominant_replicon =
      mode_value(combo_rep_types),
    
    dominant_replicon_prop =
      mode_prop(combo_rep_types),
    
    dominant_primary_cluster =
      mode_value(combo_primary_cluster_id),
    
    dominant_primary_cluster_prop =
      mode_prop(combo_primary_cluster_id),
    
    dominant_secondary_cluster =
      mode_value(combo_secondary_cluster_id),
    
    dominant_secondary_cluster_prop =
      mode_prop(combo_secondary_cluster_id),
    
    
    ## Length
    
    
    mean_plasmid_length =
      mean(autocycler_contig_length),
    
    median_plasmid_length =
      median(autocycler_contig_length),
    
    max_plasmid_length =
      max(autocycler_contig_length),
    
    
    ## Gene content
    
    
    sum_amrfinder =
      sum(n_amrfidner_hits),
    
    sum_amr =
      sum(n_amr_genes),
    
    sum_virulence =
      sum(n_virulence_genes),
    
    sum_stress =
      sum(n_stress_genes),
    
    
    ## Unique AMR classes
    
    
    amr_classes =
      paste(
        sort(
          unique(
            unlist(
              strsplit(
                na.omit(amr_classes),
                ";"
              )
            )
          )
        ),
        collapse = ";"
      ),
    
    n_unique_amr_classes =
      length(
        unique(
          unlist(
            strsplit(
              na.omit(amr_classes),
              ";"
            )
          )
        )
      ),
    
    .groups = "drop"
  )


# id multi-clustr STs
dcj_cluster_summary <- dcj_cluster_summary |>
  mutate(
    lineage_class =
      if_else(
        n_ST == 1,
        "Single ST",
        "Multiple STs"
      )
  )

# descriptive stats
dcj_cluster_summary |>
  group_by(lineage_class) |>
  summarise(
    n_clusters = n(),
    median_cluster_size =   median(cluster_size),
    mean_cluster_size =     mean(cluster_size),
    median_plasmid_length =     median(max_plasmid_length),
    median_amrfinder =    median(sum_amrfinder),
    median_amr =  median(sum_amr)
  )

dcj_cluster_summary <- dcj_cluster_summary |>
  mutate(amr_mean_classes = n_unique_amr_classes / cluster_size,
         amr_mean_amr = sum_amr / cluster_size,
         amr_mean_stress = sum_stress / cluster_size,
         amr_mean_virulence = sum_virulence / cluster_size,
  )
# cluster sizes
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Single ST",]$cluster_size)
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Multiple STs",]$cluster_size)
wilcox.test(cluster_size ~ lineage_class, data = dcj_cluster_summary)
# W = 8816, p-value = 1.886e-10

# length
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Multiple STs",]$median_plasmid_length)
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Single ST",]$median_plasmid_length)
wilcox.test(median_plasmid_length ~ lineage_class,data = dcj_cluster_summary)
# W = 3676, p-value = 1.523e-07
# so cluster length and size different but AMR burden not

# amr classes per plasmid
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Multiple STs",]$amr_mean_classes) # 0 (0-4) 0-10
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Single ST",]$amr_mean_classes) # 0 (0-9) 0-42
wilcox.test(amr_mean_classes ~ lineage_class, data = dcj_cluster_summary) # <0.0001


# amr genes per plasmid
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Multiple STs",]$amr_mean_amr) # 0 (0-4) 0-10
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Single ST",]$amr_mean_amr) # 0 (0-9) 0-42
wilcox.test(amr_mean_amr ~ lineage_class, data = dcj_cluster_summary) # # W = 5700.5, p-value = 0.2357
# p=0.0516

# viruence gene sper plasmid
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Multiple STs",]$amr_mean_virulence) # 0 (0-4) 0-10
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Single ST",]$amr_mean_virulence) # 0 (0-9) 0-42
wilcox.test(amr_mean_virulence ~ lineage_class, data = dcj_cluster_summary) # # W = 5700.5, p-value = 0.2357
# p=0.00426


# stess genes per plasmid
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Multiple STs",]$amr_mean_stress) # 0 (0-4) 0-10
summary(dcj_cluster_summary[dcj_cluster_summary$lineage_class == "Single ST",]$amr_mean_stress) # 0 (0-9) 0-42
wilcox.test(amr_mean_stress ~ lineage_class, data = dcj_cluster_summary) # # W = 5700.5, p-value = 0.2357
# p=0.00426

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# make publication df
dcj_cluster_summary_pub <- dcj_cluster_summary |>
  filter(cluster_size >= 5) |>
  arrange(desc(cluster_size)) |>
  mutate(dom_spp = paste0(dominant_species, " (", round(dominant_species_prop*100) ,"%)"),
         dom_ST = paste0(dominant_ST, " (", round(dominant_ST_prop*100) ,"%)"),
         dom_sampletype = paste0(dominant_sampletype, " (", round(dominant_sampletype_prop*100) ,"%)"),
         dom_region = paste0(dominant_region, " (", round(dominant_region_prop*100) ,"%)"),,
         dom_hcai = paste0(dominant_healthcare, " (", round(dominant_healthcare_prop*100) ,"%)"),
         MOB_cluster = paste0(dominant_primary_cluster, ":", dominant_secondary_cluster)
  ) |>
  select(dcj_0_cluster, cluster_size, lineage_class, dominant_replicon, MOB_cluster, median_plasmid_length,
         dom_spp, dom_ST, dom_sampletype, dom_region, dom_hcai, 
         amr_classes, amr_mean_classes, amr_mean_amr, amr_mean_virulence, amr_mean_stress)
#select(dcj_0_clusters, cluster_size, median_plasmid_ lineage class, n_STs, )
#View(dcj_cluster_summary_pub)
write.csv(dcj_cluster_summary_pub, "dcj_0_cluster_summary_pub_5_or_more.csv", row.names = FALSE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
dcj_cluster_summary_long <- dcj_cluster_summary |>
  mutate(amr_mean_classes = n_unique_amr_classes / cluster_size,
         amr_mean_amr = sum_amr / cluster_size,
         amr_mean_stress = sum_stress / cluster_size,
         amr_mean_virulence = sum_virulence / cluster_size,
  ) |>
  select(dcj_0_cluster, lineage_class, cluster_size, median_plasmid_length, amr_mean_classes,
         amr_mean_amr, amr_mean_virulence, amr_mean_stress) |>
  pivot_longer(cols = c(cluster_size, median_plasmid_length, amr_mean_classes,
                        amr_mean_amr, amr_mean_virulence, amr_mean_stress), 
               names_to = "var", values_to = "value") |>
  mutate(var = case_when(var == "cluster_size" ~ "Cluster size",
                         var == "median_plasmid_length" ~ "Median plasmid length",
                         var == "amr_mean_classes" ~ "Mean number of\nantibiotic class resistances",
                         var == "amr_mean_amr" ~ "Mean number of\nAMR genes per plasmids",
                         var == "amr_mean_virulence" ~ "Mean number of\nvirulence genes per plasmids",
                         var == "amr_mean_stress" ~ "Mean number of\nstress genes per plasmids",
                         TRUE ~ NA))
#View(dcj_cluster_summary_long)
dcj_cluster_summary_long$var <- factor(dcj_cluster_summary_long$var, 
                                       levels = c("Cluster size",
                                                  "Median plasmid length",
                                                  "Mean number of\nantibiotic class resistances",
                                                  "Mean number of\nAMR genes per plasmids",
                                                  "Mean number of\nvirulence genes per plasmids",
                                                  "Mean number of\nstress genes per plasmids"
                                       ))

multi_st_colours <- c(
  "Single ST" = "grey75",
  "Multiple STs"  = "#D55E00"
)
multi_STs_fill <- c(
  "Single ST"  = "grey96",
  "Multiple STs" = "#FEE8C8"
)

multi_vs_sinlge_cluster_plot <- ggplot(data = dcj_cluster_summary_long,  aes(x = lineage_class, y = value, colour = lineage_class, fill = lineage_class)) +
  geom_point(size = 0.8, alpha = 0.8, position = position_jitter(width = 0.1)) +
  geom_boxplot(fill = NA, width = 0.8) +
  facet_wrap( ~ var, scales = "free", ncol = 3) +
  scale_fill_manual(values = multi_STs_fill) +
  scale_colour_manual(values = multi_st_colours) +
  labs(x = "ST specificity of DCJ/indel-0 plasmid close cluster",
       y = "",
       fill = "ST specificity",
       colour = "ST specificity") +
  theme_minimal()

multi_vs_sinlge_cluster_plot
ggsave("multi_vs_sinlge_cluster_plot_dcj_0_close_clusters.png", multi_vs_sinlge_cluster_plot, 
       units = "in", width = 10, height = 8, dpi = 300)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


# function to test enrichment
run_cluster_fisher <- function(data,
                               variable){
  
  categories <- unique(data[[variable]])
  
  map_dfr(categories,function(cat){
    
    dat <- data |>
      
      mutate(
        
        in_category =
          .data[[variable]] == cat
      )
    
    map_dfr(
      unique(dat$dcj_0_cluster),
      
      function(cluster){
        
        dat2 <- dat |>
          
          mutate(
            
            in_cluster =
              dcj_0_cluster == cluster
          )
        
        tab <-
          
          table(
            dat2$in_category,
            dat2$in_cluster
          )
        
        ft <- fisher.test(tab)
        
        tibble(
          
          cluster,
          
          variable,
          
          category = cat,
          
          odds_ratio =
            unname(ft$estimate),
          
          conf_low =
            ft$conf.int[1],
          
          conf_high =
            ft$conf.int[2],
          
          p =
            ft$p.value
        )
        
      })
    
  })
  
}

# run for all 4 metadata variables
metadata_vars <- c(
  "region",
  "sampletype_cat",
  "healthcare_exposure",
  "species_complex"
)

# one observation per ST within each cluster
dcj_cluster_ST <- dcj_clusters_filtered |>
  distinct(
    dcj_0_cluster,
    kleborate_mlst_ST,
    region,
    sampletype_cat,
    healthcare_exposure,
    species_complex
  )

#nrow(dcj_cluster_ST) # 588
#nrow(dcj_clusters_filtered) # 737
cluster_associations <- map_dfr(metadata_vars,
                                ~run_cluster_fisher(
                                  dcj_cluster_ST,
                                  .x
                                )) |>
  
  mutate(p_fdr = p.adjust(p,method = "BH"))
View(cluster_associations)

cluster_assoc_sig <- cluster_associations |>
  filter(p_fdr < 0.05)

sig_screening_plasmids <- plasmid_subcommunities_meta |>
  filter(dcj_0_cluster == "community_0_subcommunity_1108_1") 

sig_screening_plasmid_amr <- amrfinder_metadata_autocycler_contigs |>
  mutate(plasmid_id = paste0(sample, "_", autocycler_contig_id)) |>
  filter(plasmid_id %in% sig_screening_plasmids$plasmid)
#View(amrfinder_metadata_autocycler_contigs)
#View(sig_screening_plasmid_amr)
nrow(sig_screening_plasmid_amr) # 20

# export
sig_contigs_lookup <- sig_screening_plasmid_amr |>
  select(sample, run.x, Contig.id)

write.table(sig_contigs_lookup, file = "sig_contigs_lookup.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
)

# check MGEs
#mges <- read.csv("mges_cleaned.csv")
sig_kaerogenes_plasmids <- plasmid_subcommunities_meta |>
  filter(dcj_0_cluster == "community_100_subcommunity_0_1") 

sig_kaerogenes_plasmid_amr <- amrfinder_metadata_autocycler_contigs |>
  mutate(plasmid_id = paste0(sample, "_", autocycler_contig_id)) |>
  filter(plasmid_id %in% sig_kaerogenes_plasmids$plasmid)
#View(sig_kaerogenes_plasmid_amr)
# 2 camb, ox isolates, carried no AMR
# CB182, AX309


# get AMR list for dcj_0 clusters
amr_selected <- amrfinder_metadata_autocycler_contigs |>
  mutate(plasmid_id = paste0(sample, "_", autocycler_unitig_id))

dcj_amr <- plasmid_subcommunities_meta |>
  left_join(amr_selected, by = c("plasmid" = "plasmid_id")) |>
  filter(!is.na(dcj_0_cluster),
         Type == "AMR") |>
  group_by(dcj_0_cluster, dcj_0_cluster_size) |>
  summarise(gene_list = paste(unique(Element.symbol), collapse = ","))
View(dcj_amr)
