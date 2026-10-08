# DENV plasma proteomics: processing and quantification
# This script performs precursor filtering, ComBat batch correction, MaxLFQ
# quantification, concentration estimation, contaminant filtering, and export.

# ---- Project configuration -------------------------------------------------

project_dir <- "PATH/TO/DENV_PROJECT"

input_dir <- file.path(project_dir, "data", "raw_data")
processed_dir <- file.path(project_dir, "data", "processed_data")
source_dir <- file.path(project_dir, "data", "source_data")

dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Libraries -------------------------------------------------------------

library(diann)
library(tidyverse)
library(arrow)
library(sva)
library(ggpubr)
library(readxl)

# ---- Input files -----------------------------------------------------------

raw_report_file <- file.path(
  input_dir,
  "20250923_DENV_blood7.4_uniprot_ISO_poly_report.parquet"
)

sampling_file <- file.path(source_dir, "sampling_contaminant.csv")
metadata_file <- file.path(source_dir, "cohort_metadata.xlsx")
proteome_file <- file.path(
  source_dir,
  "231122_uniprotkb_proteome_UP000005640_2023_11_16.csv"
)
reference_protein_file <- file.path(source_dir, "ref-prot-list.csv")
plate_info_file <- file.path(source_dir, "sample_table_plates_info.csv")

fulldf <- read_parquet(raw_report_file)
sampling <- read_csv(sampling_file, show_col_types = FALSE)
metadata <- read_excel(metadata_file)
MM <- read_csv(proteome_file, show_col_types = FALSE)
Ref_prot <- read_csv(reference_protein_file, show_col_types = FALSE)
plate_info <- read_csv(plate_info_file, show_col_types = FALSE)

# ---- Prepare the DIA-NN report ---------------------------------------------

fulldf <- fulldf %>%
  mutate(Run = str_replace(Run, "_S\\d+.*", ""))

# ---- Filter low-depth runs -------------------------------------------------

fulldf <- fulldf %>%
  group_by(Run) %>%
  mutate(Precursor.Count = n_distinct(Precursor.Id)) %>%
  ungroup()

threshold <- median(fulldf$Precursor.Count, na.rm = TRUE) * 0.6

fulldf <- fulldf %>%
  filter(Precursor.Count >= threshold)

# ---- Filter high-confidence precursors -------------------------------------

X2 <- fulldf %>%
  filter(
    Q.Value <= 0.01,
    Lib.Q.Value <= 0.01,
    Lib.PG.Q.Value <= 0.01
  ) %>%
  group_by(Genes) %>%
  mutate(Total.Peptide.Count = n_distinct(Stripped.Sequence)) %>%
  ungroup()

X2_HighQual <- X2 %>%
  filter(Total.Peptide.Count >= 2)

oxidation_peptides <- X2_HighQual %>%
  filter(str_detect(Precursor.Id, "\\(UniMod:35\\)")) %>%
  pull(Precursor.Id) %>%
  unique()

X2_HighQual <- X2_HighQual %>%
  filter(!Precursor.Id %in% oxidation_peptides)

# ---- Add sample and plate information --------------------------------------

X2_HighQual <- X2_HighQual %>%
  rename(File.Name = Run) %>%
  mutate(
    Group = case_when(
      str_detect(File.Name, "_A\\d+") ~ "AP",
      str_detect(File.Name, "_C\\d+") ~ "CP",
      str_detect(File.Name, "_E\\d+") ~ "ERP",
      str_detect(File.Name, "_R\\d+") ~ "RP",
      str_detect(File.Name, "_HD\\d+") ~ "HD",
      str_detect(File.Name, "_AMSBIO_\\d+") ~ "AMSBIO",
      str_detect(File.Name, "_QC_\\d+") ~ "QC",
      TRUE ~ NA_character_
    ),
    Prob.Nr. = case_when(
      str_detect(File.Name, "_A(\\d+)") ~ str_extract(File.Name, "_A(\\d+)") %>% str_extract("\\d+"),
      str_detect(File.Name, "_C(\\d+)") ~ str_extract(File.Name, "_C(\\d+)") %>% str_extract("\\d+"),
      str_detect(File.Name, "_E(\\d+)") ~ str_extract(File.Name, "_E(\\d+)") %>% str_extract("\\d+"),
      str_detect(File.Name, "_R(\\d+)") ~ str_extract(File.Name, "_R(\\d+)") %>% str_extract("\\d+"),
      str_detect(File.Name, "_HD(\\d+)") ~ str_extract(File.Name, "HD(\\d+)"),
      str_detect(File.Name, "_AMSBIO_\\d+") ~ "AMSBIO",
      str_detect(File.Name, "_QC_\\d+") ~ "QC",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(plate_info, by = "File.Name")

# ---- Prepare peptide matrix for ComBat -------------------------------------

input_data_removeBatchEffect_log2 <- X2_HighQual %>%
  unite(
    ProteinGroups_ModifiedSequence,
    Protein.Group,
    Modified.Sequence,
    sep = "_|_"
  ) %>%
  select(
    ProteinGroups_ModifiedSequence,
    File.Name,
    Precursor.Normalised
  ) %>%
  group_by(ProteinGroups_ModifiedSequence, File.Name) %>%
  summarise(
    Precursor.Normalised = mean(Precursor.Normalised, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = File.Name,
    values_from = Precursor.Normalised
  ) %>%
  as.data.frame()

rownames(input_data_removeBatchEffect_log2) <-
  input_data_removeBatchEffect_log2$ProteinGroups_ModifiedSequence

input_data_removeBatchEffect_log2 <-
  input_data_removeBatchEffect_log2[, -1, drop = FALSE]

input_data_removeBatchEffect_log2 <- input_data_removeBatchEffect_log2 %>%
  mutate(across(everything(), ~ {
    x <- .x
    min_val <- suppressWarnings(min(x[x > 0], na.rm = TRUE))

    if (is.infinite(min_val)) {
      return(x)
    }

    x[is.na(x)] <- min_val / 2
    x[x == 0] <- min_val / 2
    log2(x)
  }))

# ---- ComBat batch correction ------------------------------------------------

sample_metadata <- X2_HighQual %>%
  select(File.Name, Group, Plate) %>%
  distinct()

input_batch <- sample_metadata$Plate[
  match(colnames(input_data_removeBatchEffect_log2), sample_metadata$File.Name)
]

number_of_cores <- max(1, parallel::detectCores() - 2)

output_ComBat_BatchEffects <- sva::ComBat(
  dat = as.matrix(input_data_removeBatchEffect_log2),
  batch = input_batch,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = TRUE,
  mean.only = FALSE,
  ref.batch = "P02",
  BPPARAM = BiocParallel::SnowParam(number_of_cores)
)

# ---- Convert corrected peptide intensities ---------------------------------

output_ComBat_BatchEffects_tibble <- as_tibble(output_ComBat_BatchEffects) %>%
  mutate(tmp = rownames(output_ComBat_BatchEffects)) %>%
  separate(
    tmp,
    into = c("Protein.Group", "Modified.Sequence"),
    sep = "_\\|_"
  ) %>%
  relocate(Protein.Group, Modified.Sequence)

output_ComBat_BatchEffects_tibble[, -(1:2)] <-
  lapply(
    output_ComBat_BatchEffects_tibble[, -(1:2)],
    function(x) as.numeric(2^x)
  )

peptide_intensity_ComBat <- output_ComBat_BatchEffects_tibble %>%
  pivot_longer(
    cols = -(Protein.Group:Modified.Sequence),
    names_to = "File.Name",
    values_to = "batch_adjusted_peptide_intensity"
  )

write_csv(
  peptide_intensity_ComBat,
  file.path(processed_dir, "full_report_batch-corrected.csv")
)

peptide_intensity <- X2_HighQual %>%
  left_join(
    peptide_intensity_ComBat,
    by = c("Protein.Group", "Modified.Sequence", "File.Name")
  )

# ---- Peptide-level MaxLFQ --------------------------------------------------

peptides_maxlfq <- diann_maxlfq(
  peptide_intensity[
    peptide_intensity$Q.Value <= 1 &
      peptide_intensity$PG.Q.Value <= 1,
  ],
  group.header = "Stripped.Sequence",
  id.header = "Precursor.Id",
  quantity.header = "batch_adjusted_peptide_intensity"
) %>%
  as_tibble(rownames = "Precursor.Id")

peptide_lookup <- peptide_intensity %>%
  select(
    Stripped.Sequence,
    Protein.Group,
    Modified.Sequence,
    Protein.Ids,
    Protein.Names,
    Genes,
    Proteotypic,
    Precursor.Charge
  ) %>%
  distinct(Stripped.Sequence, .keep_all = TRUE)

peptides_maxlfq <- peptides_maxlfq %>%
  left_join(
    peptide_lookup,
    by = c("Precursor.Id" = "Stripped.Sequence")
  ) %>%
  select(
    Protein.Group,
    Protein.Ids,
    Protein.Names,
    Genes,
    Proteotypic,
    Stripped.Sequence = Precursor.Id,
    Modified.Sequence,
    everything()
  )

peptides_maxlfq_unique <- peptides_maxlfq %>%
  filter(Proteotypic == 1)

write_csv(
  peptides_maxlfq,
  file.path(processed_dir, "Peptides_LFQ_DENV_Combat.csv")
)

write_csv(
  peptides_maxlfq_unique,
  file.path(processed_dir, "Peptides_unique_LFQ_DENV_Combat.csv")
)

# ---- Protein-level MaxLFQ --------------------------------------------------

LFQ_ComBat <- diann_maxlfq(
  peptide_intensity,
  sample.header = "File.Name",
  group.header = "Genes",
  id.header = "Precursor.Id",
  quantity.header = "batch_adjusted_peptide_intensity"
) %>%
  as_tibble(rownames = "Genes")

LFQ_ComBat_long <- LFQ_ComBat %>%
  pivot_longer(
    cols = -Genes,
    names_to = "Run",
    values_to = "LFQ"
  ) %>%
  mutate(
    Group = case_when(
      str_detect(Run, "_A\\d+") ~ "AP",
      str_detect(Run, "_C\\d+") ~ "CP",
      str_detect(Run, "_E\\d+") ~ "ERP",
      str_detect(Run, "_R\\d+") ~ "RP",
      str_detect(Run, "_HD\\d+") ~ "HD",
      str_detect(Run, "_AMSBIO_\\d+") ~ "AMSBIO",
      str_detect(Run, "_QC_\\d+") ~ "QC",
      TRUE ~ NA_character_
    ),
    Prob.Nr. = case_when(
      str_detect(Run, "_A(\\d+)") ~ str_extract(Run, "_A(\\d+)") %>% str_extract("\\d+"),
      str_detect(Run, "_C(\\d+)") ~ str_extract(Run, "_C(\\d+)") %>% str_extract("\\d+"),
      str_detect(Run, "_E(\\d+)") ~ str_extract(Run, "_E(\\d+)") %>% str_extract("\\d+"),
      str_detect(Run, "_R(\\d+)") ~ str_extract(Run, "_R(\\d+)") %>% str_extract("\\d+"),
      str_detect(Run, "_HD(\\d+)") ~ str_extract(Run, "HD(\\d+)"),
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Age = metadata$Age[match(Prob.Nr., metadata$Prob.Nr.)],
    Sex = metadata$Sex[match(Prob.Nr., metadata$Prob.Nr.)],
    Classification = metadata$Classification[
      match(Prob.Nr., metadata$Prob.Nr.)
    ]
  )

# ---- Reference-protein calibration -----------------------------------------

LFQ_ComBat_long <- LFQ_ComBat_long %>%
  mutate(
    logI = log10(LFQ),
    Mass = MM$Mass[match(Genes, MM$Gene.Names)]
  )

amsbio_reference <- LFQ_ComBat_long %>%
  filter(Group == "AMSBIO") %>%
  mutate(
    LogC = Ref_prot$LogC[match(Genes, Ref_prot$Genes)]
  ) %>%
  filter(!is.na(LogC))

calibration_fit <- lm(logI ~ LogC, data = amsbio_reference)

slope <- coef(calibration_fit)[["LogC"]]
intercept <- coef(calibration_fit)[["(Intercept)"]]

calibration_plot <- ggplot(
  amsbio_reference,
  aes(x = LogC, y = logI)
) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Linear Regression",
    x = "logC",
    y = "logI"
  ) +
  stat_regline_equation(
    label.x = 3,
    label.y = 2.5,
    formula = y ~ x
  ) +
  stat_cor(
    label.x = 3,
    label.y = 3,
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))
  ) +
  theme_grey() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  file.path(processed_dir, "linear_regression_postComBat.png"),
  plot = calibration_plot,
  dpi = 300
)

LFQ_ComBat_long <- LFQ_ComBat_long %>%
  mutate(Cexp.nM = 10^((logI - intercept) / slope) * 1000)

# ---- Annotate contaminants and remove variable immunoglobulins -------------

LFQ_ComBat_long <- LFQ_ComBat_long %>%
  separate_rows(Genes, sep = ";") %>%
  left_join(
    sampling %>% select(Genes, Source),
    by = "Genes"
  ) %>%
  group_by(across(-c(Genes, Source))) %>%
  summarise(
    Genes = paste(unique(Genes), collapse = ";"),
    Source = paste(na.omit(unique(Source)), collapse = ";"),
    .groups = "drop"
  ) %>%
  select(Genes, Source, everything()) %>%
  filter(
    !str_detect(Genes, "IGLV"),
    !str_detect(Genes, "IGKV"),
    !str_detect(Genes, "IGHV")
  )

LFQ_ComBat_noCont_long <- LFQ_ComBat_long %>%
  filter(
    !str_detect(
      Source,
      regex("keratin|erythrocyte|platelet|coagulation", ignore_case = TRUE)
    )
  )

# ---- Build wide LFQ tables -------------------------------------------------

run_metadata <- LFQ_ComBat_long %>%
  select(Run, Prob.Nr., Sex, Age, Group, Classification) %>%
  distinct()

LFQ_ComBat_noCont_wide <- LFQ_ComBat_noCont_long %>%
  select(Genes, LFQ, Run) %>%
  pivot_wider(names_from = Genes, values_from = LFQ) %>%
  left_join(run_metadata, by = "Run") %>%
  mutate(
    Subject = Prob.Nr.,
    Donor = paste(Subject, Classification, Sex, Age, Group, sep = "|")
  ) %>%
  select(
    Run, Donor, Classification, Subject, Sex, Age, Group,
    everything()
  ) %>%
  select(
    all_of(c("Run", "Donor", "Classification", "Subject", "Sex", "Age", "Group")),
    sort(setdiff(names(.), c(
      "Run", "Donor", "Classification", "Subject", "Sex", "Age", "Group"
    )))
  )

LFQ_ComBat_Cont_wide <- LFQ_ComBat_long %>%
  select(Genes, LFQ, Run) %>%
  pivot_wider(names_from = Genes, values_from = LFQ) %>%
  left_join(run_metadata, by = "Run") %>%
  mutate(Subject = Prob.Nr.) %>%
  select(
    Run, Subject, Sex, Age, Group, Classification,
    everything()
  )

# ---- Build concentration tables --------------------------------------------

Conc_ComBat_noCont_wide <- LFQ_ComBat_noCont_long %>%
  select(Genes, Cexp.nM, Run) %>%
  pivot_wider(names_from = Genes, values_from = Cexp.nM) %>%
  left_join(run_metadata, by = "Run") %>%
  mutate(
    Subject = Prob.Nr.,
    Donor = paste(Subject, Classification, Sex, Age, Group, sep = "|")
  ) %>%
  select(
    Run, Donor, Classification, Subject, Sex, Age, Group,
    everything()
  ) %>%
  select(
    all_of(c("Run", "Donor", "Classification", "Subject", "Sex", "Age", "Group")),
    sort(setdiff(names(.), c(
      "Run", "Donor", "Classification", "Subject", "Sex", "Age", "Group"
    )))
  )

# The original script referenced LFQ_ComBat_Contaminant_long, which was never
# created. Here the contaminant-containing concentration table is derived
# from the complete calibrated LFQ table.
Conc_ComBat_Cont_wide <- LFQ_ComBat_long %>%
  select(Genes, Cexp.nM, Run) %>%
  pivot_wider(names_from = Genes, values_from = Cexp.nM) %>%
  left_join(run_metadata, by = "Run") %>%
  mutate(Subject = Prob.Nr.) %>%
  select(
    Run, Subject, Sex, Age, Group, Classification,
    everything()
  )

# ---- Export processed cohort data ------------------------------------------

LFQ_ComBat_noCont_wide_export <- LFQ_ComBat_noCont_wide %>%
  filter(!Group %in% c("AMSBIO", "QC"))

Conc_ComBat_noCont_wide_export <- Conc_ComBat_noCont_wide %>%
  filter(!Group %in% c("AMSBIO", "QC"))

LFQ_ComBat_noCont_wide_AMSBIO <- LFQ_ComBat_noCont_wide %>%
  filter(Group == "AMSBIO")

LFQ_ComBat_Cont_wide_woAMSBIO <- LFQ_ComBat_Cont_wide %>%
  filter(!Group %in% c("AMSBIO", "QC")) %>%
  select(where(~ !all(is.na(.x))))

LFQ_ComBat_Cont_wide_AMSBIO <- LFQ_ComBat_Cont_wide %>%
  filter(Group %in% c("AMSBIO", "QC"))

write_csv(
  LFQ_ComBat_noCont_wide_export,
  file.path(processed_dir, "LFQ_DENV.csv")
)

write_csv(
  Conc_ComBat_noCont_wide_export,
  file.path(processed_dir, "Conc_DENV.csv")
)

write_csv(
  LFQ_ComBat_noCont_wide_AMSBIO,
  file.path(processed_dir, "LFQ_AMSBIO.csv")
)

write_csv(
  LFQ_ComBat_Cont_wide_woAMSBIO,
  file.path(processed_dir, "LFQ_DENV_wContaminant.csv")
)

write_csv(
  LFQ_ComBat_Cont_wide_AMSBIO,
  file.path(processed_dir, "LFQ_AMSBIO_wContaminant.csv")
)

# ---- Create Supplementary Data 1 -------------------------------------------

DENV_df <- LFQ_ComBat_noCont_wide_export %>%
  pivot_longer(
    cols = -(1:7),
    names_to = "Genes",
    values_to = "LFQ"
  ) %>%
  arrange(Genes) %>%
  filter(!str_detect(Genes, fixed(".")))

DENV_df <- DENV_df %>%
  mutate(
    Classification = case_when(
      Classification %in% c("DF", "DSS", "DHF") ~ "hospitalized",
      Classification %in% c("ASD", "VMD", "MD") ~ "subclinical",
      Classification == "HD" ~ "healthy",
      TRUE ~ Classification
    )
  )

genes_to_keep <- DENV_df %>%
  group_by(Genes, Group) %>%
  summarise(
    n_missing = sum(is.na(LFQ)),
    n_subjects_present = n_distinct(Subject),
    prop_missing = n_missing / n_subjects_present,
    .groups = "drop"
  ) %>%
  select(Genes, Group, prop_missing) %>%
  pivot_wider(names_from = Group, values_from = prop_missing) %>%
  filter(
    !if_all(
      any_of(c("AP", "CP", "ERP", "RP", "HD")),
      ~ .x > 0.7
    )
  ) %>%
  pull(Genes)

supplementary_data_s1 <- DENV_df %>%
  filter(Genes %in% genes_to_keep) %>%
  pivot_wider(names_from = Genes, values_from = LFQ)

age_breaks <- seq(0, 70, 5)

age_labels <- paste0(
  head(age_breaks, -1),
  "-",
  head(age_breaks + 4, -1)
)

supplementary_data_s1 <- supplementary_data_s1 %>%
  mutate(
    Age_range = cut(
      Age,
      breaks = age_breaks,
      labels = age_labels,
      right = FALSE,
      include.lowest = TRUE
    )
  ) %>%
  mutate(
    across(
      where(is.character),
      identity
    )
  )

gene_columns <- setdiff(
  names(supplementary_data_s1),
  c(
    "Run", "Donor", "Classification", "Subject",
    "Sex", "Age", "Group", "Age_range"
  )
)

supplementary_data_s1 <- supplementary_data_s1 %>%
  mutate(across(all_of(gene_columns), as.numeric))

write_csv(
  supplementary_data_s1,
  file.path(processed_dir, "Supplementary_Data_S1.csv")
)
