# DENV Plasma Proteomics: Statistical Analysis and Visualization
# Downstream analysis of processed LFQ data


library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggsignif)
library(readxl)
library(dplyr)
library(png)
library(arrow)
library(ggplot2)
library(reshape2)
library(psych)
library(sva)
library(FSA)
library(purrr)
library(vegan)
library(ggrepel)
library(multcomp)
library(ggnewscale)
library(yarrr)
library(egg)

project_dir <- "PATH/TO/DENV_PROJECT"

out_folder <- file.path(project_dir, "data", "processed_data")
source_data <- file.path(project_dir, "data", "source_data")
figure_folder <- file.path(project_dir, "figures")

dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load processed LFQ data and metadata
# ============================================================

DENV_df <- read.csv(
  file.path(out_folder, "Supplementary_Data_S1.csv"),
  check.names = FALSE
)

metadata <- read_excel(file.path(source_data, "cohort_metadata.xlsx"))

DENV_df <- DENV_df[, -1] %>%
  select(Run, Donor, Classification, Subject, Sex, Age, Group, everything()) %>%
  pivot_longer(cols = -(1:7), names_to = "Genes", values_to = "LFQ") %>%
  arrange(Genes)

DENV_df$Classification_02 <- metadata$Classification[
  match(DENV_df$Subject, metadata$Prob.Nr.)
]

DENV_df <- DENV_df %>%
  mutate(
    Classification_02 = case_when(
      Classification_02 %in% c("DSS", "DHF") ~ "DHF/DSS",
      TRUE ~ Classification_02
    )
  )

# ============================================================
# Log2 transformation and missing-value imputation
# ============================================================

DENV_imp <- DENV_df %>%
  mutate(LFQ = log2(LFQ))

set.seed(123)

DENV_imp <- DENV_imp %>%
  complete(Genes, Subject, Group, Donor, fill = list(LFQ = NA)) %>%
  filter(!is.na(Run)) %>%
  group_by(Genes) %>%
  mutate(
    min_val = min(LFQ[LFQ > 0], na.rm = TRUE),
    LFQ = if_else(
      is.na(LFQ),
      runif(n(), min = 0.5 * min_val, max = min_val),
      LFQ
    )
  ) %>%
  ungroup() %>%
  select(-min_val)

# ============================================================
# AP versus RP analysis
# ============================================================

DENV_imp_phase <- DENV_imp %>%
  filter(Group != "HD") %>%
  select(Subject, Group, Classification, Donor, Genes, LFQ) %>%
  mutate(
    Group = factor(Group, levels = c("AP", "CP", "ERP", "RP")),
    Classification = factor(Classification)
  )

fit_phase_model <- function(subdf) {
  gene <- unique(subdf$Genes)

  model <- lmer(
    LFQ ~ Group * Classification + (1 | Subject),
    data = subdf,
    REML = FALSE
  )

  emmeans_results <- emmeans(model, ~ Group | Classification)

  contrasts <- contrast(
    emmeans_results,
    method = "pairwise"
  ) %>%
    summary(infer = TRUE)

  contrasts$Genes <- gene

  list(gene = gene, contrasts = contrasts, model = model)
}

phase_results <- DENV_imp_phase %>%
  group_split(Genes) %>%
  map(fit_phase_model)

phase_contrasts <- bind_rows(map(phase_results, "contrasts")) %>%
  group_by(contrast) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

sig_AP_RP <- phase_contrasts %>%
  filter(contrast == "AP - RP", padj < 0.05)

dir.create(file.path(out_folder, "Final_tables"), recursive = TRUE, showWarnings = FALSE)

write.csv(
  sig_AP_RP,
  file.path(out_folder, "Final_tables", "significant_AP_RP_contrasts.csv"),
  row.names = FALSE
)

# ============================================================
# AP versus RP visualization
# ============================================================

genes_to_plot <- c(
  "VCAM1", "CTSD", "FUCA1", "HSP90B1", "SERPINA3", "FGL1",
  "SAA1", "ORM1", "SERPINA10", "GOLM1", "LRG1", "B2M",
  "POSTN", "APOH", "ITIH1", "C2", "CRP", "VWF"
)

sig_labels_AP_RP <- sig_AP_RP %>%
  filter(Genes %in% genes_to_plot) %>%
  mutate(
    annotation = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  group_by(Genes) %>%
  slice(1) %>%
  ungroup()

plot_df_AP_RP <- DENV_imp %>%
  filter(Group %in% c("AP", "RP"), Genes %in% genes_to_plot) %>%
  group_by(Genes) %>%
  mutate(
    zscore = (LFQ - mean(LFQ, na.rm = TRUE)) / sd(LFQ, na.rm = TRUE)
  ) %>%
  ungroup()

gene_order_AP_RP <- plot_df_AP_RP %>%
  group_by(Genes, Group) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = mean_z) %>%
  mutate(diff = AP - RP) %>%
  arrange(desc(diff))

plot_df_AP_RP <- plot_df_AP_RP %>%
  mutate(Genes = factor(Genes, levels = gene_order_AP_RP$Genes))

y_pos_AP_RP <- plot_df_AP_RP %>%
  group_by(Genes) %>%
  summarise(y_position = max(zscore, na.rm = TRUE) + 0.6, .groups = "drop")

signif_df_AP_RP <- sig_labels_AP_RP %>%
  left_join(y_pos_AP_RP, by = "Genes")

plot_AP_RP <- ggplot(
  plot_df_AP_RP,
  aes(x = Group, y = zscore, fill = Group)
) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  facet_wrap(~ Genes, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c("AP" = "#531253", "RP" = "#9cf6f6")) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  labs(x = NULL, y = "Standardized abundance (z-score)") +
  theme_light() +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(angle = 90, vjust = 0.5, size = 7),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    legend.key.size = unit(3, "mm")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

if (nrow(signif_df_AP_RP) > 0) {
  plot_AP_RP <- plot_AP_RP +
    geom_signif(
      data = signif_df_AP_RP,
      aes(
        xmin = "RP",
        xmax = "AP",
        annotations = annotation,
        y_position = y_position
      ),
      manual = TRUE,
      inherit.aes = FALSE,
      linewidth = 0.3,
      textsize = 2.5,
      tip_length = 0.01,
      vjust = 0.3
    )
}

print(plot_AP_RP)

ggsave(
  file.path(figure_folder, "violinplot_AP_RP.pdf"),
  plot_AP_RP,
  width = 19,
  height = 11,
  units = "cm"
)

# ============================================================
# Acute-phase DF versus DHF/DSS analysis
# ============================================================

DENV_imp_AP <- DENV_imp %>%
  filter(Group == "AP", Classification_02 %in% c("DF", "DHF/DSS")) %>%
  select(Subject, Group, Classification_02, Donor, Genes, LFQ) %>%
  mutate(
    Classification_02 = factor(
      Classification_02,
      levels = c("DF", "DHF/DSS")
    )
  )

fit_AP_model <- function(subdf) {
  gene <- unique(subdf$Genes)

  model <- lm(LFQ ~ Classification_02, data = subdf)

  emmeans_results <- emmeans(model, ~ Classification_02)

  contrasts <- contrast(
    emmeans_results,
    method = "pairwise"
  ) %>%
    summary(infer = TRUE)

  contrasts$Genes <- gene

  list(gene = gene, contrasts = contrasts, model = model)
}

AP_results <- DENV_imp_AP %>%
  group_split(Genes) %>%
  map(fit_AP_model)

AP_contrasts <- bind_rows(map(AP_results, "contrasts")) %>%
  group_by(contrast) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

sig_AP <- AP_contrasts %>%
  filter(p.value < 0.01) %>%
  pull(Genes)

# ============================================================
# Acute-phase DF versus DHF/DSS visualization
# ============================================================

df_ap <- DENV_imp %>%
  filter(
    Group == "AP",
    Genes %in% sig_AP,
    Classification_02 %in% c("DF", "DHF/DSS")
  ) %>%
  group_by(Genes) %>%
  mutate(
    zscore = (LFQ - mean(LFQ, na.rm = TRUE)) / sd(LFQ, na.rm = TRUE)
  ) %>%
  ungroup()

fc_df <- df_ap %>%
  group_by(Genes, Classification_02) %>%
  summarise(mean_LFQ = mean(LFQ, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Classification_02, values_from = mean_LFQ) %>%
  mutate(
    log2FC = `DHF/DSS` - DF,
    absFC = abs(log2FC)
  )

df_plot_DF_DHF <- df_ap %>%
  left_join(fc_df %>% select(Genes, log2FC, absFC), by = "Genes") %>%
  mutate(
    Genes = factor(Genes, levels = fc_df %>% arrange(absFC) %>% pull(Genes))
  )

plot_DF_DHF <- ggplot(
  df_plot_DF_DHF,
  aes(x = Genes, y = zscore, fill = Classification_02)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.9,
    width = 0.7,
    color = "black",
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    aes(color = Classification_02),
    size = 2,
    alpha = 0.6,
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)
  ) +
  scale_fill_manual(values = c("DHF/DSS" = "#DB3A34", "DF" = "#2274A5")) +
  scale_color_manual(values = c("DHF/DSS" = "#DB3A34", "DF" = "#2274A5")) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  theme_light(base_size = 13) +
  labs(
    x = NULL,
    y = "z-score",
    title = "Differentiated proteins between DF and DHF/DSS",
    subtitle = "p-value < 0.01 | ranked on increasing fold change"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank()
  )

print(plot_DF_DHF)

ggsave(
  file.path(figure_folder, "boxplot_AP_lmr_DF_DHF.pdf"),
  plot_DF_DHF,
  width = 15,
  height = 7,
  dpi = 300
)

sig_df_dhf <- AP_contrasts %>%
  filter(p.value < 0.01)

write.csv(
  sig_df_dhf,
  file.path(
    out_folder,
    "Final_tables",
    "significant_DF_DHF_DSS_contrasts.csv"
  ),
  row.names = FALSE
)




#------------------- subclinical vs hospitalized lmm model ---------------------

DENV_imp <- DENV_imp %>% filter(Group != c("HD"))

DENV_imp <- DENV_imp %>%
  dplyr::select(Subject,Group,Classification,Donor,Genes,LFQ) %>%
  mutate(
    Group = factor(Group, levels = c("AP","CP","ERP","RP")),
    Classification = factor(Classification)
  )

group_colors2 <- c(
  "AP" = "#531253", 
  "CP" = "#DAA49A",   
  "ERP" = "#5588bb",  
  "RP" = "#9cf6f6",
  "HD" = "#66A61E"  
)

############## Fit model in each phase separately ##############

DENV_imp_AP <- DENV_imp %>% filter(Group == c("AP"))
DENV_imp_AP <- DENV_imp_AP %>%
  dplyr::select(Subject,Group,Classification,Donor,Genes,LFQ) %>%
  mutate(
    Classification = factor(Classification)
  )
set.seed(123)
fit_one <- function(subdf) {
  gene <- unique(subdf$Genes)
  m <- lm(LFQ ~ Classification, data = subdf)
  a <- anova(m)
  em <- emmeans(m, ~  Classification ) 
  contr <- contrast(em, method = "pairwise") %>% summary(infer = TRUE)
  contr$Genes <- gene
  list(gene = gene, a = a, contrasts = contr, model = m)
}

# Run for all genes (use group_split to preserve Genes column)
res_list_AP <- DENV_imp_AP %>% group_split(Genes) %>% map(fit_one)

# Bind all contrasts into a dataframe
contrasts_df_AP <- bind_rows(map(res_list_AP, "contrasts"))

# Adjust p-values per contrast (or globally)
contrasts_df_AP <- contrasts_df_AP %>%
  group_by(contrast) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

sig_AP <- contrasts_df_AP %>%
  filter(padj < 0.05) %>%
  pull(Genes)

sig_labels_AP <- contrasts_df_AP %>%
  filter(Genes %in% sig_AP) %>%
  mutate(
    annotation = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "ns"
    )
  ) %>%
  group_by(Genes) %>%
  slice(1) %>%      # only one contrast per gene
  ungroup()

df_AP <- DENV_imp %>%
  filter(Group == "AP", Genes %in% sig_AP) %>%
  mutate(
  ) %>%
  group_by(Genes) %>%
  mutate(
    zscore = (LFQ - mean(LFQ, na.rm = TRUE)) / sd(LFQ, na.rm = TRUE)
  ) %>%
  ungroup()

gene_order <- df_AP %>%
  group_by(Genes, Classification) %>%
  dplyr::summarize(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Classification, values_from = mean_z) %>%
  mutate(diff = hospitalized - subclinical) %>%
  arrange(desc(diff))

df_AP <- df_AP %>%
  mutate(Genes = factor(Genes, levels = gene_order$Genes))

y_pos_AP <- df_AP %>%
  group_by(Genes) %>%
  summarise(y_position = max(zscore, na.rm = TRUE) + 0.6)

signif_df_AP <- left_join(sig_labels_AP, y_pos_AP, by = "Genes")

# violin plot ###
plot_AP <- ggplot(df_AP, aes(x = Classification, y = zscore, fill = Classification)) +
  geom_violin(trim = FALSE,alpha = 0.9) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.9) +
  facet_wrap(~ Genes, scales = "free_y", ncol = 5) +
  scale_fill_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
  theme_light() +
  labs(
    x = "Classification",
    y = "Standardized abundance (z-score)",
    title = ""
  ) +
  theme(
    text            = element_text(size = 6),
    axis.text.y     = element_text(size = 6,angle = 0),
    axis.title.y    = element_text(size = 7,vjust = 0.5,angle = 90),
    strip.text      = element_text(size = 7, face = "bold"),
    legend.text     = element_blank(),
    axis.text.x     = element_blank(),
    axis.title.x    = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    legend.title     = element_blank()
  ) +
  geom_signif(
    data = signif_df_AP,
    aes(
      xmin = "subclinical",
      xmax = "hospitalized",
      annotations = annotation,
      y_position = y_position
    ),
    manual = TRUE,         
    inherit.aes = FALSE,
    size = 0.25,
    textsize = 2,
    tip_length = 0.01,
    vjust = 0.3
  )+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.13)))

print(plot_AP)

ggsave(file.path(figure_folder, "violinplot_AP_lm_zscore.pdf"),
  plot_AP, 
  width = 11, 
  height = 11,
  units = "cm",
  dpi = 300)



########################### Fit model in critical phase #################################

DENV_imp_CP <- DENV_imp %>% filter(Group == c("CP"))
DENV_imp_CP <- DENV_imp_CP %>%
  dplyr::select(Subject,Group,Classification,Donor,Genes,LFQ) %>%
  mutate(
    Classification = factor(Classification)
  )
set.seed(123)
fit_one <- function(subdf) {
  gene <- unique(subdf$Genes)
  m <- lm(LFQ ~ Classification, data = subdf)
  a <- anova(m)
  em <- emmeans(m, ~  Classification ) 
  contr <- contrast(em, method = "pairwise") %>% summary(infer = TRUE)
  contr$Genes <- gene
  list(gene = gene, a = a, contrasts = contr, model = m)
}

res_list_CP <- DENV_imp_CP %>% group_split(Genes) %>% map(fit_one)

contrasts_df_CP <- bind_rows(map(res_list_CP, "contrasts"))

contrasts_df_CP <- contrasts_df_CP %>%
  group_by(contrast) %>%
  mutate(padj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

sig_CP <- contrasts_df_CP %>%
  filter(padj < 0.05) %>%
  pull(Genes)

sig_labels_CP <- contrasts_df_CP %>%
  filter(Genes %in% sig_CP) %>%
  mutate(
    annotation = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "ns"
    )
  ) %>%
  group_by(Genes) %>%
  slice(1) %>%     
  ungroup()


df_CP <- DENV_imp %>%
  filter(Group == "CP", Genes %in% sig_CP) %>%
  mutate(
  ) %>%
  group_by(Genes) %>%
  mutate(
    zscore = (LFQ - mean(LFQ, na.rm = TRUE)) / sd(LFQ, na.rm = TRUE)
  ) %>%
  ungroup()

gene_diff <- df_CP %>%
  group_by(Genes, Classification) %>%
  summarise(mean_z = mean(zscore, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Classification, values_from = mean_z) %>%
  mutate(diff = hospitalized - subclinical)

df_CP <- df_CP %>%
  left_join(gene_diff %>% dplyr::select(Genes, diff), by = "Genes") %>%
  mutate(
    Genes = forcats::fct_reorder(Genes, diff, .desc = TRUE)
  )

y_pos_CP <- df_CP %>%
  group_by(Genes) %>%
  summarise(y_position = max(zscore, na.rm = TRUE) + 0.6)

signif_df_CP <- left_join(sig_labels_CP, y_pos_CP, by = "Genes")

# violin plot
plot_CP <- ggplot(df_CP, aes(x = Classification, y = zscore, fill = Classification)) +
  geom_violin(trim = FALSE,alpha = 0.9) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 0.9) +
  facet_wrap(~ Genes,scales = "free_y",ncol=5) +
  scale_fill_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", size = 0.4) +
  theme_light() +
  labs(
    x = "Classification",
    y = "Standardized abundance (z-score)",
    title = ""
  ) +
  theme(
    text            = element_text(size = 6),
    axis.text.y     = element_text(size = 6,angle = 0),
    axis.title.y    = element_text(size = 7,vjust = 0.5,angle = 90),
    strip.text      = element_text(size = 7, face = "bold"),
    legend.text     = element_blank(),
    axis.text.x     = element_blank(),
    axis.title.x    = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position  = "none",
    legend.title     = element_blank()
  )+
  geom_signif(
    data = signif_df_CP,
    aes(
      xmin = "subclinical",
      xmax = "hospitalized",
      annotations = annotation,
      y_position = y_position
    ),
    manual = TRUE,         
    inherit.aes = FALSE,
    size = 0.25,
    textsize = 2,
    tip_length = 0.01,
    vjust = 0.3
  )+ scale_y_continuous(expand = expansion(mult = c(0.05, 0.13)))
print(plot_CP)

ggsave(file.path(figure_folder, "violinplot_CP_lm_zscore.png"), 
       plot_CP, 
       width = 11, 
       height = 10,
       units = "cm", 
       dpi = 300)


#------------------- FC plot of AP and CP ---------------------------------

df_all <- bind_rows(
  df_AP %>% mutate(Group = "AP"),
  df_CP %>% mutate(Group = "CP")
)

df_all <- df_all %>%
  mutate(LFQ = 2 ^ LFQ)

df_fc <- df_all %>%
  group_by(Group, Genes, Classification) %>%
  dplyr::summarize(
    mean_LFQ = mean(LFQ, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Classification,
    values_from = mean_LFQ
  ) %>%
  mutate(
    fold_change = hospitalized / subclinical,
    log2FC = log2(fold_change)
  )

df_plot <- df_fc %>%
  mutate(
    Group = factor(Group, levels = c("AP", "CP"))
  ) %>%
  group_by(Genes) %>%
  arrange(desc(Group), .by_group = TRUE) %>%
  mutate(
    total_log2FC = sum(log2FC, na.rm = TRUE),
    cumulative_log2FC = cumsum(log2FC)
  ) %>%
  ungroup()

fc_dotplot <- ggplot(
  df_plot,
  aes(
    x = reorder(Genes, total_log2FC),
    y = cumulative_log2FC,
    colour = Group
  )
) +
  
  geom_segment(
    aes(
      x = reorder(Genes, total_log2FC),
      xend = reorder(Genes, total_log2FC),
      y = 0,
      yend = total_log2FC
    ),
    inherit.aes = FALSE,
    colour = "grey75",
    linewidth = 0.6
  ) +
  geom_point(
    size = 2.2,
    stroke = 0
  ) +
  
  coord_flip() +
  
  scale_y_continuous(
    limits = c(-3.5, 3.5)
  ) +
  
  scale_colour_manual(
    values = c(
      "AP" = "#531253",
      "CP" = "#DAA49A"
    )
  ) +
  
  theme_light() +
  
  labs(
    title = "Fold change in p < 0.05 Genes",
    x = "Gene names",
    y = "Differential abundance [log2 FC]"
  ) +
  
  theme(
    text = element_text(family = "sans"),
    strip.text = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7, hjust = 1),
    axis.title.x = element_text(size = 8),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title.y = element_text(size = 8),
    plot.background = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    legend.title = element_blank(),
    legend.key.size = unit(3, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(3, "mm"),
    title = element_text(size = 7)
  )

print(fc_dotplot)

ggsave(file.path(figure_folder,"FC_AP_CP_dotplot.pdf"),
  fc_dotplot,
  width = 9,
  height = 22,
  units = "cm",
  dpi = 300
)

############### Longitudinal dot plots of differentiated proteins #################

AP <- c("CRP",'B2M',"C2",'VCAN','C1QA','CP','AAT','C1QB','C1QC',"SAA1","SAA2","LRG1",
        'HP','S100A8','S100A9','FGA','FGB',"FGG",'A1AT','AACT','LBP',
        'CD14',"A2M","C4A","C4B","C9",'ITIH3',"SERPINA3","SERPINA1",
        "SERPINA10","SERPING1","SERPINF2","ORM1","ORM2","CFP","TF",'TTR',
        'ITIH1','ITIH2','SERPINA4','SERPIND1',"APOA1",'APOA4','APOA2',"RBP4",
        "SELENOP","AHSG")

endoth <- c("APOE","CD14","LBP","ORM1","PTX3","PVR","SLC3A2","VCAM1","VCAN")
viral <- c("CKM","ENPP2","GOLM1","LCP1","LGALS3BP","PCSK9","WARS1","APOH","AZGP1")

# Prepare data
df_long_AP <- DENV_imp %>%
  filter(Genes %in% c(sig_AP, sig_CP))%>%
  filter(!grepl("^HD", Group)) %>%
  filter(Genes %in% AP) %>%
  mutate(
    LFQ = 2 ^ LFQ,                            
    Genes = factor(Genes, levels = AP)         
  ) %>%
  group_by(Genes, Group, Classification) %>%
  dplyr::summarize(
    mean_LFQ = mean(LFQ, na.rm = TRUE),
    sd_LFQ   = sd(LFQ, na.rm = TRUE),
    .groups = "drop"
  )

# Plot all genes together
plot_longit_AP <- ggplot(df_long_AP, aes(x = Group, y = mean_LFQ, group = Classification, color = Classification)) +
  geom_ribbon(aes(ymin = mean_LFQ - sd_LFQ, ymax = mean_LFQ + sd_LFQ, fill = Classification), alpha = 0.2, color = NA) +
  geom_line(size = 0.5) +
  geom_point(aes(fill = Classification), colour = "black", size = 1, pch = 21) +
  scale_color_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  scale_fill_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  scale_y_continuous(labels = scales::scientific) +
  facet_wrap(~ Genes, scales = "free_y", ncol = 5) +
  theme_light(base_size = 6) +
  labs(
    x = "",
    y = "Protein abundance [LFQ]",
    title = "proteins involved in inflammation",
    color = "Classification",
    fill = "Classification"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "none",
    legend.title = element_blank(),
    panel.background = element_blank(),
    axis.title   = element_text(size = 6),
    axis.text    = element_text(size = 6),
    axis.title.y = element_text(
      size = 6,
      margin = margin(r = 4)
    )
  )


df_long_endo <- DENV_imp %>%
  filter(Genes %in% c(sig_AP, sig_CP))%>%
  filter(!grepl("^HD", Group)) %>%
  filter(Genes %in% endoth) %>%
  mutate(
    LFQ = 2 ^ LFQ,                          
    Genes = factor(Genes, levels = endoth)         
  ) %>%
  group_by(Genes, Group, Classification) %>%
  dplyr::summarize(
    mean_LFQ = mean(LFQ, na.rm = TRUE),
    sd_LFQ   = sd(LFQ, na.rm = TRUE),
    .groups = "drop"
  )

# Plot all genes together
plot_longit_endo <- ggplot(df_long_endo, aes(x = Group, y = mean_LFQ, group = Classification, color = Classification)) +
  geom_ribbon(aes(ymin = mean_LFQ - sd_LFQ, ymax = mean_LFQ + sd_LFQ, fill = Classification), alpha = 0.2, color = NA) +
  geom_line(size = 0.5) +
  geom_point(aes(fill = Classification), colour = "black", size = 1, pch = 21) +
  scale_color_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  scale_fill_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  scale_y_continuous(labels = scales::scientific) +
  facet_wrap(~ Genes, scales = "free_y", ncol = 3) +
  theme_light(base_size = 6) +
  labs(
    x = "",
    y = "Protein abundance [LFQ]",
    title = "endothelial dysregulation",
    color = "Classification",
    fill = "Classification"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "none",
    legend.title = element_blank(),
    panel.background = element_blank(),
    axis.title   = element_text(size = 6),
    axis.text    = element_text(size = 6),
    axis.title.y = element_text(
      size = 6,
      margin = margin(r = 4)  
    ))

df_long_viral <- DENV_imp %>%
  filter(Genes %in% c(sig_AP, sig_CP))%>%
  filter(!grepl("^HD", Group)) %>%
  filter(Genes %in% viral) %>%
  mutate(
    LFQ = 2 ^ LFQ,                             
    Genes = factor(Genes, levels = viral)        
  ) %>%
  group_by(Genes, Group, Classification) %>%
  dplyr::summarize(
    mean_LFQ = mean(LFQ, na.rm = TRUE),
    sd_LFQ   = sd(LFQ, na.rm = TRUE),
    .groups = "drop"
  )

# Plot all genes together
plot_longit_viral <- ggplot(df_long_viral, aes(x = Group, y = mean_LFQ, group = Classification, color = Classification)) +
  geom_ribbon(aes(ymin = mean_LFQ - sd_LFQ, ymax = mean_LFQ + sd_LFQ, fill = Classification), alpha = 0.2, color = NA) +
  geom_line(size = 0.5) +
  geom_point(aes(fill = Classification), colour = "black", size = 1, pch = 21) +
  scale_color_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  scale_fill_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  scale_y_continuous(labels = scales::scientific) +
  facet_wrap(~ Genes, scales = "free_y", ncol = 3) +
  theme_light(base_size = 6) +
  labs(
    x = "",
    y = "Protein abundance [LFQ]",
    title = "viral-induced inflammation",
    color = "Classification",
    fill = "Classification"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    legend.position = "none",
    legend.title = element_blank(),
    panel.background = element_blank(),
    axis.title   = element_text(size = 6),
    axis.text    = element_text(size = 6),
    axis.title.y = element_text(
      size = 6,
      margin = margin(r = 4) 
    ))

plot_longit_AP_fixed <- set_panel_size(
  plot_longit_AP,
  width  = unit(2, "cm"),
  height = unit(1.6, "cm")
)

plot_longit_endo_fixed <- set_panel_size(
  plot_longit_endo,
  width  = unit(2.3, "cm"),
  height = unit(1.8, "cm")
)

plot_longit_viral_fixed <- set_panel_size(
  plot_longit_viral,
  width  = unit(2.3, "cm"),
  height = unit(1.8, "cm")
)

# Save plots
ggsave(file.path(figure_folder, "lineplot_ap_proteins.pdf"), plot_longit_AP_fixed, width = 15, height = 8,units="cm", dpi = 300)
ggsave(file.path(figure_folder, "lineplot_endo_proteins.pdf"), plot_longit_endo_fixed, width = 10, height = 9, units="cm",dpi = 300)
ggsave(file.path(figure_folder, "lineplot_viral_proteins.pdf"), plot_longit_viral_fixed, width = 10, height = 9, units="cm",dpi = 300)

# Save table
contrasts_df_AP$Group <- "AP"
contrasts_df_CP$Group <- "CP"
sig_AP_df <- contrasts_df_AP %>%
  filter(padj < 0.05)
sig_CP_df <- contrasts_df_CP %>%
  filter(padj < 0.05)
sig_AP_CP <- bind_rows(sig_AP_df, sig_CP_df)

write.csv(sig_AP_CP, "significant_AP_RP_contrasts.csv", row.names = FALSE)


################################ Immunoglobulin subclasses #####################################

# List of genes to plot
genes_to_plot <- c("IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM", "IGHA1", "IGHA2","IGHD")

# Prepare data
df_ig <- DENV_imp %>%
  filter(!grepl("^HD", Group)) %>%
  filter(Genes %in% genes_to_plot) %>%
  mutate(LFQ = 2 ^ LFQ) %>%
  group_by(Genes, Group, Classification) %>%
  mutate(Genes = factor(Genes, levels = genes_to_plot))%>%
  dplyr::summarize(
    mean_LFQ = mean(LFQ, na.rm = TRUE),
    sd_LFQ   = sd(LFQ, na.rm = TRUE),
    .groups = "drop"
  )

# Plot all genes together
plot_ig <- ggplot(df_ig, aes(x = Group, y = mean_LFQ, group = Classification, color = Classification)) +
  geom_ribbon(aes(ymin = mean_LFQ - sd_LFQ, ymax = mean_LFQ + sd_LFQ, fill = Classification), alpha = 0.2, color = NA) +
  geom_line(size = 0.5) +
  geom_point(aes(fill = Classification), colour = "black", size = 1, pch = 21) +
  scale_color_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  scale_fill_manual(values = c("hospitalized" = "#5A5895FF", "subclinical" = "#E5BA3AFF")) +
  scale_y_continuous(labels = scales::scientific) +
  facet_wrap(~ Genes, scales = "free_y", ncol = 4) +
  theme_light(base_size = 7) +
  labs(
    x = "",
    y = "Protein abundance [LFQ]",
    title = "immunoglobulin subclasses",
    subtitle = "longitudinal comparison in hospitalized and subclinical patients (mean ± sd)",
    color = "Classification",
    fill = "Classification"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, size = 6),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    axis.text.y = element_text(angle = 0, size = 6),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text =  element_text(angle = 0, size = 7),
    legend.key.size   = unit(3, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.key.width  = unit(3, "mm"),
    axis.title.y = element_text(
      size = 7,
      margin = margin(r = 4)
    )
  )

plot_longit_ig_fixed <- set_panel_size(
  plot_ig,
  width  = unit(3.2, "cm"),
  height = unit(2.6, "cm")
)

ggsave(file.path(figure_folder,"lineplot_IG_genes_faceted.pdf"), 
       plot_longit_ig_fixed, 
       width = 18, 
       height = 10, 
       units="cm",
       dpi = 300)



df_ig <- DENV_imp %>%
  mutate(LFQ = 2 ^ LFQ) %>%
  filter(Genes %in% genes_to_plot)

df_HD <- df_ig %>% filter(Group == "HD")

df_HD_expanded <- df_HD %>%
  rename(OriginalGroup = Group) %>%
  tidyr::crossing(Group = c("AP", "CP", "ERP", "RP")) %>%
  dplyr::select(-OriginalGroup)

df_all <- df_ig %>%
  filter(Group != "HD") %>%
  bind_rows(df_HD_expanded) %>%
  mutate(
    Genes = factor(Genes, levels = genes_to_plot),
    Classification = factor(Classification,
                            levels = c("healthy","subclinical","hospitalized"))
  )

pal <- c(
  "healthy"      = "#66A61E",
  "subclinical"  = "#E5BA3AFF",
  "hospitalized" = "#5A5895FF"
)

plot_gene_group <- function(df, gene, group) {
  
  df_sub <- df %>% filter(Group == group, Genes == gene)
  
  p <- ggbetweenstats(
    data              = df_sub,
    x                 = Classification,
    y                 = LFQ,
    type              = "parametric",
    pairwise.display  = "significant",
    p.adjust.method   = "BH",
    results.subtitle  = FALSE,
    messages          = FALSE,
    title             = paste0(gene),
    ylab              = "Protein abundance [LFQ]",
    ggtheme           = theme_light(),
    point.args        = list(alpha = 1,size = 1),
    pairwise.annotation.args = list(size = 1.8),
    ggsignif.args = list(textsize = 1.8, tip_length = 0.01),
    centrality.plotting = FALSE,
    ggplot.component = list(theme(title = element_text(size = 6)))
  )
  
  p + scale_fill_manual(values = pal) +
    scale_color_manual(values = pal) +
    scale_y_continuous(labels = scales::scientific) +
    theme(
      axis.text.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 7, face = "bold"),
      axis.text.y = element_text(size = 7),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6, angle = 90),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text  = element_text(size = 7)
    )
}

groups <- c("AP", "CP", "ERP", "RP")

for (g in groups) {
  
  message("Processing group: ", g)
  
  dir.create(file.path(output_dir, g), showWarnings = FALSE)
  
  plot_list <- list()
  
  for (gene in genes_to_plot) {
    
    p <- plot_gene_group(df_all, gene, g)
    plot_list[[gene]] <- p
    
    ggsave(
      file.path(output_dir, g, paste0(gene, "_", g, ".pdf")),
      p,
      width = 4, height = 3, units = "cm"
    )
  }
  
  
  combined <- combine_plots(
    plot_list,
    plotgrid.args = list(nrow = 2L),
    annotation.args = list(
      title = paste0("Immunoglobulin subclasses - ",g),
      caption = "comparison in hospitalized, subclinical and healthy individuals",
      theme = theme(legend.position="bottom",
                    plot.subtitle = element_text(size = 8),
                    plot.title = element_text(size = 8),
                    plot.caption = element_text(size = 7))
    ))
  ggsave(
    file.path(figure_folder, paste0("violinplots_Igs_", g,".tiff")),
    combined,
    width  = 21,   
    height = 17,  
    units  = "cm",
    dpi    = 300,
    compression = "lzw"
  )
}


# ------------------- NS1 plots -------------------------------

# non-imputed dataset
DENV_df <- DENV_df[!is.na(DENV_df$LFQ), ]

classification_order <- c("healthy","subclinical","hospitalized")

DENV_df <- DENV_df %>%
  mutate(Subject = gsub("^HD", "", Subject)) %>%
  mutate(Classification = factor(Classification, levels = classification_order)) %>%
  mutate(Group = factor(Group, levels = c("AP", "CP", "ERP", "RP", "HD")))

DENV_df <- DENV_df %>%
  mutate(
    is_HD = grepl("^0", Subject),  
    Subject_numeric = as.numeric(Subject)  
  )

ordered_subjects <- DENV_df %>%
  mutate(Classification = factor(Classification, levels = classification_order)) %>%
  arrange(Classification, Subject) %>%
  pull(Subject) %>%
  unique()

DENV_df <- DENV_df %>%
  mutate(Subject = factor(Subject, levels = ordered_subjects))

ns1_data <- DENV_df %>%
  filter(str_detect(Genes, "NS1")) %>%
  arrange(Classification, Subject)

subjects_with_ns1 <- ns1_data %>%
  filter(!is.na(LFQ)) %>%
  pull(Subject) %>%
  unique()

ns1_data <- ns1_data %>%
  filter(Subject %in% subjects_with_ns1) %>%
  mutate(Subject = droplevels(Subject))

ns1_data <- ns1_data %>%
  group_by(Classification) %>%
  mutate(
    Subject = factor(Subject, levels = unique(Subject)),
  ) %>%
  ungroup()


error_data <- ns1_data %>%
  group_by(Classification, Subject) %>%
  summarise(
    Mean_Cexp = mean(LFQ, na.rm = TRUE),
    SD_Cexp = sd(LFQ, na.rm = TRUE),
    .groups = "drop"
  )

main_plot <- ggplot(ns1_data, aes(x = Subject, y = LFQ, color = Group)) +
  geom_point(aes(fill = Group), shape = 21, size = 6,
             color = "black", stroke = 0.4) +
  geom_errorbar(
    data = error_data,
    aes(
      x = Subject,
      ymin = Mean_Cexp - SD_Cexp,
      ymax = Mean_Cexp + SD_Cexp
    ),
    inherit.aes = FALSE,
    width = 0.3,
    color = "gray57"
  )+
  scale_y_continuous(labels = scales::scientific)+
  ggh4x::facet_grid2(Genes ~ Classification, scales = "free", independent = "x") +
  ggh4x::force_panelsizes(cols = c(0.3,2.3,2))+
  labs(
    title = "NS1",
    y = "LFQ",
    color = "Phase"
  ) +
  theme_light() +
  theme( text = element_text(family = "sans"),
         plot.title = element_text(size = 22),
         plot.subtitle = element_text(size = 19),
         legend.title = element_text(size = 18),
         legend.text = element_text(size = 18),
         panel.grid.major = element_line(color = "grey90", size = 0.2),
         panel.grid.minor = element_blank(),
         #panel.border = element_rect(color = "black", fill = NA, size = 0.8),
         axis.text.x = element_text(size = 15),
         axis.ticks.x = element_line(linewidth = 0.4),
         axis.text.y = element_text(size = 15),
         axis.title.y = element_text(size = 18, face = "bold"),
         axis.title.x = element_blank(),
         plot.background = element_rect(fill = "white", color = NA),
         strip.text.y.right = element_blank(),
         strip.text = element_text(size = 20, face = "bold")
  ) +
  scale_fill_manual(values = group_colors2)


vline_data <- ns1_data %>%
  group_by(Genes, Classification) %>%
  summarise(
    n_subjects = n_distinct(Subject),
    .groups = "drop"
  ) %>%
  mutate(
    xintercept = map(n_subjects, ~ seq(1.5, .x - 0.5, by = 1))
  ) %>%
  tidyr::unnest(xintercept)

annotation_plot <- ggplot(ns1_data) +
  ggnewscale::new_scale_color() +
  geom_text(
    aes(x = Subject, y = 0, label = ifelse(Sex == "female", "F", "M")),
    size = 6,
    color = "gray57",
    vjust = 0.5
  ) +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = c(0),
    labels = c("Sex")
  )+
  ggh4x::facet_grid2(Genes ~ Classification, scales = "free", independent = "x") +
  ggh4x::force_panelsizes(cols = c(0.3,2.3))+
  geom_vline(
    data = vline_data,
    aes(xintercept = xintercept),
    color = "gray90",
    linewidth = 0.5,
    inherit.aes = FALSE
  )+
  labs(x = "Subject", y = NULL) +
  theme_light() +
  theme(
    text = element_text(family = "sans"),
    strip.text = element_blank(),
    axis.text.y = element_text(size = 14, hjust = 0.7),
    axis.ticks.x = element_line(linewidth = 0.4),
    axis.ticks.y = element_line(linewidth = 0.4),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray90"),
    axis.text.x = element_text(size = 14),
    plot.background = element_rect(fill = "white", color = NA),
    axis.title.x = element_text(size = 18, face = "bold")
  )

final_plot <- (
  main_plot + theme(legend.position = "top")
) /
  (annotation_plot + theme(legend.position = "none"))+ 
  plot_layout(heights = c(2.5, 0.3)) 

print(final_plot)
ggsave(
  file.path(figure_folder, "NS1_singleplot.pdf"),
  plot = final_plot,
  width = 22,
  height = 9,
  device = cairo_pdf
)
