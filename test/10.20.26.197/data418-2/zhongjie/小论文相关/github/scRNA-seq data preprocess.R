# Combine all data together
# Load data----
library(Seurat)
library(qs)
library(dplyr)
library(tidyr)
setwd("/data/scRNA/")
MTG<-readRDS(file = "MTG.rds")
DLPFC<-readRDS(file = "DLPFC.rds")

# Preprocess MTG DLPFC----
MTGsub<-subset(MTG,subset = ADNC%in%c("Not AD","Low"))
DLPFCsub<-subset(DLPFC,subset = ADNC%in%c("Not AD","Low"))
# Remove donors from MTG that are not present in DLPFC
MTGsub@meta.data$donor_id <- as.character(MTGsub@meta.data$donor_id)
cells_to_keep <- colnames(MTGsub)[MTGsub@meta.data$donor_id != "H20.33.043"]
MTGsub <- subset(MTGsub, cells = cells_to_keep)
## Get geneset----
library(msigdbr)
kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG")
kegg<- split(kegg$gene_symbol, kegg$gs_name)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark<- split(hallmark$gene_symbol, hallmark$gs_name)
#symbol
fatty_genes_symbol<- unique(hallmark$HALLMARK_FATTY_ACID_METABOLISM,
                            kegg$KEGG_FATTY_ACID_METABOLISM)
fatty_genes_set_symbol<-list(fatty_acid=fatty_genes_symbol)
#ensembl
library(biomaRt)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_conversion <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = fatty_genes_symbol,  # Your Gene Symbol list
  mart = mart
)
convert_dict <- setNames(gene_conversion$ensembl_gene_id, gene_conversion$hgnc_symbol)
fatty_genes_ensembl <- unname(convert_dict[fatty_genes_symbol])
fatty_genes_set_ensembl <- list(fatty_acid = fatty_genes_ensembl)
#Standardize metadata names----

DLPFCsub$Age_death<-DLPFCsub$`Age at death`
MTGsub$Age_death<-MTGsub$`Age at death`

#Calculate scores----
MTGsub<-AddModuleScore(MTGsub,features = fatty_genes_set_ensembl,name = "FAscore")
DLPFCsub<-AddModuleScore(DLPFCsub,features = fatty_genes_set_ensembl,name = "FAscore")

#Compare scores-----
seurat_objects <- list(
  "MTG" = MTGsub,
  "DLPFC" = DLPFCsub
)

# Define color scheme
color_schemes <- list(
  "MTG" = c("#d3d3d3", "#54278f"),   # From light gray to deep purple
  "DLPFC" = c("#d3d3d3", "#d62728")    # From light gray to deep red
)

# Load necessary packages
library(dplyr)
library(ggplot2)
library(Seurat)
library(patchwork)
library(scales)  # For pretty_breaks()

# Load necessary packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# Assume seurat_objects already exists in workspace
plots <- list()
all_data_combined <- list()

for (name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[name]]
  
  # 1. Get data and convert Age_death, sex, ADNC to factors
  data_to_plot <- FetchData(seurat_obj, vars = c("FAscore1", "donor_id", "sex", "Age_death", "ADNC")) %>%
    mutate(across(c(Age_death, sex, ADNC), as.factor))
  
  # 2. Calculate mean FAscore for each combination (sex, Age_death, ADNC)
  data_combined_filtered <- data_to_plot %>%
    group_by(sex, Age_death, ADNC) %>%
    summarise(mean_FAscore = mean(FAscore1, na.rm = TRUE), .groups = "drop") %>%
    group_by(sex, Age_death) %>%
    filter(n_distinct(ADNC) == 2) %>%
    mutate(Group = paste0(sex, "_", Age_death))
  
  # 3. Calculate effect size
  effect_df <- data_combined_filtered %>%
    pivot_wider(names_from = ADNC, values_from = mean_FAscore) %>%
    mutate(effect = `Low` - `Not AD`) %>%
    dplyr::select(sex, Age_death, effect)
  
  # 4. Wilcoxon test and generate significance annotations
  sig_results <- data_to_plot %>%
    group_by(sex, Age_death) %>%
    filter(n_distinct(ADNC) == 2) %>% 
    summarise(
      p_value = wilcox.test(FAscore1[ADNC == "Not AD"],
                            FAscore1[ADNC == "Low"],
                            paired = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      adjp = p.adjust(p_value, method = "BH"),
      p_label = case_when(
        adjp < 0.001 ~ "***",
        adjp < 0.01  ~ "**",
        adjp < 0.05  ~ "*",
        TRUE         ~ "ns"
      ),
      Group = paste0(sex, "_", Age_death)
    )
  
  # 5. Merge effect values with significance results
  data_combined_filtered <- left_join(data_combined_filtered, effect_df, by = c("sex", "Age_death"))
  data_combined_filtered <- left_join(
    data_combined_filtered, 
    dplyr::select(ungroup(sig_results), Group, adjp, p_label),
    by = "Group"
  )
  
  # 6. Set colors based on significance and effect direction
  data_combined_filtered <- data_combined_filtered %>%
    mutate(sig_color = case_when(
      adjp >= 0.05 ~ "grey",
      adjp < 0.05 & effect > 0 & ADNC == "Low"    ~ "lightpink",
      adjp < 0.05 & effect > 0 & ADNC == "Not AD" ~ "lightblue",
      adjp < 0.05 & effect < 0 & ADNC == "Low"    ~ "lightblue",
      adjp < 0.05 & effect < 0 & ADNC == "Not AD" ~ "lightpink"
    ))
  
  # 7. Optional: Filter unwanted groups
  data_combined_filtered <- data_combined_filtered %>%
    filter(Group != "male_65 to 77 years old")
  
  # ✅ 8. Add Group order setting for "middle two rows swap"
  desired_order <- c(
    "female_78 to 89 years old",
    "male_78 to 89 years old", 
    "female_90+ years old",         
    "male_90+ years old"
  )
  # Only keep currently existing groups
  desired_order <- intersect(desired_order, unique(data_combined_filtered$Group))
  data_combined_filtered$Group <- factor(data_combined_filtered$Group, levels = rev(desired_order))
  
  # 9. Plot
  plot <- ggplot(data_combined_filtered, aes(x = ADNC, y = Group)) +
    geom_point(aes(size = mean_FAscore, fill = sig_color),
               shape = 21, color = "lightgray") +
    geom_text(aes(label = round(mean_FAscore, 3)), color = "black", size = 2.5) +
    geom_text(
      data = subset(data_combined_filtered, ADNC == "Low"),
      aes(label = p_label),
      color = "black", size = 3, nudge_x = 0.65
    ) +
    scale_size_continuous(name = "FAscore", range = c(5, 15)) +
    scale_fill_manual(name = "Significance",
                      values = c("grey" = "grey",
                                 "lightpink" = "#7a0177",
                                 "lightblue" = "#fff7bc")) +
    theme_minimal() +
    labs(x = "Pathology stage", y = "Group (Sex_Age)", title = name) +
    coord_cartesian(clip = "off") 
  
  plots[[name]] <- plot
  all_data_combined[[name]] <- data_combined_filtered
}

# Debug info: Check number of generated plots
print(paste("Number of generated plots:", length(plots)))

# Combine plots (2 columns)
wrap_plots(plots, ncol = 2)

#Draw sunburst plot-----
library(plotly)
library(dplyr)

df <- MTGsub@meta.data %>%
  distinct(donor_id, .keep_all = TRUE) %>%   # Ensure each donor appears only once
  # Convert relevant columns to character type
  mutate(
    ADNC = as.character(ADNC),
    Age_death  = as.character(Age_death),
    sex        = as.character(sex)
  ) %>%
  group_by(ADNC, Age_death, sex) %>%
  summarise(donor_count = n(), .groups = "drop") %>%
  mutate(
    sequence = paste(ADNC, Age_death, sex, sep = "-")
  )

# First level: ADNC
level1 <- df %>% 
  group_by(ADNC) %>% 
  summarise(value = sum(donor_count)) %>% 
  mutate(id = ADNC,
         label = ADNC,
         parent = "")

# Second level: Age_death (age groups under each ADNC)
level2 <- df %>% 
  group_by(ADNC, Age_death) %>% 
  summarise(value = sum(donor_count)) %>% 
  mutate(id = paste(ADNC, Age_death, sep = "-"),
         label = Age_death,
         parent = ADNC)

# Third level: sex (gender under each ADNC and Age_death)
level3 <- df %>% 
  group_by(ADNC, Age_death, sex) %>% 
  summarise(value = sum(donor_count)) %>% 
  mutate(id = paste(ADNC, Age_death, sex, sep = "-"),
         label = sex,
         parent = paste(ADNC, Age_death, sep = "-"))

# Merge data from three levels
sunburst_df <- bind_rows(level1, level2, level3)
sunburst_df <- sunburst_df %>%
  mutate(
    # Simple judgment: parent=="" -> level=1; otherwise if parent in {"1","2"} then level=2; otherwise level=3
    level = case_when(
      parent == "" ~ 1,
      parent %in% c("Not AD","Low") ~ 2,
      TRUE ~ 3
    )
  )

# Assign colors for each level's categories
# 1) level=1: Distinguish by ADNC
color_map_l1 <- c("Not AD" = "#317cb6", "Low" = "#b52330")

# 2) level=2: Distinguish by Age_death
color_map_age <- c("65 to 77 years old" = "#6eaed1", "78 to 89 years old" = "#dd6f58",
                   "90+ years old"="#d2d9dc","Less than 65 years old"="#ebebec")

# 3) level=3: Distinguish by sex
color_map_sex <- c("female" = "#f8b193", "male" = "#b7d8e9")

# Assign colors based on level
sunburst_df$color <- NA_character_

sunburst_df$color[sunburst_df$level == 1] <- 
  color_map_l1[sunburst_df$ADNC[sunburst_df$level == 1]]

sunburst_df$color[sunburst_df$level == 2] <- 
  color_map_age[sunburst_df$Age_death[sunburst_df$level == 2]]

sunburst_df$color[sunburst_df$level == 3] <- 
  color_map_sex[sunburst_df$sex[sunburst_df$level == 3]]

# Draw plot
fig <- plot_ly(
  data         = sunburst_df,
  ids          = ~id,
  labels       = ~label,
  parents      = ~parent,
  values       = ~value,
  type         = "sunburst",
  branchvalues = "total",
  hoverinfo    = "label+value+percent parent",
  marker       = list(line = list(width = 1), colors = ~color),
  maxdepth     = 3
) %>%
  layout(margin = list(l = 0, r = 0, t = 0, b = 0))

fig
#Export results
# First save as HTML
htmlwidgets::saveWidget(fig, "sunburst_MTG.html")
# Convert to PDF
webshot::webshot("MTG.html", "sunburst_MTG.pdf")

#Draw umap plot merge donors----
#MTG
p <- DimPlot(MTGsub, reduction = "umap", group.by = "ADNC", 
             cols = c("Not AD" = "#aadaf2", "Low" = "#eb6063"))
p + labs(x = "UMAP 1", y = "UMAP 2", title = "MTG")
#DLPFC
p <- DimPlot(DLPFCsub, reduction = "umap", group.by = "ADNC", 
             cols = c("Not AD" = "#aadaf2", "Low" = "#eb6063"))
p + labs(x = "UMAP 1", y = "UMAP 2", title = "DLPFC")

#Draw umap plot merge brain regions----
#M_D
M_D<-merge(MTGsub,DLPFCsub)
M_D <- NormalizeData(M_D)
M_D <- FindVariableFeatures(M_D)
M_D <- ScaleData(M_D)
M_D <- RunPCA(M_D)
M_D <- RunUMAP(object = M_D, dims = 1:30)
p<-DimPlot(M_D,reduction = "umap",group.by = "donor_id")
M_D$sample <- ifelse(Cells(M_D) %in% Cells(MTGsub), "MTGsub",
                     ifelse(Cells(M_D) %in% Cells(DLPFCsub), "DLPFCsub", NA))
p<-DimPlot(M_D,reduction = "umap",group.by = "sample",
           cols = c("MTGsub" = "#d06553", "DLPFCsub" = "#6495c0"),
           raster=FALSE)
p + labs( title = "Dataset1")+
  theme(
    plot.title = element_text(family = "Arial", face = "plain", size = 14)  # face = "plain" removes bold
  )

#findmarkers----
library(Seurat)
library(msigdbr)
library(dplyr)
#findmarkers
fatty_genes_symbol<- unique(c(hallmark$HALLMARK_FATTY_ACID_METABOLISM,
                              kegg$KEGG_FATTY_ACID_METABOLISM,
                              kegg$KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS,
                              kegg$KEGG_PPAR_SIGNALING_PATHWAY))
gene_conversion <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = fatty_genes_symbol,  # Your Gene Symbol list
  mart = mart
)
convert_dict <- setNames(gene_conversion$ensembl_gene_id, gene_conversion$hgnc_symbol)
fatty_genes_ensembl <- unname(convert_dict[fatty_genes_symbol])

valid_fatty_genes_ensembl <- intersect(fatty_genes_ensembl, rownames(MTGsub))

#findmarkers heatmap----
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

seurat_objects <- list(
  "MTG" = MTGsub,
  "DLPFC" = DLPFCsub
)

FA_logFC_all <- list()  
FA_gene_logFC_matrix <- list()

# For each brain region
for (name in names(seurat_objects)) {
  seurat_obj <- seurat_objects[[name]]
  
  # Construct group information: sex + Age_death
  meta_data <- FetchData(seurat_obj, vars = c("sex", "Age_death", "ADNC"))
  meta_data$Group <- paste0(meta_data$sex, "_", meta_data$Age_death)
  seurat_obj <- AddMetaData(seurat_obj, metadata = meta_data$Group, col.name = "Group")
  seurat_obj <- subset(seurat_obj, subset = Group != "male_65 to 77 years old")
  meta_data_filtered <- FetchData(seurat_obj, vars = c("Group", "ADNC"))
  
  # Find valid Groups that contain both Phase1 and Phase2
  group_valid <- meta_data_filtered %>%
    group_by(Group) %>%
    summarise(n_phase = n_distinct(ADNC)) %>%
    filter(n_phase == 2) %>%
    pull(Group)
  
  region_results <- list()        # Store differential results for each Group in current brain region
  top_genes_union_region <- c()   # Union of top10 genes in current brain region
  
  for (g in group_valid) {
    # Extract cells for this Group (only keep Phase 1 and Phase 2 cells)
    seurat_sub <- subset(seurat_obj, subset = Group == g & ADNC %in% c("Not AD", "Low"))
    Idents(seurat_sub) <- "ADNC"
    
    # Differential expression analysis (only analyze fatty acid metabolism related genes)
    markers <- tryCatch({
      FindMarkers(seurat_sub,
                  ident.1 = "Low", ident.2 = "Not AD",
                  features = valid_fatty_genes_ensembl,
                  logfc.threshold = 0, min.pct = 0)
    }, error = function(e) return(NULL))
    
    if (!is.null(markers) && nrow(markers) > 0) {
      # Convert rownames to gene column and add current Group information
      markers <- markers %>% tibble::rownames_to_column("gene")
      markers$Group <- g
      
      # Save complete differential results (don't delete non-significant genes)
      region_results[[g]] <- markers
      
      # Filter significant genes with adj.p < 0.05 for selecting top10
      markers_sig <- markers %>% filter(p_val_adj < 0.05)
      
      if (nrow(markers_sig) > 0) {
        top10 <- markers_sig %>% 
          arrange(desc(abs(avg_log2FC))) %>% 
          head(30) %>% 
          pull(gene)
        # Update union of top genes for current brain region
        top_genes_union_region <- union(top_genes_union_region, top10)
      }
    }
  }
  
  # Merge differential results from all Groups in current brain region, save as long format dataframe
  if (length(region_results) > 0) {
    region_df <- bind_rows(region_results)
    FA_logFC_all[[name]] <- region_df
  }
  
  # Construct gene × group log2FC matrix for current brain region, only keep top genes for current brain region
  df <- FA_logFC_all[[name]]
  df_sub <- df %>% filter(gene %in% top_genes_union_region)
  
  logFC_mat <- df_sub %>%
    dplyr::select(Group, gene, avg_log2FC) %>%
    pivot_wider(names_from = Group, values_from = avg_log2FC) %>%
    column_to_rownames("gene")
  
  FA_gene_logFC_matrix[[name]] <- as.data.frame(logFC_mat)
}

#Convert ensemble to symbol
library(biomaRt)
library(dplyr)

# Initialize Ensembl human database
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

# Merge Ensembl IDs from both brain regions
ensembl_ids <- unique(unlist(lapply(FA_gene_logFC_matrix, rownames)))

# Query Ensembl to Symbol mapping
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Remove rows without symbol and deduplicate
gene_map <- gene_map %>%
  filter(hgnc_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Replace rownames in each table
for (name in names(FA_gene_logFC_matrix)) {
  df <- FA_gene_logFC_matrix[[name]]
  
  common_ids <- intersect(rownames(df), gene_map$ensembl_gene_id)
  symbol_names <- gene_map$hgnc_symbol[match(common_ids, gene_map$ensembl_gene_id)]
  
  # Replace Ensembl IDs in rownames with symbols
  rownames(df)[match(common_ids, rownames(df))] <- symbol_names
  
  # Update back to list
  FA_gene_logFC_matrix[[name]] <- df
}

#Heatmap
# Get DLPFC data
mat <- FA_gene_logFC_matrix[["MTG"]]
# Remove rows where no value has absolute value greater than 1
mat_filtered <- mat[apply(mat, 1, function(x) any(abs(x) > 1)), ]
# Calculate maximum absolute log2FC for each gene (row) across all groups
gene_max <- apply(mat_filtered, 1, function(x) max(abs(x)))
# Sort by maximum absolute value in descending order, take names of top 10 genes
top_genes <- names(sort(gene_max, decreasing = TRUE))[1:10]
# Build new matrix based on selected genes
mat_top10 <- mat_filtered[top_genes, ]
my_colors <- colorRampPalette(c("navy", "white", "firebrick"))(100)
max_val <- max(abs(mat_top10), na.rm = TRUE)

# Set breakpoints based on maximum absolute value (from -max_val to max_val)
my_breaks <- seq(-max_val, max_val, length.out = 101)
# Draw heatmap
pheatmap(mat_top10,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         scale = "none",
         color = my_colors,
         breaks = my_breaks,
         main = "Top 10 fatty acid-related DEGs (DLPFC)",
         fontsize_row = 8,
         fontsize_col = 9)
#write.csv(mat_filtered,file = "DEGs_DLPFC.csv")
#write.csv(mat_filtered,file = "DEGs_MTG.csv")

#astrocyte marker----
cell_marker<-read.csv(file="cell_marker.csv",header = T,sep = ",")
library(biomaRt)
mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
map <- getBM(c("hgnc_symbol", "ensembl_gene_id"), 
             filters = "hgnc_symbol", 
             values = unique(cell_marker$Cell.marker), 
             mart = mart)
cell_marker <- dplyr::left_join(cell_marker, map, by = c("Cell.marker" = "hgnc_symbol")) %>%
  dplyr::mutate(Cell.marker = ensembl_gene_id) %>%
  dplyr::filter(!is.na(Cell.marker))
#Merge two seurat objects
Seurat<-merge(MTGsub,DLPFCsub)

cell_markers <- cell_marker %>% 
  group_by(Cell.name) %>% 
  summarise(markers = list(Cell.marker))

# Create empty list to store expression for each cell type
expression_results <- list()

# Loop through markers for each cell type
for (i in 1:nrow(cell_markers)) {
  cell_type <- cell_markers$Cell.name[i]
  markers <- cell_markers$markers[[i]]
  
  # Extract expression of these markers in single cell data
  expression <- FetchData(Seurat, vars = markers,layer = "counts")
  
  # Calculate mean expression of these markers
  mean_expression <- rowMeans(expression, na.rm = TRUE)
  
  # Save results
  expression_results[[cell_type]] <- mean_expression
}

# Convert to dataframe
expression_summary <- do.call(rbind, lapply(names(expression_results), function(cell_type) {
  data.frame(Cell.name = cell_type, MeanExpression = expression_results[[cell_type]])
}))

# Calculate average expression by cell type
summary_by_celltype <- expression_summary %>%
  group_by(Cell.name) %>%
  summarise(AverageExpression = mean(MeanExpression))

# Draw summary bar plot
summary_by_celltype$Cell.name <- factor(summary_by_celltype$Cell.name, 
                                        levels = summary_by_celltype$Cell.name[order(summary_by_celltype$AverageExpression)])

ggplot(summary_by_celltype, aes(x = Cell.name, y = AverageExpression, fill = Cell.name)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Expression of Cell Type Markers", x = "Cell Type", y = "Expression Level") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  coord_flip()

ggplot(summary_by_celltype, aes(x = Cell.name, y = AverageExpression, fill = Cell.name)) +
  geom_bar(stat = "identity") +
  labs(title = "Mean Expression of Cell Type Markers", x = "Cell Type", y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set3")