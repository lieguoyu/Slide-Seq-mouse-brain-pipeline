# path
path <- getwd()
setwd(path)

# Check if the package is installed
if (!require("Seurat")) {
    # If not installed, install the package
    install.packages("Seurat")
    # Loading the package after installation
    library("Seurat")
}
if (!require("ggplot2")) {
    install.packages("ggplot2")
    library("ggplot2")}
if (!require("ggrepel")) {
    install.packages("ggrepel")
    library("ggrepel")}
if (!require("openxlsx")) {
    install.packages("openxlsx")
    library("openxlsx")}
if (!require("org.Mn.eg.db")) {
    install.packages("org.Mn.eg.db")
    library("org.Mn.eg.db")}
if (!require("scCsutomize")) {
    install.packages("scCsutomize")
    library("scCsutomize")}
if (!require("tidyverse")) {
    install.packages("tidyverse")
    library("tidyverse")}
if (!require("RColorBrewer")) {
    install.packages("RColorBrewer")
    library("RColorBrewer")}
if (!require("clusterProfiler")) {
    install.packages("clusterProfiler")
    library("clusterProfiler")}
if (!require("dplyr")) {
    install.packages("dplyr")
    library("dplyr")}

# set seed
set.seed(3)

# Loading the sample
pfc <- readRDS("Prefrontal_Cortex_Sample.rds")

dir.create("./Plots")

# Graphics
jpeg("Plots/pfc_UMAP.jpeg", units = "in", height = 10, width = 15, res = 300)
DimPlot(pfc, reduction = "umap", label = T, label.box = T)
dev.off()

jpeg("Plots/pfc_spatial.jpeg", units = "in", height = 10, width = 15, res = 300)
SpatialPlot(pfc)
dev.off()

jpeg("Plots/pfc_ncount.jpeg", units = "in", height = 10, width = 15, res = 300)
SpatialFeaturePlot(pfc, features = "log10_nCount_RNA")
dev.off()

jpeg("Plots/pfc_nfeature.jpeg", units = "in", height = 10, width = 15, res = 300)
SpatialFeaturePlot(pfc, features = "log10_nFeature_RNA")
dev.off()

jpeg("Plots/pfc_numreads.jpeg", units = "in", height = 10, width = 15, res = 300)
SpatialFeaturePlot(pfc, features = "log10_numReads")
dev.off()


# Searching for marker genes for each cluster
dir.create("./Gene_markers")
gene_marker <- FindAllMarkers(pfc, only.pos = T, logfc.threshold = 0.1)
write.xlsx(gene_marker, "Gene_markers/geneM_per_cluster.xlsx", rowNames = T)

# Saving data from each cluster in an excel sheet
wb <- createWorkbook()
for (i in 0:(length(unique(pfc$seurat_clusters)) - 1)) {
    cluster_nombre <- paste("cluster", i, sep = "")
    addWorksheet(wb, sheetName = cluster_nombre)
    writeData(wb, sheet = cluster_nombre, x = subset(list_expr, subset = list_expr$cluster == i))
}

saveWorkbook(wb, file = "Gene_markers/gene_marker_sheet.xlsx", overwrite = T)

rm(list = setdiff(ls(), "pfc", "gene_marker"))

# Moving the genes column as the first column
wb <- loadWorkbook("Gene_markers/gene_marker_sheet.xlsx")
for (i in 0:(length(unique(pfc$seurat_clusters)) - 1)) {
    cluster_nombre <- paste("cluster", i, sep = "")
    hoja <- read.xlsx(wb, sheet = cluster_nombre)
    ultima <- colnames(hoja)[ncol(hoja)]
    ind_colm <- which(colnames(hoja) == ultima)
    transf <- c(ultima, colnames(hoja)[-colm])
    hoja <- hoja[, transf]
    writeData(wb, sheet = cluster_nombre, x = hoja)
}
saveWorkbook(wb, file = "Gene_markers/gene_marker_sheet.xlsx", overwrite = T)

rm(list = setdiff(ls(), c("pfc", "gene_marker", "wb")))


# Analysis GO
dir.create("./Analysis.GO")
#EnrichGO
wb_new <- createWorkbook()
for (i in 0:(length(unique(pfc$seurat_clusters)) - 1)) {
    cluster_nombre <- paste("cluster", i, sep = "")
    hoja <- read.xlsx(wb, sheet = cluster_nombre)
    genes_list <- hoja[,1]
    ego <- enrichGO(gene          = genes_list,
                    OrgDb         = "org.Mm.eg.db",
                    ont           = "ALL",
                    keyType = "SYMBOL",
                    readable      = TRUE)
    cluster_summary <- data.frame(ego)
    addWorksheet(wb_new, sheetName = cluster_nombre)
    writeData(wb_new, sheet = cluster_nombre, x = cluster_summary)
}
saveWorkbook(wb_new, file = "Analysis.GO/enrichGO_pfc.xlsx")

rm(list = setdiff(ls(), c("pfc", "gene_marker", "wb")))

#GroupGO
wb_new <- createWorkbook()
for (i in 0:(length(unique(pfc$seurat_clusters)) - 1)) {
    cluster_nombre <- paste("cluster", i, sep = "")
    addWorksheet(wb_new, sheetName = cluster_nombre)
    hoja <- read.xlsx(wb, sheet = cluster_nombre)
    genes_list <- hoja[,1] #Using function groupGO across 3 ontology types
    
    ggo_bp <- groupGO(gene = genes_list,
                      OrgDb = "org.Mm.eg.db",
                      keyType = "SYMBOL",
                      ont = "BP",
                      level = 5,
                      readable = TRUE)
    bp_df <- data.frame(ggo_bp)
    bp_df$ontology <- "BP"
    
    ggo_mf <- groupGO(gene = genes_list,
                      OrgDb = "org.Mm.eg.db",
                      keyType = "SYMBOL",
                      ont = "MF",
                      level = 5,
                      readable = TRUE)
    mf_df <- data.frame(ggo_mf)
    mf_df$ontology <- "MF"
    
    ggo_cc <- groupGO(gene = genes_list,
                      OrgDb = "org.Mm.eg.db",
                      keyType = "SYMBOL",
                      ont = "CC",
                      level = 5,
                      readable = TRUE)
    cc_df <- data.frame(ggo_cc)
    cc_df$ontology <- "CC"
    #merge all data frames
    all_ggo <- rbind(bp_df, mf_df, cc_df)
    #move last column to first position
    n = ncol(all_ggo)
    new.order = c(n, 1:(n - 1))
    all_ggo <- all_ggo[, new.order]
    
    all_ggo <- subset(all_ggo, all_ggo$Count > 0) #Eliminate 0 count
    writeData(wb_new, sheet = cluster_nombre, x = all_ggo)  #Writing excel
}
saveWorkbook(wb_new, file = "Analysis.GO/groupGO_pfc.xlsx")

rm(list = setdiff(ls(), c("pfc", "gene_marker")))

# Laboratory classification method
dir.create("./lab_classification")

# Loading reference database
ref <- read.table("Annotation_DBbrain.txt", sep = "\t", header = T)

wb <- createWorkbook()
for (i in 0:(length(unique(pfc$seurat_clusters)) - 1)) {
    clusters <- subset(gene_marker, gene_marker$cluster == i)
    list_genes <- toupper(clusters$gene)
    re_enrich <- enricher(gene = list_genes, 
                          TERM2GENE = ref, 
                          pvalueCutoff = 1, 
                          qvalueCutoff = 1, 
                          minGSSize = 1,
                          maxGSSize = 100000)
    #Handle errors
    tryCatch({
        res <- re_enrich@result
    }, error = function(e) {
        # If there is a error, print a massage and write it in the excel sheet
        cat("Error en el cluster", i, ":", conditionMessage(e), "\n")
        note <- paste("No gene can be mapped in cluster", i)
        addWorksheet(wb, sheetName = paste("cluster", i, sep = ""))
        writeData(wb, sheet = paste("cluster", i, sep = ""), x = note, startCol = 1, startRow = 1)
        return(NULL)
    })
    # If there is not a error, continue with the loop and write down the results in the excel sheet
    if (!is.null(re_enrich)) {
        res <- re_enrich@result
        cluster_nombre <- paste("cluster", i, sep = "")
        addWorksheet(wb, sheetName = cluster_nombre)
        writeData(wb, sheet = cluster_nombre, x = res)
    }
}

saveWorkbook(wb, file = "lab_classification/pfc_lab_class.xlsx", overwrite = T)

rm(list = setdiff(ls(), "pfc"))

# Add the classification in the datasets
pfc_class <- c("Neurons","Astrocytes", "Oligodendrocytes", "Neurons", "Unassigned", "Unassigned",
                 "Unassigned", "Vascular_cells", "Vascular_cells", "Oligodendrocytes", "Neurons",
                 "Unassigned", "Oligodendrocytes", "Astrocytes", "Astrocytes",
                 "Oligodendrocytes", "Unassigned", "Unassigned", "Astrocytes", "Oligodendrocytes",
                 "Unassigned", "Astrocytes", "Astrocytes", "Unassigned", "Unassigned", "Unassigned")
names(pfc_class) <- c(0,seq(1:25))

pfc <- Rename_Clusters(seurat_object = pfc, new_idents = pfc_class, meta_col_name = "lab_class")

# Add the classificacion as metadata 
class <- as.data.frame(Idents(pfc))
pfc <- AddMetaData(pfc, class, "class_name")

saveRDS(pfc, "Prefrontal_Cortex_Sample.rds")

jpeg("Plots/DimPlot_pfc_lab_class.jpeg", units = "in", height = 10, width = 15, res = 300)
DimPlot(objeto_seu_test, label = T, label.box = T, pt.size = 0.5) 
dev.off()

rm(list = setdiff(ls(), "pfc"))

# Apply the same method to the 'Unassigned' class
dir.create("lab_classification/unassigned")
unsg <- subset(pfc, idents = "Unassigned")
unsg <- SCTransform(unsg, conserve.memory = T, verbose = T)
unsg <- RunPCA(unsg, verbose = T)
unsg <- RunUMAP(unsg, dims = 1:10, verbose = T)
unsg <- FindNeighbors(unsg, dims = 1:10, verbose = T)
unsg <- FindClusters(unsg, verbose = T)

gm <- FindAllMarkers(unsg, logfc.threshold = 0.1, only.pos = T)

ref <- read.table("Annotation_DBbrain.txt", sep = "\t", header = T)

wb <- createWorkbook()
for (i in 0:(length(unique(unsg$seurat_clusters)) - 1)) {
    clusters <- subset(gm, gm$cluster == i)
    list_genes <- toupper(clusters$gene)
    re_enrich <- enricher(gene = list_genes, 
                          TERM2GENE = ref, 
                          pvalueCutoff = 1, 
                          qvalueCutoff = 1, 
                          minGSSize = 1,
                          maxGSSize = 100000)
    #Handle errors
    tryCatch({
        res <- re_enrich@result
    }, error = function(e) {
        # If there is a error, print a massage and write it in a excel sheet
        cat("Error en el cluster", i, ":", conditionMessage(e), "\n")
        note <- paste("No gene can be mapped in cluster", i)
        addWorksheet(wb, sheetName = paste("cluster", i, sep = ""))
        writeData(wb, sheet = paste("cluster", i, sep = ""), x = note, startCol = 1, startRow = 1)
        return(NULL)
    })
    # If there is not a error, continue with the loop and write down the results in excel sheet
    if (!is.null(re_enrich)) {
        res <- re_enrich@result
        cluster_nombre <- paste("cluster", i, sep = "")
        addWorksheet(wb, sheetName = cluster_nombre)
        writeData(wb, sheet = cluster_nombre, x = res)
    }
}
saveWorkbook(wb, file = "lab_classification/unassigned/pfc_unassigned.xlsx", overwrite = T)

unsg_class <- c("Neurons_unclassified","Grial_cells_unclassified", "Unassigned", "Neurons_unclassified", "Astrocytes_unclassified",
                 "Oligodendrocytes_unclassified", "Oligodendrocytes_unclassified", "Unassigned", "Oligodendrocytes_unclassified", "Neurons_unclassified",
                 "Neurons_unclassified", "Oligodendrocytes_unclassified", "Neurons_unclassified", "Neurons_unclassified")
names(unsg_class) <- c(0,seq(1:13))
unsg <- Rename_Clusters(seurat_object = unsg, new_idents = unsg_class, meta_col_name = "unsg_class")
saveRDS(unsg, "lab_classification/unassigned/unsg.rds")

rm(list = setdiff(ls(), "pfc"))

# Apply the same method to the 'Neurons' class
dir.create("lab_classification/neurons")
neu <- subset(pfc, idents = "Neurons")
neu <- SCTransform(unsg, conserve.memory = T, verbose = T)
neu <- RunPCA(unsg, verbose = T)
neu <- RunUMAP(unsg, dims = 1:10, verbose = T)
neu <- FindNeighbors(unsg, dims = 1:10, verbose = T)
neu <- FindClusters(unsg, verbose = T)

gm <- FindAllMarkers(neu, logfc.threshold = 0.1, only.pos = T)

ref <- read.table("Annotation_DBonlyneurons.txt", sep = "\t", header = T)

wb <- createWorkbook()
for (i in 0:(length(unique(neu$seurat_clusters)) - 1)) {
    clusters <- subset(gm, gm$cluster == i)
    list_genes <- toupper(clusters$gene)
    re_enrich <- enricher(gene = list_genes, 
                          TERM2GENE = ref, 
                          pvalueCutoff = 1, 
                          qvalueCutoff = 1, 
                          minGSSize = 1,
                          maxGSSize = 100000)
    #Handle errors
    tryCatch({
        res <- re_enrich@result
    }, error = function(e) {
        # If there is a error, print a massage and write it in a excel sheet
        cat("Error en el cluster", i, ":", conditionMessage(e), "\n")
        note <- paste("No gene can be mapped in cluster", i)
        addWorksheet(wb, sheetName = paste("cluster", i, sep = ""))
        writeData(wb, sheet = paste("cluster", i, sep = ""), x = note, startCol = 1, startRow = 1)
        return(NULL)
    })
    # If there is not a error, continue with the loop and write down the results in excel sheet
    if (!is.null(re_enrich)) {
        res <- re_enrich@result
        cluster_nombre <- paste("cluster", i, sep = "")
        addWorksheet(wb, sheetName = cluster_nombre)
        writeData(wb, sheet = cluster_nombre, x = res)
    }
}
saveWorkbook(wb, file = "lab_classification/neurons/pfc_neurons.xlsx", overwrite = T)

rm(list = setdiff(ls(), "pfc"))

# Allen classification method
dir.create("./allen_classification")

# Loading file.csv with allen classification
allen <- read.csv("pfc_allen.csv", header = F)

#Filter rows and colomns
allen <- allen[5:nrow(allen),]
colnames(allen) <- allen[1,]
allen <- allen[2:nrow(allen),]
row.names(allen) <- allen[,1]
allen <- allen[,2:ncol(allen)]
allen <- allen[, 1:6]
allen <- allen[, -c(1,4)]

write.csv(allen, "allen_classification/allen_mycellmap_pfc.csv")

#Select class and subclass along with it's correlation values
allen_class <- allen[, 1, drop = F]
allen_sublcass <- allen[, 3, drop = F]
allen_correlation <- allen[, 2, drop = F]

# Add these new tables metadata 
pfc <- AddMetaData(pfc, allen_class, col.name = "allen_className")
pfc <- AddMetaData(pfc, allen_sublcass, col.name = "allen_subclassName")
pfc <- AddMetaData(pfc, allen_correlation, col.name = "allen_classCorrelation")

saveRDS(pfc, "Prefrontal_Cortex_Sample.rds")

mycolor <- c(
    "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", 
    "#FB9A99", "palegreen2","#CAB2D6", "#FDBF6F", "darkslategrey", "khaki2", "maroon", "orchid1", 
    "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown", "lightpink", "tomato", "coral3", "seagreen", "slategrey", 
    "tan", "midnightblue", "maroon4", "plum"
)

jpeg("Plots/SpatialDimplot_allen_pfc.jpeg", units = "in", height = 10, width = 15, res = 300)
SpatialDimPlot(lc, group.by = "allen_className") + scale_fill_manual(values = c34) + 
    guides(fill = guide_legend(override.aes = list(size = 5)))
dev.off()

jpeg("Plots/Dimplot_allen_pfc.jpeg", units = "in", height = 10, width = 15, res = 300)
DimPlot(lc, group.by = "allen_className", cols = c34, pt.size = 0.5)
dev.off()

rm(list = setdiff(ls(), c("pfc", "allen_class")))

# Simplifying the annotation of categories
glut <- grep("Glut", allen_class$class_name)
allen_class$class_name[glut] <- "Glutamatergic"
gaba <- grep("GABA", allen_class$class_name)
allen_class$class_name[gaba] <- "GABAergic"
inm <- grep("34 Immune", allen_class$class_name)
allen_class$class_name[inm] <- "Immune"
vas <- grep("33 Vascular", allen_class$class_name)
allen_class$class_name[vas] <- "Vascular"
oli <- grep("31 OPC-Oligo", allen_class$class_name)
allen_class$class_name[oli] <- "OPC-Oligodendrocytes"
ast <- grep("30 Astro-Epen", allen_class$class_name)
allen_class$class_name[ast] <- "Astrocytes"
oec <- grep("32 OEC", allen_class$class_name)
allen_class$class_name[oec] <- "OEC"
sero <- grep("22 MB-HB Sero", allen_class$class_name)
allen_class$class_name[sero] <- "Serotonergic"
dopa <- grep("21 MB Dopa", allen_class$class_name)
allen_class$class_name[dopa] <- "Dopaminergic"

pfc <- AddMetaData(pfc, allen_class, "allen_class_simplified")
saveRDS(pfc, "Prefrontal_Cortex_Sample.rds")

c9 <- c("tan", "green4","skyblue2", "gold1","#6A3D9A", "#FF7F00", "lightpink", "orchid1", "#E31A1C")

jpeg("Plots/SpatialDimplot_allenclass_pfc_simplified.jpeg", units = "in", height = 10, width = 15, res = 300)
SpatialDimPlot(pfc, group.by = "allen_class_simplified") + scale_fill_manual(values = c9) + 
    guides(fill = guide_legend(override.aes = list(size = 5), title = "allen_class_simplified")) 
dev.off()

jpeg("Plots/Dimplot_allenclass_pfc_simplified.jpeg", units = "in", height = 10, width = 15, res = 300)
DimPlot(pfc, group.by = "allen_class_simplified") + scale_color_manual(values = c9) + 
    guides(fill = guide_legend(override.aes = list(size = 5), title = "allen_class_simplified")) 
dev.off()

rm(list = setdiff(ls(), "pfc"))

# Contrast between laboratory method vs. allen's method
dir.create("./lab_vs_allen")
dir.create("lab_vs_allen/subset_class")
dir.create("lab_vs_allen/freq_allen")


#ast
ast <- subset(pfc, idents = "Astrocytes")
ast_lab <- as.data.frame(ast$class_name) 
colnames(ast_lab) <- "class"
ast_allen <- as.data.frame(ast$allen_className)
colnames(ast_allen) <- "class"
ast_allen <- subset(ast_allen, subset = ast_allen$class == "30 Astro-Epen")

reg <- data.frame(
    id = "astrocytes",
    celltype_lab = as.numeric(length(rownames(ast_lab))),
    cel_contras_allen = as.numeric(length(rownames(ast_lab)[rownames(ast_lab) %in% rownames(ast_allen)]))
)
tabla_allen <- as.data.frame(table(ast$allen_className))
rownames(tabla_allen) <- tabla_allen$Var1
tabla_allen <- tabla_allen[,-1, drop = F] 
write.csv(tabla_allen, "lab_vs_allen/freq_allen/freq_ast_allen_class.csv")
cells_ast <- rownames(ast_lab)[rownames(ast_lab) %in% rownames(ast_allen)]
df_ast <- data.frame(row.names = cells_ast, class_confirmed = rep("Astrocytes", length(cells_ast))) 
write.csv(df_ast, "lab_vs_allen/subset_class/df.ast.csv")
rm(list = setdiff(ls(), c("pfc", "reg")))

#Vas
vas <- subset(pfc, idents = "Vascular_cells")
vas_lab <- as.data.frame(vas$class_name) 
colnames(vas_lab) <- "class"
vas_allen <- as.data.frame(vas$allen_className)
colnames(vas_allen) <- "class"
vas_allen <- subset(vas_allen, subset = vas_allen$class == "33 Vascular")
reg[nrow(reg) + 1,] <-  c("Vascular_cells", 
                          as.numeric(length(rownames(vas_lab))), 
                          as.numeric(length(rownames(vas_lab)[rownames(vas_lab) %in% rownames(vas_allen)]))
)
tabla_allen <- as.data.frame(table(vas$allen_className))
rownames(tabla_allen) <- tabla_allen$Var1
tabla_allen <- tabla_allen[,-1, drop = F] 
write.csv(tabla_allen, "lab_vs_allen/freq_allen/freq_vas_allen_class.csv")
cells_vas <- rownames(vas_lab)[rownames(vas_lab) %in% rownames(vas_allen)]
df_vas <- data.frame(row.names = cells_vas, class_confirmed = rep("Vascular_cells", length(cells_vas))) 
write.csv(df_vas, "lab_vs_allen/subset_class/df.vas.csv")
rm(list = setdiff(ls(), c("pfc", "reg")))

#Oli
oli <- subset(pfc, idents = "Oligodendrocytes")
oli_lab <- as.data.frame(oli$class_name)
colnames(oli_lab) <- "class"
oli_allen <- as.data.frame(oli$allen_className)
colnames(oli_allen) <- "class"
oli_allen <- subset(oli_allen, subset = oli_allen$class == "31 OPC-Oligo")
reg[nrow(reg) + 1,] <-  c("Oligodendrocytes", 
                          as.numeric(length(rownames(oli_lab))), 
                          as.numeric(length(rownames(oli_lab)[rownames(oli_lab) %in% rownames(oli_allen)]))
)
tabla_allen <- as.data.frame(table(oli$allen_className))
rownames(tabla_allen) <- tabla_allen$Var1
tabla_allen <- tabla_allen[,-1, drop = F] 
write.csv(tabla_allen, "lab_vs_allen/freq_allen/freq_oli_allen_class.csv")
cells_oli <- rownames(oli_lab)[rownames(oli_lab) %in% rownames(oli_allen)]
df_oli <- data.frame(row.names = cells_oli, class_confirmed = rep("Oligodendrocytes", length(cells_oli))) 
write.csv(df_oli, "lab_vs_allen/subset_class/dfoli.csv")
rm(list = setdiff(ls(), c("pfc", "reg")))

#Mcg 
mcg <- subset(pfc, idents = "Microglia")
mcg_lab <- as.data.frame(mcg$class_name)
colnames(mcg_lab) <- "class"
mcg_allen <- as.data.frame(mcg$allen_className)
colnames(mcg_allen) <- "class"
mcg_allen <- subset(mcg_allen, subset = mcg_allen$class == "34 Immune")
reg[nrow(reg) + 1,] <-  c("Microglia", 
                          as.numeric(length(rownames(mcg_lab))), 
                          as.numeric(length(rownames(mcg_lab)[rownames(mcg_lab) %in% rownames(mcg_allen)]))
)
tabla_allen <- as.data.frame(table(mcg$allen_className))
rownames(tabla_allen) <- tabla_allen$Var1
tabla_allen <- tabla_allen[,-1, drop = F] 
write.csv(tabla_allen, "lab_vs_allen/freq_allen/freq_mcg_allen_class.csv")
cells_mcg <- rownames(mcg_lab)[rownames(mcg_lab) %in% rownames(mcg_allen)]
df_mcg <- data.frame(row.names = cells_mcg, class_confirmed = rep("Microglia", length(cells_mcg))) 
write.csv(df_mcg, "lab_vs_allen/subset_class/dfmcg.csv")
rm(list = setdiff(ls(), c("pfc", "reg")))

#Neu
neu <- subset(pfc, idents = "Neurons")
neu_lab <- as.data.frame(neu$class_name)
colnames(neu_lab) <- "class"
neu_allen <- as.data.frame(neu$allen_className)
colnames(neu_allen) <- "class"

#Eliminating not neuron classes
neu_allen <- subset(neu_allen, neu_allen$class != "30 Astro-Epen")
neu_allen <- subset(neu_allen, neu_allen$class != "31 OPC-Oligo")
neu_allen <- subset(neu_allen, neu_allen$class != "34 Immune")
neu_allen <- subset(neu_allen, neu_allen$class != "33 Vascular")
neu_allen <- subset(neu_allen, neu_allen$class != "32 OEC")
reg[nrow(reg) + 1,] <-  c("Neurons", 
                          as.numeric(length(rownames(neu_lab))), 
                          as.numeric(length(rownames(neu_lab)[rownames(neu_lab) %in% rownames(neu_allen)]))
)

write.table(reg, "lab_vs_allen/class_labvsAllenDB_pfc.txt", sep = "\t", dec = ",")

neu_class <- neu_lab$class[rownames(neu_lab) %in% row.names(neu_allen)]
neu_cells <- rownames(neu_lab)[rownames(neu_lab) %in% row.names(neu_allen)]
neu_class_allen <- neu_allen$class[row.names(neu_allen) %in% row.names(neu_lab)]

df_neu <- data.frame(row.names = neu_cells, class_allen = neu_class_allen, class_confirmed = neu_class)
write.csv(df_neu, "lab_vs_allen/subset_class/df_neu_pfc_allen.csv")

tabla_allen <- as.data.frame(table(neu$allen_className))
rownames(tabla_allen) <- tabla_allen$Var1
neu <- tabla_allen[,-1, drop = F]
write.csv(tabla_allen, "lab_classification/freq_allen/freq_neu_allen_class.csv")
rm(list = setdiff(ls(), c("pfc", "neu")))

# Classification neurons cells
neu$class_confirmed <- 0
glut <- grep("Glut", neu$class_allen)
gaba <- grep("GABA", neu$class_allen)
sero <- grep("Sero", neu$class_allen)
dopa <- grep("Dopa", neu$class_allen)
neu$class_confirmed[glut] <- "Glutamatergic"
neu$class_confirmed[gaba] <- "GABAergic"
neu$class_confirmed[sero] <- "Serotonergic"
neu$class_confirmed[dopa] <- "Dopaminergic"
neu <- neu[,2, drop = F]
write.csv(neu, "lab_vs_allen/subset_class/df_neu_pfc_allen.csv")

rm(list = setdiff(ls(), c("pfc", "neu")))
# Add confirmed cell type
ast <- read.csv("lab_vs_allen/subset_class/dfast.csv", row.names = 1)
oli <- read.csv("lab_vs_allen/subset_class/dfoli.csv", row.names = 1)
vas <- read.csv("lab_vs_allen/subset_class/dfvas.csv", row.names = 1)
mcg <- read.csv("lab_vs_allen/subset_class/dfmcg.csv", row.names = 1)
df_merged <- rbind(neu, ast, oli, vas, mcg)

pfc <- AddMetaData(pfc, df_merged, "celltype_confirmed")

saveRDS(pfc, "Prefrontal_Cortex_Sample.rds")
rm(list = setdiff(ls(), "pfc"))

# Cell type counts
dir.create("./Count")
t_final <- as.data.frame(pfc$celltype_confirmed)
colnames(t_final) <- "Classification"
t_final <- replace(t_final, is.na(t_final), "Unconfirmed")
df_statistic <- data.frame(
    row.names = "Cells",
    Count = as.numeric(length(t_final$Classification)),
    Count_class_confirmed = as.numeric(length(t_final$Classification[t_final$Classification != "Unconfirmed"]))
)
df_statistic$Percent <- round((df_statistic[1,2] * 100)/df_statistic[1,1], 2)
write.table(df_statistic, "lc_table_%_accurance.txt",sep = "\t", dec = ",")

df <- data.frame(Classification = unique(t_final$Classification), Count = NA, Percent = NA)
for (i in 1:nrow(df)) {
    num <- length(grep(df$Classification[i], t_final$Classification))
    df$Count[i] <- num
    perc <- round((df$Count[i] * 100)/length(t_final$Classification), 2)
    df$Percent[i] <- perc
}
write.table(df, "Count/pfc_table_count_per_class.txt", sep = "\t", dec = ",")
rm(list = setdiff(ls(), c("df", "df_statistic", "pfc")))

# Chart Pie plots
dir.create("./piechart_plots")

# Chart Pie plot1, percent confirmed & uncofirmed
dataf <- data.frame(
    Labels = c("Class_confirmed", "Class_unconfirmed"),
    Percent = c(df_statistic[1,3], df[1,3])
)

jpeg("piechart_plots/Piechartplot_pfc_conf_vs_unconf.jpeg", units = "in", height = 10, width = 15, res = 300)
ggplot(dataf, aes(x = "", y = Percent, fill = Labels)) + geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) + theme_void() +
    geom_text(aes(label = paste0(Percent, "%")), position = position_stack(vjust = 0.5)) +
    scale_fill_brewer(palette = "Set1") 
dev.off()

# Chart Pie plot2, percent of each confirmed celltype
# Transform dataframe
df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(Percent))), 
           pos = Percent/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Percent/2, pos))

jpeg("piechart_plots/Piechartplot_pfc_confirmed_celltype.jpeg", units = "in", height = 10, width = 15, res = 300)
ggplot(df, aes(x = "" , y = Percent, fill = fct_inorder(Classification))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Set1") +
    geom_label_repel(data = df2,
                     aes(y = pos, label = paste0(Percent, "%")),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Classification")) +
    theme_void()
dev.off()

rm(list = setdiff(ls(), "pfc"))

# Chart Pie plot3, percent of each unconfirmed celltype
# Table with uncofirmed celltype 
t_allen <- as.data.frame(pfc$celltype_confirmed)
colnames(t_allen) <- "Classification"
t_allen <- replace(t_allen, is.na(t_allen), "Unconfirmed")
t_allen_u <- subset(t_allen, t_allen$Classification == "Unconfirmed")

# Table with lab classification
t_class <- as.data.frame(lc$class_name)
colnames(t_class) <- "Classification"

# Table with cells common to both tables
b_table <- subset(t_class, rownames(t_class) %in% rownames(t_allen_u))
b_table$Classification <- as.character(b_table$Classification)
rm(list = setdiff(ls(), c("pfc","b_table", "t_allen")))

# Combine tables
t_allen_c <- subset(t_allen, t_allen$Classification != "Unconfirmed")
rm(list = setdiff(ls(), c("pfc","b_table", "t_allen_c")))

#Count of each celltype
df <- data.frame(Classification = sort(unique(b_table$Classification)), Count = NA, Percent = NA)
for (i in 1:nrow(df)) {
    num <- length(grep(paste0("^", df$Classification[i], "$"), b_table$Classification))
    df$Count[i] <- num
    perc <- round((df$Count[i] * 100)/length(b_table$Classification), 2)
    df$Percent[i] <- perc
}

write.table(df, "Count/pfc_count_unconfirmed_celltype.txt", sep = "\t", dec = ",")

df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(Percent))), 
           pos = Percent/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Percent/2, pos))

jpeg("piechart_plots/Piechartplot_pfc_unconfirmed_celltype.jpeg", units = "in", height = 10, width = 15, res = 300)
ggplot(df, aes(x = "" , y = Percent, fill = fct_inorder(Classification))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_brewer(palette = "Set1") +
    geom_label_repel(data = df2,
                     aes(y = pos, label = paste0(Percent, "%")),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Classification")) +
    theme_void()
dev.off()

rm(list = setdiff(ls(), c("pfc", "b_table", "t_allen_c")))

# Chart Pie plot4, percent of final classification
comb_t <- rbind(b_table, t_allen_c)
pfc <- AddMetaData(pfc, comb_t, "final_classification")
rm(list = setdiff(ls(), c("pfc","comb_t")))

# Count of each celltype 
df <- data.frame(Classification = unique(comb_t$Classification), Count = NA, Percent = NA)
for (i in 1:nrow(df)) {
    num <- length(grep(paste0("^", df$Classification[i], "$"), comb_t$Classification))
    df$Count[i] <- num
    perc <- round((df$Count[i] * 100)/length(comb_t$Classification), 2)
    df$Percent[i] <- perc
}

if (sum(df$Count) == length(comb_t$Classification)) {
    cat("all fine")
} else {
    cat("something was wrong")
}

write.table(df, "Count/pfc_count_final_classification.txt", sep = "\t", dec = ",")
rm(list = setdiff(ls(), c("pfc","df")))

df2 <- df %>% 
    mutate(csum = rev(cumsum(rev(Percent))), 
           pos = Percent/2 + lead(csum, 1),
           pos = if_else(is.na(pos), Percent/2, pos))

mycolors = c(brewer.pal(name = "Dark2", n = 8), brewer.pal(name = "Paired", n = 6))

jpeg("piechart_plots/Piechartplot_pfc_finalcellstypes.jpeg", units = "in", height = 10, width = 15, res = 300)
ggplot(df, aes(x = "" , y = Percent, fill = fct_inorder(Classification))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    scale_fill_manual(values = mycolors) +
    geom_label_repel(data = df2,
                     aes(y = pos, label = paste0(Percent, "%")),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    guides(fill = guide_legend(title = "Classification")) +
    theme_void()
dev.off()

rm(list = setdiff(ls(), "pfc"))

# Cortical layer analysis
dir.create("./cortical_layer")

# Selecting cells expressing marker genes
ctip2_cell <- WhichCells(pfc, expression = Bcl11b > 0) 
brn1_cell <- WhichCells(pfc, expression = Pou3f3 > 0) 
tbr1_cell <- WhichCells(pfc, expression = Tbr1 > 0) 
satb2_cell <- WhichCells(pfc, expression = Satb2 > 0) 
cux1_cell <- WhichCells(pfc, expression = Cux1 > 0) 

tb <- data.frame(row.names = tbr1_cell, layer_marker = rep("Tbr1", length(tbr1_cell)))
ct2 <- data.frame(row.names = ctip2_cell, layer_marker = rep("Ctip2", length(ctip2_cell)))
stb2 <- data.frame(row.names = satb2_cell, layer_marker = rep("Satb2", length(satb2_cell)))
br1 <- data.frame(row.names = brn1_cell, layer_marker = rep("Brn1", length(brn1_cell)))
cx1 <- data.frame(row.names = cux1_cell, layer_marker = rep("Cux1", length(cux1_cell)))

# Removing duplicates
#tb
tbrVSctip <- rownames(tb)[rownames(tb) %in% rownames(ct2)]
tb <- tb[!(rownames(tb) %in% tbrVSctip), ,drop = F]
trbVSstb <- rownames(tb)[rownames(tb) %in% rownames(stb2)]
tb <- tb[!(rownames(tb) %in% trbVSstb), ,drop = F]
tbrVScx1 <- rownames(tb)[rownames(tb) %in% rownames(cx1)]
tb <- tb[!(rownames(tb) %in% tbrVScx1), ,drop = F]
trbVSbr <- rownames(tb)[rownames(tb) %in% rownames(br1)]
tb <- tb[!(rownames(tb) %in% trbVSbr), ,drop = F]
#ctip
ctipVSsatb2 <- rownames(ct2)[rownames(ct2) %in% rownames(stb2)]
ct2 <- ct2[!(rownames(ct2) %in% ctipVSsatb2), ,drop = F]
ctipVScx1 <- rownames(ct2)[rownames(ct2) %in% rownames(cx1)]
ct2 <- ct2[!(rownames(ct2) %in% ctipVScx1), ,drop = F]
ctipVSbr1 <- rownames(ct2)[rownames(ct2) %in% rownames(br1)]
ct2 <- ct2[!(rownames(ct2) %in% ctipVSbr1), ,drop = F]
#satb2
satb2VScux1 <- rownames(stb2)[rownames(stb2) %in% rownames(cx1)]
satb2VSbr1 <- rownames(stb2)[rownames(stb2) %in% rownames(br1)]
#cux1
cux1VSbrn1 <- rownames(cx1)[rownames(cx1) %in% rownames(br1)]
cx1 <- cx1[!(rownames(cx1) %in% cux1VSbrn1), ,drop = F]

df <- rbind(tb, ct2, stb2, br1, cx1)
write.csv(df, "cortical_layer/cells_cortical_layer.csv")

pfc <- AddMetaData(pfc, metadata = df, col.name = "layer_marker")
col <- c("Tbr1" = "red", "Ctip2" = "blue", "Satb2" = "green", "Cux1" = "yellow", "Brn1" = "orange")
jpeg("Dimplot_marker_layer.jpeg", units = "in", height = 10, width = 15, res = 300)
SpatialDimPlot(pfc, group.by = "layer_marker", cols = col) + 
    scale_fill_manual(na.value = "white", values = col)
dev.off()

saveRDS(pfc, "Prefrontal_Cortex_Sample.rds")

rm(list = ls())
