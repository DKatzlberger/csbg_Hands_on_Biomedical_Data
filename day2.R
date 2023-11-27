#' ---
#' title: "Day 2 Hands-on-Biomedical-Data"
#' author: Daniel Katzlberger
#' output: github_document
#' ---
#'# Introduction in differential expression analysis
#'Load packages
#+ message=FALSE
#+ 
require(tidyverse)
require(limma)
require(pheatmap)
require(ComplexHeatmap)
require(enrichR)

#'Load data
data <- readRDS("data.RDS") 
design <- readRDS("design.RDS") 
gmap <- readRDS("gmap.RDS") 

#'### Exercise 2.1
#' The `data` object is a matrix, `design` and `gmap` are dfs <br>
#' Checking rows and columns
dim(data)
#' 22055 rows and 230 cols
dim(design)
#' 230 rows and 3 cols
dim(gmap)
#' 21923 rows and 2 cols <br>
#'What is contained in the object `data`
str(data)
#+ results = "hide"
head(data)
#' Samples in columns, gene name in rows and the table cells contain 
#' gene expression (RNA counts) <br>

#' What is contained in the object `design`
str(design)
#+ results = "hide"
head(design)
#' Cell_type, stimulus, organ columns and rows contain samples <br>
#' What is contained in the object `gmap`
str(gmap)
#+ results = "hide"
head(gmap)
#' Different gene names to map 
#+
#'## Subsetting the data
#'Only working with liver_fibroblasts

sub.design <- design[grepl("Liver_Fibroblasts", row.names(design)),]
sub.design <- sub.design |> 
    filter(cell_type == "GP38posCD31neg" & stimulus == "IFNa" | stimulus == "PBS")

sub.data <- data[,row.names(sub.design)]
stopifnot(colnames(sub.data) == row.names(sub.design)) # checking if columns and rows align

#' Renaming the columns to make them easier to interpret
colnames(sub.data) <- gsub("^Liver_Fibroblasts_(.+)_RNA_(\\d)$", "\\1_\\2", colnames(sub.data))
row.names(sub.design) <- gsub("^Liver_Fibroblasts_(.+)_RNA_(\\d)$", "\\1_\\2", row.names(sub.design))

#'## Making correlation Heatmap
#'### Exercise 2.2
#' Correlating the data
corMT <- cor(sub.data, method="spearman")
diag(corMT) <- NA
heatmap(corMT, main="Correlation between liver_fibroblasts with different stimuli", symm = TRUE)

#' Looking at the heatmap the IFNa_1 highly correlates with IFNa_2 and IFNa_3. 
#' That can be expected as the stimulus was the same. However, it also highly correlates 
#' with PBS_1 and also shows noticeable correlation with PBS_2. <br>

#'## MDS projection
#'
data.frame(cmdscale(dist(2-corMT),eig=TRUE, k=2)$points) |> 
    add_column(stimulus = sub.design$stimulus) |> 
    rownames_to_column("sample") |> 
    mutate(sn = gsub("^.+?_(\\d)$", "\\1", sample)) |> 
    ggplot(aes(x=X1,y=X2)) + 
    geom_point(aes(color=stimulus)) + 
    geom_text(aes(label=sn), hjust = 1.5) +
    theme_bw()

#'# Differential expression
#'### Exercise 2.3
#' Adjust the design matrix to have PBS as intercept 
sub.design <- sub.design |> 
    mutate(stimulus = factor(stimulus, levels=c("PBS", "IFNa")))

model <- model.matrix(~stimulus, data=sub.design)

Heatmap(model)

#'## Use limma voom to normalize the data

dataVoom <- voom(sub.data, design=model, plot = TRUE) 
norm.data <- dataVoom$E

#'Types of objects
str(dataVoom)
#'`dataVoom`is a class from `limma` called `EList`
str(norm.data)
#'`norm.data`is a `matrix` <br>
#+
#'## Plotting the data before and after normalization
boxplot(sub.data[1:30,])
#' First code shows the gene expression per stimulus
boxplot(t(sub.data[1:30,]))
#' Second code shows the gene expression per gene <br>
#' Doing the same for normalized data
boxplot(norm.data[1:30,])
boxplot(t(norm.data[1:30,]))
#' Variance can be better explained/seen <br>
#' Looking at density
plot(density(data[8,]))
plot(density(dataVoom$E[8,]))
plot(density(log2(data[8,])))
#' Distribution for the normalized data is better visible. Maybe two states of 
#' the gene are visible
#+
#'## Differential expression
limmaFit <- lmFit(dataVoom, design=model)
limmaFit <- eBayes(limmaFit)

#' Look at coeffs and store them in a list
head(coef(limmaFit)) 
limmaRes <- list()
for(coefx in colnames(coef(limmaFit))){ 
    limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) 
}
limmaRes <- bind_rows(limmaRes, .id = "coef")
limmaRes <- filter(limmaRes, coef != "(Intercept)")

#'# Data interpretation
#'## Vulcano plot
#'### Exercise 2.4
threshold <- abs(limmaRes$logFC) > 2 & limmaRes$P.Value < 0.05

ggplot(limmaRes, aes(x = logFC, y= -log10(P.Value), color = threshold)) +
    geom_point(alpha = 0.5) +
    scale_color_brewer(palette = "Set2") +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(0.05)) + 
    geom_vline(xintercept = c(-2, 2)) 

ggplot(limmaRes, aes(x = logFC, y= -log10(P.Value))) +
    geom_hex()

#'## P-value distribution
#'### Exercise 2.5
ggplot(limmaRes, aes(x=P.Value, fill=factor(floor(AveExpr)))) + 
    geom_histogram()

#' There are a lot of not significant genes with low expression, that might be houshold genes but 
#' have nothing to do with the treatment of the cells. We dont trust them.
#' 
#'## Number of hits
#'### Exercise 2.6
#' Genes tested:
limmaRes |> 
    rownames_to_column("Genes") |> 
    count(unique("Genes"))

#' Significant genes:
limmaRes |> 
    rownames_to_column("Gene") |> 
    filter(adj.P.Val < 0.01 & AveExpr > -4 ) |>  
    count()


limmaResSig <- limmaRes |> 
    rownames_to_column("Gene") |> 
    filter(adj.P.Val < 0.01 & AveExpr > -4 ) # i am very strict on what is significant


#'## Visualizing one gene
#'### Exercise 2.7
#' Get the Gene that is the most upregulated and downregulated
upregGenname <- limmaResSig |> 
    filter(logFC == max(logFC)) |> 
    pull(Gene)

downregGenname <- limmaResSig |> 
    filter(logFC == min(logFC)) |> 
    pull(Gene)

#'Adding the Gene to the design table and visualizing
sub.design <- sub.design |> 
    mutate( upregGen = norm.data[upregGenname,], downregGen = norm.data[downregGenname,])

ggplot(sub.design, aes(x=stimulus, y=upregGen)) +
    geom_point() +
    ylab("E") +
    ggtitle("Upregulatetd", upregGenname)
#' The `logFC` seems plausible as the gene in `IFNa` is higher expressed than in the control `PBS`

ggplot(sub.design, aes(x=stimulus, y=downregGen)) +
    geom_point() +
    ylab("E") +
    ggtitle("Downregulated", downregGenname)

#'## Visualizing multiple genes
#'### Exercise 2.8
#' Getting top 30 absolute differently expressed genes 
goi <- limmaResSig |> 
    slice_max(abs(logFC), n=30) |> 
    pull(Gene)
 

scale(norm.data[goi,]) |> 
     Heatmap(row_split=ifelse(limmaRes[goi,]$logFC > 0, "up", "down"), column_split=sub.design$stimulus)

#'## Enrichment analysis
goi <- limmaResSig |> 
    filter(logFC > 0) |> 
    pull(Gene)

#' Mapping the gene names
goi <- gmap[goi,]$external_gene_name |>  unique()

databases <-  c("MSigDB_Hallmark_2020", "GO_Biological_Process_2021")

enr.res <- enrichr(goi, databases)

#' Hallmark
enr.res$MSigDB_Hallmark_2020 |> 
    ggplot( aes(x=log(Odds.Ratio), y=fct_reorder(Term, Odds.Ratio), size = -log10(P.value))) +
    geom_point() 

#' Bioprocess
#+ plot1, fig.width=10
enr.res$GO_Biological_Process_2021 |> 
    filter(P.value < 0.01) |> 
    ggplot( aes(x=log(Odds.Ratio), y=fct_reorder(Term, Odds.Ratio), size = -log10(P.value))) +
    geom_point() 

#'## Final questions
#'### Exercise 2.10
#'Some genes in PBS_2 are  low expressed, that are also low expressed in IFNa. 
#'Looking at the MDS also this sample appears to be different from the other PBS samples 
#'and somewhat similar to the IFNa samples. Moreover, IFNa_2 expression levels seem
#' to be similar to the PBS samples according to the MDS plot. <br>
#' In this case filtering lowly expressed genes might not be necessary as none 
#' of them are significant according to the p-value distributions. However, the 
#' applied cut off filter them anyways. <br>
#' I would trust the results more would I have a little bit more insight in the
#' applied functions. <br>
#' The results of the enrichment analysis suggests an increase in INFa response and 
#' immune response, which likely is true, because the cells were treated with IFNa.