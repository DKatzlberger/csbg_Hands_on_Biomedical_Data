#' ---
#' title: "Day 2 Hands-on-Biomedical-Data"
#' author: Daniel Katzlberger
#' output: github_document
#' ---
#'# My own Dataset
#' Load packages
#+ message=FALSE
require(tidyverse)
require(limma)
require(patchwork)
require(pheatmap)
require(ComplexHeatmap)
require(enrichR)
require(ggrepel)
require(gemma.R)

#'# Getting Data

gse <- "GSE12172"

get_datasets(gse) |> 
    select(taxon.Name, experiment.Accession, experiment.SampleCount, experiment.Description)

design <- get_dataset_design(gse) # get design table
e <- get_dataset_processed_expression(gse) # get expression
e <- as.data.frame(e) # make it to data frame

data <- as.matrix(e[,row.names(design)]) 
row.names(data) <- e$GeneSymbol # make it expression data
data <- as.data.frame(data)

dim(design)
dim(data)
#' Data contains all 90 samples
#+
#'## Further inside in the data
design |> 
    count(disease.staging)

design |> 
    count(phenotype)

#' Most samples are from disease stage III with 57 samples. There are only 7 samples from disease stage II. 
#' The phenotype contains the information from which tumor the sample origins.The description of the 
#' data says that 30 samples are derived from LMP tumors and 60 from Serous tumors found in the variable `phenotypes`. <br>
#' For the first part the differential expression between disease stages is analysed. 

# new.names.df <- design |> 
#    rownames_to_column("sample") |> 
#    unite(new_name, sample, disease.staging) |> 
#    select(new_name)

# colnames(data) <- new.names.df$new_name

#'## Filter samples with `<NA>` and removing it from design table

NA.column <- names(which(sapply(data, function(x)all(is.na(x))))) # sample that is removed

data <- data |> select(where(function(x) all(!is.na(x)))) # removing sample from data
design <- design[rownames(design) %in% colnames(data),] # removing sample from design table

stopifnot(colnames(data) == row.names(design)) # checking if columns and rows align

#'# Correlation Analysis
#'## Correlation heatmap

corMT <- cor(data, method="spearman") 
diag(corMT) <- NA
#+ plot1, fig.width=10
Heatmap(corMT)

#'## MDS projection
#+ plot2, fig.width=10
data.frame(cmdscale(dist(2-corMT),eig=TRUE, k=2)$points) |> 
    rownames_to_column("sample") |> 
    mutate(sample = str_replace_all(sample, "^[^_]*_","")) |> 
    ggplot(aes(x=X1,y=X2)) + 
    geom_point() +
    geom_text_repel(aes(label=sample)) +
    theme_bw()

#' The heatmap and the MDS projection suggest that there might be two clusters in our data.
#' Try to figure out what the clusters might be. 
#+
#'## Clustering
tree <- hclust(dist(t(data)))
#+ plot3, fig.width=10
plot(tree)

cutree <- as.data.frame(cutree(tree, k = 2)) |> 
    mutate(hcluster = cutree(tree, k = 2)) |> 
    select(hcluster)

cutree |> 
    count(hcluster)

#' It seems the clustering  differentiates the two phenotypes. However, it is not able
#' to distinguish all samples correctly. 
#+
#'# Differential expression and data normalization
#'## Model matrix
#' Disease stage I is used as intercept

design <- design |> 
    mutate(disease.staging = factor(disease.staging, ordered = FALSE)) |> 
    mutate(disease.staging = relevel(disease.staging, "stage_I")) # relevel to make stage I control


model <- model.matrix(~disease.staging, data=design)
#+ plot4, fig.width=10
Heatmap(model)

#'## Distribution of data
#' Looking at distribution to decide if normalization is needed.
#+ plot5, fig.width=10
boxplot(data[1:30,])
#+ plot6, fig.width=10
boxplot(t(data[1:30,]))

#' Data seems already normalized. No normalization needed.
#+
#'## Perform differential expression

limmaFit <- lmFit(data, design=model)
limmaFit <- eBayes(limmaFit)


limmaRes <- list() 
for(coefx in colnames(coef(limmaFit))){
    limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) |> 
        rownames_to_column("gene")
}

limmaRes <- bind_rows(limmaRes, .id = "coef") 
limmaRes <- filter(limmaRes, coef != "(Intercept)") 

#'# Data interpretation
#'## Vulcano plot

threshold <- abs(limmaRes$logFC) > 2 & limmaRes$P.Value < 0.05 # set threshold
#+ plot7, fig.width=10
limmaRes |> 
    mutate(coef = str_replace(coef, "disease.staging", "")) |> 
    ggplot(aes(x = logFC, y= -log10(P.Value), color = threshold)) +
    geom_point(alpha = 0.2) +
    scale_color_brewer(palette = "Set2") +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(0.05)) + 
    geom_vline(xintercept = c(-2, 2)) +
    facet_wrap(coef~.)

#' The vulcano plot reveals that there are no big difference between stage I and 
#' stage II. However, with later disease stages the difference increases.

#'## P-value distribution
#+ plot8, fig.width=10
ggplot(limmaRes, aes(x=P.Value, fill=factor(floor(AveExpr)))) + 
    geom_histogram() +
    facet_wrap(coef~.)
#' The p.value distribution seems good.
#+
#'## Number of hits  
#'Genes tested per disease stage and total number of tests performed

limmaRes |> dplyr::count(coef)
dim(limmaRes) 

#'## Filter significant genes
#'Because the `p.value` distribution is good no cut off is set.

limmaRes |> 
    filter(adj.P.Val < 0.01) |> # strict cut off
    count(coef)

limmaResSig <- limmaRes |> 
    filter(adj.P.Val < 0.01)

#'# Visualizing results
#'## Visualizing one gene

Genename <- limmaResSig |> 
    filter(logFC == max(logFC)) |> # getting gene with the  highest logFC
    pull(gene)

data <- as.matrix(data) # convert data to matrix for later functions to work properly
#+ plot9, fig.width=10
design |> 
    mutate(upregGene = data[Genename,]) |> 
    ggplot(aes(x=disease.staging, y=upregGene)) +
    geom_point() +
    ggtitle(Genename) +
    ylab("E")

limmaRes |> 
    filter(gene == Genename) |> 
    pull(logFC, coef)

#' `logFC` are not that high for stage II and stage III
#+
#'## Visualizing multiple genes
#' Get the top 5 significant genes by logFC for each disease stage


goi.all <- limmaResSig |> 
    group_by(coef) |> 
    slice_max(logFC, n=5) |> # looking just at the upregulated
    pull(gene)

#+ plot10, fig.width=10
(p.coef <- limmaRes |> 
        filter(gene %in% goi.all) |> 
        ggplot(aes(y=gene, x=str_remove(coef, "disease.staging"), color=logFC, size=-log10(adj.P.Val))) + 
        geom_point() +
        scale_color_gradient2(high="red", low="blue") +
        xlab("disease.staging") +
        theme_bw())

dat.list <- list()
for(gg in goi.all){
    dat.list[[gg]] <- design |> 
        mutate(E=scale(data[gg,])) |> 
        rownames_to_column("sample") |> 
        remove_rownames()
}
#+ plot11, fig.width=10
(p.vals <- bind_rows(dat.list, .id="gene") |> 
        ggplot(aes(x=sample, y=gene, fill=E)) + 
        geom_tile() +
        facet_grid(. ~ disease.staging, space ="free", scales = "free") +
        scale_fill_gradient2(low="blue", high="red") +
        theme(axis.text.x=element_text(angle=90,hjust=1)))

#' Final plot. 
#+ plot12, fig.width=15
p.vals + p.coef
#' Interesting that for some samples particular genes are extremly upregulated in stage II.
#'# Enrichment analysis
#' Getting databases
dbs <- listEnrichrDbs()
databases <-  c("MSigDB_Hallmark_2020")
#'## Upregulated genes


enr.res.list <- list()
for(coefx in unique(limmaResSig$coef)){
    
    # Extract genes of interests (GOI) for a given coefficient (see yesterday's example)
    goi.E <-  limmaResSig |> 
        filter(coef == coefx & logFC > 1) |> 
        pull(gene) |> 
        unique()
    
    
    # Add code here to perform enrichment analysis (see yesterday's example)
    enr.res <- enrichr(goi.E, databases)
    
    # The results will be a list, where each entry is one database. We will combine those into one long table
    enr.res <- bind_rows(enr.res, .id="db")
    
    # Store results in the list
    enr.res.list[[coefx]] <- enr.res
}

enr.res.all <- bind_rows(enr.res.list, .id="coef")

#'### Plotting the results of the enrichment

filterTerm <- enr.res.all |> 
    filter(Adjusted.P.value < 0.05) |> 
    pull(Term)
#+ plot13, fig.width=10
enr.res.all |> 
    filter(Term %in% filterTerm) |> 
    ggplot(aes(x=str_remove(coef, "disease.staging"), y=Term, size = -log10(Adjusted.P.value), color = log(Odds.Ratio))) +
    geom_point() +
    scale_colour_gradient2(low='blue', mid = "white", high='red') +
    xlab("disease stage")




#'# Discussion
#' Particular genes get upregulated with disease stage.
#' The enrichment analysis proposed an enrichment in the G2-M Checkpoint and E2F-Targets 
#' and others including immun activity. However, the most differences was 
#' observed for the phenotypes, that almost cluster in the MDS projection and the hierachical clustering. 
#+
#'# Comparing disease stages in different phenotypes.
#' Just comparing disease stage I with disease stage IV
#'## Subsetting the data

sub.design <- design |> 
    filter(disease.staging %in% c("stage_I", "stage_IV")) |> 
    mutate(disease.staging = factor(disease.staging, ordered = FALSE)) |> 
    mutate(disease.staging = relevel(disease.staging, "stage_I")) |> 
    mutate(phenotype = factor(phenotype, ordered = FALSE)) |> 
    mutate(phenotype = relevel(phenotype, "low_malignant_potential")) # relevel to make low_mal control

sub.data <- data[,row.names(sub.design)] 

stopifnot(colnames(sub.data) == row.names(sub.design)) # checking if columns and rows align

dim(sub.data)
dim(sub.design)


sub.design |> 
    count(phenotype, disease.staging)
    
#' The availability of the samples, might cause problems calculating the model. There is only one sample 
#' in the coefficient: `low_malignant_potential + stage_I`.
#+
#'## Model matrix
model <- model.matrix(~disease.staging*phenotype, data=sub.design)
#+ plot14, fig.width=10
Heatmap(model)

#'## Perform differential expression

limmaFit <- lmFit(sub.data, design=model)
limmaFit <- eBayes(limmaFit)

limmaRes <- list() 
for(coefx in colnames(coef(limmaFit))){
    limmaRes[[coefx]] <- topTable(limmaFit, coef=coefx,number = Inf) |> 
        rownames_to_column("gene")
}

limmaRes <- bind_rows(limmaRes, .id = "coef") 
limmaRes <- filter(limmaRes, coef != "(Intercept)") 

#'## Fit contrast
#' creating contrats matrix
contrast.mt <- cbind(malignent_stage_IV = c(0,1,0,1))
row.names(contrast.mt) <- colnames(coef(limmaFit))

#  fitting contrast
limmaFit.contrast <- contrasts.fit(limmaFit,contrast.mt)
limmaFit.contrast <- eBayes(limmaFit.contrast)

# extracting results
limmaRes.contrast <- topTable(limmaFit.contrast, coef=colnames(contrast.mt),number = Inf) |> 
    rownames_to_column("gene") |> 
    mutate(coef=colnames(contrast.mt))

# add to rest of table
limmaRes <- rbind(limmaRes.contrast, limmaRes)
table(limmaRes$coef)

#' A little bit of clean up

limmaRes <- limmaRes |> 
    mutate(coef = str_replace(coef, "disease.stagingstage_IV:phenotypemalignant$", "interaction")) |> 
    mutate(coef = str_replace(coef, "disease.staging", "")) |> 
    mutate(coef = str_replace(coef, "phenotype", "")) 

#'# Data interpretation

#'## Vulcano plot
#+ plot15, fig.width=10
threshold <- abs(limmaRes$logFC) > 2 & limmaRes$P.Value < 0.05 # set threshold
ggplot(limmaRes, aes(x = logFC, y= -log10(P.Value), color = threshold)) +
    geom_point(alpha = 0.2) +
    scale_color_brewer(palette = "Set2") +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(0.05)) + 
    geom_vline(xintercept = c(-2, 2)) +
    facet_wrap(coef~.)

#'## P-value distribution
#+ plot16, fig.width=10
ggplot(limmaRes, aes(x=P.Value, fill=factor(floor(AveExpr)))) + 
    geom_histogram() +
    facet_wrap(coef~.)

#' p.value distribution for interaction and low_malignant_potential_stage_IV (in the graph only stage_IV) 
#' term is off, maybe duo to the sample size. There is only one sample for low_malignant_potential_stage_IV.
#' This also might effect the interaction term. Furthermore, there are no significant genes. This might also be 
#' a problem of the low sample size.
#+
#'## Number of hits

limmaRes |> dplyr::count(coef)


limmaRes |> 
    filter(adj.P.Val < 0.01) |> 
    count(coef)


#' No significant genes. This might also be duo to the sample distribution. 
