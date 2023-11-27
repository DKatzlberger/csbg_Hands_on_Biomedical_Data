#' ---
#' title: "Day 1 Hands-on-Biomedical-Data"
#' author: Daniel Katzlberger
#' output: github_document
#' ---
#' The output is suppressed for all the not so interesting calls
#' 
#'# Basics
#' Load required packages
#' 
require(tidyverse)
require(pheatmap)
require(gemma.R)

#'## Variables
#+ eval=FALSE
a <- 1
b <- 3
typeof(a)
typeof(b)
str(a)
str(b)
a
b

#'## Functions
#+ eval=FALSE
sum(a, b)
a + b
sum(5,6)
c <- sum(a, b)
(c <- sum(a, b))
str(c)
sum(a,c)

#+ error=TRUE
a <- "abc"
b <- "hello"
sum(a, b)
#+ eval=FALSE
typeof(a)
typeof(b)
str(a)
str(b)

a <- TRUE
b <- FALSE
sum(b, b)
sum(a, a)
typeof(a)
typeof(b)
str(a)
str(b)

#'## Vectors
#+ eval=FALSE
vec <- c(1,2,3,4,1,2,3)
names(vec) <- LETTERS[1:4]
str(vec)
sum(vec)
mean(vec)
vec[2]
vec[3]
vec[6] # exceeding the range of the names, there is no value names
vec["A"]
vec[c("A", "D")]
unique(vec)
unique(names(vec)) # each index with no name is assigned <NA>
1:10
-5:5
vec <- 10:1
vec[2:3]


#' Creating my own vector with the letters `B`, `C`, `D`
#+ eval=TRUE

?LETTERS
your.vector <- LETTERS[2:4]
print(your.vector)
stopifnot(your.vector == c("B", "C", "D"))

#'## Lists
#+ eval=FALSE
list_x <- list("a", 1, "b", "xyz", TRUE)
str(list_x)
list_x[1]
str(list_x[[1]])
str(list_x[1])
list_x[[2]]
list_x[2:4]
str(list_x[2:4])
vec_x <- c("a", 1, "b", "xyz", TRUE)
str(vec_x)

#'## Loops and conditions
#+ eval=FALSE
for(x in 1:6){
    print(paste("x =", x)) # this is like f string in python
    if(x > 3){
        print("...x is greater than three")
    }
    if(x == 2){
        print("...x equals two")
    }
    if(x != 4){
        print("...x is not four")
    }
    if(x %% 2 == 0){
        print("...x is even")
    } else {
        print("...x is uneven")
    }
    if(x %in% c(3,5)){
        print("...x is three or five")
    }
    print("-------")
}
#'# Visualizations
#'## Matrices
#' Load the data into R
m <- readRDS("data.RDS")

#' Summarize the matrix
#+ eval=FALSE
str(m)
dim(m)
head(m)
row.names(m)
colnames(m)
1:20
dim(m)
dim(m[1:20,])

#' Get the colnames with `"Liver_Fibroblasts"` in the name
#+ eval=FALSE
colnames(m)
grepl("Liver_Fibroblasts", colnames(m))
dim(m)
dim(m[, grepl("Liver_Fibroblasts", colnames(m))])



#' Subset first 20 rows of `m` containing `"Liver_Fibroblasts"` and rename cols and rows
m <- m[1:20, grepl("Liver_Fibroblasts", colnames(m))]
colnames(m) <- gsub("^Liver_Fibroblasts_(.+)_RNA_(\\d)$", "\\1_\\2", colnames(m))
print(colnames(m))

row.names(m) <- paste0("g", 1:nrow(m))
print(row.names(m))

#' `m`looks now like this
print(m)
#' It is possible to transpose the matrix `t()`
#+ eval=FALSE
t(m)

#' Now we correlate the genes
#+ eval=TRUE
cor(m, method="spearman")

p <- pheatmap(cor(m, method="spearman"))
p

#' Here the diagonal is set to `<NA>`
cMT <- cor(m, method="spearman")
diag(cMT) <- NA

p <- pheatmap(cMT)
p

#' Subsetting the matrix for 30 rows and first 10 columns
m1 <- readRDS("data.RDS")
#+ eval=FALSE
m1 <- m1[1:40, grepl("Liver_Fibroblasts", colnames(m1))] # i take 40 rows to later subset 30
colnames(m1) <- gsub("^Liver_Fibroblasts_(.+)_RNA_(\\d)$", "\\1_\\2", colnames(m1))
row.names(m1) <- paste0("g", 1:nrow(m1))
#+ eval=TRUE
dim(m1)

matrix.check <- m1[1:30,1:10]
stopifnot(all(dim(matrix.check) == c(30,10)))
dim(matrix.check)

#'## Data frames and dplyr
#+ eval=FALSE
starwars
str(starwars$name)
#' Check uniqueness of names in the data
stopifnot(length(unique(starwars$name))==nrow(starwars))


sw <- starwars |> 
    select(where(function(x) !is.list(x))) |> 
    as.data.frame()

sw |> count(homeworld) 

#' Add column with firstname
sw <- sw |> 
    mutate(firstname = str_remove(name, " .+$")) 

#'### Exercise 1.1
sw |> count(gender)

sw |> count(skin_color)

sw |> 
    pull("skin_color") |> 
    str_split(", ") |> 
    unlist() |> 
    table() |> 
    sort()

#' The difference is that the strings, that contain more than one color are split
#' into their individual colors. Therefore, not counted as an extra color.
#+
#'### Exercise 1.2

tallcharacters <- sw |> 
    filter(height > 200) |> 
    pull(name)

#'### Exercise 1.3
#'Pipe 1
start.time <- Sys.time()
starwars |> 
    filter(hair_color == "blond") |> 
    filter(sex == "male")
end.time <- Sys.time()
time.taken1 <- end.time - start.time

#' The function filter is called twice in this pipe. <br>
#' Pipe 2
start.time <- Sys.time()
starwars |> 
    filter(hair_color == "blond" & sex == "male")
end.time <- Sys.time()
time.taken2 <- end.time - start.time

#' Here the function filter is not called twice but within filter two conditions
#' are set and combined with `&`
#' meaning `hair_color` must be blonde and `sex` must be male.

print(c(time.taken1,time.taken2))
#'Second pipe executes faster <br>
#'Pipe 3
starwars |> 
    filter(hair_color == "blond" | sex == "male")

#' Here the conditions are combined with `|` meaning or/oder. This pipe prints all 
#' rows either with the `hair_color`
#' blond or with `sex` male.

#'### Exercise 1.4
sw <- sw |> 
    mutate(bmi = mass/(height/100)^2)
#' Look at the quartiles of the `bmi`column

quantile(sw$bmi, na.rm = TRUE)

# data.frames are lists, thus the same syntax with two squared brackets works
quantile(sw[["bmi"]], na.rm = TRUE)

# the function 'pull' can also extract columns
quantile(pull(sw, "bmi"), na.rm = TRUE)

#'## Plotting
ggplot(sw, aes(x=mass, y=height)) + geom_point(na.rm = TRUE)
px <- ggplot(sw, aes(x=mass, y=height))



#'### Exercise 1.5
px + geom_hex()
px + geom_point(color = "red", shape = 1) + 
    geom_text(data = sw |> filter( mass > 1000), aes(label = firstname))

ggplot(sw, aes(x=height)) + geom_histogram()
ggplot(sw, aes(x=height)) + geom_density() + stat_ecdf()


sw |> 
    group_by(gender) |> 
    top_n(3, bmi) |> 
    ggplot(aes(x=fct_reorder(firstname, bmi), y=bmi, fill=gender)) + 
    geom_bar(stat="identity") +
    facet_grid(rows = vars(gender), space = "free_x", scales = "free_y")

#'## Factors
ggplot(sw, aes(x=sex,y =height)) + geom_violin()

sw |> 
    filter(sex %in% c("male", "female")) |> 
    ggplot(aes(x=sex,y =height)) + geom_violin()

str(sw$sex)

# Here we just look at the column 'sex' and see that it is a vector of characters.
sw |> 
    filter(sex %in% c("male", "female")) |> 
    pull(sex) |> 
    str()

# Now it is converted to a factor with the levels 'female' and then 'male'
sw |> 
    filter(sex %in% c("male", "female")) |> 
    mutate(sex = factor(sex)) |> 
    pull(sex) |> 
    str()

# You can choose the order of levels, which has effects on plots and other 
# functions (design matrices)
sw |> 
    filter(sex %in% c("male", "female")) |> 
    mutate(sex = factor(sex, levels=c("male", "female"))) |> 
    pull(sex) |> 
    str()

#'### Exercise 1.6
sw |> 
    filter(sex %in% c("male", "female")) |> 
    mutate(sex = factor(sex, levels=c("male", "female"))) |> 
    ggplot(aes(x=sex,y =height)) + geom_violin()


#'## GEMMA
#'### Exercise 1.7
#+ eval=FALSE
get_datasets("HIV1", limit = 100, taxa = "human") |> 
    filter(geeq.batchCorrected == TRUE) |> 
    select(taxon.Name, taxon.ID, experiment.Accession, experiment.SampleCount)

gse <- "GSE21589"

get_datasets(gse) |> 
    select(experiment.ShortName, experiment.Name, experiment.ID, experiment.Description)

d <- get_dataset_design(gse)
str(d)
head(d)
with(d, table(batch))
row.names(d)

e <- get_dataset_processed_expression(gse)
str(e)
colnames(e)
e <- as.data.frame(e)

dataMT <- as.matrix(e[,row.names(d)])
str(dataMT)
row.names(dataMT) <- e$GeneSymbol

boxplot(dataMT, las=2)



