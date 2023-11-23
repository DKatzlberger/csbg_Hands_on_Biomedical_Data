day1
================
Daniel Katzlberger
2023-11-23

# Basics

Load required packages

``` r
require(tidyverse)
require(pheatmap)
```

Variables

``` r
a <- 1
b <- 3
typeof(a)
```

    ## [1] "double"

``` r
typeof(b)
```

    ## [1] "double"

``` r
str(a)
```

    ##  num 1

``` r
str(b)
```

    ##  num 3

``` r
a
```

    ## [1] 1

``` r
b
```

    ## [1] 3

Functions

``` r
sum(a, b)
```

    ## [1] 4

``` r
a + b
```

    ## [1] 4

``` r
sum(5,6)
```

    ## [1] 11

``` r
c <- sum(a, b)
(c <- sum(a, b))
```

    ## [1] 4

``` r
str(c)
```

    ##  num 4

``` r
sum(a,c)
```

    ## [1] 5

``` r
a <- "abc"
b <- "hello"
sum(a, b)
```

    ## Error in sum(a, b): invalid 'type' (character) of argument

``` r
typeof(a)
```

    ## [1] "character"

``` r
typeof(b)
```

    ## [1] "character"

``` r
str(a)
```

    ##  chr "abc"

``` r
str(b)
```

    ##  chr "hello"

``` r
a <- TRUE
b <- FALSE
sum(b, b)
```

    ## [1] 0

``` r
sum(a, a)
```

    ## [1] 2

``` r
typeof(a)
```

    ## [1] "logical"

``` r
typeof(b)
```

    ## [1] "logical"

``` r
str(a)
```

    ##  logi TRUE

``` r
str(b)
```

    ##  logi FALSE

Vectors

``` r
vec <- c(1,2,3,4,1,2,3)
names(vec) <- LETTERS[1:4]
str(vec)
```

    ##  Named num [1:7] 1 2 3 4 1 2 3
    ##  - attr(*, "names")= chr [1:7] "A" "B" "C" "D" ...

``` r
sum(vec)
```

    ## [1] 16

``` r
mean(vec)
```

    ## [1] 2.285714

``` r
vec[2]
```

    ## B 
    ## 2

``` r
vec[3]
```

    ## C 
    ## 3

``` r
vec[6] # exceeding the range of the letters, there is no value names
```

    ## <NA> 
    ##    2

``` r
vec["A"]
```

    ## A 
    ## 1

``` r
vec[c("A", "D")]
```

    ## A D 
    ## 1 4

``` r
unique(vec)
```

    ## [1] 1 2 3 4

``` r
unique(names(vec)) # each index with no letter is assigned <NA>
```

    ## [1] "A" "B" "C" "D" NA

``` r
1:10
```

    ##  [1]  1  2  3  4  5  6  7  8  9 10

``` r
-5:5
```

    ##  [1] -5 -4 -3 -2 -1  0  1  2  3  4  5

``` r
vec <- 10:1
vec[2:3]
```

    ## [1] 9 8

Load the data into R

``` r
m <- readRDS("data.RDS")
```

Subset first 20 rows of `m` containing `"Liver_Fibroblasts"`

``` r
m <- m[1:20, grepl("Liver_Fibroblasts", colnames(m))]
print(dim(m))
```

    ## [1] 20 20

Get the first 30 rows and the first 10 columns of `m`

``` r
print(gsub("^Liver_Fibroblasts_(.+)_RNA_(\\d)$", "\\1_\\2", colnames(m)))
```

    ##  [1] "IL6_2"  "TGFb_3" "TNF_1"  "IFNa_2" "PBS_1"  "IFNa_3" "TGFb_1" "IL3_1"  "IL3_3" 
    ## [10] "IL6_1"  "PBS_3"  "IL3_2"  "IL6_3"  "IFNa_1" "TGFb_2" "TNF_2"  "PBS_2"  "IFNg_2"
    ## [19] "IFNg_3" "IFNg_1"
