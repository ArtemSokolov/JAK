library(tidyverse)
library(synapser)
library(synExtra)
library(pheatmap)

suppressMessages(synLogin())

syn <- synExtra::synDownloader("/data/JAK/DGE", ifcollision="overwrite.local")
read_csv2 <- partial(read_csv, col_types=cols())

## Download transcript counts and metadata
x <- syn("syn20820769") %>% read_csv2()
y <- syn("syn20820771") %>% read_csv2()

## Cluster all the stressors and controls
vStress <- c( "LPS", "dsRNA +lipo", "naked dsRNA", "Lipo control", "Drug control" )
    
## Replaces 0s with NAs
zero2na <- function( v ) { ifelse(v == 0, NA, v) }
    
## Identify the wells of interest
## Isolate genes that were detected in all wells of interest
y1 <- y %>% filter(Drug %in% vStress) %>% select(Well, Drug)
x1 <- x %>% select(HUGO, y1$Well) %>% mutate_all( zero2na ) %>% na.omit

## Compute pairwise correlations
rr <- x1 %>% as.data.frame %>% column_to_rownames("HUGO") %>%
    cor( method="sp", use="pairwise.complete.obs" )
ry <- y1 %>% as.data.frame %>% #mutate( Drug = as.factor(Drug) ) %>%
    column_to_rownames("Well")

## Set up a palette
pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2") %>%
    setNames(vStress) %>% list(Drug = .)
gg <- pheatmap(rr, annotation_row = ry, annotation_col = ry,
               cellwidth = 12, cellheight = 12,
               annotation_colors = pal, silent=TRUE)

ggsave("dsRNA-cluster.pdf", gg, width=9, height=9)
