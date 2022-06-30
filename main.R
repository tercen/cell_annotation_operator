suppressPackageStartupMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
})

ctx = tercenCtx()

if(length(ctx$rnames) != 1) stop("Only one row factor must be projected.")

annotation_level <- ctx$op.value("annotation_level", as.character, "1")

rvals <- ctx$rselect()[[1]]

if(any(grepl("Comp|::", rvals))) {
  rvals <- unlist(lapply(
    strsplit(rvals, ":|,| "),
    tail,
    1
  ))
}

annot <- read.csv("./default_annotation.csv", header = TRUE)

annot <- annot %>% filter(level == annotation_level)
pos_list <- strsplit(gsub(" ", "", annot$high_markers), "_")
neg_list <- strsplit(gsub(" ", "", annot$low_markers), "_")

# splitted <- strsplit(annot$markers, "(?<=[+-])", perl = TRUE)
# pos_list <- lapply(splitted, function(x) gsub("[+]", "", x[grep("[+]", x)]))
# neg_list <- lapply(splitted, function(x) gsub("[-]", "", x[grep("[-]", x)]))

rvals_in <- rvals[rvals %in% unlist(c(pos_list, neg_list))]

if(length(rvals_in) == 0) stop("No markers found in the annotation. Please check the marker names.")
msg <- paste0("The following markers have been found in the annotation: ", paste0(rvals_in, collapse = ", "))
ctx$log(msg)

mat <- ctx$as.matrix()
rownames(mat) <- rvals

probas <- apply(mat, 2, function(m) {
  sum_z_pos <- unlist(lapply(pos_list, function (x) { sum(m[x], na.rm = TRUE) }))
  sum_z_neg <- unlist(lapply(neg_list, function (x) { sum(-m[x], na.rm = TRUE) }))
  sum_z <- sum_z_pos + sum_z_neg
  sum_z[abs(sum_z) < 0] <- 0
  pv <- pnorm(sum_z, lower.tail = TRUE)
  prob <- round(pv / sum(pv), 4)
  prob
})

rownames(probas) <- annot$population
colnames(probas) <- as.character(1:ncol(probas))

df_out <- probas %>% 
  as_tibble() %>%
  mutate(population = rownames(probas)) %>%
  pivot_longer(names_to = ".ci", values_to = "prob", -population) %>%
  mutate(.ci = as.integer(.ci) - 1)

df_out %>%
  ctx$addNamespace() %>%
  ctx$save()
