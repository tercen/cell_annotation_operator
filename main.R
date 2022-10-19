suppressPackageStartupMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
})

#options("tercen.workflowId" = "75cc2bda2d854b35362e4f45150136b9")
#options("tercen.stepId"     = "5ae79b70-5195-4578-8439-eaf21b9a250b")

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
### Extract the population table with the documentId
doc.id.tmp<-as_tibble(ctx$select())
doc.id<-doc.id.tmp[[grep("documentId" , colnames(doc.id.tmp))]][1]

annot_tbl<-ctx$client$tableSchemaService$select(doc.id) %>%
  as_tibble()

annot_df <- as.data.frame(annot_tbl)
#annot_df <- read.csv("./default_annotation.csv", header = TRUE)

annot <- annot_df #%>% filter(level == annotation_level)
pos_list <- strsplit(gsub(" ", "", annot$high_markers), "_")
neg_list <- strsplit(gsub(" ", "", annot$low_markers), "_")

# splitted <- strsplit(annot$markers, "(?<=[+-])", perl = TRUE)
# pos_list <- lapply(splitted, function(x) gsub("[+]", "", x[grep("[+]", x)]))
# neg_list <- lapply(splitted, function(x) gsub("[-]", "", x[grep("[-]", x)]))

rvals_in <- rvals[rvals %in% unlist(c(pos_list, neg_list))]

if(length(rvals_in) == 0) warning("No markers found in the annotation. Please check the marker names.")
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
  mutate(.ci = as.integer(.ci) - 1L) %>%
  ctx$addNamespace() 

df_out %>%
  ctx$save()
