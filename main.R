suppressPackageStartupMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
})

# http://127.0.0.1:5400/test/w/0d1335bc08a28da4ef5f89545f004c3d/ds/44c1b2c1-8b13-478a-bd98-a9ce43be196f
# options("tercen.workflowId" = "0d1335bc08a28da4ef5f89545f004c3d")
# options("tercen.stepId"     = "44c1b2c1-8b13-478a-bd98-a9ce43be196f")

ctx = tercenCtx()

if(length(ctx$rnames) != 1) stop("Only one row factor must be projected.")

pthres <- ctx$op.value('P-value threshold', as.double, 0.05)

rvals <- ctx$rselect()[[1]]

### Extract the population table with the documentId
annot <- ctx$client$tableSchemaService$select(
  ctx$select(ctx$labels[[1]], nr = 1)[[1]]
) %>%
  as_tibble()

pos_list <- strsplit(gsub(" ", "", annot$pos_markers), "_")
neg_list <- strsplit(gsub(" ", "", annot$neg_markers), "_")

# marker specificity
full_list <- lapply(seq_along(pos_list), function(x) c(pos_list[[x]], neg_list[[x]]))
unique_markers <- unique(unlist(c(pos_list, neg_list)))
ct <- sapply(unique_markers, function(m) {
  ct_p <- sum(unlist(lapply(pos_list, function(x) m %in% x)))
  ct_n <- sum(unlist(lapply(neg_list, function(x) m %in% x)))
  c(ct_p, ct_n)
})
specs <- apply(ct, 1, function(x) {
  1 - (x - min(x)) / (length(pos_list) - min(x))
})

rvals_in <- rvals[rvals %in% unlist(c(pos_list, neg_list))]

if(length(rvals_in) == 0) {
  stop("No markers found in the annotation. Please check the marker names.")
}
ctx$log(paste0(
  "The following markers have been found in the annotation: ",
  paste0(rvals_in, collapse = ", ")
))

### Probabilities computation
mat <- ctx$as.matrix()
rownames(mat) <- rvals

mat[mat == -20] <- min(mat[mat > -20])
mat[mat == 20] <- max(mat[mat < 20])

probs <- apply(mat, 2, function(m) {
  sum_z_pos <- unlist(lapply(pos_list, function (x) { sum(m[x] * specs[x, 1] / sqrt(length(x)), na.rm = TRUE) }))
  sum_z_neg <- unlist(lapply(neg_list, function (x) { sum(m[x] * specs[x, 2] / sqrt(length(x)), na.rm = TRUE) }))
  sum_z <- sum_z_pos - sum_z_neg
  sum_z[abs(sum_z) < 0] <- 0
  pnorm(sum_z, lower.tail = FALSE, log.p = TRUE)
})

rownames(probs) <- annot$population
colnames(probs) <- as.character(1:ncol(probs))

p_values <- exp(probs)

probas <- apply(-probs, 2, function(n) {
  prob <- round(n / sum(n), 4)
  prob
})

## Output
# 1 Probability table
df_prob <- probas %>%
  as_tibble(rownames = "pop.name") %>%
  tidyr::pivot_longer(cols = !matches("pop.name"), names_to = ".ci", values_to = "best.pop.probability") %>%
  mutate(.ci = as.integer(.ci) - 1L)

df_probs <- p_values %>%
  as_tibble(rownames = "pop.name") %>%
  tidyr::pivot_longer(cols = !matches("pop.name"), names_to = ".ci", values_to = "best.pop.pvalue") %>%
  mutate(.ci = as.integer(.ci) - 1L) %>%
  merge(df_prob) %>% 
  mutate("best.pop.-log(pvalue)" = -log10(best.pop.pvalue))

# 1 Cluster annotation table

probable_pop <- df_probs %>% 
  group_by(.ci, .drop = FALSE) %>%
  filter(best.pop.pvalue < pthres) %>%
  summarise(across(pop.name, paste, collapse = ", ")) %>%
  mutate(.ci = as.integer(.ci)) %>%
  arrange(.ci) %>%
  rename(pop.list = pop.name)

df_clust_annot <- df_probs %>% 
  group_by(.ci) %>%
  # filter(neglog_p_value == max(neglog_p_value)) %>%
  filter(`best.pop.-log(pvalue)` == max(`best.pop.-log(pvalue)`)) %>%
  arrange(.ci) %>%
  rename(max_pop = pop.name) %>%
  ungroup() %>%
  left_join(probable_pop, ".ci") %>%
  mutate(probable_pop = if_else(pop.list == "", "Unknown", pop.list)) %>%
  ctx$addNamespace() 
  

### Return results to tercen
df_probs_out <- df_probs %>%
  # rename(per_pop_pval = p_value, per_pop_pob = best.pop.probability) %>%
  select(-`best.pop.-log(pvalue)`) %>%
  ctx$addNamespace()

ctx$save(list(df_clust_annot, df_probs_out))
