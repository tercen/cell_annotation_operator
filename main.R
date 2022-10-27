suppressPackageStartupMessages({
  library(tercen)
  library(dplyr, warn.conflicts = FALSE)
  library(tidyr)
})

#options("tercen.workflowId" = "75cc2bda2d854b35362e4f45150136b9")
#options("tercen.stepId"     = "5ae79b70-5195-4578-8439-eaf21b9a250b")

ctx = tercenCtx()

if(length(ctx$rnames) != 1) stop("Only one row factor must be projected.")

pthres <- ctx$op.value('P-value threshold', as.double, 0.05)

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

annot <- as.data.frame(annot_tbl)
pos_list <- strsplit(gsub(" ", "", annot$pos_markers), "_")
neg_list <- strsplit(gsub(" ", "", annot$neg_markers), "_")

rvals_in <- rvals[rvals %in% unlist(c(pos_list, neg_list))]

if(length(rvals_in) == 0) warning("No markers found in the annotation. Please check the marker names.")
msg <- paste0("The following markers have been found in the annotation: ", paste0(rvals_in, collapse = ", "))
ctx$log(msg)

mat <- ctx$as.matrix()
rownames(mat) <- rvals

p_val <- apply(mat, 2, function(m) {
  sum_z_pos <- unlist(lapply(pos_list, function (x) { sum(m[x], na.rm = TRUE) }))
  sum_z_neg <- unlist(lapply(neg_list, function (x) { sum(-m[x], na.rm = TRUE) }))
  sum_z <- sum_z_pos + sum_z_neg
  sum_z[abs(sum_z) < 0] <- 0
  pv <- 1-pnorm(sum_z, lower.tail = TRUE)
  pv
})
rownames(p_val) <- annot$population
colnames(p_val) <- as.character(1:ncol(p_val))

probas <- apply(p_val, 2, function(n) {
  prob <- round(n / sum(n), 4)
  prob
})
rownames(probas) <- annot$population
colnames(probas) <- as.character(1:ncol(probas))


apply(t(p_val), 2, function(x) which.max(x))
tpv<-t(p_val)

table_pv<-cbind(rownames(tpv),"Unknow")
for (cluster in rownames(tpv)){
  
  res<-which( tpv[cluster,] == min(tpv[cluster,]))
  if (min(tpv[cluster,])[1]<pthres){
    table_pv[[as.integer(cluster),2]]<-paste(colnames(tpv)[res], collapse = '_')
  }
}
colnames(table_pv)<-c(".ci","max_population")

####output
pval_out <- p_val %>% 
  as_tibble() %>%
  mutate(population = rownames(p_val)) %>%
  pivot_longer(names_to = ".ci", values_to = "pv", -population) %>%
  mutate(.ci = as.integer(.ci) - 1L) %>%
  ctx$addNamespace() 

prob_out <- probas %>% 
  as_tibble() %>%
  mutate(population = rownames(probas)) %>%
  pivot_longer(names_to = ".ci", values_to = "prob", -population) %>%
  mutate(.ci = as.integer(.ci) - 1L) %>%
  ctx$addNamespace()

tbl_pv_out <- table_pv %>% 
  as_tibble() %>%
  mutate(.ci = as.integer(.ci) - 1L) %>%
  ctx$addNamespace() %>%
  as_relation() 

df_out<-merge(prob_out,pval_out)

join_res <- df_out %>%
  left_join_relation(ctx$crelation, ".ci", ctx$crelation$rids) %>%
  left_join_relation(tbl_pv_out, list(), list()) %>%
  as_join_operator(ctx$cnames, ctx$cnames)%>%
  save_relation(ctx)

join_res %>%
  ctx$save()
