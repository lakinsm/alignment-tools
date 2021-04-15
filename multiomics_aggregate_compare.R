#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(omicsPrint)
library(vegan)
library(grid)
library(tsne)

set.seed(2718)
options(echo=FALSE)

# 1. kraken_analytic_matrix.csv, 2. host_vcf_table.tsv, 3. viral_vcf_table.tsv, 4. target_sample_name, 5. out_dir
args = commandArgs(trailingOnly = TRUE)
kraken_fpath = args[1]
host_var_fpath = args[2]
viral_var_fpath = args[3]
target = args[4]
output_dir = args[5]

create_working_directories = function(root_dir) 
{
  # If subdirs for stats and exploratory variables don't exist, create them
  ifelse(!dir.exists(file.path(root_dir)),
         dir.create(file.path(root_dir), mode='777'), FALSE)
  ifelse(!dir.exists(file.path(root_dir, 'microbiome_matrices')),
         dir.create(file.path(root_dir, 'microbiome_matrices'), mode='777'), FALSE)
  ifelse(!dir.exists(file.path(root_dir, 'microbiome_simplex_matrices')),
         dir.create(file.path(root_dir, 'microbiome_simplex_matrices'), mode='777'), FALSE)
  ifelse(!dir.exists(file.path(root_dir, 'reports')),
         dir.create(file.path(root_dir, 'reports'), mode='777'), FALSE)
  ifelse(!dir.exists(file.path(root_dir, 'graphs')),
         dir.create(file.path(root_dir, 'graphs'), mode='777'), FALSE)
}

create_working_directories(output_dir)

eps = 1

kraken = read.csv(kraken_fpath, header=T)
kmat = data.table(kraken)

kmat[, id := rownames(kraken)]
kraken_taxonomy = data.table(id=rownames(kraken))
setDT(kraken_taxonomy)[, c('Domain',
                           'Kingdom',
                           'Phylum',
                           'Class',
                           'Order',
                           'Family',
                           'Genus',
                           'Species') := tstrsplit(id, '|', type.convert=T, fixed=T)]
setkey(kraken_taxonomy, id)
setkey(kmat, id)

kall = kraken_taxonomy[kmat]
kall[, c('id', 'Kingdom') := NULL]
kall_rsums = rowSums(kall[, .SD, .SDcols=!1:7])
kall = kall[kall_rsums >= quantile(kall_rsums, 0.05)]

kdomain = kall[!is.na(Domain) & Domain != 'NA', lapply(.SD, sum), by='Domain', .SDcols=!1:7]
kphylum = kall[!is.na(Phylum) & Phylum != 'NA', lapply(.SD, sum), by='Phylum', .SDcols=!1:7]
kclass = kall[!is.na(Class) & Class != 'NA', lapply(.SD, sum), by='Class', .SDcols=!1:7]
korder = kall[!is.na(Order) & Order != 'NA', lapply(.SD, sum), by='Order', .SDcols=!1:7]
kfamily = kall[!is.na(Family) & Family != 'NA', lapply(.SD, sum), by='Family', .SDcols=!1:7]
kgenus = kall[!is.na(Genus) & Genus != 'NA', lapply(.SD, sum), by='Genus', .SDcols=!1:7]
kspecies = kall[!is.na(Species) & Species != 'NA', lapply(.SD, sum), by='Species', .SDcols=!1:7]


### naive Bayes with Lakin & Abdo LOR IC

# Domain
kdomain_target = matrix(kdomain[[`target`]], nrow=1)
kdomain_simplex = copy(kdomain)
kdomain_simplex[, Domain := NULL]
kdomain_simplex[, eval(target) := NULL]
kdomain_simplex = as.matrix(kdomain_simplex)
kdomain_col_counts = colSums(kdomain_simplex)
kdomain_simplex = apply(kdomain_simplex, 2, function(x) (x+eps) / sum(x+eps))
kdomain_nb_max = rep(0L, ncol(kdomain_simplex))
for(i in 1:ncol(kdomain_simplex)) {
  kdomain_nb_max[i] = kdomain_col_counts[i] * (matrix(kdomain_simplex[, i], nrow=1) %*% log(matrix(kdomain_simplex[, i], ncol=1)))
}
kdomain_nb = kdomain_target %*% log(kdomain_simplex)

kdomain_nb_sort = sort(kdomain_nb, decreasing=TRUE)
kdomain_nb1_idx = which(kdomain_nb == kdomain_nb_sort[1])
kdomain_nb1 = kdomain_nb_sort[1]
kdomain_nb2_idx = which(kdomain_nb == kdomain_nb_sort[2])
kdomain_nb2 = kdomain_nb_sort[2]

kdomain_lor_ic = (kdomain_nb1 - kdomain_nb_max[kdomain_nb1_idx]) - (kdomain_nb2 - kdomain_nb_max[kdomain_nb2_idx])
kdomain_top = colnames(kdomain_simplex)[kdomain_nb1_idx]
kdomain_top_ic = colnames(kdomain_simplex)[ifelse(kdomain_lor_ic >= 0, kdomain_nb1_idx, kdomain_nb2_idx)]

# Phylum
kphylum_target = matrix(kphylum[[`target`]], nrow=1)
kphylum_simplex = copy(kphylum)
kphylum_simplex[, Phylum := NULL]
kphylum_simplex[, eval(target) := NULL]
kphylum_simplex = as.matrix(kphylum_simplex)
kphylum_col_counts = colSums(kphylum_simplex)
kphylum_simplex = apply(kphylum_simplex, 2, function(x) (x+eps) / sum(x+eps))
kphylum_nb_max = rep(0L, ncol(kphylum_simplex))
for(i in 1:ncol(kphylum_simplex)) {
  kphylum_nb_max[i] = kphylum_col_counts[i] * (matrix(kphylum_simplex[, i], nrow=1) %*% log(matrix(kphylum_simplex[, i], ncol=1)))
}
kphylum_nb = kphylum_target %*% log(kphylum_simplex)

kphylum_nb_sort = sort(kphylum_nb, decreasing=TRUE)
kphylum_nb1_idx = which(kphylum_nb == kphylum_nb_sort[1])
kphylum_nb1 = kphylum_nb_sort[1]
kphylum_nb2_idx = which(kphylum_nb == kphylum_nb_sort[2])
kphylum_nb2 = kphylum_nb_sort[2]

kphylum_lor_ic = (kphylum_nb1 - kphylum_nb_max[kphylum_nb1_idx]) - (kphylum_nb2 - kphylum_nb_max[kphylum_nb2_idx])
kphylum_top = colnames(kphylum_simplex)[kphylum_nb1_idx]
kphylum_top_ic = colnames(kphylum_simplex)[ifelse(kphylum_lor_ic >= 0, kphylum_nb1_idx, kphylum_nb2_idx)]

# Class
kclass_target = matrix(kclass[[`target`]], nrow=1)
kclass_simplex = copy(kclass)
kclass_simplex[, Class := NULL]
kclass_simplex[, eval(target) := NULL]
kclass_simplex = as.matrix(kclass_simplex)
kclass_col_counts = colSums(kclass_simplex)
kclass_simplex = apply(kclass_simplex, 2, function(x) (x+eps) / sum(x+eps))
kclass_nb_max = rep(0L, ncol(kclass_simplex))
for(i in 1:ncol(kclass_simplex)) {
  kclass_nb_max[i] = kclass_col_counts[i] * (matrix(kclass_simplex[, i], nrow=1) %*% log(matrix(kclass_simplex[, i], ncol=1)))
}
kclass_nb = kclass_target %*% log(kclass_simplex)

kclass_nb_sort = sort(kclass_nb, decreasing=TRUE)
kclass_nb1_idx = which(kclass_nb == kclass_nb_sort[1])
kclass_nb1 = kclass_nb_sort[1]
kclass_nb2_idx = which(kclass_nb == kclass_nb_sort[2])
kclass_nb2 = kclass_nb_sort[2]

kclass_lor_ic = (kclass_nb1 - kclass_nb_max[kclass_nb1_idx]) - (kclass_nb2 - kclass_nb_max[kclass_nb2_idx])
kclass_top = colnames(kclass_simplex)[kclass_nb1_idx]
kclass_top_ic = colnames(kclass_simplex)[ifelse(kclass_lor_ic >= 0, kclass_nb1_idx, kclass_nb2_idx)]

# Order
korder_target = matrix(korder[[`target`]], nrow=1)
korder_simplex = copy(korder)
korder_simplex[, Order := NULL]
korder_simplex[, eval(target) := NULL]
korder_simplex = as.matrix(korder_simplex)
korder_col_counts = colSums(korder_simplex)
korder_simplex = apply(korder_simplex, 2, function(x) (x+eps) / sum(x+eps))
korder_nb_max = rep(0L, ncol(korder_simplex))
for(i in 1:ncol(korder_simplex)) {
  korder_nb_max[i] = korder_col_counts[i] * (matrix(korder_simplex[, i], nrow=1) %*% log(matrix(korder_simplex[, i], ncol=1)))
}
korder_nb = korder_target %*% log(korder_simplex)

korder_nb_sort = sort(korder_nb, decreasing=TRUE)
korder_nb1_idx = which(korder_nb == korder_nb_sort[1])
korder_nb1 = korder_nb_sort[1]
korder_nb2_idx = which(korder_nb == korder_nb_sort[2])
korder_nb2 = korder_nb_sort[2]

korder_lor_ic = (korder_nb1 - korder_nb_max[korder_nb1_idx]) - (korder_nb2 - korder_nb_max[korder_nb2_idx])
korder_top = colnames(korder_simplex)[korder_nb1_idx]
korder_top_ic = colnames(korder_simplex)[ifelse(korder_lor_ic >= 0, korder_nb1_idx, korder_nb2_idx)]

# Family
kfamily_target = matrix(kfamily[[`target`]], nrow=1)
kfamily_simplex = copy(kfamily)
kfamily_simplex[, Family := NULL]
kfamily_simplex[, eval(target) := NULL]
kfamily_simplex = as.matrix(kfamily_simplex)
kfamily_col_counts = colSums(kfamily_simplex)
kfamily_simplex = apply(kfamily_simplex, 2, function(x) (x+eps) / sum(x+eps))
kfamily_nb_max = rep(0L, ncol(kfamily_simplex))
for(i in 1:ncol(kfamily_simplex)) {
  kfamily_nb_max[i] = kfamily_col_counts[i] * (matrix(kfamily_simplex[, i], nrow=1) %*% log(matrix(kfamily_simplex[, i], ncol=1)))
}
kfamily_nb = kfamily_target %*% log(kfamily_simplex)

kfamily_nb_sort = sort(kfamily_nb, decreasing=TRUE)
kfamily_nb1_idx = which(kfamily_nb == kfamily_nb_sort[1])
kfamily_nb1 = kfamily_nb_sort[1]
kfamily_nb2_idx = which(kfamily_nb == kfamily_nb_sort[2])
kfamily_nb2 = kfamily_nb_sort[2]

kfamily_lor_ic = (kfamily_nb1 - kfamily_nb_max[kfamily_nb1_idx]) - (kfamily_nb2 - kfamily_nb_max[kfamily_nb2_idx])
kfamily_top = colnames(kfamily_simplex)[kfamily_nb1_idx]
kfamily_top_ic = colnames(kfamily_simplex)[ifelse(kfamily_lor_ic >= 0, kfamily_nb1_idx, kfamily_nb2_idx)]

# Genus
kgenus_target = matrix(kgenus[[`target`]], nrow=1)
kgenus_simplex = copy(kgenus)
kgenus_simplex[, Genus := NULL]
kgenus_simplex[, eval(target) := NULL]
kgenus_simplex = as.matrix(kgenus_simplex)
kgenus_col_counts = colSums(kgenus_simplex)
kgenus_simplex = apply(kgenus_simplex, 2, function(x) (x+eps) / sum(x+eps))
kgenus_nb_max = rep(0L, ncol(kgenus_simplex))
for(i in 1:ncol(kgenus_simplex)) {
  kgenus_nb_max[i] = kgenus_col_counts[i] * (matrix(kgenus_simplex[, i], nrow=1) %*% log(matrix(kgenus_simplex[, i], ncol=1)))
}
kgenus_nb = kgenus_target %*% log(kgenus_simplex)

kgenus_nb_sort = sort(kgenus_nb, decreasing=TRUE)
kgenus_nb1_idx = which(kgenus_nb == kgenus_nb_sort[1])
kgenus_nb1 = kgenus_nb_sort[1]
kgenus_nb2_idx = which(kgenus_nb == kgenus_nb_sort[2])
kgenus_nb2 = kgenus_nb_sort[2]

kgenus_lor_ic = (kgenus_nb1 - kgenus_nb_max[kgenus_nb1_idx]) - (kgenus_nb2 - kgenus_nb_max[kgenus_nb2_idx])
kgenus_top = colnames(kgenus_simplex)[kgenus_nb1_idx]
kgenus_top_ic = colnames(kgenus_simplex)[ifelse(kgenus_lor_ic >= 0, kgenus_nb1_idx, kgenus_nb2_idx)]

# Species
kspecies_target = matrix(kspecies[[`target`]], nrow=1)
kspecies_simplex = copy(kspecies)
kspecies_simplex[, Species := NULL]
kspecies_simplex[, eval(target) := NULL]
kspecies_simplex = as.matrix(kspecies_simplex)
kspecies_col_counts = colSums(kspecies_simplex)
kspecies_simplex = apply(kspecies_simplex, 2, function(x) (x+eps) / sum(x+eps))
kspecies_nb_max = rep(0L, ncol(kspecies_simplex))
for(i in 1:ncol(kspecies_simplex)) {
  kspecies_nb_max[i] = kspecies_col_counts[i] * (matrix(kspecies_simplex[, i], nrow=1) %*% log(matrix(kspecies_simplex[, i], ncol=1)))
}
kspecies_nb = kspecies_target %*% log(kspecies_simplex)

kspecies_nb_sort = sort(kspecies_nb, decreasing=TRUE)
kspecies_nb1_idx = which(kspecies_nb == kspecies_nb_sort[1])
kspecies_nb1 = kspecies_nb_sort[1]
kspecies_nb2_idx = which(kspecies_nb == kspecies_nb_sort[2])
kspecies_nb2 = kspecies_nb_sort[2]

kspecies_lor_ic = (kspecies_nb1 - kspecies_nb_max[kspecies_nb1_idx]) - (kspecies_nb2 - kspecies_nb_max[kspecies_nb2_idx])
kspecies_top = colnames(kspecies_simplex)[kspecies_nb1_idx]
kspecies_top_ic = colnames(kspecies_simplex)[ifelse(kspecies_lor_ic >= 0, kspecies_nb1_idx, kspecies_nb2_idx)]

# Aggregate
nb_results = data.table(
  taxon_level=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
  top_match=c(kdomain_top, kphylum_top, kclass_top, korder_top, kfamily_top, kgenus_top, kspecies_top),
  top_match_ic=c(kdomain_top_ic, kphylum_top_ic, kclass_top_ic, korder_top_ic, kfamily_top_ic, kgenus_top_ic, kspecies_top_ic),
  nb1=c(kdomain_nb1, kphylum_nb1, kclass_nb1, korder_nb1, kfamily_nb1, kgenus_nb1, kspecies_nb1),
  nb2=c(kdomain_nb2, kphylum_nb2, kclass_nb2, korder_nb2, kfamily_nb2, kgenus_nb2, kspecies_nb2),
  lor_ic=c(kdomain_lor_ic, kphylum_lor_ic, kclass_lor_ic, korder_lor_ic, kfamily_lor_ic, kgenus_lor_ic, kspecies_lor_ic)
)
# print(kdomain_nb)
# print(kphylum_nb)
# print(kclass_nb)
# print(korder_nb)
# print(kfamily_nb)
# print(kgenus_nb)
# print(kspecies_nb)


### Multinomial Likelihood Classification
exact_log_limit = function(theta, n) 
{
  if(n > 0) {
    return(sum(log(theta + seq(0, n))))
  }
  else {
    return(0)
  }
}

dm_marginal_log_likelihood = function(qvec, X)
{
  # this_X = X + eps
  this_X = X
  log_marginal_probs = c()
  for(c in 1:ncol(X)) {
    tvec = this_X[, c]
    qsum_d1 = round(sum(qvec))
    sum_theta = sum(tvec)
    lprob = 0
    lprob = lprob + exact_log_limit(1, qsum_d1 - 1)
    lprob = lprob - exact_log_limit(sum_theta, qsum_d1 - 1)
    for(d in 1:length(qvec)) {
      if(qvec[d] > 0) {
        lprob = lprob + exact_log_limit(tvec[d], round(qvec[d] - 1)) - exact_log_limit(1, round(qvec[d] - 1))
      }
    }
    log_marginal_probs[c] = lprob
  }
  return(log_marginal_probs)
}

kdomain_mlprobs = dm_marginal_log_likelihood(as.vector(kdomain_target), kdomain_simplex)
kdomain_mlprobs_sort = sort(kdomain_mlprobs, decreasing = TRUE)
kdomain_ml1_idx = which.max(kdomain_mlprobs)
kdomain_ml1_prob = kdomain_mlprobs[kdomain_ml1_idx]
kdomain_ml_top = colnames(kdomain_simplex)[kdomain_ml1_idx]
kdomain_ml2_prob = kdomain_mlprobs_sort[2]
kdomain_ml2_idx = which(kdomain_mlprobs == kdomain_ml2_prob)
kdomain_llr = kdomain_ml1_prob - kdomain_ml2_prob

kphylum_mlprobs = dm_marginal_log_likelihood(as.vector(kphylum_target), kphylum_simplex)
kphylum_mlprobs_sort = sort(kphylum_mlprobs, decreasing = TRUE)
kphylum_ml1_idx = which.max(kphylum_mlprobs)
kphylum_ml1_prob = kphylum_mlprobs[kphylum_ml1_idx]
kphylum_ml_top = colnames(kphylum_simplex)[kphylum_ml1_idx]
kphylum_ml2_prob = kphylum_mlprobs_sort[2]
kphylum_ml2_idx = which(kphylum_mlprobs == kphylum_ml2_prob)
kphylum_llr = kphylum_ml1_prob - kphylum_ml2_prob

kclass_mlprobs = dm_marginal_log_likelihood(as.vector(kclass_target), kclass_simplex)
kclass_mlprobs_sort = sort(kclass_mlprobs, decreasing = TRUE)
kclass_ml1_idx = which.max(kclass_mlprobs)
kclass_ml1_prob = kclass_mlprobs[kclass_ml1_idx]
kclass_ml_top = colnames(kclass_simplex)[kclass_ml1_idx]
kclass_ml2_prob = kclass_mlprobs_sort[2]
kclass_ml2_idx = which(kclass_mlprobs == kclass_ml2_prob)
kclass_llr = kclass_ml1_prob - kclass_ml2_prob

korder_mlprobs = dm_marginal_log_likelihood(as.vector(korder_target), korder_simplex)
korder_mlprobs_sort = sort(korder_mlprobs, decreasing = TRUE)
korder_ml1_idx = which.max(korder_mlprobs)
korder_ml1_prob = korder_mlprobs[korder_ml1_idx]
korder_ml_top = colnames(korder_simplex)[korder_ml1_idx]
korder_ml2_prob = korder_mlprobs_sort[2]
korder_ml2_idx = which(korder_mlprobs == korder_ml2_prob)
korder_llr = korder_ml1_prob - korder_ml2_prob

kfamily_mlprobs = dm_marginal_log_likelihood(as.vector(kfamily_target), kfamily_simplex)
kfamily_mlprobs_sort = sort(kfamily_mlprobs, decreasing = TRUE)
kfamily_ml1_idx = which.max(kfamily_mlprobs)
kfamily_ml1_prob = kfamily_mlprobs[kfamily_ml1_idx]
kfamily_ml_top = colnames(kfamily_simplex)[kfamily_ml1_idx]
kfamily_ml2_prob = kfamily_mlprobs_sort[2]
kfamily_ml2_idx = which(kfamily_mlprobs == kfamily_ml2_prob)
kfamily_llr = kfamily_ml1_prob - kfamily_ml2_prob

kgenus_mlprobs = dm_marginal_log_likelihood(as.vector(kgenus_target), kgenus_simplex)
kgenus_mlprobs_sort = sort(kgenus_mlprobs, decreasing = TRUE)
kgenus_ml1_idx = which.max(kgenus_mlprobs)
kgenus_ml1_prob = kgenus_mlprobs[kgenus_ml1_idx]
kgenus_ml_top = colnames(kgenus_simplex)[kgenus_ml1_idx]
kgenus_ml2_prob = kgenus_mlprobs_sort[2]
kgenus_ml2_idx = which(kgenus_mlprobs == kgenus_ml2_prob)
kgenus_llr = kgenus_ml1_prob - kgenus_ml2_prob

kspecies_mlprobs = dm_marginal_log_likelihood(as.vector(kspecies_target), kspecies_simplex)
kspecies_mlprobs_sort = sort(kspecies_mlprobs, decreasing = TRUE)
kspecies_ml1_idx = which.max(kspecies_mlprobs)
kspecies_ml1_prob = kspecies_mlprobs[kspecies_ml1_idx]
kspecies_ml_top = colnames(kspecies_simplex)[kspecies_ml1_idx]
kspecies_ml2_prob = kspecies_mlprobs_sort[2]
kspecies_ml2_idx = which(kspecies_mlprobs == kspecies_ml2_prob)
kspecies_llr = kspecies_ml1_prob - kspecies_ml2_prob

# Aggregate
custom_aic2 = function(ll, K, n)
{
  return((-2 * ll) + (2 * K) + (2 * K * (K + 1) / (n - K - 1)))
}


llr_vec = c(kdomain_llr, kphylum_llr, kclass_llr, korder_llr, kfamily_llr, kgenus_llr, kspecies_llr)
top_prob_vec = c(kdomain_ml1_prob, kphylum_ml1_prob, kclass_ml1_prob, korder_ml1_prob, 
                 kfamily_ml1_prob, kgenus_ml1_prob, kspecies_ml1_prob)
param_counts = c(nrow(kdomain), nrow(kphylum), nrow(kclass), nrow(korder), nrow(kfamily), nrow(kgenus), nrow(kspecies))
mlc_results = data.table(
  taxon_level=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
  top_match=c(kdomain_ml_top, kphylum_ml_top, kclass_ml_top, korder_ml_top, kfamily_ml_top, kgenus_ml_top, kspecies_ml_top),
  top_prob=top_prob_vec,
  likelihood_ratio=llr_vec,
  aic=custom_aic2(top_prob_vec, param_counts, 1)
)

# cat('\n\n')
# print(kdomain_mlprobs)
# print(kphylum_mlprobs)
# print(kclass_mlprobs)
# print(korder_mlprobs)
# print(kfamily_mlprobs)
# print(kgenus_mlprobs)
# print(kspecies_mlprobs)


## Ordination on simplex

# Domain
# kdomain_ord = kdomain[, .SD, .SDcols=!'Domain']
# kdomain_ord = t(apply(kdomain_ord, 2, function(x) x / sum(x)))
# 
# domain.pca.res = prcomp(kdomain_ord)
# domain.pca.proj = data.table(scale(kdomain_ord, domain.pca.res$center, domain.pca.res$scale) %*% domain.pca.res$rotation[, 1:2])
# domain.pca.proj[, samples := colnames(kdomain)[2:ncol(kdomain)]]
# domain.pca.proj[, taxonomy_level := rep('Domain', nrow(domain.pca.proj))]
# 
# domain.nmds.res = metaMDS(kdomain_ord, k=2, try=50, trymax=100)
# domain.nmds.proj = data.table(domain.nmds.res$points)
# domain.nmds.proj[, samples := colnames(kdomain)[2:ncol(kdomain)]]
# domain.nmds.proj[, taxonomy_level := rep('Domain', nrow(domain.nmds.proj))]
# 
# # domain.tsne.res = tsne(kdomain_ord, k=2, perplexity=30, max_iter=1000, epoch=500)
# # domain.tsne.proj = data.table(domain.tsne.res)
# # domain.tsne.proj[, samples := colnames(kdomain)[2:ncol(kdomain)]]
# # domain.tsne.proj[, taxonomy_level := rep('Domain', nrow(domain.tsne.proj))]
# 
# 
# # Phylum
# kphylum_ord = kphylum[, .SD, .SDcols=!'Phylum']
# kphylum_ord = t(apply(kphylum_ord, 2, function(x) x / sum(x)))
# 
# phylum.pca.res = prcomp(kphylum_ord)
# phylum.pca.proj = data.table(scale(kphylum_ord, phylum.pca.res$center, phylum.pca.res$scale) %*% phylum.pca.res$rotation[, 1:2])
# phylum.pca.proj[, samples := colnames(kphylum)[2:ncol(kphylum)]]
# phylum.pca.proj[, taxonomy_level := rep('Phylum', nrow(phylum.pca.proj))]
# 
# phylum.nmds.res = metaMDS(kphylum_ord, k=2, try=50, trymax=100)
# phylum.nmds.proj = data.table(phylum.nmds.res$points)
# phylum.nmds.proj[, samples := colnames(kphylum)[2:ncol(kphylum)]]
# phylum.nmds.proj[, taxonomy_level := rep('Phylum', nrow(phylum.nmds.proj))]
# 
# # phylum.tsne.res = tsne(kphylum_ord, k=2, perplexity=30, max_iter=1000, epoch=500)
# # phylum.tsne.proj = data.table(phylum.tsne.res)
# # phylum.tsne.proj[, samples := colnames(kphylum)[2:ncol(kphylum)]]
# # phylum.tsne.proj[, taxonomy_level := rep('Phylum', nrow(phylum.tsne.proj))]
# 
# 
# # Class
# kclass_ord = kclass[, .SD, .SDcols=!'Class']
# kclass_ord = t(apply(kclass_ord, 2, function(x) x / sum(x)))
# 
# class.pca.res = prcomp(kclass_ord)
# class.pca.proj = data.table(scale(kclass_ord, class.pca.res$center, class.pca.res$scale) %*% class.pca.res$rotation[, 1:2])
# class.pca.proj[, samples := colnames(kclass)[2:ncol(kclass)]]
# class.pca.proj[, taxonomy_level := rep('Class', nrow(class.pca.proj))]
# 
# class.nmds.res = metaMDS(kclass_ord, k=2, try=50, trymax=100)
# class.nmds.proj = data.table(class.nmds.res$points)
# class.nmds.proj[, samples := colnames(kclass)[2:ncol(kclass)]]
# class.nmds.proj[, taxonomy_level := rep('Class', nrow(class.nmds.proj))]
# 
# # class.tsne.res = tsne(kclass_ord, k=2, perplexity=30, max_iter=1000, epoch=500)
# # class.tsne.proj = data.table(class.tsne.res)
# # class.tsne.proj[, samples := colnames(kclass)[2:ncol(kclass)]]
# # class.tsne.proj[, taxonomy_level := rep('Class', nrow(class.tsne.proj))]
# 
# 
# # Order
# korder_ord = korder[, .SD, .SDcols=!'Order']
# korder_ord = t(apply(korder_ord, 2, function(x) x / sum(x)))
# 
# order.pca.res = prcomp(korder_ord)
# order.pca.proj = data.table(scale(korder_ord, order.pca.res$center, order.pca.res$scale) %*% order.pca.res$rotation[, 1:2])
# order.pca.proj[, samples := colnames(korder)[2:ncol(korder)]]
# order.pca.proj[, taxonomy_level := rep('Order', nrow(order.pca.proj))]
# 
# order.nmds.res = metaMDS(korder_ord, k=2, try=50, trymax=100)
# order.nmds.proj = data.table(order.nmds.res$points)
# order.nmds.proj[, samples := colnames(korder)[2:ncol(korder)]]
# order.nmds.proj[, taxonomy_level := rep('Order', nrow(order.nmds.proj))]
# 
# # order.tsne.res = tsne(korder_ord, k=2, perplexity=30, max_iter=1000, epoch=500)
# # order.tsne.proj = data.table(order.tsne.res)
# # order.tsne.proj[, samples := colnames(korder)[2:ncol(korder)]]
# # order.tsne.proj[, taxonomy_level := rep('Order', nrow(order.tsne.proj))]
# # 
# 
# # Family
# kfamily_ord = kfamily[, .SD, .SDcols=!'Family']
# kfamily_ord = t(apply(kfamily_ord, 2, function(x) x / sum(x)))
# 
# family.pca.res = prcomp(kfamily_ord)
# family.pca.proj = data.table(scale(kfamily_ord, family.pca.res$center, family.pca.res$scale) %*% family.pca.res$rotation[, 1:2])
# family.pca.proj[, samples := colnames(kfamily)[2:ncol(kfamily)]]
# family.pca.proj[, taxonomy_level := rep('Family', nrow(family.pca.proj))]
# 
# family.nmds.res = metaMDS(kfamily_ord, k=2, try=50, trymax=100)
# family.nmds.proj = data.table(family.nmds.res$points)
# family.nmds.proj[, samples := colnames(kfamily)[2:ncol(kfamily)]]
# family.nmds.proj[, taxonomy_level := rep('Family', nrow(family.nmds.proj))]
# # 
# # family.tsne.res = tsne(kfamily_ord, k=2, perplexity=30, max_iter=1000, epoch=500)
# # family.tsne.proj = data.table(family.tsne.res)
# # family.tsne.proj[, samples := colnames(kfamily)[2:ncol(kfamily)]]
# # family.tsne.proj[, taxonomy_level := rep('Family', nrow(family.tsne.proj))]
# 
# 
# # Genus
# kgenus_ord = kgenus[, .SD, .SDcols=!'Genus']
# kgenus_ord = t(apply(kgenus_ord, 2, function(x) x / sum(x)))
# 
# genus.pca.res = prcomp(kgenus_ord)
# genus.pca.proj = data.table(scale(kgenus_ord, genus.pca.res$center, genus.pca.res$scale) %*% genus.pca.res$rotation[, 1:2])
# genus.pca.proj[, samples := colnames(kgenus)[2:ncol(kgenus)]]
# genus.pca.proj[, taxonomy_level := rep('Genus', nrow(genus.pca.proj))]
# 
# genus.nmds.res = metaMDS(kgenus_ord, k=2, try=50, trymax=100)
# genus.nmds.proj = data.table(genus.nmds.res$points)
# genus.nmds.proj[, samples := colnames(kgenus)[2:ncol(kgenus)]]
# genus.nmds.proj[, taxonomy_level := rep('Genus', nrow(genus.nmds.proj))]
# 
# # genus.tsne.res = tsne(kgenus_ord, k=2, perplexity=30, max_iter=1000, epoch=500)
# # genus.tsne.proj = data.table(genus.tsne.res)
# # genus.tsne.proj[, samples := colnames(kgenus)[2:ncol(kgenus)]]
# # genus.tsne.proj[, taxonomy_level := rep('Genus', nrow(genus.tsne.proj))]
# 
# 
# # Species
# kspecies_ord = kspecies[, .SD, .SDcols=!'Species']
# kspecies_ord = t(apply(kspecies_ord, 2, function(x) x / sum(x)))
# 
# species.pca.res = prcomp(kspecies_ord)
# species.pca.proj = data.table(scale(kspecies_ord, species.pca.res$center, species.pca.res$scale) %*% species.pca.res$rotation[, 1:2])
# species.pca.proj[, samples := colnames(kspecies)[2:ncol(kspecies)]]
# species.pca.proj[, taxonomy_level := rep('Species', nrow(species.pca.proj))]
# 
# species.nmds.res = metaMDS(kspecies_ord, k=2, try=50, trymax=100)
# species.nmds.proj = data.table(species.nmds.res$points)
# species.nmds.proj[, samples := colnames(kspecies)[2:ncol(kspecies)]]
# species.nmds.proj[, taxonomy_level := rep('Species', nrow(species.nmds.proj))]
# # 
# # species.tsne.res = tsne(kspecies_ord, k=2, perplexity=30, max_iter=1000, epoch=500)
# # species.tsne.proj = data.table(species.tsne.res)
# # species.tsne.proj[, samples := colnames(kspecies)[2:ncol(kspecies)]]
# # species.tsne.proj[, taxonomy_level := rep('Species', nrow(species.tsne.proj))]
# 
# 
# # Aggregate
# pca.proj = rbind(domain.pca.proj, phylum.pca.proj, class.pca.proj, order.pca.proj, family.pca.proj,
#                  genus.pca.proj, species.pca.proj)
# nmds.proj = rbind(domain.nmds.proj, phylum.nmds.proj, class.nmds.proj, order.nmds.proj, family.nmds.proj,
#                   genus.nmds.proj, species.nmds.proj)
# tsne.proj = rbind(domain.tsne.proj, phylum.tsne.proj, class.tsne.proj, order.tsne.proj, family.tsne.proj,
#                   genus.tsne.proj, species.tsne.proj)

# Plot
# basesize = 32
# ptsize = 4
# 
# g_pca = ggplot(pca.proj, aes(x=PC1, y=PC2, color=samples)) +
#   geom_point(size=ptsize) +
#   facet_wrap(~taxonomy_level, nrow=2) +
#   custom_theme_legend(base_size = basesize)
# 
# g_nmds = ggplot(nmds.proj, aes(x=MDS1, y=MDS2, color=samples)) +
#   geom_point(size=ptsize) +
#   facet_wrap(~taxonomy_level, nrow=2) +
#   custom_theme_legend(base_size = basesize)
# 
# png(file.path(output_dir, 'graphs', 'microbiome_pca_multi_ordination.png'), width=3200, height=2000)
# print(g_pca)
# dev.off()
# 
# png(file.path(output_dir, 'graphs', 'microbiome_nmds_multi_ordination.png'), width=3200, height=2000)
# print(g_nmds)
# dev.off()


### Write results from microbiome analysis
write.csv(kdomain, file=file.path(output_dir, 'microbiome_matrices', 'domain.csv'), row.names=FALSE)
write.csv(kphylum, file=file.path(output_dir, 'microbiome_matrices', 'phylum.csv'), row.names=FALSE)
write.csv(kclass, file=file.path(output_dir, 'microbiome_matrices', 'class.csv'), row.names=FALSE)
write.csv(korder, file=file.path(output_dir, 'microbiome_matrices', 'order.csv'), row.names=FALSE)
write.csv(kfamily, file=file.path(output_dir, 'microbiome_matrices', 'family.csv'), row.names=FALSE)
write.csv(kgenus, file=file.path(output_dir, 'microbiome_matrices', 'genus.csv'), row.names=FALSE)
write.csv(kspecies, file=file.path(output_dir, 'microbiome_matrices', 'species.csv'), row.names=FALSE)

write.csv(cbind(as.vector(kdomain_target), kdomain_simplex), 
          file=file.path(output_dir, 'microbiome_simplex_matrices', 'domain.csv'), row.names=FALSE)
write.csv(cbind(as.vector(kphylum_target), kphylum_simplex), 
          file=file.path(output_dir, 'microbiome_simplex_matrices', 'phylum.csv'), row.names=FALSE)
write.csv(cbind(as.vector(kclass_target), kclass_simplex), 
          file=file.path(output_dir, 'microbiome_simplex_matrices', 'class.csv'), row.names=FALSE)
write.csv(cbind(as.vector(korder_target), korder_simplex), 
          file=file.path(output_dir, 'microbiome_simplex_matrices', 'order.csv'), row.names=FALSE)
write.csv(cbind(as.vector(kfamily_target), kfamily_simplex), 
          file=file.path(output_dir, 'microbiome_simplex_matrices', 'family.csv'), row.names=FALSE)
write.csv(cbind(as.vector(kgenus_target), kgenus_simplex), 
          file=file.path(output_dir, 'microbiome_simplex_matrices', 'genus.csv'), row.names=FALSE)
write.csv(cbind(as.vector(kspecies_target), kspecies_simplex), 
          file=file.path(output_dir, 'microbiome_simplex_matrices', 'species.csv'), row.names=FALSE)

# write.csv(nb_results, file=file.path(output_dir, 'reports', 'naive_bayes_results.csv'), row.names=FALSE)
write.csv(mlc_results, file=file.path(output_dir, 'reports', 'multinomial_likelihood_results.csv'), row.names=FALSE)
# write.csv(pca.proj, file=file.path(output_dir, 'graphs', 'microbiome_pca_multiordination.csv'), row.names=FALSE)
# write.csv(nmds.proj, file=file.path(output_dir, 'graphs', 'microbiome_nmds_multiordination.csv'), row.names=FALSE)




### Fingerprinting analysis
virus = data.table(read.csv(viral_var_fpath, stringsAsFactors = FALSE, header=TRUE))
host = data.table(read.csv(host_var_fpath, stringsAsFactors = FALSE, header=TRUE))

fingerprint_res = data.table(
  sample=character(),
  data_type=character(),
  analysis=character(),
  top_match=character(),
  value=numeric()
)

virus_compare = copy(virus)
virus_compare[, VariantID := NULL]
virus_compare[, eval(target) := NULL]
virus_target = virus[[`target`]]
virus_match_perc = c()
virus_mismatch_perc = c()
for(i in 1:ncol(virus_compare)) {
  virus_match_perc[i] = 100 * sum(virus_target == virus_compare[, ..i]) / nrow(virus)
  virus_mismatch_perc[i] = 100 - virus_match_perc[i]
}

for(i in 1:length(virus_match_perc)) {
  fingerprint_res = rbind(fingerprint_res,
                          data.table(
                            sample=colnames(virus_compare)[i],
                            data_type='Virus',
                            analysis='% Concordance',
                            top_match=colnames(virus_compare)[which.max(virus_match_perc)],
                            value=max(virus_match_perc)
                          ))
}



host_compare = copy(host)
host_compare[, VariantID := NULL]
host_compare[, eval(target) := NULL]
host_target = host[[`target`]]
host_match_perc = c()
host_mismatch_perc = c()
for(i in 1:ncol(host_compare)) {
  host_match_perc[i] = 100 * sum(host_target == host_compare[, ..i]) / nrow(host)
  host_mismatch_perc[i] = 100 - host_match_perc[i]
}

for(i in 1:length(host_match_perc)) {
  fingerprint_res = rbind(fingerprint_res,
                          data.table(
                            sample=colnames(host_compare)[i],
                            data_type='Host',
                            analysis='% Concordance',
                            top_match=colnames(host_compare)[which.max(host_match_perc)],
                            value=max(host_match_perc)
                          ))
}

sig_thresh = 0

virus_omics_x = as.matrix(as.numeric(virus_target))
colnames(virus_omics_x) = target
rownames(virus_omics_x) = virus[['VariantID']]

virus_omics_y = as.matrix(sapply(virus_compare, as.numeric))
colnames(virus_omics_y) = colnames(virus_compare)
rownames(virus_omics_y) = virus[['VariantID']]

virus_omics_all = cbind(virus_omics_x, virus_omics_y)

print(apply(virus_omics_all, 2, class))

virus_omics_res = alleleSharing(virus_omics_all, alpha=sig_thresh)
print(virus_omics_res)

cat('\n\n')

print(host)
print(host[['VariantID']][duplicated(host[['VariantID']])])
print(apply(host, 2, class))

host_omics_x = as.matrix(as.numeric(host_target))
colnames(host_omics_x) = target
rownames(host_omics_x) = host[['VariantID']]

host_omics_y = as.matrix(sapply(host_compare, as.numeric))
colnames(host_omics_y) = colnames(host_compare)
rownames(host_omics_y) = host[['VariantID']]
host_omics_y = data.frame(apply(host_omics_y, 2, function(x) as.numeric(as.character(x))))

host_omics_all = cbind(host_omics_x, host_omics_y)

print(apply(host_omics_all, 2, class))

host_omics_res = alleleSharing(host_omics_all, alpha=sig_thresh)
print(host_omics_res)



