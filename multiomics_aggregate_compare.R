#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(omicsPrint)

options(echo=FALSE)

# 1. kraken_analytic_matrix.csv, 2. host_vcf_table.tsv, 3. viral_vcf_table.tsv, 4. target_sample_name, 5. out_dir
args = commandArgs(trailingOnly = TRUE)
kraken_fpath = args[1]
host_vcf_fpath = args[2]
viral_vcf_fpath = args[3]
target = args[4]
output_dir = args[5]

create_working_directories = function(root_dir) 
{
  # If subdirs for stats and exploratory variables don't exist, create them
  ifelse(!dir.exists(file.path(root_dir)),
         dir.create(file.path(root_dir), mode='777'), FALSE)
  ifelse(!dir.exists(file.path(root_dir, 'microbiome_matrices')),
         dir.create(file.path(root_dir, 'microbiome_matrices'), mode='777'), FALSE)
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


### Write results from microbiome analysis
write.csv(kdomain, file=file.path(output_dir, 'microbiome_matrices', 'domain.csv'), row.names=FALSE)
write.csv(kphylum, file=file.path(output_dir, 'microbiome_matrices', 'phylum.csv'), row.names=FALSE)
write.csv(kclass, file=file.path(output_dir, 'microbiome_matrices', 'class.csv'), row.names=FALSE)
write.csv(korder, file=file.path(output_dir, 'microbiome_matrices', 'order.csv'), row.names=FALSE)
write.csv(kfamily, file=file.path(output_dir, 'microbiome_matrices', 'family.csv'), row.names=FALSE)
write.csv(kgenus, file=file.path(output_dir, 'microbiome_matrices', 'genus.csv'), row.names=FALSE)
write.csv(kspecies, file=file.path(output_dir, 'microbiome_matrices', 'species.csv'), row.names=FALSE)

# write.csv(nb_results, file=file.path(output_dir, 'reports', 'naive_bayes_results.csv'), row.names=FALSE)
write.csv(mlc_results, file=file.path(output_dir, 'reports', 'multinomial_likelihood_results.csv'), row.names=FALSE)

print(cbind(as.vector(korder_target), korder_simplex))
