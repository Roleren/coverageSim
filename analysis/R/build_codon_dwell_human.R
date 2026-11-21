# Build inst/extdata/codon_dwell_human.csv (S12, explainability loop; run once on the shared server).
# Sources (read-only):
#  - H1 codon-dwell GLM (Codex explainer_h1_codon_dwell_20261118/model.rds): 16 offsets -10..+5 relative to the A site
#    x 60 codon contrasts, reference = model$sense[1] (AAA, beta 0); beta index (k + 10) * 60 + a - 1 for codon sense[a].
#  - Human dwell atlas per-library betas (claude-megafst m6_study_fingerprint/beta_matrix_w10_5.rds: 2,960 libraries,
#    367 studies, 16 offsets x 61 codons): consensus = median over studies of the study means; study_sd = per offset
#    root mean square of study mean - consensus; library_sd = per offset root mean square of library - study mean
#    (includes sampling noise).
#  - Striated (muscle + heart) delta vs consensus at E / P / A (claude-megafst dwell_by_group/muscle_detail.csv, M1-a).
library(data.table)
h1 <- readRDS("/home/rstudio/shared_workspace/hakon_playground/agent_sync/codex_tasks/explainer_h1_codon_dwell_20261118/model.rds")
atlas <- readRDS("/home/rstudio/shared_workspace/open_workspace/claude-megafst/m6_study_fingerprint/beta_matrix_w10_5.rds")
muscle <- fread("/home/rstudio/shared_workspace/open_workspace/claude-megafst/dwell_by_group/muscle_detail.csv")
muscle[, `:=`(offset = as.integer(sub(":.*", "", feature)), codon = sub(".*:", "", feature))]

sense <- h1$sense
out <- CJ(offset = -10:5, codon = sense, sorted = FALSE)
a <- match(out$codon, sense)
beta_index <- (out$offset + 10L) * 60L + a - 1L  # 0 for the reference codon of offset -10
out[, h1_beta := 0]
out[a > 1L, h1_beta := h1$beta[beta_index[a > 1L]]]

X <- atlas$X; st <- atlas$m$BioProject
smean <- rowsum(X, st) / as.vector(table(st)[rownames(rowsum(X, st))])
consensus <- apply(smean, 2, stats::median)
sdev <- sweep(smean, 2, consensus); ldev <- X - smean[st, ]
off <- as.integer(sub(":.*", "", colnames(X))); cod <- sub(".*:", "", colnames(X))
key <- paste(out$offset, out$codon, sep = ":")
out[, atlas_consensus := unname(consensus[key])]
per <- data.table(offset = sort(unique(off)))
per[, study_sd := vapply(offset, function(o) sqrt(mean(sdev[, off == o]^2)), numeric(1))]
per[, library_sd := vapply(offset, function(o) sqrt(mean(ldev[, off == o]^2)), numeric(1))]
out <- merge(out, per, by = "offset", sort = FALSE)
out[, striated_delta := 0]
out[muscle, striated_delta := i.delta_muscle, on = c(codon = "codon", offset = "offset")][]
stopifnot(nrow(out) == 976L, !anyNA(out))
setorder(out, offset, codon)
fwrite(out[, .(offset, codon, h1_beta = signif(h1_beta, 6), atlas_consensus = signif(atlas_consensus, 6),
               study_sd = signif(study_sd, 6), library_sd = signif(library_sd, 6),
               striated_delta = signif(striated_delta, 6))],
       "inst/extdata/codon_dwell_human.csv")
