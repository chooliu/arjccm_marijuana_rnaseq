# ==============================================================================
# 09_examine_batcheffect_and_RUV.R
# check for expression diffs between Batch 1 and Batch 2
# determine number of RUV components (RUVr, based on first-pass resids)
# time-intensive to run, can be skipped
# ==============================================================================



# RLE plots
# ------------------------------------------------------------------------------
edgeR_obj$counts %>%
  `+`(., 0.01) %>%
  log() %>%
  apply(., 1, (function(x) { x - median(x) })) %>%
  cbind(., sample = metadata_filtered$ID, Batch = metadata_filtered$Batch) %>%
  as_tibble() %>%
  gather(key, value, -c("sample", "Batch")) %>%
  mutate(value = as.numeric(value)) %>%
  arrange(Batch) %>%
  mutate(sample = as.factor(sample) %>% fct_inorder) %>%
  ggplot(data = ., aes(x = sample, y = value, fill = Batch)) +
  geom_hline(color = "#555555", yintercept = 0) +
  geom_boxplot(outlier.shape = NA, notch = T, alpha = 0.5) +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(limits = c(-2, 2)) +
  theme_few() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ef8a62", "#67a9cf")) +
  xlab("Sample") + ylab("Relative Log Expression")
ggsave("./FiguresTables/FigS1A.png", width = 4, height = 2.5, dpi = 300)




# obtain first-pass residuals (pre-RUV)
# ------------------------------------------------------------------------------

# no batch correction
dm_firstpass <-
  model.matrix(~ Group + Age + Sex + Obese,
               data = metadata_filtered[filter_subjects_to_keep, ])
edgeR_obj_firstpass <- estimateDisp(edgeR_obj, dm_firstpass)
edgeR_obj_firstpass <- calcNormFactors(edgeR_obj_firstpass)
edgeR_fit_firstpass <- glmQLFit(edgeR_obj_firstpass, dm_firstpass)
resid_firstpass <- residuals(edgeR_fit_firstpass, type = "deviance")

# explicit correction for batch
dm_firstpass_batch <-
  model.matrix(~ Group + Age + Sex + Obese + Batch,
               data = metadata_filtered[filter_subjects_to_keep, ])
edgeR_obj_firstpass_batch <- estimateDisp(edgeR_obj, dm_firstpass_batch)
edgeR_obj_firstpass_batch <- calcNormFactors(edgeR_obj_firstpass_batch)
edgeR_fit_batch <- glmQLFit(edgeR_obj_firstpass_batch, dm_firstpass_batch)
resid_firstpass_batch <- residuals(edgeR_fit_batch, type = "deviance")




# check RLE, by k (only needs to be run once)
# k = 3 seems to yield a plateau in btwn-subject variability
# ------------------------------------------------------------------------------
RLE_by_RUVk <-
  function (k, resid_in, designmat) {

    ruv_results <-
      RUVr(
        edgeR_obj$counts,
        T,
        k = k,
        resid_in)

    dm_add_RUV <- cbind(designmat, ruv_results$W)
    edgeR_obj_ruv <- estimateDisp(edgeR_obj, dm_add_RUV)
    edgeR_obj_ruv <- calcNormFactors(edgeR_obj_ruv)
    edgeR_fit_ruv <- glmQLFit(edgeR_obj_ruv, dm_add_RUV)

    resid_with_ruv <-
      residuals(edgeR_fit_ruv, type = "deviance")

    resid_with_ruv %>%
      apply(., 1, function(x) { x - median(x) })



  }







# plot RLE first pass
# ------------------------------------------------------------------------------
resid_firstpass %>%
  apply(., 1, (function(x) { x - median(x) })) %>%
  cbind(., sample = metadata_filtered$ID, Batch = metadata_filtered$Batch) %>%
  as_tibble() %>%
  gather(key, value, -c("sample", "Batch")) %>%
  mutate(value = as.numeric(value)) %>%
  arrange(Batch) %>%
  mutate(sample = as.factor(sample) %>% fct_inorder) %>%
  ggplot(data = ., aes(x = sample, y = value, fill = Batch)) +
  geom_hline(color = "#555555", yintercept = 0) +
  geom_boxplot(outlier.shape = NA, notch = T, alpha = 0.5) +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(limits = c(-4, 4)) +
  theme_few() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ef8a62", "#67a9cf")) +
  xlab("Sample") + ylab("Residuals, No Batch Adjustment")



# firstpass + batch correct
# ------------------------------------------------------------------------------
resid_firstpass_batch %>%
  apply(., 1, (function(x) { x - median(x) })) %>%
  cbind(., sample = metadata_filtered$ID, Batch = metadata_filtered$Batch) %>%
  as_tibble() %>%
  gather(key, value, -c("sample", "Batch")) %>%
  mutate(value = as.numeric(value)) %>%
  arrange(Batch) %>%
  mutate(sample = as.factor(sample) %>% fct_inorder) %>%
  ggplot(data = ., aes(x = sample, y = value, fill = Batch)) +
  geom_hline(color = "#555555", yintercept = 0) +
  geom_boxplot(outlier.shape = NA, notch = T, alpha = 0.5) +
  scale_x_discrete(breaks = NULL) +
  scale_y_continuous(limits = c(-4, 4)) +
  theme_few() +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("#ef8a62", "#67a9cf")) +
  xlab("Sample") + ylab("Residuals, Adjusting for Batch")









# plot RLEs at varying k
# ------------------------------------------------------------------------------



make_RUV_boxplot <- function(k, label = NULL,
                             resid_in = resid_firstpass, distmat = dm_firstpass) {

  RLE_by_RUVk(k, resid_in, distmat) %>%
    as_tibble() %>%
    cbind(sample = metadata_filtered$ID, Batch = metadata_filtered$Batch) %>%
    gather(key, value, -c("sample", "Batch")) %>%
    mutate(value = as.numeric(value)) %>%
    arrange(Batch) %>%
    mutate(sample = as.factor(sample) %>% fct_inorder,
           Batch = as.factor(Batch)) %>%
    ggplot(data = ., aes(x = sample, y = value, fill = Batch)) +
    geom_hline(color = "#555555", yintercept = 0) +
    geom_boxplot(outlier.shape = NA, notch = T, alpha = 0.5) +
    scale_x_discrete(breaks = NULL) +
    scale_y_continuous(limits = c(-3.5, 3.5)) +
    theme_few() +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("#ef8a62", "#67a9cf")) +
    xlab("Sample") + ylab(if_else(is.null(label), paste0("Residuals, RUV (k = ", k, ")"), label))

}

make_RUV_boxplot(1)
make_RUV_boxplot(2)
make_RUV_boxplot(3)

make_RUV_boxplot(1, label = "Deviance Residuals,\nRUV (k = 1) + Batch",
                 resid_in = resid_firstpass_batch, distmat = dm_firstpass_batch)
make_RUV_boxplot(2, label = "Deviance Residuals,\nRUV (k = 2) + Batch",
                 resid_in = resid_firstpass_batch, distmat = dm_firstpass_batch)
make_RUV_boxplot(3, label = "Deviance Residuals,\nRUV (k = 3) + Batch",
                 resid_in = resid_firstpass_batch, distmat = dm_firstpass_batch)
ggsave("./FiguresTables/FigS2A.png", width = 4, height = 2.5, dpi = 300)






# check k
# ------------------------------------------------------------------------------

evaluate_btwnparticipant_by_k <-
  function (k, variance_function = iqr,
            resid_in = resid_firstpass, dm_in = dm_firstpass) {

    ruv_results <-
      RUVr(
        edgeR_obj$counts,
        T,
        k = k,
        resid_in)

    ruv_results$normalizedCounts %>% # 17602 x 41
      apply(., 1, variance_function) %>% # within-subject variability
      (function(x) { max(x) - min(x) }) # between-subjects

  }

plot(1:15,
     sapply(1:15, evaluate_btwnparticipant_by_k,
            variance_function = iqr,
            resid_in = resid_firstpass_batch,
            dm_in = dm_firstpass_batch))

cv_filter <-
  apply(edgeR_obj$counts, 1,
        function(x) { sd(x) / mean(x) }) %>%
  order(decreasing = T) %>%
  .[1:2500]

dist_by_sample <-
  log(edgeR_obj$counts[cv_filter, ] + 1) %>% t %>% dist

evaluate_ruv_k_adonis <- function(k_to_try) {
  ruv_results <-
    RUVr(
      edgeR_obj$counts,
      T,
      k = k_to_try,
      resid_firstpass_batch,
    )
  
  adonis(dist_by_sample ~ ruv_results$W,
                data = NULL, permutations = 1000) %>%
    .$aov.tab %>% .$R2 %>% .[1]
  
}


plot(1:15, sapply(1:15, evaluate_ruv_k_adonis))













