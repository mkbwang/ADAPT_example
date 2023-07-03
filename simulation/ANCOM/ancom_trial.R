source('simulation/ANCOM/ancombc_prep.R')
source('simulation/ANCOM/ancom_utils.R')


data(atlas1006, package = "microbiome")

out <- ancom(data = atlas1006, assay_name = "counts",
                        tax_level = NULL, phyloseq = NULL,
                       p_adj_method = "BY", prv_cut = 0, lib_cut = 0,
                main_var = "bmi_group", adj_formula = NULL,
                      rand_formula = NULL, lme_control = NULL,
                      struc_zero = TRUE, neg_lb = FALSE, alpha = 0.05, n_cl = 1)


