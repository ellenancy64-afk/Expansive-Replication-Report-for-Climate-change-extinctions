library(tidyverse)
library(rstan)
library(bayesplot)
options(mc.cores = parallel::detectCores())


####################################################################
dataP <- read.table(".MetaRisk2 aggthres 5.txt", header = T)
sjPlot::view_df(dataP)

data.use <- dataP %>%
  filter(is.finite(Pre.Ind.Rise)) %>%
  mutate(percent2 = case_when(
    adj.percent == 0 ~ 0.001,
    adj.percent == 1 ~ 0.999,
    .default = adj.percent
  )) %>%
  mutate(habs = case_when(
    Region == "Marine" ~ "Marine",
    Fresh  == "Y"      ~ "Fresh",
    .default           = "Terrestrial"
  ))


data.use %>%
  ggplot(aes(x = Pre.Ind.Rise, y = percent2, size = log(Total.N))) +
  geom_point(alpha = 0.6, shape = 20, color = "#5485A0") +
  scale_x_continuous(breaks = seq(0, 5, 1), limits = c(0, 5.5)) +
  scale_y_continuous(breaks = seq(0, .35, 0.05), limits = c(0, .35)) +
  scale_fill_continuous(low = "white", high = "#5485A0") +
  xlab("Pre-industrial rise in temperature (C)") +
  ylab("Predicted extinction risk") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18),
    title      = element_text(size = 20),
    axis.text  = element_text(size = 16),
    legend.position = "none"
  )
####################################################################






####################################################################
# create model matrix for coefficients
betamat <- model.matrix(~ 0 + Pre.Ind.Rise:Region, data = data.use)

stan.data <- list(
  N       = nrow(data.use),
  percent = data.use$percent2,
  betamat = betamat,
  phi     = data.use$Total.N,
  S       = length(unique(data.use$Study)),
  P       = ncol(betamat),
  Study   = as.integer(unclass(factor(data.use$Study)))
)

mod <- stan(
  file    = ".MetaRisk2 RSTAN quad.stan",
  data    = stan.data,
  chains  = 4,
  iter    = 8000,
  warmup  = 4000,
  control = list(adapt_delta = 0.9, max_treedepth = 15)
)
#######################################################




#######################################################
tree.rich <- c(10441, 7035, 6680, 7035, 8646, 27186)
vert.rich <- c(4646, 6597, 1973, 550, 4085, 5620)
cont.rich <- tree.rich + vert.rich
rel.rich  <- cont.rich / sum(cont.rich)
#######################################################




#######################################################
pred.reg.df <- mod %>% 
  posterior::as_draws_df() %>% 
  select(mu, 
         beta.1 = "beta[1]",
         beta.2 = "beta[2]",
         beta.3 = "beta[3]",
         beta.4 = "beta[4]",
         # beta.5 = "beta[5]", # oceans not treated  here
         beta.6 = "beta[6]",
         beta.7 = "beta[7]"
  ) %>% 
  mutate(beta.cont.w = 
           beta.1 * rel.rich[1] + 
           beta.2 * rel.rich[2] + 
           beta.3 * rel.rich[3] + 
           beta.4 * rel.rich[4] + 
           beta.6 * rel.rich[5] + 
           beta.7 * rel.rich[6]
  ) %>%
  select(mu, beta.cont.w) %>% 
  expand(
    nesting(mu, beta.cont.w),
    P.Ind = seq(from = 0, to = 5.5, by = .1)
  ) %>%
  mutate(pred.reg.quantw = brms::inv_logit_scaled(mu + beta.cont.w * P.Ind)) %>%
  group_by(P.Ind) %>%
  tidybayes::mean_qi(pred.reg.quantw, .width = 0.95)
####################################################################




####################################################################
pred.reg.df %>%
  ggplot() +
  geom_point(
    data = data.use,
    aes(x = Pre.Ind.Rise, y = percent2, size = log(Total.N)),
    alpha = 0.6, shape = 20, color = "#5485A0"
  ) +
  geom_ribbon(
    aes(x = P.Ind, ymin = .lower, ymax = .upper),
    alpha = .2, fill = "darkred"
  ) +
  geom_line(
    aes(x = P.Ind, y = pred.reg.quantw),
    linewidth = 3, color = "darkred"
  ) +
  xlab("Pre-industrial rise in temperature (C)") +
  ylab("Predicted extinction risk") +
  theme_classic() +
  scale_x_continuous(
    breaks = seq(0, 5, 1),
    limits = c(0, 5.5)
  ) +
  scale_y_continuous(
    breaks = seq(0, .35, 0.05),
    limits = c(0, .4)
  ) +
  theme(
    axis.title      = element_text(size = 18),
    title           = element_text(size = 20),
    axis.text       = element_text(size = 16),
    legend.position = "none"
  )
####################################################################
library(tidyverse)
library(tidybayes)
library(brms)
library(posterior)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
draws_df <- mod %>% 
  as_draws_df() %>% 
  as_tibble() %>%  # ???????????????????????? tibble
  select(mu, 
         beta.1 = `beta[1]`, beta.2 = `beta[2]`, beta.3 = `beta[3]`,
         beta.4 = `beta[4]`, beta.6 = `beta[6]`, beta.7 = `beta[7]`)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
long_draws <- draws_df %>%
  mutate(
    `1_Original`  = beta.1*w_orig[1] + beta.2*w_orig[2] + beta.3*w_orig[3] + 
      beta.4*w_orig[4] + beta.6*w_orig[5] + beta.7*w_orig[6],
    
    `2_Equal`     = beta.1*w_equal[1] + beta.2*w_equal[2] + beta.3*w_equal[3] + 
      beta.4*w_equal[4] + beta.6*w_equal[5] + beta.7*w_equal[6],
    
    `3_Tree_Only` = beta.1*w_tree[1] + beta.2*w_tree[2] + beta.3*w_tree[3] + 
      beta.4*w_tree[4] + beta.6*w_tree[5] + beta.7*w_tree[6],
    
    `4_Vert_Only` = beta.1*w_vert[1] + beta.2*w_vert[2] + beta.3*w_vert[3] + 
      beta.4*w_vert[4] + beta.6*w_vert[5] + beta.7*w_vert[6]
  ) %>%
  select(mu, `1_Original`, `2_Equal`, `3_Tree_Only`, `4_Vert_Only`) %>%
  pivot_longer(cols = starts_with(c("1", "2", "3", "4")), 
               names_to = "Scheme", 
               values_to = "beta_w")
# print(colnames(long_draws))

# -------------------------------------------------------------------
# -------------------------------------------------------------------
temp_range <- seq(from = 0, to = 5.5, by = 0.1)

pred.reg.df_expanded <- long_draws %>%
  # ?????? crossing ?????? Scheme, mu, beta_w ??? temp_range ??????????????????
  crossing(P.Ind = temp_range) %>%
  # ????????????????????? plogis(mu + beta_w * x) ????????? inv_logit
  mutate(pred_risk = plogis(mu + beta_w * P.Ind)) %>%
  # ??????????????????
  group_by(Scheme, P.Ind) %>%
  mean_qi(pred_risk, .width = 0.95)

# -------------------------------------------------------------------
# -------------------------------------------------------------------
ggplot(pred.reg.df_expanded, aes(x = P.Ind, y = pred_risk, color = Scheme, fill = Scheme)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.1, color = NA) +
  geom_line(size = 1.2) +
  scale_color_viridis_d(option = "magma", end = 0.8) +
  scale_fill_viridis_d(option = "magma", end = 0.8) +
  labs(x = "Temperature Rise (C)", y = "Extinction Risk", title = "Weighting Sensitivity Analysis") +
  theme_classic()

regions <- c("Region1", "Region2", "Region3", "Region4", "Region6", "Region7")
leave_one_out_results <- list()

# -------------------------------------------------------------------
# Africa, Asia, Australia, Europe, N.America, S.America
regions_names <- c("Africa", "Asia", "Australia", "Europe", "N.America", "S.America")

# -------------------------------------------------------------------
loo_results <- list()

# -------------------------------------------------------------------
for(i in 1:6) {

  w_temp <- cont.rich 
  w_temp[i] <- 0           
  w_temp <- w_temp / sum(w_temp) 
  
  loo_results[[i]] <- draws_df %>%
    as_tibble() %>%
    mutate(beta_w = beta.1*w_temp[1] + beta.2*w_temp[2] + beta.3*w_temp[3] + 
             beta.4*w_temp[4] + beta.6*w_temp[5] + beta.7*w_temp[6]) %>%
    crossing(P.Ind = seq(0, 5, 0.1)) %>%
    mutate(pred_risk = plogis(mu + beta_w * P.Ind)) %>%
    group_by(P.Ind) %>%
    mean_qi(pred_risk) %>%
    # ??????????????????????????????????????????????????????
    mutate(Excluded_Region = regions_names[i])
}

# 4. ????????????
loo_data <- bind_rows(loo_results)

# 5. ?????????????????????????????????????????????
ggplot(loo_data, aes(x = P.Ind, y = pred_risk, color = Excluded_Region)) +
  geom_line(linewidth = 1) +
  # ????????????????????????????????????????????????????????????
  geom_line(data = pred.reg.df_expanded %>% filter(Scheme == "1_Original"),
            aes(x = P.Ind, y = pred_risk), 
            color = "black", linetype = "dashed", linewidth = 1.2) +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Regional Sensitivity (Leave-One-Out)",
    subtitle = "Dashed black line = All regions; Colored lines = One region excluded",
    x = "Temperature Rise (C)",
    y = "Global Predicted Risk",
    color = "Excluded Region"
  ) +
  theme_classic()

# -------------------------------------------------------------------

???
target_temps <- c(1.5, 3.0)


violin_data <- long_draws %>%
  crossing(P.Ind = target_temps) %>%
  mutate(pred_risk = plogis(mu + beta_w * P.Ind)) %>%
  mutate(Temp_Label = paste0(P.Ind, "??C Rise"))


library(ggplot2)

p_violin <- ggplot(violin_data, aes(x = Scheme, y = pred_risk, fill = Scheme)) +
  geom_violin(alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.1, fill = "white", color = "gray30", outlier.shape = NA) +
  facet_wrap(~Temp_Label, scales = "free_y") + 
  scale_fill_viridis_d(option = "mako", end = 0.8) +
  labs(
    title = "Distribution of Extinction Risk at Key Temperature Thresholds",
    subtitle = "Comparing posterior uncertainty across weighting schemes",
    y = "Predicted Extinction Risk",
    x = "Weighting Scheme"
  ) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10), 
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none" # X??????????????????????????????????????????
  )

print(p_violin)
