# Arkansas data processed at Colorado Mesa University
# Mostly by Nia Taubr
# Other help from:
# Andy Stack, Noah Enoch, Carly Scheck

# mle tidy

# `sizeSpectra` is not hosted on CRAN
# to install the package directly from github, you need to have the `devtools` package installed. Run the following once if you don't already have it downloaded
# install.packages("devtools") 

#To install the latest version of sizeSpectra, run the following:
# devtools::install_github("andrew-edwards/sizeSpectra")

library(sizeSpectra)
library(tidyverse)


# custom function 
MLE_tidy <- function(df, rsp_var){
  # define variables
  x <- df[[rsp_var]]
  xmin = min(x)
  xmax = max(x)
  log.x = log(x)
  sum.log.x = sum(log.x)
  
  # initial starting point for parameter estimate
  PL.bMLE = 1/(log(min(x)) - sum.log.x/length(x)) - 1
  
  # non-linear minimization  
  PLB.minLL = nlm(negLL.PLB, 
                  p = PL.bMLE, x = x, n = length(x), 
                  xmin = xmin, xmax = xmax,
                  sumlogx = sum.log.x)
  
  # estimate for b
  PLB.bMLE = PLB.minLL$estimate
  # minimum estimate of b
  PLB.minNegLL = PLB.minLL$minimum
  
  ## 95% CI calculation
  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 1e-05)
  PLB.LLvals = vector(length = length(bvec))
  for (i in 1:length(bvec)) {
    PLB.LLvals[i] = negLL.PLB(bvec[i],
                              x = x,
                              n = length(x), 
                              xmin = xmin,
                              xmax = xmax,
                              sumlogx = sum.log.x)
  }
  critVal = PLB.minNegLL + qchisq(0.95, 1)/2
  bIn95 = bvec[PLB.LLvals < critVal]
  # confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  if (PLB.MLE.bConf[1] == min(bvec) | 
      PLB.MLE.bConf[2] == max(bvec)) {
    dev.new()
    plot(bvec, PLB.LLvals)
    abline(h = critVal, col = "red")
    stop("Need to make bvec larger - see R window")
  }
  # return b estimate and min/max 95% CI
  return(data.frame(b = PLB.bMLE,
                    minCI = min(bIn95),
                    maxCI = max(bIn95)))
}

# read in dry weight data
dat <- readRDS("cmu/derived_data/cmu_ark_dw.RDS")
dat %>%
  distinct(year) 

# how many individuals by site and year?
dat %>%
  group_by(site, year) %>%
  count() %>%
  arrange(year)

# how many individuals per sample?
dat %>%
  group_by(site, rep, year) %>%
  count() 

# plot of the count of individuals by sample
dat %>%
  group_by(site, rep, year) %>%
  count() %>%
  ggplot(aes(x = year,
             y = n, 
             color = site,
             shape = rep)) +
  geom_point(
    position = position_jitter(
      width = 0.1,
      height = 0
    ))

# actually estimate the isd exponent
mle_lambda_rep <- dat %>%
  group_by(year, site, rep) %>% # site_number # NOT rep/sample
  nest() %>%
  mutate(lambda = map(data, # lambda = exponent/slope
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup() %>%
  mutate(year_0 = year - min(year)) # get rid 


mle_lambda_rep %>%
  ggplot(aes(x = year, 
             y = b,
             ymin = minCI, 
             ymax = maxCI,
             color = site, 
             group = rep)) +
  geom_pointrange(
    size = 1,
    position = position_dodge(width = 0.75)
  ) +
  scale_color_manual(values = c("#019AFF", "#FF1984")) +
  scale_fill_manual(values = c("#019AFF", "#FF1984")) +
  stat_smooth(method = "lm",
              inherit.aes = FALSE,
              aes(x = year, y = b,
                  color = site,
                  fill = site),
              alpha = 0.15) +
  theme_bw() +
  labs(y = "Estimated \U03BB") +
  NULL
ggsave("plots/cmu_showcase_mle.png",
       units = "in",
       width = 8,
       height = 6)


# old ordinary least squares regression
# don't need to do this
ols <- lm(b~site*year_0, dat = mle_lambda_rep)
summary(ols)
ggplot(ols,
       aes(x = .fitted,
           y = .resid)) +
  geom_point()

# weighted regression 
# gamma = SE
# gamma = (b_high - b_low) / (2 * 1.96)

mle_lambda_rep %>%
  ggplot(aes(sample = b)) +
  stat_qq()+ 
  stat_qq_line()

mle_lambda_rep <- mle_lambda_rep %>% 
  mutate(se = (maxCI - minCI) / 2 * 1.96,
         var = se**2)

summary(lm(b~site*year_0, dat = mle_lambda_rep, weights = 1 / var))

ggplot(mle_lambda_rep,
       aes(x = year_0,
           y = b,
           ymin = minCI,
           ymax = maxCI,
           color = site)) +
  geom_pointrange() +
  geom_smooth(method = "lm")


mle_lambda_rep %>%
  group_by(site, year) %>%
  summarize(mean_b = mean(b)) %>%
  ggplot(aes(x = year,
             y = mean_b,
             color = site)) +
  geom_point() +
  geom_smooth(method = "lm")

mle_lambda_rep %>%
  #group_by(site, year) %>%
  #summarize(mean_b = mean(b)) %>%
  mutate(year_c = year - mean(year)) %>%
  lm(b ~ site*year_c, data = .) %>%
  summary()


# combine reps ------------------------------------------------------------



# estimate lambda exponent of a power law
# N ~ M^b
# N = number of individuals
# M = individual body size
# b = lambda, exponent describing power law
## This is unknown parameter that we are estimating 
mle_lambda <- dat %>%
  #filter(dw >0.0026) %>%
  group_by(site, year) %>%
  nest() %>%
  mutate(lambda = map(data,
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup()

mle_lambda <- mle_lambda %>% 
  mutate(se = (maxCI - minCI) / 2 * 1.96,
         var = se**2,
         year_0 = year - min(year))

weighted_ols <- lm(b~site*year_0,
                   dat = mle_lambda,
                   weights = 1 / var) # rerun with new weights
summary(weighted_ols)

ggplot(mle_lambda,
       aes(x = year,
           y = b,
           ymin = minCI,
           ymax = maxCI,
           color = site)) +
  geom_pointrange() +
  scale_color_manual(values = c("#019AFF", "#FF1984")) +
  geom_smooth(method = "lm", 
              aes(weight = 1/var),
              alpha = 0.15) +
  labs(y = "Estimated \U03BB",
       x = "Year") +
  theme_bw()
ggsave("plots/mle_weighted_ols_Apr_2024.png",
       units = "in",
       width = 7,
       height = 6)