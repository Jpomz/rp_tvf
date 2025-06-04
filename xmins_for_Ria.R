# Example of "trimming" undersampled data for Ria

library(sizeSpectra)
library(poweRlaw)
library(tidyverse)

# the example data set is from one of my students  
# It looks at macroinvertebrate body sizes at two sites
# one site is going to be restored in the future, and the other site is an upstream control  
# The restoration has not yet happened, so this is just the "before" data for a BACI study design. 

# Load the data (csv was included in the email I sent you)
### Change the file path in the quotes as needed
dat <- read_csv("dat_dw.csv") #### CHANGE THIS ####

head(dat)

# the dw column is the estimated dry weight of each individual. 
#site_number == 01 is the upstream reference
#site_number == 03 is the restoration site

# clean up the data and remove NA's in the dw column

dat <- dat |>
  filter(!is.na(dw))

# custom function to perform  MLE on "tidy" data
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

# estimate the isd exponent
raw_mle_lambda <- dat %>%
  group_by(site_number) %>% 
  nest() %>% # this creates "list-columns"
  # the nest() function makes it easier to plug this data set into the custom function written above
  mutate(lambda = map(data, # this estimates lambda w/ MLE
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup() 

# this will probably throw an error about:
# nlm(); NA/Inf replaced by max positive value  
# you can ignore this  

# What are our estimates?

raw_mle_lambda

# both sites are basically the same, and much "shallower" than would be expected  
# when looking at individual-level data, lambda should be ~-1.95 to -2 (Andersen and Beyer 2006) 

# plot the "raw" results:
raw_mle_lambda %>%
  ggplot(aes(x = site_number, 
             y = b, # "lmabda"
             ymin = minCI, # CI's for lambda
             ymax = maxCI,
             color = site_number)) +
  geom_pointrange(
    size = 1) +
  scale_color_manual(values = c("#019AFF", "#FF1984")) +
  theme_bw() +
  labs(y = "Estimated \U03BB") +
  NULL

### Trim data ####
# the following code uses the "xmin" method of Clauset et al. 2009. 
# Basically, it estimates the tail of the data which best fits a power law  
# this uses functions from the poweRlaw package

# 1) first, make a list of data sets  
dat_list <- dat %>% 
  group_by(site_number) %>% 
  group_split()

# 2) create empty list to fill with xmins
xmin_list = list() 
xmin_list

# 3) get list of xmins for each sample
# use set seed for reproducibility
# there is some slight variation if you use different seeds due to random sampling  
# I haven't tested this extenisively, but it seems to generally estimate similar values of xmin
set.seed(202002)
for(i in 1:length(dat_list)){
  powerlaw = conpl$new(dat_list[[i]]$dw) # get power law estimate from poweRlaw package
  xmin_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin, # extract the xmin from the poweRlaw package
                          site_number = unique(dat_list[[i]]$site_number))
}
# note that this estimates power law exponents, but they are in a different form than the exponenets discussed in Edwards et al. 2017, and, to be honest, I don't really know how to interpret them or how to "convert" them to estimates as in Edwards  


# get the xmin for each site
xmins_clauset = bind_rows(xmin_list)

# so this is saying that all the body sizes < 0.0746 in site 01, and all the body sizes < 0.196 in site 03 are undersampled. 

# So now, we "trim" the data and only keep the body sizes which are greater than those values. 

dat <- left_join(dat, xmins_clauset)

# filter all dw that are greater than the new `xmin_clauset` column

dat_trimmed <- dat |>
  filter(dw > xmin_clauset)

# re-estimate lambda  
# actually estimate the isd exponent
trimmed_mle_lambda <- dat_trimmed %>%
  group_by(site_number) %>% 
  nest() %>% # this creates "list-columns"
  mutate(lambda = map(data, # this estimates lambda w/ MLE
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup() 

trimmed_mle_lambda

# now the estimates are not the same for each site, and, tellingly, they are closer to the "expected" value of ~-2

# plot the "raw" and "trimmed" results for comparison:
trimmed_mle_lambda %>%
  mutate(data = "trimmed") |>
  bind_rows(raw_mle_lambda |>
              mutate(data = "raw")) |>
  ggplot(aes(x = data, 
             y = b, # "lmabda"
             ymin = minCI, # CI's for lambda
             ymax = maxCI,
             color = site_number)) +
  geom_pointrange(
    size = 1,
    position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("#019AFF", "#FF1984")) +
  theme_bw() +
  labs(y = "Estimated \U03BB")  +
  NULL

# now, the question is, how do we know if this is "right"?
# well, we don't
# however, I'm currently working on a simulation for a new manuscript which takes data from a known power law and artificially undersamples observations below some value  
# this is still a work in progress, but it generally seems like undersampled data is "shallower" and closer to -1 instead of -2 (like in this example)
# It also seems like the undersampled data washes out known differences in lambda  
# in other words, when I sample data from a lambda = -1 and a different data set with lambda = -2 and then undersample the small body sizes, the lmabda_estimates are similar to each other, and the differences between the groups is no longer detectable, or not as extreme. 
# anyways, we are now recommending trimming the data. Hopefully I'll have a draft of that new paper this summer and can send it along when it's ready if you're interested. 