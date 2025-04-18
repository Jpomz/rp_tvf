# ruth powell TVF start
# Nia Taubr student grant awardee

library(sizeSpectra)
library(poweRlaw)
library(tidyverse)

dat <- read_csv("dat_dw.csv")

dat_list <- dat %>% 
  filter(!is.na(dw)) |>
  group_by(site_number) %>% 
  group_split()

# 3) create empty list to fill with xmins
xmin_list = list() 

# 4) get list of xmins for each sample
set.seed(202002)
for(i in 1:length(dat_list)){
  powerlaw = conpl$new(dat_list[[i]]$dw) # get power law estimate from poweRlaw package
  xmin_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin, # extract the xmin from the poweRlaw package
                          site_number = unique(dat_list[[i]]$site_number))
}

xmins_clauset = bind_rows(xmin_list)
