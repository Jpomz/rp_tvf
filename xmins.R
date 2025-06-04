# ruth powell TVF start
# Nia Taubr student grant awardee

library(sizeSpectra)
# install.packages("poweRlaw")
library(poweRlaw)
library(tidyverse)

# custom function
# custom function 
MLE_tidy <- function(df, rsp_var){
  # define variables
  x <- df[[rsp_var]]
  x_n <- length(x)
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
                    maxCI = max(bIn95),
                    x_n = x_n))
}

# make sure path to your data is correct
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

(xmins_clauset = bind_rows(xmin_list))

dat <- left_join(dat, xmins_clauset) |>
  filter(!is.na(dw))


# actually estimate the isd exponent
mle_lambda <- dat %>%
  filter( dw >= xmin_clauset) |>
  group_by(site_number) %>% # site_number # NOT rep/sample
  nest() %>%
  mutate(lambda = map(data, # lambda = exponent/slope
                      MLE_tidy,
                      "dw")) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup() 


mle_lambda %>%
  ggplot(aes(x = site_number, 
             y = b,
             ymin = minCI, 
             ymax = maxCI,
             color = site_number)) +
  geom_pointrange(
    size = 1,
    position = position_dodge(width = 0.75)
  ) +
  scale_color_manual(values = c("#019AFF", "#FF1984")) +
  scale_fill_manual(values = c("#019AFF", "#FF1984")) +
  theme_bw() +
  labs(y = "Estimated \U03BB") +
  NULL


mle_dat <- dat %>%
  filter( dw >= xmin_clauset)

mle_dat |>
  group_by(site_number) |>
  count()

xmins_clauset

dat |>
  group_by(site_number) |>
  summarize(min = min(dw), 
            max = max(dw))

dat |>
  mutate(fill = dw >= xmin_clauset) |>
  #filter(dw < 7) |>
  ggplot(aes(x = (dw),
             fill = fill)) +
  geom_histogram(binwidth = 0.025) +
  facet_wrap(~site_number,
             scales = "free_x") +
  #scale_y_log10() +
  scale_x_log10() +
  NULL



# waterfall plot
mle_dat |>
  group_by(site_number) |>
  summarize(min = min(dw), 
            max = max(dw)) |>
  left_join(mle_lambda)
xmins_clauset
mle_lambda



# n_line_01 <- 1000
# #n_line_02 <- 156
# line_dat <- data.frame(
#   line_01 = (1 - pPLB(x = seq(from = 0.196,
#                          to = 76.7,
#                          length.out = n_line_01),
#                  b = -1.83, xmin = 0.196, xmax = 76.7)) * n_line_01,
#   line_03 = (1 - pPLB(x = seq(from = 0.196,
#                              to = 76.7,
#                              length.out = n_line_01),
#                      b = -2.14, xmin = 0.196, xmax = 76.7)) * n_line_01,
#   x = seq(from = 0.196,
#               to = 76.7,
#           length.out = n_line_01)) |>
#   pivot_longer(1:2, names_to = "line", values_to = "y") |>
#   group_by(line) %>%
#   arrange(line, -x) %>%
#   mutate(order = row_number())


x_01 <- seq(from = 0.196,
            to = 76.7,
            length.out = 806)
line_01 = (1 - pPLB(x = x_01,
                    b = -1.83, xmin = 0.196, xmax = 76.7)) * 806

x_02 <- seq(from = 0.196,
            to = 6.86,
            length.out = 806)
line_02 = (1 - pPLB(x = x_02,
                    b = -2.14, xmin = 0.196, xmax = 6.86)) * 806

line_01_df <- data.frame(x = x_01,
                         y = line_01,
                         line = "01")
line_02_df <- data.frame(x = x_02,
                         y = line_02,
                         line = "03")
line_dat <- bind_rows(line_01_df, line_02_df)
  
mle_dat %>% 
  filter(dw >= 0.196) |>
  group_by(site_number) %>%
  sample_n(1000, replace = TRUE) |>
  arrange(site_number, -dw) %>%
  mutate(order = row_number()) %>%
  ggplot(aes(x = dw,
             y = order,
             #size = body_mass,
             color = site_number)) +
  geom_point(shape = 1, alpha = 0.5) +
  #theme_dark() +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  labs(y = "Number of values \u2265 x",
       x = "Individual body mass") +
  geom_line(data = line_dat, 
            aes(x = x, 
                y = y, 
                color = line),
            inherit.aes = FALSE,
            linewidth = 3, 
            alpha = 0.5) +
  ### add CI numbers here ###
  annotate(geom = "text", x = 10, y = 0.8, label = "CI") +
  ### add another annotate with other CI numbers here ###
  NULL


# figure out log-log plots here ####
bin_and_center <- function(data, var, breaks, ...){
  # data is a data frame
  # var is a string, and is the name of a column in data which you want to bin
  # breaks controls the number of bins as defined in hist() 
  # See ?hist for details
  
  # bin values using hist()
  binned_hist = hist(data[[var]], 
                     breaks = breaks, # need to predefine breaks
                     # e.g. Log2 breaks = 2^seq(min, max) 
                     # Log10 breaks = 10^seq(min, max)
                     include.lowest = TRUE, plot = FALSE)
  # calculate "left" and "right" edge of bins
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]
  # total bin width = right edge - left edge
  break_width = breaks_offset - breaks_orig
  count = binned_hist$counts 
  dataout = data.frame(
    # normalize counts =count / width (White et al 1997)
    log_count_corrected = log10(count / break_width),
    # original midpoint of bin log10 transformed
    log_mids = log10(binned_hist$mids),
    log_mids_center = NA)
  # remove bins with 0 counts
  # -Inf comes from log10(count / break_width) above
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  # recenter data at x=0
  mid_row = ceiling(nrow(dataout)/2)
  # subtract value of mid row from all mids
  dataout$log_mids_center = 
    dataout[,"log_mids"] - dataout[mid_row,"log_mids"]
  dataout
}

breaks <- 2^seq(floor(log2(min(dat$dw))),
                ceiling(log2(max(dat$dw))))
breaks

log_lambda <- dat %>%
  filter( dw >= xmin_clauset) |>
  group_by(site_number) %>% # site_number # NOT rep/sample
  nest() %>%
  mutate(lambda = map(data, # lambda = exponent/slope
                      bin_and_center,
                      var = "dw",
                      breaks = breaks)) %>%
  unnest(cols = lambda) %>%
  select(-data) %>%
  ungroup() 

log_lambda |>
  ggplot(aes(x = log_mids, 
             y = log_count_corrected,
             color = site_number,
             shape = site_number,
             fill = log_mids,
             group = as.factor(site_number))) +
  geom_point(size = 5,
             stroke = NA) +
  scale_fill_viridis_c(option = "turbo") +
  stat_smooth(aes(x = log_mids, 
                  y = log_count_corrected,
                  color = site_number),
              inherit.aes = FALSE,
              method = "lm",
              se = FALSE) +
  scale_shape_manual(values = c(21, 24)) +
  scale_color_manual(values = c("red", "black")) +
  theme_bw() + 
  labs(y = "Log(Count)",
       x = "Log(Body Size)") +
  guides(fill = "none")



summary(lm(log_count_corrected ~ log_mids_center * as.factor(site_number), data = log_lambda))
