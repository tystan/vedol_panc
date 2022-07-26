
# ---- libs ----

library("arrow") # reading in parquet data
library("dplyr") # piping and manipulation of data
library("tibble") # prettier data.frames
library("ggplot2") # plotting
library("knitr") # to print kable(.) tables
library("tidyr") # pivot_wider() function
library("pharmsignal") # disproportionality statistics 


# ---- load ----

# load set_membership data set that has set identifiers for each id

# read in line data using arrow::read_parquet()
# (data saved in parquet format for speed/storage size)
# takes ~ 2 sec on i5-8400/32GB@2133MHz(CL15)/500GB 970 EVO Plus
system.time({
  # set_membership
  set_mem <-
    read_parquet(
      file = "data/set_mem.parquet"
    )
})

# check data is in tibble format
is_tibble(set_mem)
class(set_mem)
nrow(set_mem)


# ---- make_signal_data ----

# standardise the summary data once the dataset has been filtered as required
mk_signal_data <- function(dat, comparator_str) {
  
  dat %>%
    mutate(
      comparator = comparator_str, 
      med = if_else(med == "[vedol]", med, "not [vedol]")
    ) %>%
    select(comparator, exposure = med, outcome = reac) %>%
    group_by(comparator, exposure, outcome) %>%
    summarise(n = n(), .groups = "keep") %>%
    ungroup() %>%
    arrange(comparator, exposure, outcome) 
  
}

# reports before 2015 and having [panc] indication are not required
set_mem <-
  set_mem %>% 
  dplyr::filter(ind == "not [panc]", date == ">=2015")

# these are the medicine classifications
### NOTE:
# "[IRA]" excludes vedolizumab
# "[mAb]" excludes vedolizumab and natalizumab
with(set_mem, table(med, useNA = "ifany"))
# med
# [IRA]   [mAb] [other] [vedol] 
# 35486  522715 6509008   19034 

# turn row data into summarised counts for analysis by sequentially
# filtering the data for each comparator of interest and then combine
signal_data <-
  set_mem %>%
  mk_signal_data(., "(1) All")

signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(med %in% c("[vedol]", "[IRA]", "[mAb]")) %>%
      mk_signal_data(., "(2) mAbs")
  )

signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(med %in% c("[vedol]", "[IRA]")) %>%
      mk_signal_data(., "(3) IRAs")
  )

signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(ind_ibd == "[IBD]") %>%
      mk_signal_data(., "(4) All IBD indi")
  )

signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(ind_cd == "[CD]") %>%
      mk_signal_data(., "(4a) CD indi")
  )

signal_data <-
  bind_rows(
    signal_data,
    set_mem %>%
      dplyr::filter(ind_cu == "[CU]") %>%
      mk_signal_data(., "(4b) CU indi")
  )

# free up memory of patient level data
rm(list = "set_mem")
gc()

# ---- wrangle ----


# have a look
signal_data %>% kable(.)

# create a,b,c,d cell counts as columns  
signal_data_wide <-
  signal_data %>%
  pivot_wider(
    names_from = c("exposure", "outcome"), 
    values_from = "n",
    names_sep = "_"
  )

# have a look
signal_data_wide

# edit exposure and outcome values for shorter expressions
# they now take form: "ExOy" = `Exposure <x> Outcome <y>` where
# <x> = "p" (positive) if exposure == "vedol", "n" (negative) otherwise
# <y> = "p" (positive) if outcome == "panc", "n" (negative) otherwise
colnames(signal_data_wide) <- gsub("\\[vedol\\]", "[E]", colnames(signal_data_wide))
colnames(signal_data_wide) <- gsub("\\[panc\\]", "[O]", colnames(signal_data_wide))
colnames(signal_data_wide) <- gsub("not \\[([EO])\\]", "[\\1]n", colnames(signal_data_wide))
colnames(signal_data_wide) <- gsub("\\](_|$)", "]p\\1", colnames(signal_data_wide))
# rm punctuation greedily
colnames(signal_data_wide) <- gsub("(\\[|\\]|_)", "", colnames(signal_data_wide)) 

# have a look
signal_data_wide

# ---- generate_stats ----

# calculate disproportionality statistics seen in tab 2/fig 1 in paper
signal_tab <-
  with(
    signal_data_wide, 
    bcpnn_mcmc_signal(
      a = EpOp,
      b = EpOn,
      c = EnOp,
      d = EnOn
    )
  )

# add analysis names to results
signal_tab <-
  bind_cols(
    Comparator = signal_data_wide$comparator,
    signal_tab
  )

# take from log2(ratio) scale to ratio scale
signal_tab <-
  signal_tab %>% 
  mutate(
    ci_lo     = 2^ci_lo,  
    ci_hi     = 2^ci_hi, 
    est       = 2^est, 
    est_scale = "orig scale"
  ) 

# ---- table2 -----

table2 <-
  signal_tab %>%
  rename(
    `N vedol and panc` = n11, 
    `N vedol` = `drug margin`
  ) %>%
  mutate(
    `RSIC = 2^IC` = sprintf("%3.2f (%3.2f, %3.2f)", est , ci_lo, ci_hi),
    `Significant` = (ci_lo > 1) | (ci_hi < 1),
    `RSIC = 2^IC` = paste0(`RSIC = 2^IC`, if_else((est > 2) & Significant, "*", "")),
    `N comparator and panc` = `event margin` - `N vedol and panc`,
    `N comparator` = `n..` - `N vedol`
  ) %>%
  dplyr::select(
    Comparator, 
    `N vedol and panc`, 
    `N vedol`, 
    `N comparator and panc`, 
    `N comparator`, 
    `RSIC = 2^IC`
  )  

# print table 2
table2 %>%
  kable(.)


# ---- figure1 -----

# get min and max values for the y-axis
ymi <- min(signal_tab[["ci_lo"]])
yma <- max(signal_tab[["ci_hi"]])

# create colour palette for comparators
pal_use <- "Tableau 10" 
lvls  <- sort(unique(signal_tab[["Comparator"]]))
pal_comp <- palette.colors(n = length(lvls), palette = pal_use) 
names(pal_comp) <- lvls

# use ggplot to plot disproportionality statistics
fig1 <-
  signal_tab %>%
  ggplot(aes(x = Comparator, y = est, col = Comparator)) %+%
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, alpha = 0.5) %+%
  geom_point(size = 2) %+%
  geom_hline(yintercept = 1) %+%
  geom_hline(yintercept = 2, linetype = 2) %+%
  scale_color_manual(values = pal_comp) %+%
  scale_y_continuous(
    trans = "log2", 
    limits = c(min(1, ymi), max(1, yma))
  ) %+%
  labs(
    x = "Comparator", 
    y = "RSIC = 2^IC estimate using BCPNN MCMC\n(ratio scale with 95% CI)", 
    col = "Comparator"
  ) %+%
  theme_bw() %+%
  theme(
    text = element_text(family = "serif"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) 

# have a look
fig1

# save plot as png
ggsave(fig1, filename = "fig/fig1.png", width = 5, height = 4)
