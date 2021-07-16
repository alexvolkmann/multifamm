# load libraries
library(tidyverse)
library(data.table)

# load the data
load("../snooker.RData")


# transform it to long format
snooker_flmm <- snooker %>%
  rbind(mutate(., y_raw = x_raw, y_new = x_new)) %>%
  mutate(x_raw = NULL, x_new = NULL,
         dim = rep(c("y", "x"), each = nrow(snooker)))

# downsize the data and filter it
snooker_flmm <- snooker_flmm %>%
  filter(location != "queue") %>%
  filter(!is.na(t_optim))
  #filter(coarsen_relSSE > 1e-4)

# rename variables as needed
snooker_flmm <- snooker_flmm %>%
  transmute(y_vec = y_new,
            t = t_optim,
            n_long = as.integer(ID),
            subject_long = as.integer(person),
            word_long = as.integer(run),
            dim = factor(paste(location, dim, sep = ".")),
            covariate.1 = as.numeric(leistungsstufe > 2),
            covariate.2 = as.numeric(group == "Interventionsgruppe"),
            covariate.3 = (as.integer(run)-1),
            covariate.4 = as.numeric(group == "Interventionsgruppe" &
                                       run == 2),
            relSSE = coarsen_relSSE) %>%
  group_by(word_long, subject_long, dim) %>%
  mutate(combi_long = as.integer(factor(n_long))) %>%
  ungroup()

# convert to data.table
snooker_flmm <- snooker_flmm %>% data.table

snooker <- snooker_flmm %>%
  filter(relSSE > 3e-3) %>%
  data.table() %>%
  select(-relSSE)
usethis::use_data(snooker)
