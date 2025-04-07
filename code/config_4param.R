
library(dplyr)

#list components
pb_age <- "pb(logAge_days, control = pb.control(order = 3))"
pb_ageSex <- "pb(sexMale_x_logAge, control = pb.control(order = 3))"
sex <- "sexMale"
rand_site <- "random(study_site)"
age <- "logAge_days"
ageSex <- "sexMale_x_logAge"
f_site <- "study_site"
int <- "1"


mu_list <- list(paste(pb_age, pb_ageSex, sex, rand_site, sep=" + "),
                paste(pb_age, pb_ageSex, sex, rand_site, "freesurfer", sep=" + ")
)

sig_list <- list(paste(pb_age, sex, rand_site, sep=" + "),
                 paste(pb_age, ageSex, sex, rand_site, sep=" + "),
                 paste(pb_age, pb_ageSex, sex, rand_site, sep=" + "),
                 paste(pb_age, sex, rand_site, "freesurfer", sep=" + "),
                 paste(pb_age, ageSex, sex, rand_site, "freesurfer", sep=" + "),
                 paste(pb_age, pb_ageSex, sex, rand_site, "freesurfer", sep=" + ")
                 
)

#nu - allow a bit more complexity, nonlinear terms/random effects
nu_elements <- c(pb_age, pb_ageSex, sex, rand_site, age, ageSex, f_site, "freesurfer")
tau_elements <- c(sex, age, ageSex, f_site, "freesurfer")

# Get all combinations
nu_combos <- unlist(
  lapply(1:length(nu_elements), function(k) {
    combn(nu_elements, k, simplify = FALSE)
  }),
  recursive = FALSE
)

tau_combos <- unlist(
  lapply(1:length(tau_elements), function(k) {
    combn(tau_elements, k, simplify = FALSE)
  }),
  recursive = FALSE
)

filt_duplicates <- function(x){(
    !(pb_age %in% x && age %in% x) 
    && !(pb_ageSex %in% x && ageSex %in% x)
    && !(rand_site %in% x && f_site %in% x)
  )}

# Remove nonsense configs
nu_filt_combos <- Filter(filt_duplicates, nu_combos)

nu_list <- lapply(nu_filt_combos, function(x) paste(x, collapse = " + ")) %>% unlist()
nu_list <- c("1", nu_list)

tau_list <- lapply(tau_combos, function(x) paste(x, collapse = " + ")) %>% unlist()
tau_list <- c("1", tau_list)

#add intercept-only
mod_combos <- expand.grid(mu_list, sig_list, nu_list, tau_list) %>%
  filter(!(Var3 == '1' & Var4 != '1')) %>% #keep nu more complex than tau
  rename("mu" = "Var1",
         "sigma" = "Var2",
         "nu" = "Var3",
         "tau" = "Var4")
saveRDS(mod_combos, file="four_param_model_combos.RDS")
