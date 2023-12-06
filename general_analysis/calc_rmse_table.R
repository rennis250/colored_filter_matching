library(tidyverse)

# look at our data ####
rmse <- read_csv("rmse.txt")

# rmse_fmin,rmse_rigid,rmse_affinemap,cspace,chroma_or_full,illum,body_color,light_or_dark,mask_name,ref_lev,exclude_darkpxs

rmse <- mutate(rmse,
    aav_rigid = 100 * (rmse_identity - rmse_rigid) / rmse_identity,
    aav_affinemap = 100 * (rmse_identity - rmse_affinemap) / rmse_identity,
    aav_fmin = 100 * (rmse_identity - rmse_fmin) / rmse_identity,
    aav_dzmura = 100 * (rmse_identity - rmse_dzmura) / rmse_identity
)

rmse %>%
    group_by(cspace, ref_lev) %>%
    summarise(
        fmin_avg = mean(aav_fmin), fmin_sd = sd(aav_fmin),
        dzmura_avg = mean(aav_dzmura), dzmura_sd = sd(aav_dzmura),
        rigid_avg = mean(aav_rigid), rigid_sd = sd(aav_rigid),
        affinemap_avg = mean(aav_affinemap), affinemap_sd = sd(aav_affinemap),
    ) %>%
    print(n = Inf)

rmse %>%
    group_by(cspace, ref_lev) %>%
    summarise(
        fmin_avg = mean(aav_fmin), fmin_sd = sd(aav_fmin),fmin_min = min(aav_fmin), fmin_max = max(aav_fmin),
        dzmura_avg = mean(aav_dzmura), dzmura_sd = sd(aav_dzmura), dzmura_min = min(aav_dzmura), dzmura_max = max(aav_dzmura)
    ) %>%
    print(n = Inf)

# look at it for all factors averaged together
rmse %>%
    group_by(cspace) %>%
    summarise(
        fmin_avg = mean(aav_fmin), fmin_sd = sd(aav_fmin),
        dzmura_avg = mean(aav_dzmura), dzmura_sd = sd(aav_dzmura),
    ) %>%
    print(n = Inf)

# let's break it down by many factors
# look at it for all factors averaged together
rmse %>%
    group_by(cspace, ref_lev, chroma_or_full, exclude_darkpxs) %>%
    summarise(
        fmin_avg = mean(aav_fmin), fmin_sd = sd(aav_fmin),
        dzmura_avg = mean(aav_dzmura), dzmura_sd = sd(aav_dzmura),
    ) %>%
    print(n = Inf)

# compare to d'zmura results ####
rmse_dz <- read_csv("rmse_dzmura.txt")

# rmse_fmin,rmse_rigid,rmse_affinemap,cspace,chroma_or_full,obs

rmse_dz <- mutate(rmse_dz,
    aav_rigid = 100 * (rmse_identity - rmse_rigid) / rmse_identity,
    aav_affinemap = 100 * (rmse_identity - rmse_affinemap) / rmse_identity,
    aav_fmin = 100 * (rmse_identity - rmse_fmin) / rmse_identity
)

rmse_dz %>%
    group_by(cspace) %>%
    summarise(
        fmin_avg = mean(aav_fmin), fmin_sd = sd(aav_fmin),
        rigid_avg = mean(aav_rigid), rigid_sd = sd(aav_rigid),
        affinemap_avg = mean(aav_affinemap), affinemap_sd = sd(aav_affinemap),
    )

# differences from before
