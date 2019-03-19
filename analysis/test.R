# ratios

# read in the weekly data
# df <- readRDS("data/jungle_rhythms_weekly_annotations.rds")
# df <- df[which(df$value != 0),]
# df$species_full <- paste(df$genus, df$species, sep = " ")

phenophase_selected <- "flowers"

# this is now per id
# counting unique year - id combinations for either
# the full or reduced (phenophase based) subsets

# loop over all species (full name)
out <- lapply(unique(df$species_full), function(species_selected){

    # list of all observed data regardless of phenophase
    # get then number of unique years
    ss_full <- subset(df, species_full == species_selected,
                      select = c('year','id'))
    nr_years_full <- length(unique(paste(ss_full$year, ss_full$id)))
    # subset based on phenophase of interest + species
    # get the number of unique years (with observations)
    ss_phen <- subset(df,
                      species_full == species_selected &
                      phenophase == phenophase_selected,
                      select = c('year','id'))
    nr_years_phen <- length(unique(paste(ss_phen$year, ss_phen$id)))

    # return the ratio
    ratio <- nr_years_phen/nr_years_full

    return(data.frame(
      "species_full" = species_selected,
      "phenophase"= phenophase_selected,
      "nr_years_full" = nr_years_full,
      "nr_years_phen" = nr_years_phen,
      "ratio" = ratio))
})

# bind everything row wise
out <- do.call("rbind", out)


print(head(out))
