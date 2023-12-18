#' Calculates summary statistics / characteristics at species-level
#' including number of individuals, start/end year, number of observed events,
#' total nr of observation years across individuals
#'
#' @param data junglerhythms data file
#' @param species_name list of species
#' @param pheno only one phenophase
#' @export
#' @return dataframe

overview_stats <- function(
  data = data,
  species_name = "Afzelia bipindensis",
  pheno = "leaf_turnover"){

  stats_output <- data.frame()

  for (j in 1:length(species_name)){
    data_subset <- data %>%
      filter(species_full %in% species_name[j]) %>%
      filter(phenophase %in% pheno)

    # remove all rows with NA in value -> years without observations
    data_subset <- data_subset[!(is.na(data_subset$value)),]
    # count number of observation years with each year per ID
    nr_years_full <- length(unique(paste(data_subset$year, data_subset$id)))
    nr_indiv <- length(unique(data_subset$id))
    start_year <- min(data_subset$year)
    end_year <- max(data_subset$year)
    # get the number of unique years with the specific phenophase event
    events_years <- data_subset %>%
      filter(value > 0)
    nr_years_phen <- length(unique(paste(events_years$year, events_years$id)))
    # return the ratio
    ratio <- nr_years_phen/nr_years_full

    stats_sp <- data.frame(species_full = species_name[j],
                           phenophase = pheno,
                           nr_indiv = nr_indiv,
                           start_year = start_year,
                           end_year = end_year,
                           site_years = nr_years_full,
                           site_years_with_phenophase = nr_years_phen,
                           ratio_site_years_with_phenophase = ratio)

    stats_output <- rbind(stats_output, stats_sp)
  }
  return(stats_output)
}
