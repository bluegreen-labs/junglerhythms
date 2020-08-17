#' Calculates event length and timing
#'
#' This sub-routine will take vectors from a data frame
#' for a single individual and phenophase and returns the timing of
#' phenophases and the length of the .
#'
#' df already needs to be a continuous timeline
#' call function a ID level
#'
#' @param year year vector
#' @param week week vector
#' @param value vector with observations
#' @export
#' @return timing of when an event (value) switches between no event (0)
#' and a observation of an event (1)


event_length <- function(
  data = data,
  species_name = "Afzelia bipindensis",
  pheno = "leaf_turnover"){

  el_output <- data.frame()

  for (j in 1:length(species_name)){
    data_subset <- data %>%
      filter(species_full %in% species_name[j]) %>%
      filter(phenophase %in% pheno)

    individuals <- unique(data_subset$id)
    el_sp <- data.frame()
    for (i in 1:length(individuals)){
      data_ind <-  data_subset %>%
        filter(id %in% individuals[i])

      # sort dataframe according to date
      data_ind <- data_ind %>%
        dplyr::arrange(date)

      # NA in timelines messes up eventlength if phenophase at beginning or end of split timeline
      data_ind$value <- ifelse(is.na(data_ind$value),0,data_ind$value)

      # get first differences
      diff_values <- diff(data_ind$value)

      # get matching info
      start <- data_ind[which(diff_values == 1) + 1,]
      end <- data_ind[which(diff_values == -1),]
      start$index <- as.numeric(rownames(start))
      end$index <- as.numeric(rownames(end))

      # the index is sequential, so take the difference
      # of this index, which counts weeks regardless
      # of date line transitions
      if(length(start$index) == length(end$index)){
        phenophase_length <- end$index - start$index + 1
      } else if(length(start$index) - 1 == length(end$index)){ # phenophase started but still ongoing at end of timeline
        # remove last row start
        rownumber <- length(start$index)
        start <- start[-rownumber,]
        phenophase_length <- end$index - start$index + 1
      } else if(length(start$index) == length(end$index) -1){ # phenophase already ongoing at onset of timeline, only end of phenophase registered
        # remove first row end
        end <- end[-1,]
        phenophase_length <- end$index - start$index + 1
      }


      # return data as data frame
      if(is_empty(phenophase_length)){
        el_ind <- data.frame(id = individuals[i],
                             year_start = NA,
                             week_start = NA,
                             year_end = NA,
                             week_end = NA,
                             index_start = NA,
                             index_end = NA,
                             phenophase_length = NA)
      } else {
        el_ind <- data.frame(id = individuals[i],
                             year_start = start$year,
                             week_start = start$week,
                             year_end = end$year,
                             week_end = end$week,
                             index_start = start$index,
                             index_end = end$index,
                             phenophase_length = phenophase_length)
      }

      el_sp <- rbind(el_sp, el_ind)
    }
    el_sp$species_full <- species_name[j]
    el_sp$phenophase <- pheno
    el_output <- rbind(el_output, el_sp)
  }
  return(el_output)
}


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

