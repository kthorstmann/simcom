


# functions ----------------------------------------------------------


# simulate data function ---------------------------------------------

#' Simulate Data for Components Analysis
#'
#' @param n_pop The size of the population to simulate
#' @param b_s The correlation between behavior and the situational variable on level-2, which determines the correlation of their counterparts on level-1.
#' @param b_e The correlation between the stable environment e and the behavior b.
#' @param s_p The correlation between the situational variable and the average behavior.
#' @param b_p The correlation between the average behavior and personality.
#' @param s_e The correlation between the situation and the stable environment
#' @param p_e The correlation between the personality p and the environment e.
#' @param cor_bs_l1 The average correlation of behavior and situations on level-1, within person.
#' @param l1_var The variance of the level-1 variables (the larger, the less accurate will the mean reflect the true score at fewer measurement occasions)
#' @param l1_bs_sample The number of level-1 measurement occaisons that are sampled.
#'
#' @return Returns a data frame with the following variables
#' \itemize{
#'  \item{"parameter 1"}{Stuff}
#'  \item{"parameter 2"}{Stuff}
#'  \item{"parameter 2"}{Stuff}
#'  \item{"parameter 2"}{Stuff}
#' }
#' @export
#'
#' @examples
#' data <- sim_data()
sim_data <- function(n_pop = 1000,
                     b_s = 0.26,
                     b_e = 0.19,
                     s_p = 0.14,
                     b_p = 0.27,
                     s_e = 0.25,
                     p_e = 0.23,
                     cor_bs_l1 = .2,
                     l1_var = 2,
                     l1_bs_sample = 50){
  ## setup
  # empirical values
  # b_s = 0.26,
  # b_e = 0.19,
  # s_p = 0.14,
  # b_p = 0.27,
  # s_e = 0.25,
  # p_e = 0.23,

  # b_s = 0.3
  # b_e = 0.2
  # s_p = 0.2
  # b_p = 0.4
  # s_e = 0.2
  # p_e = 0.3

  Sigma <- matrix(c(
    1, b_s, b_p, b_e,
    b_s, 1, s_p, s_e,
    b_p, s_p, 1, p_e,
    b_e, s_e, p_e, 1
  ), nrow=4, ncol=4)

  # Sigma matrix
  n_var <- ncol(Sigma)
  means_all = 0
  mu <- rep(means_all, n_var)
  rawvars <- MASS::mvrnorm(n=n_pop, mu=mu, Sigma=Sigma)
  colnames(rawvars) <- c("B", "S", "P", "E")
  raw_data <- data.frame(rawvars)

  B <- raw_data[,"B"]
  S <- raw_data[,"S"]

  id <- stringi::stri_rand_strings(length(B), 20)
  raw_data$id <- id

  simulate_l1_bs <- function(mu_S_person, mu_B_person, cor_bs_l1 = .2,
                             l1_bs_population = 1000,
                             l1_bs_sample = 50,
                             sd_within_slope = 0.05){

    # bsr <- rnorm(n=1, mean=cor_bs_l1, sd = sd_within_slope) # level 1 correlation between b and s with ind. differences in slopes
    bsr <- cor_bs_l1
    l1_bs_matrix <- matrix(c(l1_var, bsr,
                             bsr, l1_var), 2, 2, byrow=TRUE)
    l1_bs_raw <- MASS::mvrnorm(n = l1_bs_population,
                               mu = c(mu_S_person, mu_B_person),
                               Sigma = l1_bs_matrix)

    esm_occasions <- sample(1:l1_bs_population, l1_bs_sample, replace = FALSE)
    sample_occasions <- data.frame(l1_bs_raw[esm_occasions,])
    colnames(sample_occasions) <- c("s_state", "b_state")
    sample_occasions
  }

  level1_states <- purrr::map2(.x = S,
                               .y = B,
                               .f = function(x, y) {
                                 simulate_l1_bs(x, y,
                                                cor_bs_l1 = cor_bs_l1,
                                                l1_bs_sample = l1_bs_sample)
                               }
                               )


  # helper function to add the id
  cidd <- function(d, i){
    is <- rep(i, nrow(d))
    d$id <- is
    d
  }

  ## add id-variables on level-1
  level_1_id_states <- purrr::map2(.x = level1_states,
                                   .y = id,
                                   .f = function(x, y) {
                                     cidd(x, y)
                                   }
                                   )


  # bind data frames from list into one data frame
  l1_states_data <- do.call("rbind", level_1_id_states)

  # merge data on level-1 and level-2
  simulated_data <- dplyr::inner_join(raw_data, l1_states_data, by = "id")
  # return data
  return(simulated_data)

}


# components analysis  -----------------------------------------------

#' Run Components Analysis
#'
#' @param data The data frame to analyze
#' @param outcome_state The outcome or dependent variable that is predicted (the average level-1 variable)
#' @param pred_trait The predictor variables on level-2
#' @param pred_state The predictor variable on level-1
#' @param id The id/group variable of the participants
#'
#' @return Returns a list with two elements:
#' \itemize{
##'  \item{1}{The results of the components analysis}
##'  \item{2}{The explained variacen per state}
##' }
#' @export
#'
#' @examples
#' simulated_data <- sim_data()
#' data <- simulated_data
#' data <- cons_comp_analysis(simulated_data)
cons_comp_analysis <- function(data,
                               outcome_state = "b_state",
                               pred_trait = c("P", "E"),
                               pred_state = c("s_state"),
                               id = "id"){

  # max measurement occasions
  out <- split(data, f = data[, id])
  length_out <- purrr::map(out, nrow)
  max_o <- max(unlist(length_out))

  # function to make all data frames have equal rows
  make_nrow <- function(x, max_o) {
    if (nrow(x) == max_o) {
      return(x)
    } else {
      x[(nrow(x)+1):max_o,] <- NA
      return(x)
    }
  }
  # make all data frames in the list have equal length
  level_1_id_states <- purrr::map(out, ~ make_nrow(., max_o))

  # function to generate the consecutive state variables
  sample_all_o <- function(i){ # rsquare
    odf <- purrr::map(level_1_id_states, ~ .[1:i, ])

    ## aggregate within person to get the average state of the person
    aggregate_within <- function(odf){
      data.frame(id = odf[1,id],
                 out_trait = mean(odf[,outcome_state], na.rm = TRUE),
                 pre_trait = mean(odf[,pred_state], na.rm = TRUE))
    }
    # aggregate states
    aggregate_states <- purrr::map(odf, aggregate_within)

    # combine states in data frame
    aggregate_data <- do.call("rbind", aggregate_states)

    # set names of id before merging
    data.table::setnames(aggregate_data, old="id", new = id)

    # merge aggregate state and trait data
    suppressWarnings(full_data <- dplyr::inner_join(data, aggregate_data, by = id))

    # aggregate state data to trait data
    level2_data <- aggregate(full_data[c(pred_trait, "out_trait", "pre_trait")],
                             by=list(full_data[,id]), mean, na.rm=TRUE)

    # set names to initial names that were handed over to the function
    data.table::setnames(x = level2_data, old = names(level2_data),
                         c(id, pred_trait, outcome_state, pred_state))

    # do the component analysis
    apsOut = yhat::aps(level2_data, outcome_state, list(pred_trait[1],
                                                        pred_state,
                                                        pred_trait[2])
                       )


    # return parameters
    commonality <- round(yhat::commonality(apsOut)[,2], 3)
    rsquared <- sum(apsOut$APSMatrix[,2])

    return(
      list(commonality,
           rsquared)
    )
  }


  # run the above function across all state averages
  params <- purrr::map(1:max_o, sample_all_o)

  # extract all first elements of the list
  c <- lapply(params, "[[", 1)
  rsq <- unlist(lapply(params, "[[", 2))
  params_df <- do.call("rbind", c)
  params_df[params_df < 0] <- NA
  return(
    list(params_df,
         rsq)
  )
}


# Plot Stacked Areas -------------------------------------------------

#' Plot Stacked Areas
#'
#' @param data The data frame that needs to be plotted. Should only contain the variables that should appear in the plot.
#' @param poly.border.col The color of the border of the stacked areas.
#' @param poly.fill.col The fill of the stacked areas. The areas are colored from top to bottom (i.e. the highest one first).
#' @param letters The letters that should be placed in the areas. If NULL, the names of the data is used.
#' @param letter.col The colors of the letters that are plotted.
#' @param letters.right. Logical, default is FALSE. If TRUE, the letters are set to the right of the plot.
#' @param ... Further arguments passed on to the plot function plot(...)
#'
#' @return Plots areas that are scaled to 100%.
#' @export
#'
#' @examples
#' simulated_data <- sim_data()
#' comp_list <- cons_comp_analysis(simulated_data)
#' plot_data <- comp_list[[1]]
#' plot_polygons(plot_data, letters.right = TRUE)
plot_polygons <- function(data,
                          poly.border.col = "grey10",
                          poly.fill.col = c("white", "white",
                                            "white", "white", "#F8CBAE",
                                            "#FDF2CA", "#9DC3E6"),
                          letters = c("P,S,E", "S,E", "P,E", "P,S", "E", "S", "P"),
                          letter.col = "black",
                          letters.right = FALSE,
                          ...){

  is.num <-  unlist(purrr::map(data, is.numeric))
  nonNumeric <-  names(data)[!is.num]
  if (any(!is.num)) {
    warning(paste("The following items/columns are non-numeric. \n ",
               nonNumeric, collapse=" "), call.=FALSE)
  }

  xunits <- (nrow(data)) # determine the number of units, because we want to plot from 0 to n.
  rownames(data) <- 1:xunits
  # get percentages:
  total <- rowSums(data, na.rm=TRUE)
  perc.dat <- as.data.frame(apply(data, 2, function(x) x/total*100))
  t.perc.dat <- t(perc.dat)
  # make a vector over which we apply the function that computes the boundaries for the polygons
  len <- 2:ncol(data)
  # apply the function in order to get the borders of the polygons
  res <- rev(lapply(len, function(x) colSums(t.perc.dat[1:x,], na.rm=TRUE)))
  # add the last polygon to list
  res[[length(res)+1]] <- t.perc.dat[1,]
  # compute the lower border ofthe last polygon, i. e. 0, repeated as many times as data has rows
  lower <- rep(0, xunits) # +1

  # plot the polygons, define the function
  plot.poly <- function(x){
    xcoords <-  c(xunits:1, 1:xunits)
    ycoords <-  c(lower, res[[x]])
    polygon(xcoords, ycoords, col = poly.fill.col[x], border = poly.border.col)
  }

  ## plot an empty area
  plot(NA, ylim=c(1, 110), xlim = c(0, xunits), type="n",
       yaxt="n", xaxt="n", bty="n", ylab="", xlab="", ...)
  axis(1)

  # plot the polygons
  invisible(lapply(1:length(res), function(x) plot.poly(x) ))

  if (ncol(data) != length(poly.fill.col)) {
    warning("Only ", length(poly.fill.col), " colors were specified, but ", ncol(data), " areas were plotted. Colors were recycled. Please check.")
  }

  # add letters, default is to names of data
  if (is.null(letters) ) { letters <- rev(colnames(data)) }

  # add to the result list with the borders another lower border, i. e. the zeros with names
  names(lower) <- 0:(length(lower)-1)
  res[[length(res)+1]] <- lower

  # apply the function over the borders, computing the borders and plotting the polygon
  invisible(
    lapply(1:(length(res)-1), function(x){
      y = x + 1
      # make two empty lists
      list_low <- list()
      list_up <- list()
      # reduce the plotting area left and right
      list_low[[x]] <- res[[x]][ceiling(length(res[[x]])*.08):floor(length(res[[x]])*.90)]
      list_up[[x]] <- res[[y]][ceiling(length(res[[y]])*.08):floor(length(res[[y]])*.90)]

      ## add here the argument that letters can be plotted to the right
      if ( nrow(data) >= 12 & !letters.right) {
        # get the position of the letters, which is where the maximum of the difference between the borders is. if there is more than one maximum, take only the first
        posx = names(which(list_low[[x]]-list_up[[x]] == max(list_low[[x]]-list_up[[x]])))[1]
        # compute the mean difference, needs to be in the tighter/reduced lists
        posy = mean(c(list_low[[x]][posx], list_up[[x]][posx]))
        text(as.numeric(posx), posy, letters[x], col = letter.col)
      } else
      {
        posx <- (nrow(data)-1)
        # compute the mean, not in the tight lists (but in 'res', to plot the letter to the right hand side)
        posy = mean(c(res[[x]][posx], res[[y]][posx]))
        text(posx-posx*0.03, posy, letters[x], col = letter.col)
      }
    }))
}





# plot components with explained variance ----------------------------

#' Plot results of components analysis
#'
#' @param params_df The results of the function, which is a list with to entries.
#' @param ... Further arguments passed on to plot_polygons
#'
#' @return
#' Returns a plot with the components
#' @export
#'
#' @examples
#' data <- sim_data()
#' comp_list <- cons_comp_analysis(simulated_data)
#' plot_com(components = comp_list[[1]], rsquared = comp_list[[2]])
plot_com <- function(components, rsquared, ...){
  params <- components
  rsq <- rsquared

  # par (mar = c(5.1, 4.1, 4.1, 0))
  plot_polygons(params, ...)

  # plot addional explained variance in the plot
  rsq_plot <- rsq*10+100
  lines(rsq_plot, col = "black", lwd = 2, lty = 3)
  min <- round(rsq[1], 2)
  max <- round(rsq[length(rsq_plot)], 2)
  text(x = 1, y = 10 + 100, labels=min)
  text(x = length(rsq_plot), y = 10 + 100, labels=max)

  axis(2, at = c(0, 20, 40, 60, 80, 100), labels = c(0, 20, 40, 60, 80, 100))
  # add explained variance
  axis(4, at=c(100, 105, 110), labels = c(0, 50, 100))
}


# replicate sim  -----------------------------------------------------

#' Repeat Simulation
#'
#' @param REPS The number of simulations to be made of
#' @param n_pop The size of the population to simulate
#' @param b_s The correlation between behavior and the situational variable on level-2, which determines the correlation of their counterparts on level-1.
#' @param b_e The correlation between the stable environment e and the behavior b.
#' @param s_p The correlation between the situational variable and the average behavior.
#' @param b_p The correlation between the average behavior and personality.
#' @param s_e The correlation between the situation and the stable environment
#' @param p_e The correlation between the personality p and the environment e.
#' @param cor_bs_l1 The average correlation of behavior and situations on level-1, within person.
#' @param l1_var The variance of the level-1 variables (the larger, the less accurate will the mean reflect the true score at fewer measurement occasions)
#' @param l1_bs_sample The number of level-1 measurement occaisons that are sampled.
#'
#' @return Returns a list with two objects: The first is a data frame for the component analyses, the second is a vector with the explained variances.
#' @export
#'
#' @examples
#' replicate_sim(10)
replicate_sim <- function(REPS = 10,
                          n_pop = 1000,
                          b_s = 0.26,
                          b_e = 0.19,
                          s_p = 0.14,
                          b_p = 0.27,
                          s_e = 0.25,
                          p_e = 0.23,
                          cor_bs_l1 = .2,
                          l1_var = 2,
                          l1_bs_sample = 50){

  # function to replicate
  repl_f <- function(...){
    data <- sim_data(...)
    parameters <- cons_comp_analysis(data)
    parameters
  }
  # pass all parameters from above

  parameters_list <- replicate(REPS, repl_f(n_pop = n_pop, b_s = b_s, b_e = b_e,
                                            s_p = s_p, b_p = b_p, s_e = s_e,
                                            p_e = p_e, cor_bs_l1 = cor_bs_l1,
                                            l1_var = l1_var, l1_bs_sample = l1_bs_sample)
  )

  parameters_pos <- seq(1, REPS*2, by = 2)
  rsquared_pos <- seq(2, REPS*2, by = 2)

  para_list <- purrr::map(parameters_pos, ~ parameters_list[[.]])
  rsq_list <- purrr::map(rsquared_pos, ~ parameters_list[[.]])

  para_mean <-  plyr::aaply(plyr::laply(para_list, as.matrix), c(2, 3), mean, na.rm =TRUE)
  rsq_mean <- colMeans(do.call("rbind", rsq_list))
  return <- list(para_mean,
                 rsq_mean)
  return(return)
}

