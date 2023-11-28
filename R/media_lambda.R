
#' @title An auxiliary function to create visual representations of lambda-related values based on MCMC results.
#'
#' @name  media_lambda
#'
#' @description It generates spatially-oriented graphs representing lambda values using color scales.
#'
#' @param results A list containing MCMC results, including elements like "S", "burn", "lambda", and "media.lambda".
#' @param data A list containing dataset information.
#' @param grid A numeric value indicating the grid division for plotting.
#' @param coord A dataframe containing longitude and latitude coordinates.
#' @param standard_scale A logical value indicating whether to standardize color scales in the generated plots.
#' @param save_file A logical value indicating whether to save the generated plots as PDF files.
#'
#' @return This function doesn't explicitly return anything. It generates spatial graphs based on lambda-related data and saves them if specified.
#'
#' @examples
#' # Example usage:
#' # media_lambda(results, data, grid = 10, coord, standard_scale = TRUE, save_file = TRUE)
#'
#' @import fields
#'
#' @export

# Media Lambda function-------------------------------------------------
media_lambda <- function(results, data, grid, coord, standard_scale, save_file) {

  corDegrade <- colorRampPalette(c("light yellow", "red"))

  if (
    !is.list(results) ||
    !is.list(data) ||
    !is.numeric(grid) ||
    !is.data.frame(coord) ||
    !is.logical(save_file) ||
    !is.logical(standard_scale)

  ) {
    stop("Arguments do not have the correct format. Please make sure to provide results and data as lists,
         grid as numeric, coordinates as a dataframe, and save file as T or F.")

  }

  else{
    # Creates a directory named 'graphs' if it doesn't exist already
    graphs <- "graphs"
    if (!file.exists(graphs)) {
      dir.create(graphs, recursive = TRUE)
    }

    num <- length(results)
    lonvec<- data$lonvec
    latvec<- data$latvec

    # Assigns specific elements from 'results' list to variables
    S = results[["S"]]
    burn = results[["burn"]]
    lambda = results[["lambda"]]
    media.lambda = results[["media.lambda"]]

    # Determines the number of combinations from the number of columns in 'media.lambda'
    num_combinacoes <- ncol(media.lambda)

    # Creates a directory named 'graphs/media_lambda' if it doesn't exist already
    graphs_media_lambda <- "graphs/media_lambda"
    if (!file.exists(graphs_media_lambda)) {
      dir.create(graphs_media_lambda, recursive = TRUE)
    }

    # Plot media_lambda graphs in the R environment
    nrows <- ceiling(sqrt(num_combinacoes))
    ncols <- ceiling(num_combinacoes / nrows)

    par(mfrow = c(nrows, ncols), oma = c(1,1,1,1),mar = c(1,1,2,1))

    # Determine the overall range of values for standardizing color scale
    zlim_min <- min(media.lambda)
    zlim_max <- max(media.lambda)

    for (i in 1:num_combinacoes) {

      if (standard_scale == T) {
        # Loops through each combination to create plots
        fields::image.plot(lonvec, latvec, matrix(media.lambda[, i], nrow = grid, ncol = grid),
                   col = corDegrade(10), main = paste("Lambda ", i),
                   xlab = "", ylab = "", yaxt = "n", xaxt = "n", zlim = c(zlim_min, zlim_max))
        points(coord$lon[comb == i], coord$lat[comb == i])
        rect(4, 4, 6, 5, border = "black")
      }

      else{
        fields::image.plot(lonvec, latvec, matrix(media.lambda[, i], nrow = grid, ncol = grid),
                   col = corDegrade(10), main = paste("Lambda ", i),
                   xlab = "", ylab = "",yaxt = "n",xaxt = "n")
        points(coord$lon[comb == i], coord$lat[comb == i])
        rect(4, 4, 6, 5, border = "black")
      }
    }
    par(mfrow = c(1, 1))

    # Plot media_lambda graphs - save the graphs in the 'graphs' folder
    if (save_file == T) {

      for (i in 1:num_combinacoes) {
        if (standard_scale == T) {
          pdf_media_lambda <- file.path(graphs_media_lambda, paste0("plot_media_lambda", i, ".pdf"))
          pdf(file = pdf_media_lambda)
          fields::image.plot(lonvec, latvec, matrix(media.lambda[, i], nrow = grid, ncol = grid),
                     col = corDegrade(10), main = paste("Lambda ", i),
                     xlab = "", ylab = "", yaxt = "n", xaxt = "n", zlim = c(zlim_min, zlim_max))
          points(coord$lon[comb == i], coord$lat[comb == i])
          rect(4, 4, 6, 5, border = "black")
          dev.off()
        }

        else{
          pdf_media_lambda <- file.path(graphs_media_lambda, paste0("plot_media_lambda", i, ".pdf"))
          pdf(file = pdf_media_lambda)
          fields::image.plot(lonvec, latvec, matrix(media.lambda[, i], nrow = grid, ncol = grid),
                     col = corDegrade(10), main = paste("Lambda ", i),
                     xlab = "", ylab = "")
          points(coord$lon[comb == i], coord$lat[comb == i])
          rect(4, 4, 6, 5, border = "black")
          dev.off()
        }
      }
    }
  }
}


