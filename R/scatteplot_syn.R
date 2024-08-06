#' @title Scatter Plot for Synthetic Data
#'
#' @name plot_syn
#'
#' @param coord_syn A list containing synthetic coordinates.
#' @param n_syn An integer specifying the number of synthetic coordinate sets to plot (must be between 4 and 10).
#'
#' @return This function creates scatter plots for the synthetic data.
#'
#' @examples
#' # Example usage:
#' # dados_sem_coord <- as.data.frame(dados_sem_coord)
#' # coordenadas <- as.data.frame(coordenadas)
#' # n_syn <- 5
#' # coord_restricted_area <- matrix(c(
#' #   4, 4,
#' #   4, 5,
#' #   6, 5,
#' #   6, 4,
#' #   4, 4
#' # ), ncol = 2, byrow = TRUE)
#' # coordenadas_sinteticas <- syncoordinates(dataset = dados_sem_coord,
#' #                                          coord = coordenadas,
#' #                                          grid = 10,
#' #                                          continuous = 4,
#' #                                          list_mcmc = mcmc_dados,
#' #                                          n.syn = n_syn,
#' #                                          coord_restricted_area = coord_restricted_area)
#' # plot_syn(coordenadas_sinteticas, n_syn)
#'
#' @export

## Scatterplot for Synthetic Data
plot_syn <- function(coord_syn, n_syn) {
  if (n_syn < 4 || n_syn > 10) {
    stop("n_syn must be between 4 and 10.")
  }

  coords <- coord_syn[[1]]  # Extracting synthetic coordinates

  nrows <- ceiling(sqrt(n_syn))
  ncols <- ceiling(n_syn / nrows)

  # Creating a list of ggplot graphs
  plots <- lapply(1:n_syn, function(i) {
    coord_set <- coords[, , i]
    if (!is.data.frame(coord_set)) {
      coord_set <- as.data.frame(coord_set)
    }
    names(coord_set) <- c("Longitude", "Latitude")

    ggplot(coord_set, aes(x = Longitude, y = Latitude)) +
      geom_point(alpha = 0.5, color = "black", shape = 16) +
      labs(x = "Longitude", y = "Latitude") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = "white", color = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_blank(),
        axis.title = element_text(size = 12, face = "bold")
      ) +
      coord_cartesian(expand = TRUE)
  })

  grid.arrange(grobs = plots, ncol = ncols, nrow = nrows)
}
