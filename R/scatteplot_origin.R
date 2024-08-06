#' @title Scatter Plot for Original Data
#'
#' @name plot_origin
#'
#' @param coord A data frame containing the original coordinates with columns named "lon" and "lat".
#'
#' @return This function creates a scatter plot for the original data.
#'
#' @examples
#' # Example usage:
#' # coordenadas <- data.frame(
#' #   lon = c(0.2288, 3.6682, 5.2056, 5.2007, 6.6556, 9.4929),
#' #   lat = c(0.1211, 2.8985, 4.8331, 4.7690, 6.6016, 9.6699)
#' # )
#' # plot_origin(coordenadas)
#'
#' @export

## Scatterplot for Original Data
plot_origin <- function(coord) {
  ggplot(coord, aes(x = lon, y = lat)) +
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
}
