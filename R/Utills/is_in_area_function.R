#function to verify if the points are inside the interval
#passed by the user
is_in_area <- function(x, y, xmin, xmax, ymin, ymax){
  is_x_valid = (x >= xmin & x <= xmax)
  is_y_valid = (y >= ymin & y <= ymax)
  return((is_x_valid & is_y_valid))
}