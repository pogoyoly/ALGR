

#' A function for generating artificial roads
#'
#' @param road_length A number
#' @param slope A lope raster
#'
#' @return A raster
#' @export
#' @import raster sp
#' @importFrom gdistance transition geoCorrection shortestPath
#'
#' @examples
#'test<-generate_pn(200,200,1,2,3,0.005,FALSE, "land_percentage", percentage = 50)
#'test2<-generate_roads(2000,test)
#'raster::plot(test2)
#'
generate_roads <- function(road_length,slope){

  circularState <- function(state, value) {
    nextState <- (state + value - 1) %% 4 + 1
    return(nextState)
  }

  #generate a temp raster that has origin and extent that starts from zero || check with antonia if we can do this for everything from start
  slope_raster <- slope
  slope_raster <- raster::raster(slope_raster)
  raster::origin(slope_raster) <- 0
  bb <- raster::extent(0, nrow(slope_raster), 0, ncol(slope_raster))
  raster::extent(slope_raster) <- bb

  slope_raster[is.na(slope_raster[]) == TRUE] <- 255
  slope_raster[slope_raster == 0] <- 255
  slope_raster[slope_raster < 0] <- 255
  nrows <- nrow(slope_raster)
  ncols <- ncol(slope_raster)

  # Set edge values to 100
  slope_raster[1, ] <- 100        # Top edge
  slope_raster[nrows, ] <- 100    # Bottom edge
  slope_raster[, 1] <- 100        # Left edge
  slope_raster[, ncols] <- 100    # Right edge



  #plot(slope_raster)
  rows <- nrow(slope_raster) - 1
  cols <- ncol(slope_raster) - 1

  realized <- 0

  #create a dummay spatial line and then empty it
  sl <- sp::SpatialLines(LinesList = list(sp::Lines(sp::Line(matrix(0, ncol = 2)), ID = NA)))
  sl <- sl[0]


  while (realized < road_length) {
    start <- sample(1:4, 1)
    end <- sample(1:3, 1)

    #use the circle state function to make sure start and end are not on the same panel
    end <- circularState(start,end)

    #choose start on one of the panels (all next 4 ifelse arguments)
    if ( start == 1) {
      xi <- sample(2:cols, 1)
      yi <- 2
      while (is.na(slope_raster[xi,yi]) == TRUE | slope_raster[xi,yi] == 255) {
        yi <- yi + 3
        if (yi > rows | yi < 0) {
          stop("Use different method for roads")
        }

      }
    } else if ( start == 2) {
      xi <- sample(2:cols, 1)
      yi <- rows - 3
      while (is.na(slope_raster[xi,yi]) == TRUE | slope_raster[xi,yi] == 255) {
        yi <- yi - 3
        if (yi > rows | yi < 0) {
          stop("Use different method for roads")
        }

      }
    } else if ( start == 3) {
      xi <- 3
      yi <- sample(2:rows, 1)
      while (is.na(slope_raster[xi,yi]) == TRUE | slope_raster[xi,yi] == 255) {
        xi <- xi + 3
        if (xi > cols | xi < 0) {
          stop("Use different method for roads")
        }

      }
    } else {
      xi <- cols - 3
      yi <- sample(2:rows, 1)
      while (is.na(slope_raster[xi,yi]) == TRUE | slope_raster[xi,yi] == 255) {
        xi <- xi - 3
        if (xi > cols | xi < 0) {
          stop("Use different method for roads")
        }

      }

    }
    A <- c(xi,yi)

    if ( end == 1) {
      xii <- sample(2:cols, 1)
      yii <- 3
      while (is.na(slope_raster[xii,yii]) == TRUE | slope_raster[xii,yii] == 255) {
        yii <- yii + 3
        if (yi > rows | yi < 0) {
          stop("Use different method for roads")
        }

      }

    } else if ( end  == 2) {
      xii <- sample(2:cols, 1)
      yii <- rows - 3
      while (is.na(slope_raster[xii,yii]) == TRUE | slope_raster[xii,yii] == 255) {
        yii <- yii - 3
        if (yi > rows | yi < 0) {
          stop("Use different method for roads")
        }

      }

    } else if ( end == 3) {
      xii <- 2
      yii <- sample(2:rows, 1)
      while (is.na(slope_raster[xii,yii]) == TRUE | slope_raster[xii,yii] == 255) {
        xii <- xii + 1
        if (xi > cols | xi < 0) {
          stop("Use different method for roads")
        }

      }

    } else {
      xii <- cols - 1
      yii <- sample(2:rows, 1)
      while (is.na(slope_raster[xii,yii]) == TRUE | slope_raster[xii,yii] == 255) {
        xii <- xii - 1
        if (xi > cols | xi < 0) {
          stop("Use different method for roads")
        }

      }

    }

    B <- c(xii,yii)

    if (any(is.na(slope_raster[])) || any(is.infinite(1/slope_raster[]))) {
      stop("The raster contains NA or infinite values.")
    }

    trans <- gdistance::transition(slope_raster, function(x) 1/mean(x), 4)
    trans_geo <- gdistance::geoCorrection(trans, type = "c")
    path <- gdistance::shortestPath(trans_geo, A, B, output = "SpatialLines")

    realized <- realized + sum(sp::SpatialLinesLengths(path, longlat = FALSE))

    sl <- rbind(sl, path)
  }

  r <- raster::raster(sl, resolution = c(1, 1))
  r_path <- raster::rasterize(sl, r, field = 1)

  #bb <- extent(slope)
  #res(r_path) <- res(slope)
  #extent(r_path) <- bb

  return(r_path)


}
