

#' A function to distribute crops between fields on a ALGR output object
#'
#' @param output_obj A ALGR placement output object
#' @param crops_matrix A crop matrix indexing each crop and its percentage in the landscape
#'
#' @return A ALGR output object with crops assigned
#' @export
#'
#' @examples
#' r<-terra::rast(matrix(1, nrow=50, ncol=50))
#' map<-establish_pac(potential_space= r,
#'                                 cell_size=1,
#'                                 includsion_value = 1,
#'                                 mean_field_size = 100,
#'                                 sd_field_size = 25,
#'                                 distribution = "norm",
#'                                 mean_shape_index = 1,
#'                                 sd_shape_index = 0.3,
#'                                 percent = 75)
#'
#'
#'crops_matrix <- data.frame(crop = c("Wheat", "Corn", "Soybean"),
#'                           percentage = c(0.4, 0.3, 0.3),
#'                           index = c(1,2,3))  # Desired percentages
#'
#' outcome<-distrubute_crops(map,crops_matrix)
#'
#' return_by_crop(outcome)


distrubute_crops <- function(output_obj, crops_matrix){

  #calculate total arable area
  total_area <- 0
  fields <- as.data.frame(matrix(ncol = 2, nrow = 0))
  for (i in 1:length(output_obj$field_list)) {
    obj <- output_obj$field_list[[i]]
    total_area <- total_area + length(obj@location[[1]])
    vec <- c(i,length(obj@location[[1]]))
    fields <- rbind(fields,vec)
  }
  colnames(fields) <- c("num","size")


  #calculate desired area of each crop
  crops_matrix$desired_area <- crops_matrix$percentage * total_area

  #sort fields and crops
  fields <- fields[order(-fields$size), ]
  crops_matrix <- crops_matrix[order(-crops_matrix$desired_area),]


  #allocate crops to fields using a greedy allocations approach:
  #  For each field, allocate the crop with the highest remaning desired area,
  #  then remove value from crops_matrix$desired_area and reorder
  for (i in 1:nrow(fields)) {
    crops_matrix <- crops_matrix[order(-crops_matrix$desired_area),]
    output_obj$field_list[[fields[i,1]]]@crop <- crops_matrix[1,3]
    crops_matrix$desired_area[1] <- crops_matrix$desired_area[1] - fields$size[i]

  }


  #allocation_df$field_size <- fields

  return(output_obj)


}







