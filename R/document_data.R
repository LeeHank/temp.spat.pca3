#' @title  Wind speed dataset obtained from Hsu, Chang, and Huang (2012)
#' 
#' @description The data consists of half-yearly east-west component of the wind 
#' speed vector at 289 locations (17x17 regular grid) from November 1992 to February
#' 1993 over latitudes $14^oS$ − $16^oN$ and longitudes $145^oE$ − $175^oE$ in the western Pacific ocean, 
#' @format a 480(time)*289(location) data matrix
"NCEP_U"

#' @title  Coordinate matrix of wind speed dataset
#' 
#' @description Coordinate matrix of wind speed dataset
#' @format a 289(location)*2(coordinates) data matrix
"NCEPlatlon"

#' @title  Simulation data
#' 
#' @description generate 
#' @format a list contains:
#' \describe{
#' \item{y_mat}{simulation data matrix}
#' \item{xsi_mat}{random coefficient matrix}
#' \item{epsilon_mat}{measurement error matrix}
#' \item{location}{location matirx for fitting}
#' \item{phi_mat}{true phi matrix}
#' \item{location_new}{location matrix for prediction}
#' \item{phi_mat_new}{basis functions evaluated at location_new}
#' }
"sim_data"

#' @title  Tuning results for sim_data
#' 
#' @description Tuning results for sim_data
#' @format a list contains:
#' \describe{
#' \item{k_best}{the best K}
#' \item{tau_best}{the best tau}
#' \item{tau_list}{a list contains tuning process of tau in each K}
#' \item{k_process}{the negativet two log-likelihood in each K}
#' \item{best_list}{the detail of best tuning results}
#' }
"cv_result_all"