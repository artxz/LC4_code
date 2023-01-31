# code for analyze LC4 neurons and downstream target neurons

# There are 3 .R files, run "LC4_startup.R" first. 
# Then "LC4_proc.R" and "LC4_target_proc.R" can be run independently.
# To search for code for all plots, try "PLOT", 
# for a particular plot, try eg. "Fig.3g" or "ED Fig.5b"


# library and cleanup ---------------------------------------------------------------------------------------------

library(tidyverse)
library(natverse)
library(alphashape3d) #alpha shape
library(alphahull) #2D alpha hull
library(RColorBrewer)
library(sf) #compute line-polygon intersections

# library(rjson)
# library(cowplot)

# clean everything up.  
rm(list=ls())


# load Buchner '71 eye map from Andrew Straw----------------------------------------------------------------------------------------------

Npt <- 699
buchner <- read.csv("data/buchner71_tp.csv", header = FALSE)
buchner <- buchner[1:Npt,]
buchner <- buchner / pi * 180 #to degree

# # plot
# dev.new()
# plot(buchner)

# azimuth and elevation range
range(buchner[buchner[,2] > 70 & buchner[,2] < 110, 1]) 
range(buchner[,2])

# set range
buchner_phi <- c(-10, 160) # for +/ 20
buchner_theta <- c(0, 180)


# define some helper functions -------------------------------------------------------------------------------------------------------

#' make polygons from an ashape
#'
#' Intended to be used with the `$edges` dataframe returned by `alphahull::ashape`. 
#' Correct indices order for easy plotting with `polygon`. If there are multiple
#' connected groups, then return multiple polygons.
#' 
#' @param xy intended to be the edge dataframe of an ashape
#'
#' @return a list of polygons' vertices. Eg. to get the first polygon, `mkpoly(ashape$edges)[[1]][,3:4]`
#' @export
#'
#' @examples
mkpoly <- function(xy) {
  xy <- as.data.frame(xy)[,1:6]
  xy_2 <- xy
  
  # remove singlet
  ind_edge <- c(as.matrix(xy[, 1:2]))
  ind_u <- unique(ind_edge)
  ind_s <- ind_u[sapply(ind_u, function(x) sum(ind_edge %in% x)) == 1]
  while (length(ind_s) > 0) {
    ind_m <- match(ind_s[1], xy[,1]) #find this index
    if (!is.na(ind_m)) {
      xy <- xy[-ind_m,]
    } else {
      ind_m <- match(ind_s[1], xy[,2]) #maybe in the second column
      xy <- xy[-ind_m,]
    }
    # search singlet again
    ind_edge <- c(as.matrix(xy[, 1:2]))
    ind_u <- unique(ind_edge)
    ind_s <- ind_u[sapply(ind_u, function(x) sum(ind_edge %in% x)) == 1]
  }
  
  xyset <- list() # vertices for a list of polygons
  if (dim(xy)[1] > 2) {
    # xy is of [x0 y0 x1 y1]
    N <- 1
    xyset[[N]] <- xy[1,]
    xy <- xy[-1, ]
    ii <- c()
    while (dim(xy)[1] >= 1) {
      ii[1] <- match(tail(xyset[[N]], 1)[2], xy[, 1])
      ii[2] <- match(tail(xyset[[N]], 1)[2], xy[, 2])
      if (!is.na(ii[1])) {
        xyset[[N]] <- rbind(xyset[[N]], xy[ii[1], ])
        xy <- xy[-ii[1], ]
      } else if (!is.na(ii[2])){
        xytmp <- xy[ii[2], c(2,1,5,6,3,4)]
        colnames(xytmp) <- colnames(xyset[[N]])
        xyset[[N]] <- rbind(xyset[[N]], xytmp)
        xy <- xy[-ii[2], ]
      } else {
        N <- N + 1
        xyset[[N]] <- xy[1, ]
        xy <- xy[-1, ]
      }
    }
  }
  return(xyset)
}



#' Polar cooridnate to Mollweide projection
#' 
#' Based on [https://mathworld.wolfram.com/MollweideProjection.html].
#' `NA` are kept in the returned matrix.
#'
#' @param thetaphi An Nx2 matrix of polar coordinate on a unit sphere (drop radius): `[0<=theta<=180, -180<=phi<=180]`
#'
#' @return An Nx2 matrix in Mollweide coordinate
#' @export
#'
#' @examples
Mollweide <- function(thetaphi){
  N0 <- nrow(thetaphi)
  xy0 <- matrix(ncol = 2, nrow = N0)
  # deal with NA
  ii_NA <- is.na(thetaphi[,1]) | is.na(thetaphi[,2])
  thetaphi <- thetaphi[!ii_NA,]
  
  N <- nrow(thetaphi)
  lambda <- thetaphi[,2]/180*pi #longitude
  phi <- (90 - thetaphi[,1])/180*pi #latitude
  xy <- matrix(ncol = 2, nrow = N)
  for (j in 1:N) {
    theta <- asin(2*phi[j]/pi) #initial guess
    if (abs(abs(theta) - pi/2) < 1e-3) {
      xy[j,] <- c(2*sqrt(2)/pi*lambda[j]*cos(theta), sqrt(2)*sin(theta))
    } else {
      dtheta <- 1
      while (dtheta > 1e-3) {
        theta_new <- theta - (2*theta + sin(2*theta) - pi*sin(phi[j])) / (2 + 2*cos(2*theta))
        dtheta <- abs(theta_new - theta)
        theta <- theta_new
      }
      xy[j,] <- c(2*sqrt(2)/pi*lambda[j]*cos(theta), sqrt(2)*sin(theta))
    }
  }
  xy0[!ii_NA,] <- xy
  return(xy0)
}



#' Compute arc length on a unit sphere
#'
#' @param p1,p2 Two 1x2 vectors, two points on a unit sphere, `[theta,phi]` in radian.
#'
#' @return A number
#' @export
#'
#' @examples
#' arcLength([pi/2, pi/2], [pi/2, pi])
arcLength <- function(p1, p2) {
  p1_xyz <- c(sin(p1[1])*cos(p1[2]), sin(p1[1])*sin(p1[2]), cos(p1[1]))
  p2_xyz <- c(sin(p2[1])*cos(p2[2]), sin(p2[1])*sin(p2[2]), cos(p2[1]))
  c2 <- sum((p1_xyz - p2_xyz)^2)
  arc_ang <- acos((2-c2)/2)
  return(arc_ang*1)
}

# load neurons from saved data ----------------------------------------------------------------------------------------------------

load('data/LC4.rda')
load('data/neu_target.rda')
load('data/TM5.rda')

load('data/conn_target.rda')

load('data/LO_msh.rda')

# coord axes --------------------------------------------------------------

# from FAFB
axis_ori <- c(350e3, 350e3, 250e3)

axis_lat <- c(-10000, 0, 0) #green
axis_dor <- c(0, -10000, 0)
axis_post <- c(0, 0, 10000)


# palette ---------------------------------------------------------------------------------------------------------

n_lvl <- 11
getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
pal_tar <- getPalette(n_lvl - 1)



# # load neurons from CATMAID server ----------------------------------------
# 
# # public server
# # 
# 
# # load LC4 neuron using annotation
# anno_LC4 <- catmaid_query_by_annotation("^putative LC4 neuron$")
# neu_skid <- anno_LC4[,"skid"]
# neu <-  read.neurons.catmaid(neu_skid, .progress='text')
# 
# LC4 <- neu
# 
# # # DEBUG
# # nopen3d()
# # tar <- neu[[2]]
# # plot3d(tar,  col= 'grey', soma=T, WithNodes = F)
# # points3d(xyzmatrix(tar$d[match(tar$tags$startOfDendrite, tar$d$PointNo),]), size = 10, col = 'blue')
# # points3d(xyzmatrix(tar$d[match(tar$tags$startOfLobula, tar$d$PointNo),]), size = 10, col = 'black')
# # points3d(xyzmatrix(tar$d[match(tar$tags$`LO4-proximal`, tar$d$PointNo),]), size = 10, col = 'pink')
# # points3d(xyzmatrix(tar$d[match(tar$tags$'LO4-distal', tar$d$PointNo),]), size = 10, col = 'red')
# 
# 
# # - load target neuron
# target_01 <-  read.neurons.catmaid(skid = 4947529, .progress='text') #RHS GF
# target_02 <-  read.neurons.catmaid(catmaid_skids("annotation:^putative DNp02$"), .progress='text')
# target_11 <-  read.neurons.catmaid(catmaid_skids("annotation:^putative DNp11$"), .progress='text')
# target_04 <-  read.neurons.catmaid(catmaid_skids("annotation:^putative DNp04$"), .progress='text')
# neu_target <- c(target_01, target_02, target_11, target_04)
# 
# # # TM5
# # neu_JSON <- fromJSON(file = "data/Tm5_LC6 mapping.json")
# # neu_skid <- c()
# # for (j in 1:length(neu_JSON)) {
# #   neu_skid[j] <- neu_JSON[[j]]$skeleton_id
# # }
# # TM5 = read.neurons.catmaid(neu_skid, .progress='text')
# 
# # new center
# TM5 = read.neurons.catmaid(c(13622778, 13622662), .progress='text')
# 
# # nopen3d()
# # points3d(xyz_layer, color = "blue", alpha = 0.7, size = 1)
# # 
# # tar <- TM5[[2]]
# # ii <- identify3d(xyzmatrix(tar$d))
# # points3d(xyzmatrix(tar$d[ii,]), size = 20)
# # 
# # catmaid_get_labels(treenodes = tar$d[ii, 'PointNo'])
# # 
# # ii = match(tar$tags$"TM5 LO col", tar$d$PointNo)
# # catmaid_remove_labels(node = tar$d[ii, 'PointNo'], labels = 'TM5 LO col')
# # 
# # catmaid_set_labels(node = tar$d[ii, 'PointNo'], labels = 'TM5 LO col')
# 
# 
# # conn ------------------------------------------------------------------------------------------------------------
# 
# a1 = -0.95; b1 = 0.74; c1 = 2.1; d1 = -180000 #separate LO and glomerulus
# 
# # nopen3d()
# # plot3d(neu, col = 'grey')
# # planes3d(a1, b1, c1, d1)
# # plot3d(target_01, col = 'red', lwd = 2)
# # plot3d(target_02, col = 'red', lwd = 2)
# # plot3d(target_11, col = 'brown', lwd = 2)
# # plot3d(target_04, col= 'blue', lwd = 2)
# 
# 
# conn_target <- list()
# for (j in 1:length(neu_target)) {
#   tb_conn <- matrix(ncol = 7)
#   for (k in 1:length(neu)) {
#     conn_fromto <- catmaid_get_connectors_between(pre_skids = neu_target[[j]]$skid, post_skids = neu[[k]]$skid)
#     if (is.null(conn_fromto)) {
#       fromto_LO <- 0
#       fromto_glu <- 0
#     } else {
#       # mat_conn <- data.matrix(conn_fromto[, c("connector_x", "connector_y", "connector_z")])
#       mat_conn <- data.matrix(conn_fromto[, c("post_node_x", "post_node_y", "post_node_z")])
#       fromto_LO <- sum(mat_conn %*% c(a1, b1, c1) + d1 > 0)
#       fromto_glu <- dim(conn_fromto)[1] - fromto_LO
#     }
#     conn_tofrom <- catmaid_get_connectors_between(pre_skids = neu[[k]]$skid, post_skids = neu_target[[j]]$skid)
#     if (is.null(conn_tofrom)) {
#       tofrom_LO <- 0
#       tofrom_glu <- 0
#     } else {
#       mat_conn <- data.matrix(conn_tofrom[, c("post_node_x", "post_node_y", "post_node_z")])
#       tofrom_LO <- sum(mat_conn %*% c(a1, b1, c1) + d1 > 0)
#       tofrom_glu <- dim(conn_tofrom)[1] - tofrom_LO
#     }
#     Nconn_glu <- fromto_glu + tofrom_glu
#     tb_conn_tmp <- matrix(c(j, k, fromto_LO, fromto_glu, tofrom_LO, tofrom_glu, Nconn_glu), nrow = 1)
#     tb_conn <- rbind(tb_conn, tb_conn_tmp)
#   }
#   conn_target[[j]] <- as.data.frame(tb_conn[-1,])
#   colnames(conn_target[[j]]) <- c("target","LC4","fromto_LO","fromto_glu","tofrom_LO","tofrom_glu", "Nconn_glu")
# }
# 
# conn_tt <- matrix(ncol = length(neu_target), nrow = length(neu_target)) #target to target
# for (j in 1:length(neu_target)) {
#   for (k in 1:length(neu_target)) {
#     conn_fromto <- catmaid_get_connectors_between(pre_skids = neu_target[[j]]$skid, post_skids = neu_target[[k]]$skid)
#     if (!is.null(conn_fromto)) {
#       mat_conn <- data.matrix(conn_fromto[, c("post_node_x", "post_node_y", "post_node_z")])
#       fromto <- dim(conn_fromto)[1]
#       conn_tt[j,k] <- fromto
#     }
#   }
# }
