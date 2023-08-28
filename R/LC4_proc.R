
# make a layer grid using LC4 dendrite (ED Fig.5a) -----------------------------
# Fit a surface to LC4 dendrites,  
# and generate grid on the surface for later analysis

# - get all nodes downstream of a dendrite tag
LC4_dend <- list()
dend_v <- matrix(ncol = 3)
for (j in 1:length(LC4)) {
  tar <- LC4[[j]]
  targ <- as.ngraph(tar)
  # use a depth first search
  distal_points <- igraph::graph.dfs(targ, root= match(tar$tags$startOfDendrite, tar$d$PointNo), unreachable=FALSE, neimode='out')$order
  distal_tree <- subset(tar, distal_points)
  
  LC4_dend[[j]] <- xyzmatrix(distal_tree$d)
  dend_v <- rbind(dend_v, LC4_dend[[j]])
}
dend_v <- dend_v[-1,]

# - make alpha mesh
msh.a <- ashape3d(dend_v, alpha = 20000, pert = T) # 20000 look ok
msh <- as.mesh3d(msh.a)
LC4_lo <- nlapply(LC4, subset, function(x) pointsinside(x, msh))

# - fit a 2nd order polynomial surface
# this takes some time!
polyfitorder <- 2
gridMargin <- 30000 #add margin on the edge
X <- dend_v[, 1];  Y <- dend_v[, 2]; Z <- dend_v[, 3]
dx2 <- 500
dy2 <- 500
xx <- seq(range(X)[1] - gridMargin, range(X)[2] + gridMargin, by = dx2)
yy <- seq(range(Y)[1] - gridMargin, range(Y)[2] + gridMargin, by = dy2)
fitlm <- lm(Z ~ poly(X, Y, degree = polyfitorder, raw = TRUE)) #linear model fit
xygrid <- expand.grid(xx, yy)
# xygrid <- with(df_xyz, expand.grid(xx, yy))
xygrid <- setNames(data.frame(xygrid), c('X', 'Y'));
valfit <- predict(fitlm, xygrid) #generate values from the fit
xyz_lm <- cbind(xygrid, valfit)

# keep points inside the lobula mesh
ii <- pointsinside(xyz_lm, LO_msh,rval='distance') > -3000 #for plotting
xyz_layer_plot <- xyz_lm[ii,] # pts

ii <- pointsinside(xyz_lm, LO_msh,rval='distance') > -3000 #tighter for computing
xyz_layer <- xyz_lm[ii,] # pts


# PLOT, ED Fig.5a
nopen3d()
par3d('windowRect' = c(10,10,1710,1710))
points3d(xyz_layer, color = "#fdae61", alpha = 0.9, size = 2)
# points3d(dend_v, size = 2, col = 'gray90')
plot3d(LC4[-landmk], col='gray90', soma = T, lwd=1, WithNodes = F)
plot3d(LC4[[landmk[1]]],  col= "#d7191c", lwd = 5, soma=T, WithNodes = F)
plot3d(LC4[[landmk[2]]],  col= "#2c7bb6", lwd = 5, soma=T, WithNodes = F)
plot3d(nlapply(TM5[1], subset, function(x) pointsinside(x, msh,rval='distance')>-1e4), col = 'gold4', lwd = 4) #TM5
plot3d(nlapply(TM5[2], subset, function(x) pointsinside(x, msh,rval='distance')>-1.2e4), col = 'gold2', lwd = 4) 
arrow3d(axis_ori, axis_ori + axis_lat, theta = pi/6,n = 4, col="green", type = "rotation")
arrow3d(axis_ori, axis_ori + axis_dor, theta = pi/6,n = 4, col="magenta", type = "rotation")
arrow3d(axis_ori, axis_ori + axis_post, theta = pi/6,n = 4, col="blue", type = "rotation")
rgl.viewpoint(fov=0,zoom=0.8,userMatrix=rotationMatrix(170/180*pi,1,0,0) %*% rotationMatrix(30/180*pi,0,0,1) %*% rotationMatrix(-65/180*pi,0,1,0))
# plot3d(neu_target, lwd=5)

# # save
# rgl.snapshot(filename = "LC4_3d.png",fmt = "png")


# project LC4 dendrite and center-of-mass onto the layer grid (ED Fig.5b) -------------
# Also project the layer grid onto a 2D plane and use TM5 to align it.

# - project onto the layer  ****** THIS IS SLOW ***** 
row.names(xyz_layer) <-  seq(1, dim(xyz_layer)[1])
ind_pj <- list() #index of projected grid points
xyz_com <- list() # center-of-mass
xyz_pj_com <- list() # xyz of com projecting on grid
for (j in 1:length(LC4)) {
  xyz_com[[j]] <- colMeans(LC4_dend[[j]]) #center of mass
  # project dendrite end points and com to the fitted grid by shortest distance
  xyz_dend_pj <- rbind(xyz_com[[j]], LC4_dend[[j]]) #append com at the beginning
  Nsigma <- 5 #exclude > 5 sigma
  com_dist <- as.matrix(dist(xyz_dend_pj))
  # colnames(com_dist[,1]) <- 'D'
  thrhd_dist <- sd(com_dist[,1])*Nsigma
  xyz_dend_pj <- cbind(xyz_dend_pj, com_dist[,1]) %>%
    as_tibble() %>%
    filter(V4 < thrhd_dist)%>%
    select(X,Y,Z) %>%
    data.matrix()
  
  ind_min <- c()
  ind_min <-apply(xyz_dend_pj, 1, function(pt) {which.min(rowSums(sweep(xyz_layer, 2, pt) ^ 2))}) # index in xyz_layer with min distance
  # ind_min <- apply(xyz_dend, 1, function(pt){which.min(rowSums(sweep(xyz_layer,2,pt)^2))}) # index in xyz_layer with min distance
  ind_com <- ind_min[1] #index of grid point that's closest to com
  ind_min <- unique(ind_min[-1]) #index of grid points
  xyz_pj_com[[j]] <- xyz_layer[ind_com,]
  ind_pj[[j]] <- row.names(xyz_layer[ind_min,])
}

# project onto a plane def by pc3 of xyz_layer
layer_pca <- prcomp(xyz_layer)
xyz_layer_rot <- t(layer_pca$rotation) %*% t(sweep(xyz_layer, 2, layer_pca$center,))
xyz_layer_rot <- t(xyz_layer_rot)
xy_layer_rot <- xyz_layer_rot[,c(1,2)]


# TM5 positions
TM5_c_xyz <- xyzmatrix(TM5[[1]]$d[match(TM5[[1]]$tags$"TM5 LO col", TM5[[1]]$d$PointNo),])
TM5_u_xyz <- xyzmatrix(TM5[[2]]$d[match(TM5[[2]]$tags$"TM5 LO col", TM5[[2]]$d$PointNo),])

# grid points near 2 TM5s
grid_c <- which.min(rowSums((sweep(xyz_layer, 2, c(TM5_c_xyz), "-"))^2))
grid_u <- which.min(rowSums((sweep(xyz_layer, 2, c(TM5_u_xyz), "-"))^2))


# - map equator and coord system 
center_new <- xyz_layer_rot[grid_c,1:2]
x_med_new <- as.numeric(center_new[1])
y_eq_new <- as.numeric(center_new[2])
angR_new <- acos((x_med_new - xyz_layer_rot[grid_u,1])/sqrt(sum((xyz_layer_rot[grid_c,1:2] - xyz_layer_rot[grid_u,1:2])^2)))

# 2d projection with dendrites, centered by 2 TM5
ang_2 <- pi/2 - angR_new
rot_2 <- matrix(c(cos(ang_2), sin(ang_2), -sin(ang_2), cos(ang_2)), ncol = 2)
xy_layer_align <- sweep(xy_layer_rot, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
xy_layer_align <- t(rot_2 %*% t(xy_layer_align))

# reconstruct (x,y) projection - just keep x,y, then align
xy_pj <- list()
xy_pj_com <- list()
for (j in 1:length(ind_pj)){
  xy_tmp <- list()
  xy_ashape <- list()
  ii <- sort(as.integer(ind_pj[[j]]))
  xy_tmp <- xy_layer_align[ii, ]
  xy_ashape <- ashape(xy_tmp, alpha = 6000)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_pj[[j]] <- list(xy=xy_tmp, ashape=xy_ashape, edge=xy_edge)
  
  xy_com_tmp <- unlist(xyz_pj_com[[j]]-layer_pca$center) %*% layer_pca$rotation
  xy_com_tmp <- xy_com_tmp[c(1,2)] - c(x_med_new, y_eq_new)
  xy_com_tmp <- t(rot_2 %*% xy_com_tmp)
  xy_pj_com[[j]] <- xy_com_tmp
}


# # use alpha-hull for the grid projection, -> obs. 
# xy_ashape_grid <- ashape(xy_layer_align + matrix(runif(dim(xy_layer_align)[1]*2, 1e-9, 2e-9), ncol = 2), alpha = 6000)
# xy_grid_ahull <- mkpoly(xy_ashape_grid$edges)[[1]][,3:4]
# xy_edge_grid <- xy_grid_ahull # hull edge points


# # check results
# nopen3d()
# points3d(xyz_layer, size = 5)
# points3d(xyz_layer[grid_c,], size = 30, col = 'gold4')
# points3d(xyz_layer[grid_u,], size = 30, col = 'gold2')
# points3d(xyz_pj_com[[landmk[1]]], size = 30, col = 'gold')
# points3d(xyz_pj_com[[landmk[2]]], size = 30, col = 'green')
# 
# dev.new()
# plot(xy_layer_align[seq(1,dim(xy_layer_align)[1], by = 7),])
# # points(center_new[1],center_new[2], col = 'red', cex = 3, pch = 16)
# points(xy_layer_align[grid_c,1],xy_layer_align[grid_c,2], col = 'red', cex = 4, pch = 19)
# points(xy_layer_align[grid_u,1],xy_layer_align[grid_u,2], col = 'brown', cex = 4, pch = 19)
# points(xy_pj_com[[landmk[1]]][1,1], xy_pj_com[[landmk[1]]][1,2], col = 'gold', cex = 4, pch = 19)
# points(xy_pj_com[[landmk[2]]][1,1], xy_pj_com[[landmk[2]]][1,2], col = 'green', cex = 4, pch = 19)
# polygon(xy_edge_grid)



# PLOT ED Fig.5b
dend_v_align <- sweep(dend_v, 2, layer_pca$center) %*% layer_pca$rotation
dend_v_align <- sweep(dend_v_align[,c(1,2)], MARGIN = 2, STATS = c(x_med_new, y_eq_new))
dend_v_align <- t(rot_2 %*% t(dend_v_align))

windows(record = F, width = 8, height = 8)
# pdf(file = "LC4_2d.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plot(xy_layer_align, col="#fdae61", cex = 1, pch = ".",
     ylim = rev(range(xy_layer_align[,2])),
     xlim = (range(xy_layer_align[,1])),
     asp = 1, main = "2D proj aligned")
points(dend_v_align, col = "grey80", pch = ".", cex = 1)
twig <- LC4_lo[[landmk[1]]]
pp <- as.matrix(sweep(twig$d[,c("X","Y","Z")], 2, layer_pca$center)) %*% layer_pca$rotation
pp <- sweep(pp[,1:2], 2, STATS = c(x_med_new, y_eq_new))
twig$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new))
plot(twig,  col= "#d7191c", lwd = 2, soma=T, WithNodes = F, add = T)
twig <- LC4_lo[[landmk[2]]]
pp <- as.matrix(sweep(twig$d[,c("X","Y","Z")], 2, layer_pca$center)) %*% layer_pca$rotation
pp <- sweep(pp[,1:2], 2, STATS = c(x_med_new, y_eq_new))
twig$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new))
plot(twig,  col= "#2c7bb6", lwd = 2, soma=T, WithNodes = F, add = T)
# points(xy_pj_com[[landmk[1]]], col = 'gold', cex = 4, pch = 19)
# points(xy_pj_com[[landmk[2]]], col = 'green', cex = 4, pch = 19)
lines(rbind(-xy_layer_align[grid_u,]*1.8, xy_layer_align[grid_u,]*1.8), lwd = 3, col = 'cyan')
points(matrix(c(x_med_new,y_eq_new), ncol =2), pch = 18, col = 'gold4', cex = 2)
points(matrix(xy_layer_align[grid_u,], ncol=2), pch = 18, col = 'gold2', cex = 2)
lines(c(50000,60000), c(-70000, -70000), col = 'black', lwd = 3)
points(matrix(unlist(xy_pj_com), ncol = 2, byrow = T), pch = 20, col = "blue", cex = 1.5)
text(x = 55000, -65000, labels = "10 µm")

# dev.off()


#  jet color plot LC4 synapse counts -----------------------------------------------------

getPalette <- colorRampPalette(brewer.pal(9, "Spectral"))

n <- 3 #choose a target, 1:4

ns <- conn_target[[n]]$Nconn_glu
LC4col <- rev(getPalette(max(ns)+1))

nopen3d()
plot3d(LC4, col = LC4col[ns+1])
plot3d(neu_target[[n]], col = 'grey', lwd = 1, soma = T)


# wrap the lobula layer onto a hemisphere -------------------------------------------------------------------------
# Model the field of view as a hemisphere and assume lobula maps onto it in a uniform fashion
# Imagine the layer grid as an elastic sheet, we wrap it onto a hemisphere.
# For each point on the layer grid, compute a radial distance (from the
# pole of the hemisphere to its boundary/equator) proportional to the distance from
# the center of the layer grid, through this point, to the grid boundary. 

# Each LC4's 2D projection has a center-of-mass and a polygon representing its dendrite extent.
# First get edge points on the polygon.
xy_poly <- list()
xy_bdpt <- matrix(ncol = 2, nrow = 0) #collect all LC4 polygon edge points to make a LO boundary
for (j in 1:length(ind_pj)) {
  ls_poly <- mkpoly(xy_pj[[j]]$ashape$edges)
  xy_poly[[j]] <- ls_poly[[1]][,3:4]
  xy_bdpt <- rbind(xy_bdpt, xy_pj[[j]]$edge[,3:4])
}

# use alpha-hull for the grid projection, -> new
xy_ashape_grid <- ashape(unique(xy_bdpt), alpha = 6000)
xy_grid_ahull <- mkpoly(xy_ashape_grid$edges)[[1]][,3:4]
xy_edge_grid <- xy_grid_ahull # hull edge points

grid_bdpt <- xy_edge_grid
grid_bdpt <- rbind(grid_bdpt, grid_bdpt[1,])

ymax <- max(xy_edge_grid[,2])
ymin <- min(xy_edge_grid[,2])


# LC4 polygons and centers, in eye coordinates
poly_st <- st_polygon(list(data.matrix(grid_bdpt)))
xy_ori <- c(0,0)
R <- (ymax-ymin)*2
xy_com <- list() # com of the projected LC4 dendrite in eye coord
for (j in 1:length(xy_poly)) {
  pj_com <- xy_pj_com[[j]] %>% as.matrix()
  colnames(pj_com) <- c("x1", "y1")
  xy_poly[[j]] <- rbind(pj_com, xy_poly[[j]])
  xy_poly[[j]] %<>% 
    as_tibble() %>%
    mutate(phiC = acos((x1-x_med_new)/sqrt((x1-x_med_new)^2+(y1-y_eq_new)^2))*(-1)^(y1<y_eq_new) + 2*pi*(y1<y_eq_new)) %>%  #angle
    mutate(thetaC = NA) %>% #radius
    transmute(x1, y1, phiC, thetaC) %>%
    data.matrix()
  for (k in 1:dim(xy_poly[[j]])[1]) {
    alpha <- xy_poly[[j]][k,'phiC']
    line_st <- st_linestring(rbind(xy_ori, c(xy_ori[1] + R*cos(alpha), xy_ori[2] + R*sin(alpha))))
    int <- data.matrix(st_intersection(line_st, poly_st))[2,]
    xy_poly[[j]][k,'thetaC'] <- pi/2 * dist(rbind(xy_poly[[j]][k,1:2], xy_ori)) / dist(rbind(int,xy_ori))
  }
  # now turn inside out and a 90-rotation about x-axis to have the front edge on the x-z plane, 
  # angle wrt looking from behiind the eye
  xy_poly[[j]] %<>%
    as_tibble() %>%
    # mutate(thetaC = pi - thetaC) %>% # turn inside out
    mutate(phiC = phiC ) %>% # align front edge to x-z plane
    mutate(x = sin(thetaC)*cos(phiC), y = sin(thetaC)*sin(phiC), z = cos(thetaC)) %>% #(x,y,z) coord
    mutate(xr = x, yr = z, zr = -y) %>% # +90 rotation about x-axis
    mutate(xrr = xr, yrr = yr, zrr = zr) %>% 
    mutate(theta = acos(zrr/sqrt(xrr^2+yrr^2+zrr^2)), phi = acos(xrr/sqrt(xrr^2+yrr^2))) %>%
    # mutate(theta_deg = 180 - theta/pi*180, phi_deg = phi/pi*180) %>% #the convention for x-y plot below
    mutate(theta_deg = theta/pi*180/180*diff(buchner_theta)+buchner_theta[1], 
           phi_deg = phi/pi*180/180*diff(buchner_phi)+buchner_phi[1]) %>% # longitude use buchner_phi
    # select(x1, y1, theta_deg, phi_deg) %>%
    data.matrix()
  
  xyMollweide <- Mollweide(xy_poly[[j]][,c("theta_deg", "phi_deg")])
  colnames(xyMollweide) <- c('xM','yM')
  xy_poly[[j]] <- cbind(xy_poly[[j]], xyMollweide)
  
  xy_com[[j]] <- xy_poly[[j]][1,]
  xy_poly[[j]] <- xy_poly[[j]][-1,]
}

# boundary in eye coord
grid_bdpt_tp <- grid_bdpt %>%
  as_tibble() %>%
  mutate(phiC = acos((x1-x_med_new)/sqrt((x1-x_med_new)^2+(y1-y_eq_new)^2))*(-1)^(y1<y_eq_new) + 2*pi*(y1<y_eq_new)) %>%  #angle
  mutate(thetaC = pi/2) %>% #radius
  mutate(x = sin(thetaC)*cos(phiC), y = sin(thetaC)*sin(phiC), z = cos(thetaC)) %>% #(x,y,z) coord
  mutate(xr = x, yr = z, zr = -y) %>% # +90 rotation about x-axis
  mutate(xrr = xr, yrr = yr, zrr = zr) %>% # do nothing
  mutate(theta = acos(zrr/sqrt(xrr^2+yrr^2+zrr^2)), phi = acos(xrr/sqrt(xrr^2+yrr^2))) %>%
  mutate(theta_deg = theta/pi*180/180*diff(buchner_theta)+buchner_theta[1], 
         phi_deg = phi/pi*180/180*diff(buchner_phi)+buchner_phi[1]) %>% # longitude use buchner_phi
  select(phi_deg, theta_deg) %>%
  data.matrix()


# see "stim area on sphere" section for plotting.


# Stimuli ---------------------------------------------------------------------------------------------------------

# -- stim disk
disk_theta <- 15 # radiuCs in deg

disk_phi <- seq(0, 359, by = 10)
disk_grid <- expand.grid(disk_theta, disk_phi)
colnames(disk_grid) <- c('theta', 'phi')
disk_grid %<>% as_tibble() %>%
  mutate(x = sin(theta/180*pi)*cos(phi/180*pi), y = sin(theta/180*pi)*sin(phi/180*pi), z = cos(theta/180*pi)) %>%
  data.matrix()

# stim positions
# stim_azim <- c(32.5, 45, 57.5, 70) # radius in deg
stim_elev <- c(-25, -12.5, 0, 12.5, 25) + 90
stim_azim <- seq(7.5, 157.5, by = 12.5) # radius in deg
# stim_elev <- seq(-87.5, 87.5, by = 12.5) + 90

stim_3 <- c(-30,0,30) 
stim_grid_azim <- expand.grid(stim_3 + 90, stim_azim)
colnames(stim_grid_azim) <- c('theta', 'phi')
stim_grid_elev <- expand.grid(stim_3 + 45, stim_elev)
stim_grid_elev <- stim_grid_elev[, c(2,1)]
colnames(stim_grid_elev) <- c('theta', 'phi')
stim_grid <- rbind(stim_grid_azim, stim_grid_elev)
stim_grid %<>% as_tibble() %>%
  mutate(x = sin(theta/180*pi)*cos(phi/180*pi), y = sin(theta/180*pi)*sin(phi/180*pi), z = cos(theta/180*pi)) %>%
  data.matrix()

# add 20-deg pitch
stim_grid_orig <- stim_grid
ang <- 20/180*pi
Ry <- matrix(c(cos(ang), 0 , sin(ang),
               0, 1, 0,
               -sin(ang), 0, cos(ang)), ncol = 3, byrow = T)
stim_grid[, c('x','y','z')] <- t(Ry %*% t(stim_grid[, c('x','y','z')]))
stim_grid %<>% as_tibble() %>%
  mutate(theta = acos(z)/pi*180, phi = acos(x / sqrt(x^2+y^2))/pi*180) %>%
  data.matrix()

stim_poly <- list()
stim_com <- list()
for (j in 1:dim(stim_grid)[1]) {
  t <- stim_grid[j, 'theta']/180*pi
  p <- stim_grid[j, 'phi']/180*pi
  
  Ry <- matrix(c(cos(t), 0, sin(t),
                 0, 1, 0,
                 -sin(t), 0, cos(t)),
               ncol = 3, byrow = T)
  
  Rz <- matrix(c(cos(p), -sin(p), 0,
                 sin(p), cos(p), 0,
                 0, 0, 1),
               ncol = 3, byrow = T)
  
  R <- Rz %*% Ry
  
  stim_poly[[j]] <- rbind(stim_grid[j,c('x','y','z')], disk_grid[,c('x','y','z')] %*% t(R))
  colnames(stim_poly[[j]]) <- c('x','y','z')
  
  stim_poly[[j]] %<>% as_tibble() %>%
    mutate(theta = acos(z)) %>%
    mutate(phi = 2*pi*(y < 0) + (-1)^(y < 0)*acos(x/sin(theta))) %>%
    mutate(theta = theta/pi*180, phi = phi/pi*180) %>%
    # mutate(theta = replace(phi, phi > 180, phi - 360)) %>% # not working
    mutate(phi = if_else(phi > 180, phi - 360, phi)) %>% #move to [-pi pi]
    data.matrix()
  stim_com[[j]] <- stim_poly[[j]][1,]
  stim_poly[[j]] <- stim_poly[[j]][-1,] #remove center
  xyMollweide <- Mollweide(stim_poly[[j]][,c('theta', 'phi')])
  colnames(xyMollweide) <- c("xM","yM")
  stim_poly[[j]] <- cbind(stim_poly[[j]], xyMollweide)
}

# 3d plot
nopen3d()
spheres3d(0,0,0,1, col='grey', alpha=0.2)
planes3d(0,0,1, 0, alpha = 0.2)
planes3d(0,1,0, 0, alpha = 0.2)
points3d(stim_grid[,3:5])
points3d(disk_grid[,3:5], col = 'blue')
points3d(t(Ry %*% t(disk_grid[,3:5])), col = 'blue')
points3d(disk_grid[,3:5] %*% t(Ry), col = 'blue')


# stim area on sphere (Fig.3g) -------------------------------------------------

# sphere grid for computing areas
S2grid <- matrix(c(0,0), ncol = 2) # [theta phi] in rad
darc <- pi/180*1
for (ti in 1:179) {
  theta <- ti * darc / 1
  ddeg <- darc / sin(theta) / pi *180
  phiN <- floor(180 / ddeg)
  tpmat <- cbind(rep(theta/pi*180, 2*phiN+1), seq(-phiN, phiN, by = 1)*ddeg )
  S2grid <- rbind(S2grid, tpmat)
}
S2grid <- rbind(S2grid, c(180,0)) #south pole
darea <- 4*pi/dim(S2grid)[1]
S2grid <- cbind(S2grid, Mollweide(S2grid))
colnames(S2grid) <- c('t','p','xM','yM')

# background grid
bkgd_grid <- S2grid[, c('xM','yM')] # use S2 grid 

# calculate stim LC4 overlap
stim_LC4_ol <- list()
for (j in 1:length(stim_poly)) {
  in_stim <- sp::point.in.polygon(bkgd_grid[,1], bkgd_grid[,2], stim_poly[[j]][,'xM'], stim_poly[[j]][,'yM'])
  ol_tmp <- c()
  for (k in 1:length(LC4)) {
    in_neu <- sp::point.in.polygon(bkgd_grid[,1], bkgd_grid[,2], xy_poly[[k]][,'xM'], xy_poly[[k]][,'yM'])
    ol_tmp <- c(ol_tmp, sum(in_stim & in_neu))
  }
  stim_LC4_ol[[j]] <- ol_tmp
}

# LC4 area
LC4_grid_N <- c()
for (k in 1:length(LC4)) {
  in_neu <- sp::point.in.polygon(bkgd_grid[,1], bkgd_grid[,2], xy_poly[[k]][,'xM'], xy_poly[[k]][,'yM'])
  LC4_grid_N <- c(LC4_grid_N, sum(in_neu))
}

# stim area
stim_grid_N <- c()
for (k in 1:length(stim_poly)) {
  in_neu <- sp::point.in.polygon(bkgd_grid[,1], bkgd_grid[,2], stim_poly[[k]][,'xM'], stim_poly[[k]][,'yM'])
  stim_grid_N <- c(stim_grid_N, sum(in_neu))
}


# set up background guidelines
bkgd_eq <- Mollweide(cbind(rep(90,37), seq(-180, 180, by = 10)))
bkgd_eq_p45 <- Mollweide(cbind(rep(45,37), seq(-180, 180, by = 10)))
bkgd_eq_m45 <- Mollweide(cbind(rep(135,37), seq(-180, 180, by = 10)))
bkgd_mer <- Mollweide(cbind(seq(0, 180, by = 10), rep(0,19)))
bkgd_mer_e <- Mollweide(cbind(seq(0, 180, by = 10), rep(90,19)))
bkgd_mer_ee <- Mollweide(cbind(seq(0, 180, by = 10), rep(180,19)))
bkgd_mer_w <- Mollweide(cbind(seq(0, 180, by = 10), rep(-90,19)))
bkgd_mer_ww <- Mollweide(cbind(seq(0, 180, by = 10), rep(-180,19)))

# dev.new()
# # plot(bkgd_grid, pch = '.', cex = 3)
# plot(bkgd_grid, cex = 0.6, pch='')
# lines(bkgd_eq); lines(bkgd_eq_m45); lines(bkgd_eq_p45)
# lines(bkgd_mer); lines(bkgd_mer_e); lines(bkgd_mer_w); lines(bkgd_mer_ee); lines(bkgd_mer_ww)
# for (j in 1:(length(xy_poly))  ) {
#   polygon(xy_poly[[j]][,c("xM","yM")], border = j, lwd = 2)
#   text(x = xy_com[[j]]['xM'], y = xy_com[[j]]['yM'], label = LC4_grid_N[j], pos = 2, offset = 0.2)
# }
# title("LC4_Mollweide")


# PLOT Fig.3g, LC4 polygons on background grid Mollweide
dev.new()
# pdf(file = "LC4 on eye Mollweide.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plot(bkgd_grid, cex = 0.6, pch='', xlim = c(-0.5, pi), ylim = c(1.5,-1.5))
lines(bkgd_eq); lines(bkgd_eq_m45); lines(bkgd_eq_p45)
lines(bkgd_mer); lines(bkgd_mer_e); lines(bkgd_mer_ee)
for (j in 1:length(xy_poly)) {
  if (j == landmk[1]) {
    polygon(xy_poly[[j]][,c("xM","yM")], col = landmk_col[1], density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]]['xM'], xy_com[[j]]['yM'], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
  else if (j == landmk[2]) {
    polygon(xy_poly[[j]][,c("xM","yM")], col = landmk_col[2], density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]]['xM'], xy_com[[j]]['yM'], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
  else {
    polygon(xy_poly[[j]][,c("xM","yM")], col = "gray", border = 'black', density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]]['xM'], xy_com[[j]]['yM'], col="blue", cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
  }
}
title("LC4_Mollweide")
# dev.off()

#  FIG, LC4 and stim ----------------------------------------------------------------------------------------------

# ### ###one eye equirectangular
bd_phi <- seq(buchner_phi[1], buchner_phi[2], by = 1)
# bd_theta <- seq(0, 180, by = 1)
bd_theta <- seq(buchner_theta[1], buchner_theta[2], by = 1)
xy_bd <- matrix(ncol = 2)
bd_grid <- expand.grid(bd_phi, bd_theta) 

windows(record = F, width = 8, height = 8)
# pdf(file = "stim_LC4.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
# postscript(file="stim_LC4_equi.eps",paper="special",width=8,height=8,horizontal=F)
plot(bd_grid, ylim = rev(range(bd_grid)), type = "n", axes = T, ann = F, asp = 1)
for (j in 1: (length(xy_poly)) ) {
  xy_bd <- rbind(xy_bd, xy_poly[[j]][,c('phi_deg', 'theta_deg')])
  if (j == landmk[[1]]) {
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "#d7191c", density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 3, pch = 3)
    # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 2, pch = 20) #pch=1 circle, 32+j ASCII
    # text(x = xy_com[[j]]['phi_deg'], y = xy_com[[j]]['theta_deg'], label = LC4_grid_N[j], pos = 2, offset = 0.2)
  }
  else if (j == landmk[2]) {
    polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "#2c7bb6", density = 20, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 3, pch = 3)
    # points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 2, pch = 20) #pch=1 circle, 32+j ASCII
    # text(x = xy_com[[j]]['phi_deg'], y = xy_com[[j]]['theta_deg'], label = LC4_grid_N[j], pos = 2, offset = 0.2)
  }
  else {
    # polygon(xy_poly[[j]][,c('phi_deg', 'theta_deg')], col = "grey", border = 'black', density = 10, angle = j*2, lwd = 2)
    points(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], col="blue", cex = 2, pch = 20) #pch=1 circle, 32+j ASCII
    text(xy_com[[j]][c('phi_deg')], xy_com[[j]][c('theta_deg')], labels = paste(j), pos = 2, offset = 0.2)
    # text(x = xy_com[[j]]['phi_deg'], y = xy_com[[j]]['theta_deg'], label = conn_target[[1]][j,"tofrom_glu"], pos = 2, offset = 0.2)
    # text(x = xy_com[[j]]['phi_deg'], y = xy_com[[j]]['theta_deg'], label = LC4_grid_N[j], pos = 2, offset = 0.2)
  }
}
# for (j in 1:length(stim_poly)) {
for (j in 1:(3*length(stim_azim))) { #azim sweep
  # for (j in 1:39) { #azim sweep
  phithe <- stim_poly[[j]][,c('phi', 'theta')]
  phithe[,1] <- if_else(phithe[,1]< -90, phithe[,1]+360, phithe[,1])
  polygon(phithe, col = "green", density = 0, angle = j*2, lwd = 2)
  # text(x = stim_com[[j]]['phi'], y = stim_com[[j]]['theta'], label = stim_grid_N[j], pos = 2, offset = 0.2)
}
for (j in (3*length(stim_azim)+1):length(stim_poly)) { #elev sweep
  polygon(stim_poly[[j]][,c('phi', 'theta')], col = "gold", density = 0, angle = j*2, lwd = 2)
}
# add boundary
xy_bd <- xy_bd[-1,]
hpts_2 <- chull(xy_bd)
hpts_2 <- c(hpts_2, hpts_2[1])
xy_bd_chull <- xy_bd[hpts_2,] # hull edge points
xy_bd_chull <- cbind(xy_bd_chull, Mollweide(xy_bd_chull[, c('theta_deg', 'phi_deg')]))
colnames(xy_bd_chull) <- c('phi_deg', 'theta_deg','xM', 'yM')
polygon(xy_bd_chull[, c('phi_deg', 'theta_deg')])
# add guide lines
lines(rbind(c(-11,180), c(-2,180)), lwd = 3)
text(-5, 180, labels = expression(paste("9",degree)), pos = 1, offset = 0.3)
lines(rbind(c(-11,180), c(-11,171)), lwd = 3)
text(-11, 175, labels = expression(paste("9",degree)), pos = 2, offset = 0.2)
lines(rbind(c(0,0), c(0,180)), lwd = 3) 
text(0, -5, labels = "front", pos = 1, offset = 0)
lines(rbind(c(90,0), c(90,180)), lwd = 3)
text(90, -5, labels = expression(paste("side 90",degree)), pos = 1, offset = 0)
lines(rbind(c(-12,90), c(162,90)), lwd = 3)
text(-17, 90, labels = "equator", pos = 1, offset = 0, srt = 90)
title("stim and LC4, equirectangular, Mollweide grid num, stim grid num ~ 696 +/-2.5")

# dev.off()

# lines(rbind(c(10,180), c(19,180)), lwd = 3)
# text(20, 180, labels = expression(paste("9",degree)), pos = 1, offset = 0.3)
# lines(rbind(c(10,180), c(10,171)), lwd = 3)
# text(10, 170, labels = expression(paste("9",degree)), pos = 2, offset = 0.2)
# lines(rbind(c(0,3), c(0,177)), lwd = 3) 
# text(0, -3, labels = "front", pos = 1, offset = 0)
# lines(rbind(c(90,3), c(90,177)), lwd = 3)
# text(90, -3, labels = expression(paste("side 90",degree)), pos = 1, offset = 0)
# lines(rbind(c(-8,90), c(160,90)), lwd = 3)
# text(-13, 90, labels = "equator", pos = 1, offset = 0, srt = 90)


# ### ### both eye, Mollweide
bd_grid <- bkgd_grid 
xy_bd_M <- matrix(ncol = 2)

windows(record = F, width = 8, height = 8)
# pdf(file = "stim_LC4_Mollweid.pdf", width=8, height=8,pointsize=12,family="Helvetica", useDingbats = F)
# postscript(file="stim_LC4_Mollweid.eps",paper="special",width=8,height=8,horizontal=F)
plot(bd_grid, ylim = (range(bd_grid[,'yM'])), type = "n", axes = T, ann = F, asp = 1)
for (j in 1:(length(xy_poly)) ) {
  xy_bd_M <- rbind(xy_bd_M, xy_poly[[j]][,c('xM', 'yM')])
  if (j == landmk[[1]]) {
    # polygon(xy_poly[[j]][,c('xM', 'yM')], col = "#d7191c", density = 20, angle = j*2, lwd = 2)
    polygon(xy_poly[[j]][,c('xM', 'yM')], border = "#d7191c", density=NULL, lwd = 2)
    points(xy_com[[j]][c('xM')], xy_com[[j]][c('yM')], col="#d7191c", cex = 2, pch = 16)
    # text(x = xy_com[[j]]['xM'], y = xy_com[[j]]['yM'], label = LC4_grid_N[j], pos = 2, offset = 0.2)
  }
  else if (j == landmk[2]) {
    # polygon(xy_poly[[j]][,c('xM', 'yM')], col = "#2c7bb6", density = 20, angle = j*2, lwd = 2)
    polygon(xy_poly[[j]][,c('xM', 'yM')], border = "#2c7bb6",density=NULL, lwd = 2)
    # points(xy_com[[j]][c('xM')], xy_com[[j]][c('yM')], col="blue", cex = 2, pch = 3)
    points(xy_com[[j]][c('xM')], xy_com[[j]][c('yM')], col="#2c7bb6", cex = 2, pch = 16)
    # text(x = xy_com[[j]]['xM'], y = xy_com[[j]]['yM'], label = LC4_grid_N[j], pos = 2, offset = 0.2)
  }
  else {
    points(xy_com[[j]][c('xM')], xy_com[[j]][c('yM')], col="black", cex = 1, pch = 16) #pch=1 circle, 32+j ASCII
    # text(x = xy_com[[j]]['xM'], y = xy_com[[j]]['yM'], label = LC4_grid_N[j], pos = 2, offset = 0.2)
  }
}
# for (j in 1:length(stim_poly)) {
for (j in 1:(3*length(stim_azim))) { #azim sweep
  polygon(stim_poly[[j]][,c('xM', 'yM')], col = "green", density = 0, angle = j*2, lwd = 2)
  # text(x = stim_com[[j]]['phi'], y = stim_com[[j]]['theta'], label = stim_grid_N[j], pos = 2, offset = 0.2)
}
for (j in (3*length(stim_azim)+1):length(stim_poly)) { #elev sweep
  polygon(stim_poly[[j]][,c('xM', 'yM')], col = "gold", density = 0, angle = j*2, lwd = 2)
}
lines(bkgd_eq); lines(bkgd_eq_m45); lines(bkgd_eq_p45)
lines(bkgd_mer); lines(bkgd_mer_e); lines(bkgd_mer_w); lines(bkgd_mer_ee); lines(bkgd_mer_ww)
# polygon(xy_bd_chull[, c('xM', 'yM')]) # need re-make after different projecton
title("stim and LC4, equirectangular, Mollweide grid num, stim grid num ~ 696 +/-2.5")

# add boundary for Mollweide (this is slightly different than converting xy_bd to Mollweide)
xy_bd_M <- xy_bd_M[-1,]
hpts_M <- chull(xy_bd_M)
hpts_M <- c(hpts_M, hpts_M[1])
xy_bd_chull_M <- xy_bd_M[hpts_M,] # hull edge points
colnames(xy_bd_chull_M) <- c('xM', 'yM')
polygon(xy_bd_chull_M[, c('xM', 'yM')])

# dev.off()

#  Eyal plot ------------------------------------------------------------------------------------------------------
# This shows the levels of overlap between LC4 dendrites in lobula

bd_phi2 <- seq(buchner_phi[1], buchner_phi[2], by = 2)
# bd_theta2 <- seq(1, 180, by = 2)
bd_theta2 <- seq(buchner_theta[1], buchner_theta[2], by = 2)
pt_grid <- expand.grid(bd_phi2, bd_theta2) 
pt_grid <- cbind(pt_grid, rep(0, dim(pt_grid)[1]))
colnames(pt_grid) <- c('x','y','paint')
for (j in 1:length(xy_poly)) {
  ii <- sp::point.in.polygon(pt_grid[,"x"], pt_grid[,"y"], xy_poly[[j]][,"phi_deg"], xy_poly[[j]][,"theta_deg"])
  pt_grid[,"paint"] <- pt_grid[,"paint"] + ii
}
pt_grid <- pt_grid[pt_grid[,"paint"] != 0, ] 
pt_grid[,"paint"] <- factor(pt_grid[,"paint"], labels = sort(unique(pt_grid[,"paint"])) )

# plot
dev.new()
ggplot(pt_grid, aes(x, y, z = paint)) +
  geom_raster(aes(fill = paint), interpolate = F) +
  scale_y_reverse() +
  scale_fill_manual(values = pal_tar, guide_legend("overlap count"), labels = paste(seq(1,8,length.out = 8))) +
  geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
  annotate("text", x = -5, y = 185, label = "9?") +
  geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
  annotate("text", x = -15, y = 175.5, label = "9?") +
  geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 1) +
  geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 1) +
  geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 1) +
  annotate("text", x = 5, y = 175.5, label = "9?") +
  theme_void() +
  labs(title = "overlap count in LO") +
  coord_fixed(ratio = 1)


