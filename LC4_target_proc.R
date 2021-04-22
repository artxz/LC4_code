# target neuron anatomical RF


# bd_grid <- expand.grid(bd_phi, bd_theta) # equirectangular
# bd_grid <- bkgd_grid # Mollweide


# Mollweide  -----------------------------------------------------------------------------

bd_grid <- S2grid

bkgd_eq <- Mollweide(cbind(rep(90,91), seq(0, 180, by = 2)))
colnames(bkgd_eq) <- c('xM','yM')
# ii_inpoly <- sp::point.in.polygon(bkgd_eq[,'xM'], bkgd_eq[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])
# bkgd_eq <- bkgd_eq[as.logical(ii_inpoly),]
bkgd_eq_p45 <- Mollweide(cbind(rep(45,91), seq(0, 180, by = 2)))
colnames(bkgd_eq_p45) <- c('xM','yM')
# ii_inpoly <- sp::point.in.polygon(bkgd_eq_p45[,'xM'], bkgd_eq_p45[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])
# bkgd_eq_p45 <- bkgd_eq_p45[as.logical(ii_inpoly),]
bkgd_eq_m45 <- Mollweide(cbind(rep(135,91), seq(0, 180, by = 2)))
colnames(bkgd_eq_m45) <- c('xM','yM')
# ii_inpoly <- sp::point.in.polygon(bkgd_eq_m45[,'xM'], bkgd_eq_m45[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])
# bkgd_eq_m45 <- bkgd_eq_m45[as.logical(ii_inpoly),]
bkgd_mer <- Mollweide(cbind(seq(0, 180, by = 10), rep(0,19)))
bkgd_mer_e <- Mollweide(cbind(seq(180, 0, by = -10), rep(90,19)))
bkgd_mer_ee <- Mollweide(cbind(seq(0, 180, by = 10), rep(180,19)))
bkgd_mer <- rbind(bkgd_mer,bkgd_mer_e,bkgd_mer_ee)
colnames(bkgd_mer) <- c('xM','yM')

# PLOT,  Gaussian around each com with synap# as height,  cp data via binning
ii_inpoly <- sp::point.in.polygon(bd_grid[,'xM'], bd_grid[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])
plvl <- list()
# simdata <- list()
simdata_df <- list()
# LC4_tar_median <- list()
# for (j in 1:length(conn_target)) {
#   LC4_tar_median[[j]] <- (quantile(conn_target[[j]]$tofrom_glu, c(0.0)))
# }
mat_names <- c('DNp01', 'DNp02', 'DNp11', 'DNp04')


for (j in 1:length(neu_target)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  # grid_Gaussian <- as.data.frame(bd_grid[, c('t', 'p')])
  grid_Gaussian$arcL = 0
  grid_Gaussian$Z = 0
  # colnames(grid_Gaussian) <- c("X","Y","Z")
  for (k in 1:length(neu)) {
    t0 <- xy_com[[k]]['theta_deg'] / 180 * pi
    p0 <- xy_com[[k]]["phi_deg"] / 180 * pi
    tp0 <- c(t0,p0)
    A <- conn_target[[j]][k,"tofrom_glu"]
    # DEFINE denominator as half-width at 10% ==> if -x^2/s^2, then s^2 = area/pi/log(10), area = N*darc^2
    # need arc length here
    grid_Gaussian$arcL <- apply(grid_Gaussian, 1, function(x) arcLength(tp0, x[1:2]/180*pi) ) 
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[6] + 1*A*exp(-x[5]^2 / (darc^2*LC4_grid_N[[j]]/pi/log(10)*2)) ) 
  }
  
  simdata_df[[j]] <- as.data.frame(grid_Gaussian[, c('xM','yM','Z')])
  # simdata_df[[j]] <- simdata[[j]]
  colnames(simdata_df[[j]]) <- c("x","y","z")
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, breaks_3)
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  
  plvl[[j]] <-
    ggplot(simdata_df[[j]], aes(x, y)) + 
    # geom_tile(aes(fill = z)) +
    # geom_point()
    # geom_raster(data = simdata_df[[j]], aes(x, y, fill = equalSpace), interpolate = F) +
    geom_tile(aes(fill = equalSpace), height = 0.025, width = 0.025) +
    # geom_tile(data = simdata_df[[j]][simdata_df[[j]]$z > 0.8, ], aes(x, y, fill = equalSpace), size = 3) +
    scale_y_reverse() +
    scale_fill_manual(values = pal_tar, guide_legend("synp den"),
                      labels = paste(seq(0.1,1,length.out = 10))) +
    # geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
    # annotate("text", x = -5, y = 185, label = "9?") +
    # geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
    # annotate("text", x = -15, y = 175.5, label = "9?") +
    # geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
    # geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
    # geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
    # annotate("text", x = 5, y = 175.5, label = "9?") +
    # # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
    # # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
    # # labs(title = paste("Target bi")) +
  geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
    # geom_polygon(data = as.data.frame(xy_bd_chull_M), aes(xM, yM))
    geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
    theme_void() +
    labs(title = mat_names[j]) +
    coord_fixed(ratio = 1)
  
  # for (k in 1:length(neu)) {
  #   plvl[[j]] <- plvl[[j]] +
  #     annotate("text", x = xy_com[[k]]['xM'], y = xy_com[[k]]['yM'], label = conn_target[[j]][k,"tofrom_glu"])
  # }
}

# FIG
windows(width = 8, height = 8)
# postscript(file="RF_DNp02.eps",paper="special",width=8,height=8,horizontal=F)
# pdf(file = "RF_DNp02.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plvl[[1]] 
# dev.off()
# ggsave("RF_DNp02.eps")
windows(width = 8, height = 8)
# postscript(file="RF_DNp11.eps",paper="special",width=8,height=8,horizontal=F)
plvl[[2]] 
# dev.off()
windows(width = 8, height = 8)
# postscript(file="RF_DNp04.eps",paper="special",width=8,height=8,horizontal=F)
plvl[[3]] 
# dev.off()
windows(width = 8, height = 8)
# postscript(file="RF_DNp04.eps",paper="special",width=8,height=8,horizontal=F)
plvl[[4]] 

# # plot syn num at com
# for (j in 1:2) {
#   plvl[[j]] <-
#     ggplot(as.data.frame(xy_bd_chull), aes(x = phi_deg, y = theta_deg )) +
#     geom_polygon(fill = rgb(1,1,1)) +
#     scale_y_reverse() +
#     # scale_fill_manual(values = pal_tar, guide_legend("synp den"), labels = paste(seq(0.1,1,length.out = 10))) +
#     geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
#     annotate("text", x = -5, y = 185, label = "9?") +
#     geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
#     annotate("text", x = -15, y = 175.5, label = "9?") +
#     geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
#     geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
#     geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
#     # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
#     # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
#     # labs(title = paste("Target bi")) +
#     labs(title = mat_names[1]) +
#     coord_fixed(ratio = 1)
#   
#   for (k in 1:length(neu)) {
#     plvl[[j]] <- plvl[[j]] +
#       annotate("text", x = xy_com[[k]]['phi_deg'], y = xy_com[[k]]['theta_deg'], label = conn_target[[j]][k,"tofrom_glu"])
#   }
# }




#  NEW, stim RF, (synap# x stim overlap) as height, sqrt(area) as width ---------------------------

bd_grid <- S2grid

bkgd_eq <- Mollweide(cbind(rep(90,181), seq(-180, 180, by = 2)))
colnames(bkgd_eq) <- c('xM','yM')
ii_inpoly <- sp::point.in.polygon(bkgd_eq[,'xM'], bkgd_eq[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])
bkgd_eq <- bkgd_eq[as.logical(ii_inpoly),]
bkgd_eq_p45 <- Mollweide(cbind(rep(45,181), seq(-180, 180, by = 2)))
colnames(bkgd_eq_p45) <- c('xM','yM')
ii_inpoly <- sp::point.in.polygon(bkgd_eq_p45[,'xM'], bkgd_eq_p45[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])
bkgd_eq_p45 <- bkgd_eq_p45[as.logical(ii_inpoly),]
bkgd_eq_m45 <- Mollweide(cbind(rep(135,181), seq(-180, 180, by = 2)))
colnames(bkgd_eq_m45) <- c('xM','yM')
ii_inpoly <- sp::point.in.polygon(bkgd_eq_m45[,'xM'], bkgd_eq_m45[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])
bkgd_eq_m45 <- bkgd_eq_m45[as.logical(ii_inpoly),]
bkgd_mer <- Mollweide(cbind(seq(0, 180, by = 10), rep(0,19)))
bkgd_mer_e <- Mollweide(cbind(seq(180, 0, by = -10), rep(90,19)))
bkgd_mer_ee <- Mollweide(cbind(seq(0, 180, by = 10), rep(180,19)))
bkgd_mer <- rbind(bkgd_mer,bkgd_mer_e,bkgd_mer_ee)
colnames(bkgd_mer) <- c('xM','yM')

# PLOT,  Gaussian around each com with synap# as height,  cp data via binning
ii_inpoly <- sp::point.in.polygon(bd_grid[,'xM'], bd_grid[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])
target_names <- c('DNp01', 'DNp02', 'DNp11', 'DNp04')
stim_names <- seq(1,length(stim_poly))

n_lvl <- 11
getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
pal_tar <- getPalette(n_lvl - 1)

plvl_stim <- list()
simdata_df_stim <- list()
stim_max_r <- c()
Avalue <- list()
for (m in 1:length(neu_target)) {
  simdata_df <- list()
  AA <- c()
  
  for (j in 1:length(stim_poly)) {
    
    nn <- c()
    grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
    grid_Gaussian$arcL = 0
    grid_Gaussian$Z = 0
    
    ## ###
    for (k in 1:length(neu)) {
      t0 <- xy_com[[k]]['theta_deg'] / 180 * pi
      p0 <- xy_com[[k]]["phi_deg"] / 180 * pi
      tp0 <- c(t0,p0)
      # A <- conn_target[[m]][k,"tofrom_glu"]
      A <- conn_target[[m]][k,"tofrom_glu"] * stim_LC4_ol[[j]][[k]] / LC4_grid_N[k]
      # DEFINE denominator as half-width at 10% ==> if -x^2/s^2, then s^2 = area/pi/log(10), area = N*darc^2
      # need arc length here
      if (A > 0) {
        nn <- c(nn, k)
        grid_Gaussian$arcL <- apply(grid_Gaussian, 1, function(x) arcLength(tp0, x[1:2]/180*pi) )
        grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[6] + 1*A*exp(-x[5]^2 / (darc^2*LC4_grid_N[[k]]/pi/log(10)*2)) )
      }
    }
    
    
    # ### ###
    # # - winner takes all
    # k <- which.max(stim_LC4_ol[[j]])
    # t0 <- xy_com[[k]]['theta_deg'] / 180 * pi
    # p0 <- xy_com[[k]]["phi_deg"] / 180 * pi
    # tp0 <- c(t0,p0)
    # A <- conn_target[[m]][k,"tofrom_glu"]
    # # A <- conn_target[[m]][k,"tofrom_glu"] * stim_LC4_ol[[j]][[k]] / LC4_grid_N[k]
    # AA <- c(AA, A)
    # # DEFINE denominator as half-width at 10% ==> if -x^2/s^2, then s^2 = area/pi/log(10), area = N*darc^2
    # # need arc length here
    # if (A > 0) {
    #   nn <- c(nn, k)
    #   grid_Gaussian$arcL <- apply(grid_Gaussian, 1, function(x) arcLength(tp0, x[1:2]/180*pi) )
    #   grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[6] + 1*A*exp(-x[5]^2 / (darc^2*LC4_grid_N[[k]]/pi/log(10)*2)) )
    # }
    
    
    simdata_df[[j]] <- as.data.frame(grid_Gaussian[, c('xM','yM','Z')])
    AA <- c(AA, sum(simdata_df[[j]]$Z))
    colnames(simdata_df[[j]]) <- c("x","y","z")
  }
  simdata_df_stim[[m]] <- simdata_df
  Avalue[[m]] <- AA
  # stim_max_r <- c(stim_max_r, max( unlist((lapply(simdata_df_stim[[m]], function(x) max(x$z)))) ) )
}


for (m in 1:length(neu_target)) {
  stim_max <- max( unlist((lapply(simdata_df_stim[[m]], function(x) max(x$z)))) )
  plvl <- list()
  
  for (j in 1:length(stim_poly)) { 
    simdata_df_stim[[m]][[j]]$z <- simdata_df_stim[[m]][[j]]$z / stim_max
    simdata_df_stim[[m]][[j]]$equalSpace <- cut(simdata_df_stim[[m]][[j]]$z, seq(0,1,length.out = n_lvl))
    
    plvl[[j]] <- 
      ggplot(simdata_df_stim[[m]][[j]], aes(x, y)) + 
      geom_tile(aes(fill = equalSpace), height = 0.025, width = 0.025) +
      scale_y_reverse() +
      scale_fill_manual(values = pal_tar, guide_legend("synp den"),
                        labels = paste(seq(0.1,1,length.out = 10))) +
      geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
      geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
      # geom_polygon(data = as.data.frame(xy_bd_chull_M), aes(xM, yM))
      geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
      geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
      theme_void() +
      theme(legend.position="none") +
      labs(title = stim_names[j]) +
      coord_fixed(ratio = 1)
    
    # # add ind of LC4
    # for (k in 1:length(neu)) {
    #   plvl[[j]] <- plvl[[j]] +
    #     annotate("text", x = xy_com[[k]]['phi_deg'], y = xy_com[[k]]['theta_deg'], label = conn_target[[m]][k,"tofrom_glu"])
    # }
  }
  plvl[[1]] <- plvl[[1]] + labs(title = paste(target_names[m], "azim sweep", sep = ','))
  plvl[[13]] <- plvl[[13]] + labs(title = paste(target_names[m], "elev sweep", sep = ','))
  plvl_stim[[m]] <- plvl
}


# # FIG
# dev.new()
# plvl[[1]]
# dev.new()
# plvl[[2]] 

m <- 3
windows(record = F, width = 12, height = 9)
plot_grid(plvl_stim[[m]][[1]],
          plvl_stim[[m]][[4]],
          plvl_stim[[m]][[7]],
          plvl_stim[[m]][[10]],
          plvl_stim[[m]][[2]],
          plvl_stim[[m]][[5]],
          plvl_stim[[m]][[8]],
          plvl_stim[[m]][[11]],
          plvl_stim[[m]][[3]],
          plvl_stim[[m]][[6]],
          plvl_stim[[m]][[9]],          
          plvl_stim[[m]][[12]],
          nrow = 3)
windows(record = F, width = 9, height = 15)
plot_grid(plvl_stim[[m]][[13]],
          plvl_stim[[m]][[14]],
          plvl_stim[[m]][[15]],
          plvl_stim[[m]][[16]],
          plvl_stim[[m]][[17]],
          plvl_stim[[m]][[18]],
          plvl_stim[[m]][[19]],
          plvl_stim[[m]][[20]],
          plvl_stim[[m]][[21]],
          plvl_stim[[m]][[22]],
          plvl_stim[[m]][[23]],          
          plvl_stim[[m]][[24]],
          plvl_stim[[m]][[25]],
          plvl_stim[[m]][[26]],
          plvl_stim[[m]][[27]],
          nrow = 5)





# - sum up stim to comp with exp
# azim
x <- c(32.5, 45, 57.5, 70)
pla <- list()
for (m in 1:3) {
  y <- c(sum(Avalue[[m]][1:3]),
         sum(Avalue[[m]][4:6]),
         sum(Avalue[[m]][7:9]),
         sum(Avalue[[m]][10:12]) )
  y <- y / y[1]
  # yrange <- range(y)
  df <- as.data.frame(cbind(x,y))
  pla[[m]] <- ggplot(df, aes(x,y)) +
    geom_line(colour = m) +
    # ylim(0, 4) +
    # labs(title = target_names[m])
    ylim(-0.5, 3) +
    theme_minimal()+
    labs(title = paste("azim,", target_names[m], sep = '') )
}


# elev
x <- rev(c(-25, -12.5, 0, 12.5, 25)) #reversed elev def
ple <- list()
for (m in 1:3) {
  y <- c(sum(Avalue[[m]][13:15]),
         sum(Avalue[[m]][16:18]),
         sum(Avalue[[m]][19:21]),
         sum(Avalue[[m]][22:24]),
         sum(Avalue[[m]][25:27]))
  y <- y / tail(y,1)
  df <- as.data.frame(cbind(x,y))
  ple[[m]] <- ggplot(df, aes(x,y)) +
    geom_line(colour = m) +
    ylim(-0.5, 3) +
    theme_minimal()+
    labs(title = paste("corr elev, ", target_names[m], sep = ''))
}


pl <- c(pla, ple)

windows(record = F, width = 9, height = 6)
# postscript(file="DNp_resp.eps",paper="special",width=8,height=8,horizontal=F)
plot_grid(pl[[1]], pl[[2]], pl[[3]],
          pl[[4]], pl[[5]], pl[[6]],
          nrow = 2)

dev.off()



# old equirectangular ---------------------------------------------------------------------------------------------

bd_grid <- expand.grid(bd_phi, bd_theta) # equirectangular

# PLOT,  Gaussian around each com with synap# as height,  cp data via binning
ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,'phi_deg'], xy_bd_chull[,'theta_deg'])
plvl <- list()
simdata <- list()
simdata_df <- list()
# LC4_tar_median <- list()
# for (j in 1:length(conn_target)) {
#   LC4_tar_median[[j]] <- (quantile(conn_target[[j]]$tofrom_glu, c(0.0)))
# }
mat_names <- c('DNp02', 'DNp11', 'DNp04')

avg_size_xy <- c()
for (j in 1:length(neu)) {
  tmp <- as.data.frame(xy_poly[[j]])
  avg_size_xy <- c(avg_size_xy, (max(tmp$theta_deg)-min(tmp$theta_deg)+max(tmp$phi_deg)-min(tmp$phi_deg))/4)
}
r_xy <- mean(avg_size_xy)
# half-width sqrt(-log(0.5))*r_xy*2


# ii_inpoly <- sp::point.in.polygon(bd_grid[,'xM'], bd_grid[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])


n_lvl <- 11
# breaks_3 <- seq(0,1,length.out = n_lvl)
# breaks_3 <- seq(0,70,length.out = n_lvl)
getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
pal_tar <- getPalette(n_lvl - 1)

for (j in 1:length(neu_target)) {
  # grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  # grid_Gaussian$Z = 0
  # for (k in 1:length(neu)) {
  #   x0 <- (xy_com[[k]]["phi_deg"])
  #   y0 <- (xy_com[[k]]["theta_deg"])
  #   A <- conn_target[[j]][k,"tofrom_glu"]
  #   # if (A >= LC4_tar_median[[j]]) { # selected neuron
  #   #   grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  #   # }
  #   grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
  # }
  
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$arcL = 0
  grid_Gaussian$Z = 0
  for (k in 1:length(neu)) {
    t0 <- xy_com[[k]]['theta_deg'] / 180 * pi
    p0 <- xy_com[[k]]["phi_deg"] / 180 * pi
    tp0 <- c(t0,p0)
    A <- conn_target[[j]][k,"tofrom_glu"]
    # DEFINE denominator as half-width at 10% ==> if -x^2/s^2, then s^2 = area/pi/log(10), area = N*darc^2
    # need arc length here
    grid_Gaussian$arcL <- apply(grid_Gaussian, 1, function(x) arcLength(tp0, x[c(2,1)]/180*pi) ) 
    grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[4] + 1*A*exp(-x[3]^2 / (darc^2*LC4_grid_N[[j]]/pi/log(10)*2)) ) 
  }
  
  grid_Gaussian_cut <- grid_Gaussian %>%
    # as_tibble() %>%
    # filter(X <= 117) %>% #9 x 11 + 2.25 +15 = 116.25
    # filter(Y >= 50 & Y <= 113) %>% #4 x 9 + 4.5 = 40.5, 2 x 9 +4.5 = 22.5
    # mutate(Z = Z*0.9) %>%
    as.data.frame()
  
  # # -- whole range no binning
  # simdata[[j]] <- grid_Gaussian_cut
  # simdata_df[[j]] <- simdata[[j]]
  # colnames(simdata_df[[j]]) <- c("x","y","z")
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
  # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # # simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, breaks_3)
  # simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  # -- whole range no binning
  simdata_df[[j]] <- as.data.frame(grid_Gaussian[, c(1,2,4)])
  # simdata_df[[j]] <- simdata[[j]]
  colnames(simdata_df[[j]]) <- c("x","y","z")
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  # simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, breaks_3)
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  
  plvl[[j]] <-
    ggplot(simdata_df[[j]], aes(x, y, z = z)) +
    geom_raster(aes(fill = equalSpace), interpolate = F) +
    scale_y_reverse() +
    scale_fill_manual(values = pal_tar, guide_legend("synp den"), labels = paste(seq(0.1,1,length.out = 10))) +
    geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
    annotate("text", x = -5, y = 185, label = "9?") +
    geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
    annotate("text", x = -15, y = 175.5, label = "9?") +
    geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
    geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
    geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
    annotate("text", x = 5, y = 175.5, label = "9?") +
    # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
    # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
    # labs(title = paste("Target bi")) +
    theme_void() +
    labs(title = mat_names[j]) +
    coord_fixed(ratio = 1)
  
  for (k in 1:length(neu)) {
    plvl[[j]] <- plvl[[j]] +
      annotate("text", x = xy_com[[k]]['phi_deg'], y = xy_com[[k]]['theta_deg'], label = conn_target[[j]][k,"tofrom_glu"])
  }
}

# FIG
dev.new()
plvl[[1]] 
dev.new()
plvl[[2]] 
dev.new()
plvl[[3]] 

# plot syn num at com
for (j in 1:2) {
  plvl[[j]] <-
    ggplot(as.data.frame(xy_bd_chull), aes(x = phi_deg, y = theta_deg )) +
    geom_polygon(fill = rgb(1,1,1)) +
    scale_y_reverse() +
    # scale_fill_manual(values = pal_tar, guide_legend("synp den"), labels = paste(seq(0.1,1,length.out = 10))) +
    geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
    annotate("text", x = -5, y = 185, label = "9?") +
    geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
    annotate("text", x = -15, y = 175.5, label = "9?") +
    geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
    geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
    geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
    # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
    # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
    # labs(title = paste("Target bi")) +
    labs(title = mat_names[1]) +
    coord_fixed(ratio = 1)
  
  for (k in 1:length(neu)) {
    plvl[[j]] <- plvl[[j]] +
      annotate("text", x = xy_com[[k]]['phi_deg'], y = xy_com[[k]]['theta_deg'], label = conn_target[[j]][k,"tofrom_glu"])
  }
}




#  stim RF, Gaussian around each com with (synap# x stim overlap) as height,   ----------------------------------------


ii_inpoly <- sp::point.in.polygon(bd_grid[,1], bd_grid[,2], xy_bd_chull[,'phi_deg'], xy_bd_chull[,'theta_deg'])
target_names <- c('DNp02', 'DNp11', 'DNp04')
stim_names <- seq(1,20)

n_lvl <- 11
getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
pal_tar <- getPalette(n_lvl - 1)

plvl_stim <- list()
simdata_df_stim <- list()
stim_max_r <- c()
for (m in 1:length(neu_target)) {
  simdata_df <- list()
  
  for (j in 1:length(stim_poly)) {
    grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
    grid_Gaussian$Z = 0
    colnames(grid_Gaussian) <- c("X","Y","Z")
    for (k in 1:length(neu)) {
      x0 <- (xy_com[[k]]["phi_deg"])
      y0 <- (xy_com[[k]]["theta_deg"])
      A <- conn_target[[m]][k,"tofrom_glu"] * stim_LC4_ol[[j]][[k]] / LC4_grid_N[k]
      grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x[3] + 1*A*exp(-(x[1]-x0)^2/r_xy^2 - (x[2]-y0)^2/r_xy^2))
    }
    
    # normalize
    simdata_df[[j]] <- grid_Gaussian
    colnames(simdata_df[[j]]) <- c("x","y","z")
    # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
    # simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
    # simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
    # simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  }
  
  simdata_df_stim[[m]] <- simdata_df
  stim_max_r <- c(stim_max_r, max( unlist((lapply(simdata_df_stim[[m]], function(x) max(x$z)))) ) )
  
}


for (m in 1:length(neu_target)) {
  stim_max <- max( unlist((lapply(simdata_df_stim[[m]], function(x) max(x$z)))) )
  
  plvl <- list()
  
  for (j in 1:length(stim_poly)) { 
    simdata_df_stim[[m]][[j]]$z <- simdata_df_stim[[m]][[j]]$z / stim_max
    simdata_df_stim[[m]][[j]]$equalSpace <- cut(simdata_df_stim[[m]][[j]]$z, seq(0,1,length.out = n_lvl))
    
    plvl[[j]] <-
      ggplot(simdata_df_stim[[m]][[j]], aes(x, y, z = z)) +
      geom_raster(aes(fill = equalSpace), interpolate = F) +
      scale_y_reverse() +
      # scale_fill_manual(values = pal_tar, guide_legend("synp den"), labels = paste(seq(0.1,1,length.out = 10))) +
      scale_fill_manual(values = pal_tar, guide = FALSE) +
      # scale_fill_manual(values = pal_tar) +
      geom_segment(aes(x = -11, y = 180, xend = -2, yend = 180), size=2, lineend = "round") +
      annotate("text", x = -5, y = 185, label = "9?") +
      geom_segment(aes(x = -11, y = 180, xend = -11, yend = 171), size=2, lineend = "round") +
      annotate("text", x = -15, y = 175.5, label = "9?") +
      geom_segment(aes(x = 0, y = 3, xend = 0, yend = 177), size = 2) +
      geom_segment(aes(x = 90, y = 3, xend = 90, yend = 177), size = 2) +
      geom_segment(aes(x = -12, y = 90, xend = 162, yend = 90), size = 2) +
      # annotate("text", x = 5, y = 175.5, label = "9?") +
      # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
      # theme(axis.title = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
      # labs(title = paste("Target bi")) +
      theme_void() +
      labs(title = stim_names[j]) +
      coord_fixed(ratio = 1)
    
    # # add ind of LC4
    # for (k in 1:length(neu)) {
    #   plvl[[j]] <- plvl[[j]] +
    #     annotate("text", x = xy_com[[k]]['phi_deg'], y = xy_com[[k]]['theta_deg'], label = conn_target[[m]][k,"tofrom_glu"])
    # }
  }
  plvl_stim[[m]] <- plvl
}



# # FIG
# dev.new()
# plvl[[1]] 
# dev.new()
# plvl[[2]] 

m <- 3
windows(record = F, width = 15, height = 15)

plot_grid(plvl_stim[[m]][[1]],
          plvl_stim[[m]][[2]],
          plvl_stim[[m]][[3]],
          plvl_stim[[m]][[4]],
          plvl_stim[[m]][[5]],
          plvl_stim[[m]][[6]],
          plvl_stim[[m]][[7]],
          plvl_stim[[m]][[8]],
          plvl_stim[[m]][[9]],
          plvl_stim[[m]][[10]],
          plvl_stim[[m]][[11]],          
          plvl_stim[[m]][[12]],          
          plvl_stim[[m]][[13]],          
          plvl_stim[[m]][[14]],          
          plvl_stim[[m]][[15]],          
          plvl_stim[[m]][[16]],
          plvl_stim[[m]][[17]],
          plvl_stim[[m]][[18]],
          plvl_stim[[m]][[19]],
          plvl_stim[[m]][[20]],
          nrow = 5)      

