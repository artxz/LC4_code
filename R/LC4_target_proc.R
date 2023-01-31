
load('data/S2grid.rda')
load('data/xy_bd_chull_M.rda')
load('data/xy_com.rda')
load('data/LC4_grid_N.rda')
load('data/stim_poly.rda')
load('data/stim_LC4_ol.rda')

darc <- pi/180*1

# Anatomically predicted receptive fields (Fig.3h) ---------------------------------------

bd_grid <- S2grid

# guidelines
bkgd_eq <- Mollweide(cbind(rep(90,91), seq(0, 180, by = 2)))
colnames(bkgd_eq) <- c('xM','yM')
bkgd_eq_p45 <- Mollweide(cbind(rep(45,91), seq(0, 180, by = 2)))
colnames(bkgd_eq_p45) <- c('xM','yM')
bkgd_eq_m45 <- Mollweide(cbind(rep(135,91), seq(0, 180, by = 2)))
colnames(bkgd_eq_m45) <- c('xM','yM')
bkgd_mer <- Mollweide(cbind(seq(0, 180, by = 10), rep(0,19)))
bkgd_mer_e <- Mollweide(cbind(seq(180, 0, by = -10), rep(90,19)))
bkgd_mer_ee <- Mollweide(cbind(seq(0, 180, by = 10), rep(180,19)))
bkgd_mer <- rbind(bkgd_mer,bkgd_mer_e,bkgd_mer_ee)
colnames(bkgd_mer) <- c('xM','yM')


# Gaussian density around each center-of-mass with synap# as height
ii_inpoly <- sp::point.in.polygon(bd_grid[,'xM'], bd_grid[,'yM'], xy_bd_chull_M[,'xM'], xy_bd_chull_M[,'yM'])
plvl <- list()
simdata_df <- list()
mat_names <- c('DNp01', 'DNp02', 'DNp11', 'DNp04')


for (j in 1:length(neu_target)) {
  grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
  grid_Gaussian$arcL = 0
  grid_Gaussian$Z = 0
  for (k in 1:length(LC4)) {
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
  colnames(simdata_df[[j]]) <- c("x","y","z")
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$z <- simdata_df[[j]]$z + 0.001
  simdata_df[[j]]$z <- simdata_df[[j]]$z / max(simdata_df[[j]]$z)
  simdata_df[[j]]$equalSpace <- cut(simdata_df[[j]]$z, seq(0,max(simdata_df[[j]]$z),length.out = n_lvl))
  
  
  plvl[[j]] <- ggplot(simdata_df[[j]], aes(x, y)) + 
    geom_tile(aes(fill = equalSpace), height = 0.025, width = 0.025) +
    scale_y_reverse() +
    scale_fill_manual(values = pal_tar, guide_legend("synp den"),labels = paste(seq(0.1,1,length.out = 10))) +
  geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50') +
    geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50') +
    theme_void() +
    labs(title = mat_names[j]) +
    coord_fixed(ratio = 1)
  
  # add conn num
  for (k in 1:length(LC4)) {
    plvl[[j]] <- plvl[[j]] +
      annotate("text", x = xy_com[[k]]['xM'], y = xy_com[[k]]['yM'], label = conn_target[[j]][k,"tofrom_glu"])
  }
}


# PLOT, Fig.3h plus 2 other types
windows(width = 8, height = 8)
# pdf(file = "RF_DNp01.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plvl[[1]] 
# dev.off()

windows(width = 8, height = 8)
# pdf(file = "RF_DNp02.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plvl[[2]] 
# dev.off()

windows(width = 8, height = 8)
# pdf(file = "RF_DNp11.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plvl[[3]] 
# dev.off()

windows(width = 8, height = 8)
# pdf(file = "RF_DNp04.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plvl[[4]] 
# dev.off()


#  kernel density estimate of reponse to stim (ED Fig.5d) ---------------------
# (synap# x stim overlap) as height, sqrt(area) as width

bd_grid <- S2grid

# guidelines
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

# Gaussian around each com with synapse No. as height,  
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
  AA <- c() #integrated response over eye
  
  for (j in 1:length(stim_poly)) {
    nn <- c()
    grid_Gaussian <- as.data.frame(bd_grid[ii_inpoly == 1,])
    grid_Gaussian$arcL = 0
    grid_Gaussian$Z = 0
    
    for (k in 1:length(LC4)) {
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
        grid_Gaussian$Z <- apply(grid_Gaussian, 1, function(x) x['Z'] + 1*A*exp(-1/2 * x['arcL']^2 / (darc^2*LC4_grid_N[[k]]/pi/log(10))) )
      }
    }
    
    simdata_df[[j]] <- as.data.frame(grid_Gaussian[, c('xM','yM','Z')])
    AA <- c(AA, sum(simdata_df[[j]]$Z))
    colnames(simdata_df[[j]]) <- c("x","y","z")
  }
  simdata_df_stim[[m]] <- simdata_df
  Avalue[[m]] <- AA
}

# plot
stim_elev <- c(-25, -12.5, 0, 12.5, 25) + 90
stim_azim <- seq(7.5, 157.5, by = 12.5) # radius in deg

for (m in 1:length(neu_target)) {
  stim_max <- max( unlist((lapply(simdata_df_stim[[m]], function(x) max(x$z)))) )
  plvl <- list()
  
  for (j in 1:length(stim_poly)) { 
    simdata_df_stim[[m]][[j]]$z <- simdata_df_stim[[m]][[j]]$z / stim_max
    simdata_df_stim[[m]][[j]]$equalSpace <- cut(simdata_df_stim[[m]][[j]]$z, seq(0,1,length.out = n_lvl))
    
    plvl[[j]] <- ggplot() +
      geom_tile(data=simdata_df_stim[[m]][[j]], aes(x, y, fill = equalSpace), height = 0.025, width = 0.025) +
      scale_y_reverse() +
      scale_fill_manual(values = pal_tar, guide_legend("synp den"),
                        labels = paste(seq(0.1,1,length.out = 10))) +
      geom_polygon(data=as.data.frame(stim_poly[[j]]), aes(xM,yM), colour ="gray30", fill=NA, lwd =0.5) +
      geom_path(data = as.data.frame(bkgd_mer), aes(x=xM, y=yM), colour = 'grey50', alpha=0.5) +
      geom_path(data = as.data.frame(bkgd_eq_m45), aes(x=xM, y=yM), colour = 'grey50', alpha=0.5) +
      geom_path(data = as.data.frame(bkgd_eq_p45), aes(x=xM, y=yM), colour = 'grey50', alpha=0.5) +
      geom_path(data = as.data.frame(bkgd_eq), aes(x=xM, y=yM), colour = 'grey50', alpha=0.5) +
      theme_void() +
      theme(legend.position="none") +
      labs(title = stim_names[j]) +
      coord_fixed(ratio = 1)
  }
  plvl[[1]] <- plvl[[1]] + labs(title = paste(target_names[m], "azim sweep", sep = ','))
  plvl[[3*length(stim_azim)+1]] <- plvl[[3*length(stim_azim)+1]] + labs(title = paste(target_names[m], "elev sweep", sep = ','))
  plvl_stim[[m]] <- plvl
}


# PLOT ED Fig.5d
m <- 4 # choose a target cell, {1,2,3,4}
# azimuth sweep
windows(width = 12, height = 9)
# pdf(file = paste("RF_", mat_names[m], "_azim.pdf", sep=''), width = 12, height = 9,pointsize=12,family="Helvetica", useDingbats = F)
plot_grid(plotlist = plvl_stim[[m]][1:(3*length(stim_azim))], nrow= 3, byrow= F)
# dev.off()

# elevation sweep
windows(width = 9, height =15)
# pdf(file = paste("RF_", mat_names[m], "_elev.pdf", sep=''), width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plot_grid(plotlist = plvl_stim[[m]][(3*length(stim_azim)+1):length(stim_poly)], nrow= 5)
# dev.off()


# sum up stim responses to compare with exp (Fig.3i and ED Fig.5f) ---------------

pal_ephy <- c(rgb(254,95,224,maxColorValue = 255),
              rgb(255,0,0,maxColorValue = 255),
              rgb(0,0,255,maxColorValue = 255),
              rgb(143,57,229,maxColorValue = 255))

# azim
x <- stim_azim
pla <- list()
for (m in 1:length(neu_target)) {
  y <- c()
  for (j in 1:length(stim_azim)) {
    y <- c(y, sum(Avalue[[m]][(3*j-2):(3*j)]))
  }
  y <- y / y[match(32.5, stim_azim)]
  # yrange <- range(y)
  df <- as.data.frame(cbind(x,y))
  pla[[m]] <- ggplot(df, aes(x,y)) +
    geom_line(colour = pal_ephy[m], lwd=2) +
    # ylim(0, 4) +
    # labs(title = target_names[m])
    coord_cartesian(ylim = c(0, max(3, ceiling(df$y)))) +
    scale_x_continuous(name = "angle [deg]", breaks = x, labels = x) +
    # coord_cartesian(ylim = c(0, 3)) +
    theme(panel.background = element_blank(), 
          axis.line = element_line(color='black'))+
          # axis.text.x = element_text(angle=45) ) +
    labs(title = paste("azim,", target_names[m], sep = '') )
}


# elev
x <- rev(c(-25, -12.5, 0, 12.5, 25)) #reversed elev def
# x <- stim_elev
ple <- list()
for (m in 1:length(neu_target)) {
  y <- c()
  for (j in (length(stim_azim)+1):length(c(stim_azim,stim_elev))) {
    y <- c(y, sum(Avalue[[m]][(3*j-2):(3*j)]))
  }
  y <- y / y[5]
  df <- as.data.frame(cbind(x,y))
  ple[[m]] <- ggplot(df, aes(x,y)) +
    geom_line(colour = pal_ephy[m], lwd=2) +
    theme(panel.background = element_blank(), 
          axis.line = element_line(color='black')) +
    coord_cartesian(ylim = c(0, 3)) +
    scale_x_continuous(name = "angle [deg]", breaks = x, labels = x) +
    labs(title = paste("elev, ", target_names[m], sep = ''))
}

pl <- c(pla, ple)


# PLOT, Fig.3i and ED Fig.5f
windows(record = F, width = 15, height = 6)
# postscript(file="DNp_resp.eps",paper="special",width=8,height=8,horizontal=F)
# pdf(file= "DNp_resp.pdf", width = 15, height =6,pointsize=12,family="Helvetica", useDingbats = F)
plot_grid(pl[[1]], pl[[2]], pl[[3]], pl[[4]],
          pl[[5]], pl[[6]], pl[[7]], pl[[8]],
          nrow = 2)
# dev.off()


# - overlay plot azim sweep
windows(width = 12, height =8)
# pdf(file= "DNp_azim_extended.pdf", width = 12, height =8,pointsize=12,family="Helvetica", useDingbats = F)
x <- stim_azim
for (m in 1:length(neu_target)) {
  y <- c()
  for (j in 1:length(stim_azim)) {
    y <- c(y, sum(Avalue[[m]][(3*j-2):(3*j)]))
  }
  y <- y / y[match(32.5, stim_azim)]
  # yrange <- range(y)
  if (m == 1 ) {
    plot(x, y, type = 'l', col=pal_ephy[m], xlim = range(x), ylim = c(0,6.5), xaxt="n", yaxt="n")
  } else {
    points(x, y, type = 'l', col=pal_ephy[m], xaxt="n", yaxt="n")
  }
  axis(side=1, at=x, labels = paste(x, "°", sep = ''))
  axis(side=2, at=seq(0,6), labels = seq(0,6))
}
dev.off()

# # SAVE data as csv
# x <- stim_azim
# df <- x
# for (m in 1:length(neu_target)) {
#   y <- c()
#   for (j in 1:length(stim_azim)) {
#     y <- c(y, sum(Avalue[[m]][(3*j-2):(3*j)]))
#   }
#   y <- y / y[match(32.5, stim_azim)]
#   # yrange <- range(y)
#   df <- cbind(df,y)
# }
# colnames(df) <- c('azim','DNp01','DNp02','DNp11','DNp04')
# # write.csv(df, "azim_sweep.csv", row.names=F)
# 
# x <- rev(c(-25, -12.5, 0, 12.5, 25)) #reversed elev def
# df <- x
# for (m in 1:length(neu_target)) {
#   y <- c()
#   for (j in (length(stim_azim)+1):length(c(stim_azim,stim_elev))) {
#     y <- c(y, sum(Avalue[[m]][(3*j-2):(3*j)]))
#   }
#   y <- y / y[5]
#   df <- cbind(df,y)
# }
# colnames(df) <- c('elev','DNp01','DNp02','DNp11','DNp04')
# # write.csv(df, "elev_sweep.csv", row.names=F)

