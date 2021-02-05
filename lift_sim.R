## 2D lifting
lift2D = function(x, y, f, nleft, stage, keepnbrs, Lmat, rw, method){
  require(deldir)
  require(tripack)
  require(rgeos)
  
  n = length(x)
  nlift = n - nleft           # number of points to be lifted
  poi = numeric(nlift)        # point of interest
  area.poi = numeric(nlift)   # area of poi
  dtl = numeric(nlift)        # detail coefficient calculated for poi
  id.poi = numeric(nlift)     # index of point lifted
  plot.stage = list(data.frame()) # optional vector of stages wanted to be plotted 
  nbrs = list(nlift)          # list of neighbours
  z = 1:n                     # bookkeeping
  L = diag(n)                 # lifting matrix initialized as an identity matrix
  
  lift.df = data.frame(x, y, z, f) # working data frame
  lift.df.full = lift.df           # back up the full data
  
  
  if (sum(method %in% c('convex', 'rectangle', 'double', 'base', 'augm')) == 0){
    warning('Weights are incorrectly specified')
  }
  
  if (method == 'rectangle'){
    cat('Starting lifting for rectangle boundary', date(), '\n')
    for (i in 1:nlift) {
      if (i %% floor(n/10) == 0) cat('Starting lifting stage', i, date(), '\n')
      ## Step 1 -- Select the cell with smallest area
      
      tes = deldir(lift.df, rw = rw, digits = 10) # first tessellation
      if(sum(tes$summary$dir.area == min(tes$summary$dir.area)) == 1){
        poi[i] = which.min(tes$summary$dir.area)
        #cat('Same area detected at stage', i, '\n')
      }else{
        poi[i] = sample(which(tes$summary$dir.area == min(tes$summary$dir.area)), 1)
      }
      area.poi[i] =  min(tes$summary$dir.area)
      id.poi[i] = tes$summary$z[poi[i]]
      nbrs1 = unique(unlist(tes$dirsgs[(tes$dirsgs[,5] == poi[i])|
                                         (tes$dirsgs[,6] == poi[i]),][,5:6]))
      nbrs1 = sort(nbrs1[nbrs1 != poi[i]]) # set of neighbours of poi
      
      new.tes = deldir(lift.df[-poi[i], ], rw = rw, digits = 10)
      nbrs2 = which(new.tes$summary$z %in% tes$summary$z[nbrs1]) # determine the positions of neighbours
      # after removing poi
      
      if(keepnbrs) nbrs[[i]] = tes$summary$z[nbrs1]
      
      area.nbrs1 = tes$summary$dir.area[nbrs1]    # area of neighbours before lifting the poi
      area.nbrs2 = new.tes$summary$dir.area[nbrs2]# area of neighbours after lifting the poi
      dif = area.nbrs2 - area.nbrs1 # area change of neighbours 
      
      ## Step 2 -- Predict the function value of poi
      pred.wts = c(dif/area.poi[i])  # Prediction weights
      dtl[i] = lift.df$f[poi[i]] - sum(pred.wts*lift.df$f[nbrs1]) # detail coefficient
      
      ## Step 3 -- Update the function values of neighbours
      b = (area.poi[i]*area.nbrs2)/(sum(area.nbrs2^2))
      lift.df$f[nbrs1] = lift.df$f[nbrs1] + dtl[i]*b
      
      ## Generate the lifting matrix
      if (Lmat){
        neig = tes$summary$z[nbrs1] # neighbours 
        
        L.new_star = diag(n) 
        L.new_star[id.poi[i], neig] = - pred.wts # pred. weights are placed into the identity matrix (predict)
        
        L.new = diag(n)
        for (k in 1:length(neig)) {L.new[neig[k], id.poi[i]] = b[k]} # b coefficients assigned assigned (update)
        
        L = L.new%*%L.new_star%*%L 
      }
      
      ## Step 4 -- Remove poi and update the dataset
      lift.df = lift.df[-poi[i], ]
      tes = new.tes
      
      if (i %in% stage){
        plot.stage[[which(stage == i)]] = lift.df
      }
    }## End of loop
  }## end of if statement for rectangle
  
  if (method == 'double'){
    cat('Starting lifting for area doubling', date(), '\n')
    
    tes = deldir(lift.df, rw = rw, digits = 10) # first tessellation
    edge = which(tes$summary$nbpt>0)
    area1 = tes$summary$dir.area
    area11 = area1
    area11[edge] = area1[edge]*2
    
    for (i in 1:nlift) {
      if (i %% floor(n/10) == 0) cat('Starting lifting stage', i, date(), '\n')
      ## Step 1 -- Select the cell with smallest area
      if(sum(tes$summary$dir.area == min(tes$summary$dir.area))>1){
        poi[i] = sample(which(tes$summary$dir.area == min(tes$summary$dir.area)), 1)
      }else{
        poi[i] = which.min(tes$summary$dir.area)
      }
      id.poi[i] = tes$summary$z[poi[i]]
      area.poi[i] =  area11[poi[i]]
      nbrs1 = unique(unlist(tes$dirsgs[(tes$dirsgs[,5] == poi[i])|(tes$dirsgs[,6] == poi[i]),][,5:6]))
      nbrs1 = sort(nbrs1[nbrs1 != poi[i]]) # set of neighbours of poi
      
      new.tes = deldir(lift.df[-poi[i], ], rw = rw, digits = 10)
      nbrs2 = which(new.tes$summary$z %in% tes$summary$z[nbrs1]) # determine the positions of neighbours
      # after removing poi
      
      if(keepnbrs) nbrs[[i]] = tes$summary$z[nbrs1]
      
      area.nbrs1 = area1[nbrs1]     # area of neighbours before lifting the poi
      area2  = new.tes$summary$dir.area
      area22 = area11[-poi[i]]
      
      if(id.poi[i] %in% edge) {
        area22[which(new.tes$summary$z %in% edge)] = area11[which(tes$summary$z %in% edge[-which(edge == id.poi[i])])]}
      if((id.poi[i] %in% edge) == 0) {
        area22[which(new.tes$summary$z %in% edge)] = area11[which(tes$summary$z %in% edge)]}
      
      area.nbrs2 = area2[nbrs2]
      dif = area.nbrs2 - area.nbrs1 # area change of neighbours 
      
      
      ## Step 2 -- Predict the function value of poi
      pred.wts = c(dif/tes$summary$dir.area[poi[i]])  # Prediction weights
      area.nbrs22 = area22[nbrs2] + pred.wts*area.poi[i] # area of neighbours after lifting the poi
      area22[nbrs2] = area.nbrs22
      
      dtl[i] = lift.df$f[poi[i]] - sum(pred.wts*lift.df$f[nbrs1]) # detail coefficient
      
      ## Step 3 -- Update the function values of neighbours
      b = (area.poi[i]*area.nbrs22)/(sum(area.nbrs22^2))
      lift.df$f[nbrs1] = lift.df$f[nbrs1] + dtl[i]*b
      
      
      
      ## Generate the lifting matrix
      if (Lmat){
        neig = tes$summary$z[nbrs1] # neighbours 
        
        L.new_star = diag(n) 
        L.new_star[id.poi[i], neig] = - pred.wts # pred. weights are placed into the identity matrix (predict)
        
        L.new = diag(n)
        for (k in 1:length(neig)) {L.new[neig[k], id.poi[i]] = b[k]} # b coefficients assigned assigned (update)
        
        L = L.new%*%L.new_star%*%L 
      }
      
      ## Step 4 -- Remove poi and update the dataset
      lift.df = lift.df[-poi[i], ]
      tes = new.tes
      area1 = area2
      area11 = area22
      
      if (i %in% stage){
        plot.stage[[which(stage == i)]] = lift.df
      }
    }## End of loop
    
  } ## End of if statement for double
  
  if (method == 'base'){
    cat('Starting lifting for base ensemble prediction', date(), '\n')
    
    #require(mgcv)
    
    prd.area = pred.base
    
    tes = deldir(lift.df, rw = rw, digits = 10) # first tessellation
    area1 = tes$summary$dir.area
    area11 = prd.area
    
    
    for (i in 1:nlift) {
      if (i %% floor(n/10) == 0) cat('Starting lifting stage', i, date(), '\n')
      ## Step 1 -- Select the cell with smallest area
      if(sum(area11 == min(area11))>1){
        poi[i] = sample(which(area11 == min(area11)), 1)
      }else{
        poi[i] = which.min(area11)
      }
      id.poi[i] = tes$summary$z[poi[i]]
      area.poi[i] =  area11[poi[i]]
      nbrs1 = unique(unlist(tes$dirsgs[(tes$dirsgs[,5] == poi[i])|(tes$dirsgs[,6] == poi[i]),][,5:6]))
      nbrs1 = sort(nbrs1[nbrs1 != poi[i]]) # set of neighbours of poi
      
      new.tes = deldir(lift.df[-poi[i], ], rw = rw, digits = 10)
      nbrs2 = which(new.tes$summary$z %in% tes$summary$z[nbrs1]) # determine the positions of neighbours
      # after removing poi
      
      if(keepnbrs) nbrs[[i]] = tes$summary$z[nbrs1]
      
      area.nbrs1 = area1[nbrs1]     # area of neighbours before lifting the poi
      area2  = new.tes$summary$dir.area
      area22 = area11[-poi[i]]
      
      # if(id.poi[i] %in% ne) {
      #   area22[which(new.tes$summary$z %in% ne)] = area11[which(tes$summary$z %in% ne[-which(ne == id.poi[i])])]}
      # if((id.poi[i] %in% ne) == 0) {
      #   area22[which(new.tes$summary$z %in% ne)] = area11[which(tes$summary$z %in% ne)]}
      
      area.nbrs2 = area2[nbrs2]
      dif = area.nbrs2 - area.nbrs1 # area change of neighbours 
      
      
      ## Step 2 -- Predict the function value of poi
      pred.wts = c(dif/area.poi[i])  # Prediction weights
      area.nbrs22 = area22[nbrs2] + pred.wts*area.poi[i] # area of neighbours after lifting the poi
      area22[nbrs2] = area.nbrs22
      
      dtl[i] = lift.df$f[poi[i]] - sum(pred.wts*lift.df$f[nbrs1]) # detail coefficient
      
      ## Step 3 -- Update the function values of neighbours
      b = (area.poi[i]*area.nbrs22)/(sum(area.nbrs22^2))
      lift.df$f[nbrs1] = lift.df$f[nbrs1] + dtl[i]*b
      
      
      
      ## Generate the lifting matrix
      if (Lmat){
        neig = tes$summary$z[nbrs1] # neighbours 
        
        L.new_star = diag(n) 
        L.new_star[id.poi[i], neig] = - pred.wts # pred. weights are placed into the identity matrix (predict)
        
        L.new = diag(n)
        for (k in 1:length(neig)) {L.new[neig[k], id.poi[i]] = b[k]} # b coefficients assigned assigned (update)
        
        L = L.new%*%L.new_star%*%L 
      }
      
      ## Step 4 -- Remove poi and update the dataset
      lift.df = lift.df[-poi[i], ]
      tes = new.tes
      area1 = area2
      area11 = area22
      
      if (i %in% stage){
        plot.stage[[which(stage == i)]] = lift.df
      }
    }## End of loop
    
  }## end of if statement for base
  
  if (method == 'augm'){
    cat('Starting lifting for augmented ensemble prediction', date(), '\n')
    #require(mgcv)
    
    prd.area = pred.augm
    
    tes = deldir(lift.df, rw = rw, digits = 10) # first tessellation
    area1 = tes$summary$dir.area
    area11 = abs(prd.area)
    
    
    for (i in 1:nlift) {
      if (i %% floor(n/10) == 0) cat('Starting lifting stage', i, date(), '\n')
      ## Step 1 -- Select the cell with smallest area
      if(sum(area11 == min(area11))>1){
        poi[i] = sample(which(area11 == min(area11)), 1)
      }else{
        poi[i] = which.min(area11)
      }
      id.poi[i] = tes$summary$z[poi[i]]
      area.poi[i] =  area11[poi[i]]
      nbrs1 = unique(unlist(tes$dirsgs[(tes$dirsgs[,5] == poi[i])|(tes$dirsgs[,6] == poi[i]),][,5:6]))
      nbrs1 = sort(nbrs1[nbrs1 != poi[i]]) # set of neighbours of poi
      
      new.tes = deldir(lift.df[-poi[i], ], rw = rw, digits = 10)
      nbrs2 = which(new.tes$summary$z %in% tes$summary$z[nbrs1]) # determine the positions of neighbours
      # after removing poi
      
      if(keepnbrs) nbrs[[i]] = tes$summary$z[nbrs1]
      
      area.nbrs1 = area1[nbrs1]     # area of neighbours before lifting the poi
      area2  = new.tes$summary$dir.area
      area22 = area11[-poi[i]]
      
      # if(id.poi[i] %in% ne) {
      #   area22[which(new.tes$summary$z %in% ne)] = area11[which(tes$summary$z %in% ne[-which(ne == id.poi[i])])]}
      # if((id.poi[i] %in% ne) == 0) {
      #   area22[which(new.tes$summary$z %in% ne)] = area11[which(tes$summary$z %in% ne)]}
      
      area.nbrs2 = area2[nbrs2]
      dif = area.nbrs2 - area.nbrs1 # area change of neighbours 
      
      
      ## Step 2 -- Predict the function value of poi
      pred.wts = c(dif/area.poi[i])  # Prediction weights
      area.nbrs22 = area22[nbrs2] + pred.wts*area.poi[i] # area of neighbours after lifting the poi
      area22[nbrs2] = area.nbrs22
      
      dtl[i] = lift.df$f[poi[i]] - sum(pred.wts*lift.df$f[nbrs1]) # detail coefficient
      
      ## Step 3 -- Update the function values of neighbours
      b = (area.poi[i]*area.nbrs22)/(sum(area.nbrs22^2))
      lift.df$f[nbrs1] = lift.df$f[nbrs1] + dtl[i]*b
      
      
      ## Generate the lifting matrix
      if (Lmat){
        neig = tes$summary$z[nbrs1] # neighbours 
        
        L.new_star = diag(n) 
        L.new_star[id.poi[i], neig] = - pred.wts # pred. weights are placed into the identity matrix (predict)
        
        L.new = diag(n)
        for (k in 1:length(neig)) {L.new[neig[k], id.poi[i]] = b[k]} # b coefficients assigned assigned (update)
        
        L = L.new%*%L.new_star%*%L 
      }
      
      ## Step 4 -- Remove poi and update the dataset
      lift.df = lift.df[-poi[i], ]
      tes = new.tes
      area1 = area2
      area11 = area22
      
      if (i %in% stage){
        plot.stage[[which(stage == i)]] = lift.df
      }
    }## End of loop
    
  }## end of if statement for augm
  
  if (method == 'convex'){
    cat('Starting lifting for convex hull boundary', date(), '\n')
    
    bnd1 = rw[2] - rw[1]
    bnd2 = rw[4] - rw[3]
    rw = c(rw[1]-(bnd1/2), rw[2]+(bnd1/2), rw[3]-(bnd2/2), rw[4]+(bnd2/2))
    tes = deldir(lift.df, rw = rw, digits = 10) # first tessellation
    chull1 = convex.hull(tri.mesh(lift.df.full$x, lift.df.full$y))
    poly1 = Polygon(cbind(chull1$x,chull1$y))
    p1 = SpatialPolygons(list(Polygons(list(poly1), "p1")))
    
    for (i in 1:nlift) {
      if (i %% floor(n/10) == 0) cat('Starting lifting stage', i, date(), '\n')
      ind1 = tes$dirsgs$ind1[which(in.convex.hull(tri.mesh(lift.df$x, lift.df$y), 
                                                  tes$dirsgs[,1], tes$dirsgs[,2])==F)]
      ind2 = tes$dirsgs$ind2[which(in.convex.hull(tri.mesh(lift.df$x, lift.df$y), 
                                                  tes$dirsgs[,3], tes$dirsgs[,4])==F)]
      ind = unique(c(ind1, ind2))
      
      area = tes$summary$dir.area
      subset = list()
      
      for (k in 1:length(ind)) {
        tmp = tes$dirsgs[which((tes$dirsgs$ind1 == ind[k])|(tes$dirsgs$ind2 == ind[k])), 1:4]
        subset[[k]] = tmp
        colnames(subset[[k]]) = c('x','y','x','y')
        rownames(subset[[k]]) = NULL
        subset[[k]] = rbind(subset[[k]][,1:2], subset[[k]][,3:4])
        subset[[k]] = subset[[k]][which(!duplicated(subset[[k]])),]
        
        # cells as SP class
        chull2 = convex.hull(tri.mesh(subset[[k]][,1], subset[[k]][,2]))
        poly2 = Polygon(cbind(chull2$x, chull2$y))
        p2 = SpatialPolygons(list(Polygons(list(poly2), "p2")))
        
        # intersect the cell with convex hull
        res = gIntersection(p1, p2)
        
        # area for convex hull
        area[ind[k]] = unlist(sapply(slot(res, "polygons"), function(p) sapply(slot(p, "Polygons"), slot, "area")))
      }
      
      poi[i] = which.min(area)
      area.poi[i] =  min(area)
      
      #nbrs1 = neighbours(tri.mesh(lift.df$x,lift.df$y))[[poi[i]]]
      nbrs1 = unique(unlist(tes$dirsgs[(tes$dirsgs[,5] == poi[i])|
                                         (tes$dirsgs[,6] == poi[i]),][,5:6]))
      nbrs1 = sort(nbrs1[nbrs1 != poi[i]]) # set of neighbours of poi
      
      new.tes = deldir(lift.df[-poi[i], ], rw = rw, digits = 10)
      
      nbrs2 = which(new.tes$summary$z %in% tes$summary$z[nbrs1]) # determine the positions of neighbours
      # after removing poi
      
      area2 = as.numeric()
      for (k in 1:length(nbrs2)) {
        subset[[k]] = new.tes$dirsgs[which((new.tes$dirsgs$ind1 == nbrs2[k])|(new.tes$dirsgs$ind2 == nbrs2[k])), 1:4]
        colnames(subset[[k]]) = c('x','y','x','y')
        rownames(subset[[k]]) = NULL
        subset[[k]] = rbind(subset[[k]][,1:2], subset[[k]][,3:4])
        subset[[k]] = subset[[k]][which(!duplicated(subset[[k]])),]
        
        # cells as SP class
        chull2 = convex.hull(tri.mesh(subset[[k]][,1], subset[[k]][,2]))
        poly2 = Polygon(cbind(chull2$x, chull2$y))
        p2 = SpatialPolygons(list(Polygons(list(poly2), "p2")))
        
        # intersect the cell with convex hull
        res = gIntersection(p1, p2)
        
        # area for convex hull
        area2[k] = unlist(sapply(slot(res, "polygons"), function(p) sapply(slot(p, "Polygons"), slot, "area")))
      }
      
      
      if(keepnbrs) nbrs[[i]] = tes$summary$z[nbrs1]
      
      area.nbrs1 = area[nbrs1]    # area of neighbours before lifting the poi
      area.nbrs2 = area2# area of neighbours after lifting the poi
      dif = area.nbrs2 - area.nbrs1 # area change of neighbours 
      
      ## Step 2 -- Predict the function value of poi
      pred.wts = c(dif/area.poi[i])  # Prediction weights
      dtl[i] = lift.df$f[poi[i]] - sum(pred.wts*lift.df$f[nbrs1]) # detail coefficient
      
      ## Step 3 -- Update the function values of neighbours
      b = (area.poi[i]*area.nbrs2)/(sum(area.nbrs2^2))
      lift.df$f[nbrs1] = lift.df$f[nbrs1] + dtl[i]*b
      
      id.poi[i] = tes$summary$z[poi[i]]
      
      ## Generate the lifting matrix
      if (Lmat){
        neig = tes$summary$z[nbrs1] # neighbours 
        
        L.new_star = diag(n) 
        L.new_star[id.poi[i], neig] = - pred.wts # pred. weights are placed into the identity matrix (predict)
        
        L.new = diag(n)
        for (k in 1:length(neig)) {L.new[neig[k], id.poi[i]] = b[k]} # b coefficients assigned (update)
        
        L = L.new%*%L.new_star%*%L 
      }
      
      ## Step 4 -- Remove poi and update the dataset
      lift.df = lift.df[-poi[i], ]
      tes = new.tes
      
      if (i %in% stage){
        plot.stage[[which(stage == i)]] = lift.df
      }
    }## End of loop
  }## end of if statement for convex
  
  ## Combine and sort detail coefficients with remaining f
  d_coef = data.frame(c(dtl, lift.df$f), c(id.poi, lift.df$z))
  d_coef = d_coef[order(d_coef[,2]), ]
  d_coef = c(d_coef[,1])
  
  return(list(detail = dtl, det_ordered = d_coef, remaining = lift.df, order = id.poi, 
              area = area.poi, full = lift.df.full,
              plot.stage = plot.stage, nbrs = nbrs, Lmat = L, 
              if(sum((method == 'base')|(method == 'augm'))>0){prd.area = prd.area}))
}