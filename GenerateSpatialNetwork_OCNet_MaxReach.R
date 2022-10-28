



library(OCNet)


setwd("~/getreal-model/data_GETREAL/river_nws")

sizeX = 32
sizeY = 32

for (approx_number_nodes in c(1)){ #500
  for (contador in 1:10){
    print(c(approx_number_nodes, contador))
    file_nodes = paste0("river_", approx_number_nodes, "_nodes_", contador, "_repl_NODES.csv")
    file_edges = paste0("river_", approx_number_nodes, "_nodes_", contador, "_repl_EDGES.csv")
    
    
    flag = FALSE
    while (flag == FALSE){
      OCN <- create_OCN(sizeX,sizeY, outletSide = "S", outletPos = sample(1:sizeX,1))
      OCN <- landscape_OCN(OCN)
      
      thr <- find_area_threshold_OCN(OCN, maxReachLength = floor(sizeX/5))
      indThr <- which(abs(thr$nNodesAG - approx_number_nodes) == min(abs(thr$nNodesAG - approx_number_nodes)))
      indThr <- max(indThr) # pick the last ind_thr that satisfies the condition above
      thrA <- thr$thrValues[indThr] # corresponding threshold area
      
      
      OCN <- aggregate_OCN(landscape_OCN(OCN), thrA = thrA, maxReachLength = floor(sizeX/5))
      
      OCN <- rivergeometry_OCN(OCN)
      
      if (OCN$AG$nNodes == approx_number_nodes) {
        flag = TRUE
      }  
      
    }
    
    nodes <- 1:OCN$AG$nNodes
    posoutlet = nodes[OCN$AG$downNode==0]
    edges <- as.matrix(OCN$AG$W)
    edges <- edges + t(edges)
    Xnode <- OCN$AG$X
    Ynode <- OCN$AG$Y
    connectedto <- OCN$AG$downNode
    
    area <- OCN$AG$A
    for (node1 in seq(from=OCN$AG$nNodes, to=1, by=-1)){
      for (node2 in seq(from=OCN$AG$nNodes, to=1, by=-1)){
        if (edges[node1, node2]>0 & OCN$AG$A[node1]>=OCN$AG$A[node2]){
          area[node1] = area[node1] - OCN$AG$A[node2]
        }
      }
    }
    
    width <- OCN$AG$width
    depth <- OCN$AG$depth
    vel   <- OCN$AG$velocity
    lengs <- OCN$AG$leng
    
    
    df <- data.frame(nodes = nodes,
                     posoutlet = posoutlet,
                     Xnode = Xnode,
                     Ynode = Ynode,
                     width = width,
                     depth = depth,
                     lengs = lengs,
                     vel   = vel
    )
  
    
    write.table(df, file_nodes, append = FALSE, sep = ";", dec = ".",
                row.names = FALSE, col.names = TRUE)
    
    
    
    df_edges <- data.frame(
      edges = edges
    )
    write.table(df_edges, file_edges, append = FALSE, sep = ";", dec = ".",
                row.names = FALSE, col.names = TRUE)
    
  
  }
}




