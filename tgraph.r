require(igraph)
rectKasteleyne <- function(m,n){ # m*n must be even
    x <- rep(1:m, n)
    y <- rep(1:n, each=m)

    bw <- (x + y)%%2
    bn <- paste(x,y,sep=",")[bw==0]
    wn <- paste(x,y,sep=",")[bw==1]
    K <- matrix(0, m*n/2, m*n/2)
    row.names(K) <- bn
    colnames(K) <- wn

    # determine white vertices on the boundary

    xb <- c()
    yb <- c()
    x <- 1
    y <- 1
    while (x < m){
        if ((x + y)%%2 == 1){
            xb <- c(xb,x)
            yb <- c(yb,y)
        }
        x <- x+1
    }
    while (y < n){
        if ((x + y)%%2 == 1){
            xb <- c(xb,x)
            yb <- c(yb,y)
        }
        y <- y+1
    }
    while (x > 1){
        if ((x + y)%%2 == 1){
            xb <- c(xb,x)
            yb <- c(yb,y)
        }
        x <- x-1
    }
    while (y > 1){
        if ((x + y)%%2 == 1){
            xb <- c(xb,x)
            yb <- c(yb,y)
        }
        y <- y-1
    }

    whiteboundary <- paste(xb, yb, sep=",")
    
    lab2coord <- function(s){
        as.numeric(strsplit(s,",")[[1]])
    }
    isneighbor <- function(s,t){
        d <- lab2coord(s) - lab2coord(t)
        if (d[1] == 0 & d[2] == 1){
            return(-1)
        } else if (sum(abs(d))==1){
            return(1)
        } else {
            return(0)
        }
    }
    for (b in bn){
        for (w in wn){
            K[b,w]<-isneighbor(b,w)
        }
    }
    return(list(K=K,xboundary=xb, yboundary=yb))
}

gauge <- function(m,n){
    KK <- rectKasteleyne(m,n)
    K <- KK$K
    v <- rep(0,dim(K)[1])
    v[1] <- 1
    Kd <- K %*% diag(as.vector(solve(K) %*% matrix(v,ncol=1)))
    colnames(Kd) <- colnames(K)
    rownames(Kd) <- rownames(K)
    vw <- rep(0,dim(K)[1])
    names(vw) <- colnames(K)
    blab <- paste(KK$xboundary, KK$yboundary, sep=",")
    vw[blab] <- diff(c(1+1i,KK$xboundary + KK$yboundary*1i))
    Kt <- diag(as.vector(matrix(vw,nrow=1) %*% solve(Kd)) ) %*% Kd
    colnames(Kt) <- colnames(K)
    rownames(Kt) <- rownames(K)
    return(list(m=m, n=n,K=Kt, xboundary=KK$xboundary, yboundary=KK$yboundary))
}

Kintegral <- function(KK){
    m <- KK$m
    n <- KK$n
    K <- KK$K
    
    g <- graph.empty()
    facex <- rep(1:(m-1), n-1)
    facey <- rep(1:(n-1), each=m-1)
    facename <- paste(facex,facey,sep=",")
    g <- add.vertices(g,length(facename),attr=list(name=facename,cx=facex,cy=facey))

    #for (v in V(g)){
    #    for (w in V(g)){
    #        if (abs(V(g)[v]$cx-V(g)[w]$cx) + abs(V(g)[v]$cy-V(g)[w]$cy) ==1 & (V(g)[v]$cx + V(g)[v]$cy)%%2 == 0 ){
    #            if (V(g)[v]$cx == V(g)[w]$cx){
    #                g <- add.edges(g,c(w,v),
    #                               weight=K[paste(V(g)[w]$cx,V(g)[w]$cy,sep=","),
    #                                   paste(V(g)[v]$cx,V(g)[v]$cy,sep=",")]
    #                               )
    #            }else{
    #                g <- add.edges(g,c(v,w),
    #                               weight=K[paste(V(g)[v]$cx,V(g)[v]$cy,sep=","),
    #                                   paste(V(g)[w]$cx,V(g)[w]$cy,sep=",")]
    #                               )
    #            }
    #        }
    #    }
    #}
    for (i in 1:length(facex)){
        for (j in 1:length(facex)){
            x <- facex[i]
            y <- facey[i]
            xd <- facex[j]
            yd <- facey[j]
            if (abs(x-xd)+abs(y-yd) == 1 & (x+y)%%2==0){

                if (x < xd){
                    g <- add.edges(g,c(i,j),
                                   weight=K[paste(x+1,y+1,sep=","),paste(x+1,y,sep=",")]
                                   )
                }
                if (x > xd){
                    g <- add.edges(g,c(i,j),
                                   weight=K[paste(x,y,sep=","),paste(x,y+1,sep=",")]
                                   )
                }
                if (y < yd){
                    g <- add.edges(g,c(j,i),
                                   weight=K[paste(x+1,y+1,sep=","),paste(x,y+1,sep=",")]
                                   )
                }
                if (y> yd){
                    g <- add.edges(g,c(j,i),
                                   weight=K[paste(x,y,sep=","),paste(x+1,y,sep=",")]
                                   )
                }
            }
        }
    }

    x <- 1
    y <- 1
    numout <- 0
    g <- add.vertices(g,1,name=paste("o",numout,sep=""))
    while (x < m){
        if ((x + y)%%2 == 0){
            numout <- numout+1
            g <- add.edges(g,c(length(V(g)), V(g)[paste(x,y,sep=",")]),
                           weight=K[paste(x,y,sep=","),paste(x+1,y,sep=",")]
                           )
        }else{
            g <- add.vertices(g,1,attr=list(name=paste("o",numout,sep="")))
            g <- add.edges(g,c(V(g)[paste(x,y,sep=",")], V(g)[length(V(g))]),
                           weight=K[paste(x+1,y,sep=","), paste(x,y,sep=",")]
                            )
        }
        x <- x+1
    }
    while (y < n){
        if ((x + y)%%2 == 0){
            numout <- numout+1
            g <- add.edges(g,c(V(g)[length(V(g))],V(g)[paste(x-1,y,sep=",")]))
        }else{
            g <- add.vertices(g,1,attr=list(name=paste("o",numout,sep="")))
            g <- add.edges(g,c(V(g)[paste(x-1,y,sep=",")], V(g)[length(V(g))]))
        }
        y <- y+1
    }
    while (x > 1){
        if ((x + y)%%2 == 0){
            numout <- numout+1
            g <- add.edges(g,c(V(g)[length(V(g))],V(g)[paste(x-1,y-1,sep=",")]))
        }else{
            g <- add.vertices(g,1,attr=list(name=paste("o",numout,sep="")))
            g <- add.edges(g,c(V(g)[paste(x-1,y-1,sep=",")], V(g)[length(V(g))]))
        }
        x <- x-1
    }
    while (y > 1){
        if ((x + y)%%2 == 0){
            numout <- numout+1
            g <- add.edges(g,c(V(g)[length(V(g))],V(g)[paste(x,y-1,sep=",")]))
        }else{
            g <- add.vertices(g,1,attr=list(name=paste("o",numout,sep="")))
            g <- add.edges(g,c(V(g)[paste(x,y-1,sep=",")], V(g)[length(V(g))]))
        }
        y <- y-1
    }
        
    return(g)
    
}


