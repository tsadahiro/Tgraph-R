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

    xb <- c()
    yb <- c()
    x <- 1
    y <- 1
    while (x < m){
        if ((x + y)%%2 == 0){
            xb <- c(xb,x)
            yb <- c(yb,y)
        }
        x <- x+1
    }
    while (y < n){
        if ((x + y)%%2 == 0){
            xb <- c(xb,x)
            yb <- c(yb,y)
        }
        y <- y+1
    }
    while (x > 1){
        if ((x + y)%%2 == 0){
            xb <- c(xb,x)
            yb <- c(yb,y)
        }
        x <- x-1
    }
    while (y > 1){
        if ((x + y)%%2 == 0){
            xb <- c(xb,x)
            yb <- c(yb,y)
        }
        y <- y-1
    }
        
    return(cbind(xb,yb))
    
}
