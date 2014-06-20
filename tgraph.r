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
    xb1 <- seq(2,m, by=2)
    if (m <= 2){
        xb2 <- c()
    }else{
        xb2 <- seq(m-1-((m+n)%%2),1,-2)
    }
    yb1 <- seq(2+(m+1)%%2,n,2)
    yb2 <- seq(n-1-(n+1)%%2,1,-2)
    xb <- c(xb1,rep(m,length(yb1)), xb2, rep(1,length(yb2)))
    yb <- c(rep(1,length(xb1)), yb1, rep(n,length(xb2)),yb2)
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

m <- 4
n <- 3
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

