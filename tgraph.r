K <- matrix(c(1,0,1, 1,1,-1, 0,1,1),ncol=3, byrow=T)
Kd <- K %*% diag(as.vector(solve(K) %*% matrix(c(1,0,0),ncol=1)))
Kt <- diag(as.vector(matrix(c(1,1i,-1),nrow=1) %*% solve(Kd)) ) %*% Kd

m <- 4
n <- 4
x <- rep(1:m, n)
y <- rep(1:m, each=n)
bw <- (x + y)%%2
plot(x,y,pch=20,col=bw+1)
bn <- paste(x,y,sep=",")[bw==0]
wn <- paste(x,y,sep=",")[bw==1]
K <- matrix(0, m*n/2, m*n/2)
row.names(K) <- bn
colnames(K) <- wn

lab2coord <- function(s){
      as.numeric(strsplit(s,",")[[1]])
  }
isneighbor <- function(s,t){
      d <-sum(abs(lab2coord(s) - lab2coord(t)))
        if (d == 1){
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
