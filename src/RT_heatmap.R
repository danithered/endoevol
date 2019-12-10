library(lattice)

x=0
y=0.2

( (1-x) *(1-exp(-0.05 * n)) / 2 ) * ( y *(1-exp(-0.05 * n)) )

n=40
Y= X = seq(0,n,1)
z = expand.grid(x=X, y=Y)


Z = matrix(
  ( (1-z$x/n) *(1-exp(-0.05 * n)) / 2 ) * ( z$y/n *(1-exp(-0.05 * n)) )
  , ncol=length(X)
  )

W = matrix(
  ( (1-z$x/n) *(1-exp(-0.05 * n)) / ((1-z$x/n) *(1-exp(-0.05 * n)) + (1-z$y/n) *(1-exp(-0.05 * n)) + 0.1) ) * ( z$y/n *(1-exp(-0.05 * n)) )
  , ncol=length(X)
)

image(X, Y, W
      , xlab="R in x"
      , ylab="R in y"
      , main="x is replicated by y"
      , sub=paste("number of bases in both x and y is", n)
      , col = hcl.colors(25, "YlOrRd", rev = TRUE)
      )
contour(X, Y, W, add = TRUE, drawlabels = T)


L=1:50
R=1:50
z = expand.grid(l=L, r=R)

u=  z$r/z$l *(1-exp(-0.01 * z$l^4))
u=  z$r/z$l *(1- exp(0.001*-z$l^2) )

alfa <- 5
beta <- 0.1
c <- 1
#u=  (alfa^z$r) / (alfa^z$r + alfa^(z$l-z$r) ) *(1- exp(-beta * z$l^c) )

u=  (z$r^alfa) / (z$r^alfa + (z$l-z$r)^alfa ) *(1- exp(-beta * z$l)^c )

#u=  1- 1/ (1+ z$r^(0.05*z$l) )
#u=  z$r/z$l *(1 - 1/ exp(0.01 * z$l) )


u[z$r>z$l] <- 0
     
Z = matrix(
  #( z$r/z$l *(1-exp(-0.1 * z$l^2))  )  
  u
  , ncol=length(L)
)

image(L, R, Z
      , xlab="length"
      , ylab="R"
      #, main="x is replicated by y"
      , sub=paste("alpha", alfa, ", beta", beta, ", c", c)
      , col = hcl.colors(25, "YlOrRd", rev = TRUE)
      , asp=1
)
#abline(a=-20, b=1)
#abline(h=10)
contour(L, R, Z, add = TRUE, drawlabels = T)

points(c(40, 20, 10), c(20, 10, 5), cex=2)
