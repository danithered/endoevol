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
