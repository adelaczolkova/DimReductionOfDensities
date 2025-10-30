
## US HOUSING DATA

# Analysis of bivariate densities of log-transformed house price and size in the USA

# Loading density kernel estimations
load("housing_densities.RData")

# Loading state names, abbreviations and regions (50 US states, District of Columbia and Puerto Rico)
st.abb = read.table("states.txt", header = T)
st.abb$region = factor(st.abb$region, levels = c("North Central","Northeast","South","West","Puerto Rico"))


## Initial settings ------------------------------------------------------------

# ranges
range_x = c(7.813996, 18.945409)
range_y = c(5.062595, 11.512915)

# grids
grid_x = seq(range_x[1], range_x[2], l=30)
grid_y = seq(range_y[1], range_y[2], l=30)

# colors
pos_cols = colorRampPalette(c("lightyellow","orange","red4"), space = "rgb")(100)
cols = c(colorRampPalette(c("blue4","turquoise","white"), space = "rgb")(50),
         colorRampPalette(c("white","orange","red4"), space = "rgb")(50))


## Clr transformation, centering and orthogonal decomposition ------------------

clrd = c_clrd = c_clrx = c_clry = c_clri = array(dim = c(length(grid_x), length(grid_y), nrow(st.abb)), 
                                                 dimnames = list(NULL, NULL, st.abb$state))
# clr transformation
for(s in st.abb$state){
  clrd[,,s] = log(dens[,,s]) - sum(log(dens[,,s]))/(length(grid_x)*length(grid_y))
}

# mean clr density
clrm = apply(clrd, 1:2, mean)

for(s in st.abb$state){
  # centering
  c_clrd[,,s] = clrd[,,s] - clrm
  
  # orthogonal decomposition
  c_clrx[,,s] = matrix(apply(c_clrd[,,s], 1, mean), nrow = length(grid_x), ncol = length(grid_y))
  c_clry[,,s] = matrix(apply(c_clrd[,,s], 2, mean), nrow = length(grid_x), ncol = length(grid_y), byrow = T)
  c_clri[,,s] = c_clrd[,,s] - c_clrx[,,s] - c_clry[,,s]
}


## FPCA for multivariate densities ---------------------------------------------

# Vectorizing the clr densities
# - vectorized centered clr density in each row of the matrix, rows corresponds to the states
vc_clrd = t(apply(c_clrd, 3, as.vector))

dim(vc_clrd) # (number of states) x (number of density function values)

# Computing scores and loadings using the singular value decomposition (SVD)
res = svd(vc_clrd)

# scores
S = res$u %*% diag(res$d)  

# loadings (eigenfunctions in the columns)
V = array(as.vector(res$v), dim = c(length(grid_x), length(grid_y), ncol(res$v))) 


## Decomposition of loadings and scores ----------------------------------------

# number of components for which we want to decompose loadings
ncomp = dim(V)[3]  

# Decomposing loadings (eigenfunctions)
Vx = Vy = Vi = array(dim = c(dim(V)[1:2],ncomp))
for(j in 1:ncomp){
  Vx[,,j] = matrix(apply(V[,,j], 1, mean), nrow = length(grid_x), ncol = length(grid_y))
  Vy[,,j] = matrix(apply(V[,,j], 2, mean), nrow = length(grid_x), ncol = length(grid_y), byrow = T)
  Vi[,,j] = V[,,j] - Vx[,,j] - Vy[,,j]
}

# Decomposing scores using vectorized density parts and vectorized patrs of loadings
Sx = t(apply(c_clrx, 3, as.vector)) %*% apply(Vx, 3, as.vector)
Sy = t(apply(c_clry, 3, as.vector)) %*% apply(Vy, 3, as.vector)
Si = t(apply(c_clri, 3, as.vector)) %*% apply(Vi, 3, as.vector)

# Visualization of the scores

layout(matrix(c(1,3,2,4,5,5), 2, 3), widths = c(2,2,0.8))
par(mar = c(4,4,3,1), mgp = c(2.5,0.9,0))
plot(Sx[,1],Sx[,2], pch=NA, main="x-marginal (price)", xlab="1st component", ylab="2nd component",
     panel.first = abline(h=0,v=0,col="grey"), 
     panel.last = text(Sx[,1],Sx[,2], st.abb[,2], col=as.numeric(st.abb$region), font=2, cex=0.9))
plot(Sy[,1],Sy[,2], pch=NA, main="y-marginal (size)", xlab="1st component", ylab="2nd component",
     panel.first = abline(h=0,v=0,col="grey"),
     panel.last = text(Sy[,1],Sy[,2], st.abb[,2], col=as.numeric(st.abb$region), font=2, cex=0.9))
plot(Si[,1],Si[,2], pch=NA, main="interaction", xlab="1st component", ylab="2nd component",
     panel.first = abline(h=0,v=0,col="grey"),
     panel.last = text(Si[,1],Si[,2], st.abb[,2], col=as.numeric(st.abb$region), font=2, cex=0.9))
plot(S[,1],S[,2], pch=NA, main="whole density", xlab="1st component", ylab="2nd component",
     panel.first = abline(h=0,v=0,col="grey"),
     panel.last = text(S[,1],S[,2], st.abb[,2], col=as.numeric(st.abb$region), font=2, cex=0.9))
par(mar = rep(0,4))
plot(0,0, axes=F, type="n", xlab="", ylab="")
legend("left", legend=c("North Central","Northeast","South","West","Puerto Rico (PR)"),
       col=1:5, pch=15, bty="n", cex=1.3, ncol=1, pt.cex=1.5)


## Explained variability -------------------------------------------------------

# squares of singular values = eigenvalues
Var = res$d^2  

# percents of explained variability by each component
var_perc = 100*Var/sum(Var)

par(mar = c(5,4,1,1), mfrow = c(1,1))
plot(var_perc, type = "b", pch = 20, xlab = "Components", ylab = "% of explained variance", ylim=c(0,30),
     panel.first = abline(h=0, col="darkgrey"),
     panel.last = text(c(4.5,4.5),
                       var_perc[1:2]+c(0.5,1.5),
                       c(paste(round(var_perc[1],2),"%"),paste(round(var_perc[2],2),"%"))))

# variability explained by the first two components
sum(var_perc[1:2])


## Decomposition of squared norms and variance ---------------------------------

# squared norm
norm2 = function(f) sum(f^2)/(length(grid_x)*length(grid_y))

norm_tab = cbind(apply(c_clrd, 3, norm2),
                 apply(c_clrx, 3, norm2),
                 apply(c_clry, 3, norm2),
                 apply(c_clri, 3, norm2))

# a number of states
N = nrow(st.abb)

par(mgp = c(2, 0.3, 0), mar = c(7,4,1,2), mfrow = c(1,1))
fields::image.plot(1:N, 1:4, norm_tab[,4:1], col = pos_cols,
                   breaks = quantile(as.vector(norm_tab), probs = seq(0, 1, 0.01)),
                   axes = F, xlab = "", ylab = "", axis.args = list(cex.axis = 0.8))
axis(2, at = 1:4, labels = rev(c("density","price (x)","size (y)","interaction")), 
     tick = F, las = 2, cex.axis = 0.8)
axis(1, at = 1:N-0.1, labels = rownames(norm_tab), tick = F, las = 2, cex.axis = 0.8)


# total variance
variance = colSums(norm_tab)/N

# percents of total variance contained in the orthogonal density parts
(perc = variance[-1]/variance[1]*100)
round(perc, 2)


# relative norms
rel_norm_tab = norm_tab[,-1]/norm_tab[,1]

par(mgp = c(2, 0.3, 0), mar = c(7,4,1,2), mfrow = c(1,1))
fields::image.plot(1:N, 1:3, rel_norm_tab[,3:1], col = pos_cols,
                   breaks = quantile(as.vector(rel_norm_tab), probs = seq(0,1,0.01)),
                   axes = F, xlab = "", ylab = "", axis.args = list(cex.axis = 0.8))
axis(2, at = 1:3, labels = rev(c("price (x)","size (y)","interaction")), 
     tick = F, las = 2, cex.axis = 0.8)
axis(1, at = 1:N-0.1, labels = rownames(rel_norm_tab), tick = F, las = 2, cex.axis = 0.8)


## Visualization of the loadings with precentage of explained variability ------

n2_1comp = c(norm2(V[,,1]), norm2(Vx[,,1]), norm2(Vy[,,1]), norm2(Vi[,,1]))
n2_2comp = c(norm2(V[,,2]), norm2(Vx[,,2]), norm2(Vy[,,2]), norm2(Vi[,,2]))

(perc1 = n2_1comp[2:4]/n2_1comp[1]*100)
(perc2 = n2_2comp[2:4]/n2_2comp[1]*100)

(var_perc1 = round(perc1/100*var_perc[1], 2))
(var_perc2 = round(perc2/100*var_perc[2], 2))

comp_range = range(c(Vx[,,1:2],Vy[,,1:2],Vi[,,1:2],V[,,1:2]))
comp_breaks = c(seq(comp_range[1], 0, length = 51), seq(0, comp_range[2], length = 51)[-1])

layout(matrix(c(0,10,10,
                11,1,2,
                11,3,4,
                11,5,6,
                11,7,8,
                0,9,9),3,6), widths = c(0.8,4,4,4,4,1), heights = c(0.6,4,4))
par(mar = c(4,4,2,1), mgp = c(2.5,0.9,0), cex.main = 1)

image(grid_x, grid_y, V[,,1], xlim = range_x, ylim=range_y, col=cols, breaks = comp_breaks, 
      xlab="log(price)", ylab="log(size)", main = paste0(round(var_perc[1], 2)," % (100 %)"), font.main = 1)
image(grid_x, grid_y, V[,,2], xlim = range_x, ylim=range_y, col=cols, breaks = comp_breaks,
      xlab="log(price)", ylab="log(size)", main = paste0(round(var_perc[2], 2)," % (100 %)"), font.main = 1)
#x
image(grid_x, grid_y, Vx[,,1], xlim = range_x, ylim=range_y, col=cols, breaks = comp_breaks, 
      xlab="log(price)", ylab="log(size)", main = paste0(var_perc1[1]," % (", round(perc1[1], 2)," %)"), font.main = 1)
image(grid_x, grid_y, Vx[,,2], xlim = range_x, ylim=range_y, col=cols, breaks = comp_breaks, 
      xlab="log(price)", ylab="log(size)", main = paste0(var_perc2[1]," % (", round(perc2[1], 2)," %)"), font.main = 1)
#y
image(grid_x, grid_y, Vy[,,1], xlim = range_x, ylim=range_y, col=cols, breaks = comp_breaks, 
      xlab="log(price)", ylab="log(size)", main = paste0(var_perc1[2]," % (", round(perc1[2], 2)," %)"), font.main = 1)
image(grid_x, grid_y, Vy[,,2], xlim = range_x, ylim=range_y, col=cols, breaks = comp_breaks, 
      xlab="log(price)", ylab="log(size)", main = paste0(var_perc2[2]," % (", round(perc2[2], 2)," %)"), font.main = 1)
#int
image(grid_x, grid_y, Vi[,,1], xlim = range_x, ylim=range_y, col=cols, breaks = comp_breaks, 
      xlab="log(price)", ylab="log(size)", main = paste0(var_perc1[3]," % (", round(perc1[3], 2)," %)"), font.main = 1)
image(grid_x, grid_y, Vi[,,2], xlim = range_x, ylim=range_y, col=cols, breaks = comp_breaks, 
      xlab="log(price)", ylab="log(size)", main = paste0(var_perc2[3]," % (", round(perc2[3], 2)," %)"), font.main = 1)

par(mar = c(10,0,8,3.5))
image(y = comp_breaks, z = matrix(comp_breaks,1), breaks = comp_breaks, col = cols, axes = F)
axis(4, las = 2)

par(mar=c(0,0,0,0))
plot.new()
text(0.4, 0.81, "1st component", srt=90, cex=1.2, font = 2)
text(0.4, 0.27, "2nd component", srt=90, cex=1.2, font = 2)

plot.new()
mtext(c("whole density","x-marginal (price)","y-marginal (size)","interaction"),
      at = c(0.12, 0.39, 0.66, 0.93), line = -2, font = 2, cex = 0.9)


dev.off()
