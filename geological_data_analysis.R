
## CZECH GEOLOGICAL DATA

# Analysis of bivariate densities of Cu and Zn concentrations in the soil

# Loading district names and abbreviations
districts = read.table("districts_names.txt", header = T, na.strings = NULL)

# Loading density kernel estimations
dens = array(dim = c(60, 60, nrow(districts)), dimnames = list(NULL, NULL, districts$name2))
for(i in districts$name2){
  dens[,,i] = t(as.matrix(read.table(paste0("./KDE_values_imput/zmat_imput_",i,".txt"), header = F)))
}

## Initial settings ------------------------------------------------------------

# ranges
range_Cu = c(0.588,4.58)
range_Zn = c(1.459,5.663)

# grids
grid_Cu = seq(range_Cu[1], range_Cu[2], l=60)
grid_Zn = seq(range_Zn[1], range_Zn[2], l=60)

# colors
pos_cols = colorRampPalette(c("lightyellow","orange","red4"), space = "rgb")(100)
cols = c(colorRampPalette(c("blue4","turquoise","white"), space = "rgb")(50),
         colorRampPalette(c("white","orange","red4"), space = "rgb")(50))


## Clr transformation and orthogonal decomposition (in L^2_0) ------------------

clrd = c_clrd = c_clrx = c_clry = c_clri = array(dim = c(length(grid_Cu), length(grid_Zn), nrow(districts)), 
                                                 dimnames = list(NULL, NULL, districts$name2))
# clr transformation
for(s in districts$name2){
  clrd[,,s] = log(dens[,,s]) - sum(log(dens[,,s]))/(length(grid_Cu)*length(grid_Zn))
}

# mean clr density
clrm = apply(clrd, 1:2, mean)

for(s in districts$name2){
  # centering
  c_clrd[,,s] = clrd[,,s] - clrm
  
  # orthogonal decomposition
  c_clrx[,,s] = matrix(apply(c_clrd[,,s], 1, mean), nrow = length(grid_Cu), ncol = length(grid_Zn))
  c_clry[,,s] = matrix(apply(c_clrd[,,s], 2, mean), nrow = length(grid_Cu), ncol = length(grid_Zn), byrow = T)
  c_clri[,,s] = c_clrd[,,s] - c_clrx[,,s] - c_clry[,,s]
}


## FPCA for multivariate densities ---------------------------------------------

# Centering and vectorizing the clr densities
# - vectorized centered density in each row of the matrix, rows corresponds to the districts
vc_clrd = t(apply(c_clrd, 3, as.vector))

dim(vc_clrd) # (number of districts) x (number of density functional values)

# Computing scores and loadings using the singular value decomposition (SVD)
res = svd(vc_clrd)

# scores
S = res$u %*% diag(res$d)  

# loadings (eigenfunctions in the columns)
V = array(as.vector(res$v), dim = c(length(grid_Cu),length(grid_Zn),ncol(res$v))) 


## Decomposition of loadings and scores ----------------------------------------

# number of components for which we want to decompose loadings
ncomp = dim(V)[3]  

# Decomposing loadings (eigenfunctions)
Vx = Vy = Vi = array(dim = c(dim(V)[1:2],ncomp))
for(j in 1:ncomp){
  Vx[,,j] = matrix(apply(V[,,j], 1, mean), nrow = length(grid_Cu), ncol = length(grid_Zn))
  Vy[,,j] = matrix(apply(V[,,j], 2, mean), nrow = length(grid_Cu), ncol = length(grid_Zn), byrow = T)
  Vi[,,j] = V[,,j] - Vx[,,j] - Vy[,,j]
}

# Decomposing scores using vectorized density parts and vectorized patrs of loadings
Sx = t(apply(c_clrx, 3, as.vector)) %*% apply(Vx, 3, as.vector)
Sy = t(apply(c_clry, 3, as.vector)) %*% apply(Vy, 3, as.vector)
Si = t(apply(c_clri, 3, as.vector)) %*% apply(Vi, 3, as.vector)

# Visualization of the scores

# groups
groups = cbind(districts$abbr, rep(NA,nrow(districts)))
rownames(groups) = districts$name2

unrelated = c("CK","CR","TA","TR","OC","SO","JC")
related = c("KO","NA","PU","KH","KM","MB","KV","MO","JN","BE","UL","NB")
wine_hops = c("BV","LN","HO","RA","ME","UH","LT","ZN","VY","BI")

groups[,2] = ifelse(districts$abbr %in% unrelated, "red3",
                    ifelse(districts$abbr %in% related, "green3",
                           ifelse(districts$abbr %in% wine_hops, "royalblue2", "darkgrey")))

layout(matrix(c(1,3,2,4,5,5), 2, 3), widths = c(2,2,0.8))
par(mar = c(4,4,3,1), mgp = c(2.5,0.9,0))
plot(Sx[,1],Sx[,2], pch=NA, main="x-marginal (Cu)", xlab="1st component", ylab="2nd component",
     panel.first = abline(h=0,v=0,col="lightgrey"), 
     panel.last = text(Sx[,1],Sx[,2], districts$abbr, cex=0.7, font=2, col=groups[,2])) 
plot(Sy[,1],Sy[,2], pch=NA, main="y-marginal (Zn)", xlab="1st component", ylab="2nd component",
     panel.first = abline(h=0,v=0,col="lightgrey"),
     panel.last = text(Sy[,1],Sy[,2], districts$abbr, cex=0.7, font=2, col=groups[,2]))
plot(Si[,1],Si[,2], pch=NA, main="interaction", xlab="1st component", ylab="2nd component",
     panel.first = abline(h=0,v=0,col="lightgrey"),
     panel.last = text(Si[,1],Si[,2], districts$abbr, cex=0.7, font=2, col=groups[,2]))
plot(S[,1],S[,2], pch=NA, main="whole density", xlab="1st component", ylab="2nd component",
     panel.first = abline(h=0,v=0,col="lightgrey"),
     panel.last = text(S[,1],S[,2], districts$abbr, cex=0.7, font=2, col=groups[,2]))
par(mar = rep(0,4))
plot(0,0, axes=F, type="n", xlab="", ylab="")
legend("left",legend=c("unclassified","unrelated", "related","wine/hops"),
       col=c("darkgrey","red3","green3","royalblue2"), pch=15, bty="n", cex=1.3, ncol=1, pt.cex=1.5)


## Explained variability -------------------------------------------------------

# squares of singular values = eigenvalues
Var = res$d^2  

# percents of explained variability by each component
var_perc = 100*Var/sum(Var)

par(mar=c(5,4,1,1), mfrow = c(1,1))
plot(var_perc, type = "b", pch = 20, xlab = "Components", ylab = "% of explained variance", ylim=c(0,25),
     panel.first = abline(h=0, col="darkgrey"),
     panel.last = text(c(5.5,6.5),
                       var_perc[1:2]+c(0.5,0.5),
                       c(paste(round(var_perc[1],2),"%"),paste(round(var_perc[2],2),"%"))))

# variability explained by the first two components
sum(var_perc[1:2])


## Decomposition of squared norms and variance ---------------------------------

# squared norm
norm2 = function(f) sum(f^2)/(length(grid_Cu)*length(grid_Zn))

norm_tab = cbind(apply(c_clrd, 3, norm2),
                 apply(c_clrx, 3, norm2),
                 apply(c_clry, 3, norm2),
                 apply(c_clri, 3, norm2))

# a number of districts
N = nrow(districts)

par(mgp = c(2, 0.3, 0), mar = c(7,4,1,2))
fields::image.plot(1:N, 1:4, norm_tab[,4:1], col = pos_cols,
                   breaks = quantile(as.vector(norm_tab), probs = seq(0,1,0.01)),
                   axes = F, xlab = "", ylab = "", axis.args = list(cex.axis = 0.8))
axis(2, at = 1:4, labels = rev(c("density","Cu (x)","Zn (y)","interaction")), 
     tick = F, las = 2, cex.axis = 0.8)
axis(1, at = 1:N-0.1, labels = districts$name, tick = F, las = 2, cex.axis = 0.7)


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
axis(2, at = 1:3, labels = rev(c("Cu (x)","Zn (y)","interaction")), 
     tick = F, las = 2, cex.axis = 0.8)
axis(1, at = 1:N-0.1, labels = districts$name, tick = F, las = 2, cex.axis = 0.7)


## Visualization of the loadings with percentage of explained variability ------

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
par(mar=c(4,4,2,1), mgp=c(2.5,0.9,0), cex.main = 1)

image(grid_Cu, grid_Zn, V[,,1], xlim = range_Cu, ylim=range_Zn, col = cols, breaks = comp_breaks,
      xlab="log(Cu concentration)", ylab="log(Zn concentration)", 
      main = paste0(round(var_perc[1], 2)," % (100 %)"), font.main = 1)
image(grid_Cu, grid_Zn, V[,,2], xlim = range_Cu, ylim=range_Zn, col = cols, breaks = comp_breaks,
      xlab="log(Cu concentration)", ylab="log(Zn concentration)", 
      main = paste0(round(var_perc[2], 2)," % (100 %)"), font.main = 1)
#x
image(grid_Cu, grid_Zn, Vx[,,1], xlim = range_Cu, ylim=range_Zn, col = cols, breaks = comp_breaks,
      xlab="log(Cu concentration)", ylab="log(Zn concentration)", 
      main = paste0(var_perc1[1]," % (", round(perc1[1], 2)," %)"), font.main = 1)
image(grid_Cu, grid_Zn, Vx[,,2], xlim = range_Cu, ylim=range_Zn, col = cols, breaks = comp_breaks,
      xlab="log(Cu concentration)", ylab="log(Zn concentration)", 
      main = paste0(var_perc2[1]," % (", round(perc2[1], 2)," %)"), font.main = 1)
#y
image(grid_Cu, grid_Zn, Vy[,,1], xlim = range_Cu, ylim=range_Zn, col = cols, breaks = comp_breaks,
      xlab="log(Cu concentration)", ylab="log(Zn concentration)", 
      main = paste0(var_perc1[2]," % (", round(perc1[2], 2)," %)"), font.main = 1)
image(grid_Cu, grid_Zn, Vy[,,2], xlim = range_Cu, ylim=range_Zn, col = cols, breaks = comp_breaks,
      xlab="log(Cu concentration)", ylab="log(Zn concentration)", 
      main = paste0(var_perc2[2]," % (", round(perc2[2], 2)," %)"), font.main = 1)
#int
image(grid_Cu, grid_Zn, Vi[,,1], xlim = range_Cu, ylim=range_Zn, col = cols, breaks = comp_breaks,
      xlab="log(Cu concentration)", ylab="log(Zn concentration)", 
      main = paste0(var_perc1[3]," % (", round(perc1[3], 2)," %)"), font.main = 1)
image(grid_Cu, grid_Zn, Vi[,,2], xlim = range_Cu, ylim=range_Zn, col = cols, breaks = comp_breaks,
      xlab="log(Cu concentration)", ylab="log(Zn concentration)", 
      main = paste0(var_perc2[3]," % (", round(perc2[3], 2)," %)"), font.main = 1)

par(mar = c(10,0,8,3.5))
image(y = comp_breaks, z = matrix(comp_breaks,1), breaks = comp_breaks, col = cols, axes = F)
axis(4, las=2)

par(mar=c(0,0,0,0))

plot.new()
text(0.4, 0.81, "1st component", srt=90, cex=1.2, font = 2)
text(0.4, 0.27, "2nd component", srt=90, cex=1.2, font = 2)

plot.new()
mtext(c("whole density","x-marginal (price)","y-marginal (size)","interaction"),
      at = c(0.12, 0.39, 0.66, 0.93), line = -2, font = 2, cex = 0.9)


dev.off()
