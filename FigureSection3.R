library(circular)
set.seed(123)
# generate data mixing three von mises components
angulos_datos <- c(rvonmises(60, mu=circular(pi/2), kappa=3),
                   rvonmises(30, mu=circular(3*pi/2), kappa=1.5),
                   rvonmises(10, mu=circular(pi), kappa=1))
teta <- seq(0, 2*pi, length.out=500)

# circular kde using von mises kernel
dkde_circ <- function(grid_teta, datos, kappa) {
  datos_num <- as.numeric(datos)
  # evaluate the average of kernels at each grid point
  sapply(grid_teta, function(t) {
    mean(exp(kappa*cos(t-datos_num)))/(2*pi*besselI(kappa, 0))
  })
}

par(mfrow=c(1, 2))

# helper to draw the base circle with axis labels
dibujar_circulo <- function(xlim=c(-1.6, 1.6), ylim=c(-1.6, 1.6)) {
  plot(NA, xlim=xlim, ylim=ylim,
       asp=1, axes=FALSE, xlab="", ylab="")
  lines(cos(teta), sin(teta), lwd=1, col="black")
  points(0, 0, pch=3, cex=0.8)
  text(0,  ylim[2]-0.2, expression(pi/2),   cex=0.9)
  text(0,  ylim[1]+0.2, expression(3*pi/2), cex=0.9)
  text(xlim[1]+0.15, 0, expression(pi),     cex=0.9)
  text(xlim[2]-0.15, 0, "0",                cex=0.9)
}

# left: small kappa = oversmoothing, modes are lost
dens1 <- dkde_circ(teta, angulos_datos, kappa=0.5)
radio1 <- 1+(dens1-min(dens1))/(max(dens1)-min(dens1))*0.18
dibujar_circulo()
lines(radio1*cos(teta), radio1*sin(teta), lwd=1.5)

# right: large kappa = undersmoothing, too noisy
dens2 <- dkde_circ(teta, angulos_datos, kappa=50)
radio2 <- 1+(dens2-min(dens2))/(max(dens2)-min(dens2))*0.65
dibujar_circulo(xlim=c(-1.9, 1.9), ylim=c(-1.9, 1.9))
lines(radio2*cos(teta), radio2*sin(teta), lwd=1.5)
