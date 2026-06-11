library(ggplot2)

# Data
grados <- c(35, 280)
radianes <- grados*pi/180

# Circular mean
seno_medio <- mean(sin(radianes))
coseno_medio <- mean(cos(radianes))
media_circ_rad <- atan2(seno_medio, coseno_medio)
media_circ_grados <- media_circ_rad*180/pi
if (media_circ_grados < 0) media_circ_grados <- media_circ_grados+360

# Linear mean
media_lin_grados <- mean(grados)  # 157.5
media_lin_rad <- media_lin_grados*pi/180

# Circle points
angulo <- seq(0, 2*pi, length.out=500)
circulo_df <- data.frame(x=cos(angulo), y=sin(angulo))

# Data points on circle
puntos <- data.frame(
  angulo_grados = grados,
  x = cos(radianes),
  y = sin(radianes)
)
puntos$etiqueta <- paste0(grados, "°")
puntos$lx <- 1.18*puntos$x
puntos$ly <- 1.18*puntos$y

r <- 0.97
circ_x <- r*cos(media_circ_rad)
circ_y <- r*sin(media_circ_rad)
lin_x  <- r*cos(media_lin_rad)
lin_y  <- r*sin(media_lin_rad)

rl <- 1.22
circ_lx <- rl*cos(media_circ_rad)
circ_ly <- rl*sin(media_circ_rad)
lin_lx  <- rl*cos(media_lin_rad)
lin_ly  <- rl*sin(media_lin_rad)

marcas <- data.frame(
  angulo_grados = c(0, 90, 180, 270),
  etiqueta = c("0°", "90°", "180°", "270°")
)
marcas$xlabel <- 1.2*cos(marcas$angulo_grados*pi/180)
marcas$ylabel <- 1.2*sin(marcas$angulo_grados*pi/180)

ggplot() +
  geom_path(data=circulo_df, aes(x, y), color="gray50", linewidth=0.7) +
  geom_segment(aes(x=-1.15, xend=1.15, y=0, yend=0),
               linetype="dashed", color="gray60", linewidth=0.4) +
  geom_segment(aes(x=0, xend=0, y=-1.15, yend=1.15),
               linetype="dashed", color="gray60", linewidth=0.4) +
  geom_text(data=marcas, aes(x=xlabel, y=ylabel, label=etiqueta),
            size=3.8, color="gray20") +
  geom_point(data=puntos, aes(x, y), color="#1A56FF", size=5) +
  geom_text(data=puntos, aes(x=lx, y=ly, label=etiqueta),
            color="#1A56FF", size=4, fontface="bold") +
  geom_segment(aes(x=0, y=0, xend=circ_x, yend=circ_y),
               arrow=arrow(length=unit(0.3, "cm"), type="closed"),
               color="#FF1A1A", linewidth=1.5) +
  geom_segment(aes(x=0, y=0, xend=lin_x, yend=lin_y),
               arrow=arrow(length=unit(0.3, "cm"), type="closed"),
               color="#00CC00", linewidth=1.5) +
  annotate("text", x=circ_lx+0.1, y=circ_ly,
           label=paste0("Circular mean\n", round(media_circ_grados, 1), "°"),
           color="#FF1A1A", size=3.8, fontface="bold", hjust=0) +
  annotate("text", x=lin_lx-0.1, y=lin_ly,
           label=paste0("Linear mean\n", round(media_lin_grados, 1), "°"),
           color="#00CC00", size=3.8, fontface="bold", hjust=1) +
  coord_fixed(xlim=c(-1.6, 1.6), ylim=c(-1.6, 1.6)) +
  theme_void()



# von Mises density
par(mfrow=c(1, 2))

teta <- seq(0, 2*pi, length.out=500)
ux <- cos(teta)
uy <- sin(teta)

dibujar_circulo <- function() {
  plot(NA, xlim=c(-1.6, 1.6), ylim=c(-1.6, 1.6),
       asp=1, axes=FALSE, xlab="", ylab="")
  lines(ux, uy, lwd=1, col="black")
  points(c(1, 0, -1, 0), c(0, 1, 0, -1), pch=0, cex=0.9, lwd=1.5)
  text(0,  1.42, expression(pi/2),   cex=0.95, font=2, col="black")
  text(0, -1.42, expression(3*pi/2), cex=0.95, font=2, col="black")
  text(-1.42, 0, expression(pi),     cex=0.95, font=2, col="black")
  text( 1.42, 0, "0",                cex=0.95, font=2, col="black")
}

# Left: mu = 0, kappa = 1.5
dens1 <- dvonmises(circular(teta), mu=circular(0), kappa=1.5)
radio1 <- 1+dens1/max(dens1)*0.5
dibujar_circulo()
lines(radio1*cos(teta), radio1*sin(teta), lwd=2, col="black")

# Right: mu = pi/2, kappa = 2.5
dens2 <- dvonmises(circular(teta), mu=circular(pi/2), kappa=2.5)
radio2 <- 1+dens2/max(dens2)*0.5
dibujar_circulo()
lines(radio2*cos(teta), radio2*sin(teta), lwd=2, col="black")


# Wrapped Normal
teta <- seq(0, 2*pi, length.out=500)

dnormal_envuelta <- function(teta, mu, rho, K=50) {
  dens <- rep(0, length(teta))
  for (k in -K:K) {
    dens <- dens+dnorm(teta-mu+2*pi*k, mean=0, sd=sqrt(-2*log(rho)))
  }
  return(dens)
}

par(mfrow=c(1, 2))

dibujar_circulo <- function() {
  plot(NA, xlim=c(-1.6, 1.6), ylim=c(-1.6, 1.6),
       asp=1, axes=FALSE, xlab="", ylab="")
  lines(cos(teta), sin(teta), lwd=1, col="black")
  points(c(1, 0, -1, 0), c(0, 1, 0, -1), pch=0, cex=0.9, lwd=1.5)
  text(0,  1.42, expression(pi/2),   cex=0.95, font=2, col="black")
  text(0, -1.42, expression(3*pi/2), cex=0.95, font=2, col="black")
  text(-1.42, 0, expression(pi),     cex=0.95, font=2, col="black")
  text( 1.42, 0, "0",                cex=0.95, font=2, col="black")
}

# Left: mu = 0, rho = 0.2
dens1 <- dnormal_envuelta(teta, mu=0, rho=0.2)
radio1 <- 1+dens1/max(dens1)*0.5
dibujar_circulo()
lines(radio1*cos(teta), radio1*sin(teta), lwd=2, col="black")

# Right: mu = pi/2, rho = 0.8
dens2 <- dnormal_envuelta(teta, mu=pi/2, rho=0.8)
radio2 <- 1+dens2/max(dens2)*0.5
dibujar_circulo()
lines(radio2*cos(teta), radio2*sin(teta), lwd=2, col="black")


# Wrapped Cauchy density
teta <- seq(0, 2*pi, length.out=500)

dcauchy_envuelta <- function(teta, mu, rho) {
  (1-rho^2)/(2*pi*(1+rho^2-2*rho*cos(teta-mu)))
}

par(mfrow=c(1, 2))

dibujar_circulo <- function(xlim=c(-1.6, 1.6), ylim=c(-1.6, 1.6)) {
  plot(NA, xlim=xlim, ylim=ylim,
       asp=1, axes=FALSE, xlab="", ylab="")
  lines(cos(teta), sin(teta), lwd=1, col="black")
  points(c(1, 0, -1, 0), c(0, 1, 0, -1), pch=0, cex=0.9, lwd=1.5)
  text(0,  1.42, expression(pi/2),   cex=0.95, font=2)
  text(0, -1.42, expression(3*pi/2), cex=0.95, font=2)
  text(-1.42, 0, expression(pi),     cex=0.95, font=2)
  text( 1.42, 0, "0",                cex=0.95, font=2)
}

# Left: mu = 0, rho = 0.2
dens1 <- dcauchy_envuelta(teta, mu=0, rho=0.2)
radio1 <- 1+(dens1-min(dens1))/(max(dens1)-min(dens1))*0.15
dibujar_circulo()
lines(radio1*cos(teta), radio1*sin(teta), lwd=2, col="black")

# Right: mu = pi/2, rho = 0.8
dens2 <- dcauchy_envuelta(teta, mu=pi/2, rho=0.8)
radio2 <- 1+(dens2-min(dens2))/(max(dens2)-min(dens2))*1.2
dibujar_circulo(xlim=c(-1.8, 1.8), ylim=c(-1.8, 2.4))
lines(radio2*cos(teta), radio2*sin(teta), lwd=2, col="black")


# Wrapped skew normal density
dnormal_sesgo_env <- function(teta, mu, kappa, lambda, K=50) {
  dens <- rep(0, length(teta))
  for (k in -K:K) {
    z <- (teta-mu+2*pi*k)*kappa
    dens <- dens+2*kappa*dnorm(z)*pnorm(lambda*z)
  }
  return(dens)
}

par(mfrow=c(1, 2))

dibujar_circulo <- function(xlim=c(-1.6, 1.6), ylim=c(-1.6, 1.6)) {
  plot(NA, xlim=xlim, ylim=ylim,
       asp=1, axes=FALSE, xlab="", ylab="")
  lines(cos(teta), sin(teta), lwd=1, col="black")
  points(c(1, 0, -1, 0), c(0, 1, 0, -1), pch=0, cex=0.9, lwd=1.5)
  text(0,  1.42, expression(pi/2),   cex=0.95, font=2)
  text(0, -1.42, expression(3*pi/2), cex=0.95, font=2)
  text(-1.42, 0, expression(pi),     cex=0.95, font=2)
  text( 1.42, 0, "0",                cex=0.95, font=2)
}

# Left: mu = 0, kappa = 0.5, lambda = 3
dens1 <- dnormal_sesgo_env(teta, mu=0, kappa=0.5, lambda=3)
radio1 <- 1+(dens1-min(dens1))/(max(dens1)-min(dens1))*0.5
dibujar_circulo()
lines(radio1*cos(teta), radio1*sin(teta), lwd=2, col="black")

# Right: mu = pi/2, kappa = 3, lambda = 3
dens2 <- dnormal_sesgo_env(teta, mu=pi/2, kappa=3, lambda=3)
radio2 <- 1+(dens2-min(dens2))/(max(dens2)-min(dens2))*1.1
dibujar_circulo(xlim=c(-1.8, 1.8), ylim=c(-1.8, 2.4))
lines(radio2*cos(teta), radio2*sin(teta), lwd=2, col="black")


# Mixture of von Mises
dmezcla_vonmises <- function(teta, vec_mu, vec_kappa, vec_p) {
  dens <- rep(0, length(teta))
  for (m in seq_along(vec_mu)) {
    dens <- dens+vec_p[m]*dvonmises(circular(teta),
                                    mu=circular(vec_mu[m]),
                                    kappa=vec_kappa[m])
  }
  return(as.numeric(dens))
}

par(mfrow=c(1, 2))

dibujar_circulo <- function(xlim=c(-1.6, 1.6), ylim=c(-1.6, 1.6)) {
  plot(NA, xlim=xlim, ylim=ylim,
       asp=1, axes=FALSE, xlab="", ylab="")
  lines(cos(teta), sin(teta), lwd=1, col="black")
  points(c(1, 0, -1, 0), c(0, 1, 0, -1), pch=0, cex=0.9, lwd=1.5)
  text(0,  1.42, expression(pi/2),   cex=0.95, font=2)
  text(0, -1.42, expression(3*pi/2), cex=0.95, font=2)
  text(-1.42, 0, expression(pi),     cex=0.95, font=2)
  text( 1.42, 0, "0",                cex=0.95, font=2)
}

# Left: symmetric bimodal mu=(0, pi), kappa=(4,4), p=(0.5,0.5)
dens1 <- dmezcla_vonmises(teta, vec_mu=c(0, pi), vec_kappa=c(4, 4), vec_p=c(0.5, 0.5))
radio1 <- 1+(dens1-min(dens1))/(max(dens1)-min(dens1))*0.5
dibujar_circulo()
lines(radio1*cos(teta), radio1*sin(teta), lwd=2, col="black")

# Right: asymmetric mu=(pi/2, 0), kappa=(8,2), p=(0.7,0.3)
dens2 <- dmezcla_vonmises(teta, vec_mu=c(pi/2, 0), vec_kappa=c(8, 2), vec_p=c(0.7, 0.3))
radio2 <- 1+(dens2-min(dens2))/(max(dens2)-min(dens2))*1.0
dibujar_circulo(xlim=c(-1.8, 1.8), ylim=c(-1.8, 2.3))
lines(radio2*cos(teta), radio2*sin(teta), lwd=2, col="black")


# Projected normal density
teta <- seq(0, 2*pi, length.out=500)

dnormal_proy <- function(teta, mu, sigma) {
  dens <- numeric(length(teta))
  for(i in seq_along(teta)){
    ci <- cos(teta[i]-mu)
    si <- sin(teta[i]-mu)
    dens[i] <- exp(-(si^2)/(2*sigma^2))/(2*pi*sigma^2)*
      (sigma*dnorm(ci/sigma)*sqrt(2*pi)+
         ci*(2*pnorm(ci/sigma)-1)*sqrt(2*pi)*sigma*dnorm(0))
  }
  dens <- pmax(dens, 0)
  return(dens)
}

par(mfrow=c(1, 2))

dibujar_circulo <- function(xlim=c(-1.6, 1.6), ylim=c(-1.6, 1.6)) {
  plot(NA, xlim=xlim, ylim=ylim,
       asp=1, axes=FALSE, xlab="", ylab="")
  lines(cos(teta), sin(teta), lwd=1, col="black")
  points(c(1, 0, -1, 0), c(0, 1, 0, -1), pch=0, cex=0.9, lwd=1.5)
  text(0,  ylim[2]-0.2, expression(pi/2),   cex=0.95, font=2)
  text(0,  ylim[1]+0.2, expression(3*pi/2), cex=0.95, font=2)
  text(xlim[1]+0.2, 0,  expression(pi),     cex=0.95, font=2)
  text(xlim[2]-0.2, 0,  "0",                cex=0.95, font=2)
}

# Left: mu = 0, sigma = 3
dens1 <- dnormal_proy(teta, mu=0, sigma=3)
radio1 <- 1+(dens1-min(dens1))/(max(dens1)-min(dens1))*0.15
dibujar_circulo()
lines(radio1*cos(teta), radio1*sin(teta), lwd=2)

# Right: mu = pi/2, sigma = 0.5
dens2 <- dnormal_proy(teta, mu=pi/2, sigma=0.5)
radio2 <- 1+(dens2-min(dens2))/(max(dens2)-min(dens2))*0.7
dibujar_circulo(xlim=c(-1.8, 1.8), ylim=c(-1.8, 2.5))
lines(radio2*cos(teta), radio2*sin(teta), lwd=2)
