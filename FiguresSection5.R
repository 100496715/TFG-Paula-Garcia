# LIBRARIES
library(tidyverse)
library(data.table)
library(circular)
library(NPCirc)

# LOAD DATA
ruta_base <- "C:/Users/paula/Downloads/multilevel-monitoring-of-activity-and-sleep-in-healthy-people-1.0.0 (1)/multilevel-monitoring-of-activity-and-sleep-in-healthy-people-1.0.0/MMASH/DataPaper"

cargar_usuarios <- function(nombre_archivo) {
  map_dfr(1:22, function(uid) {
    ruta <- file.path(ruta_base, paste0("user_", uid), nombre_archivo)
    if (file.exists(ruta)) {
      tryCatch({
        df <- data.table::fread(ruta, colClasses="character", drop=1)
        df <- type_convert(df, col_types=cols())
        df$user_id <- uid
        as_tibble(df)
      }, error = function(e) {
        message("Error in user_", uid, ": ", e$message)
        NULL
      })
    }
  })
}

# load each file for all users
info_usuario    <- cargar_usuarios("user_info.csv")
sueno           <- cargar_usuarios("sleep.csv")
actigrafo       <- cargar_usuarios("Actigraph.csv")
actividad       <- cargar_usuarios("Activity.csv")
saliva          <- cargar_usuarios("saliva.csv")
cuestionario    <- cargar_usuarios("questionnaire.csv")
rr              <- cargar_usuarios("RR.csv")

# LINEAR DESCRIPTIVE STATISTICS
stats_partic <- info_usuario %>%
  summarise(
    n            = n(),
    age_mean     = round(mean(Age, na.rm=TRUE), 1),
    age_sd       = round(sd(Age, na.rm=TRUE), 1),
    age_min      = min(Age, na.rm=TRUE),
    age_max      = max(Age, na.rm=TRUE),
    weight_mean  = round(mean(Weight, na.rm=TRUE), 1),
    weight_sd    = round(sd(Weight, na.rm=TRUE), 1),
    height_mean  = round(mean(Height, na.rm=TRUE), 1),
    height_sd    = round(sd(Height, na.rm=TRUE), 1)
  )
print(stats_partic)

stats_sueno <- sueno %>%
  summarise(
    tst_mean        = round(mean(`Total Sleep Time (TST)`, na.rm=TRUE), 1),
    tst_sd          = round(sd(`Total Sleep Time (TST)`, na.rm=TRUE), 1),
    efic_mean       = round(mean(Efficiency, na.rm=TRUE), 1),
    efic_sd         = round(sd(Efficiency, na.rm=TRUE), 1),
    waso_mean       = round(mean(`Wake After Sleep Onset (WASO)`, na.rm=TRUE), 1),
    waso_sd         = round(sd(`Wake After Sleep Onset (WASO)`, na.rm=TRUE), 1),
    despert_mean    = round(mean(`Number of Awakenings`, na.rm=TRUE), 1),
    despert_sd      = round(sd(`Number of Awakenings`, na.rm=TRUE), 1)
  )
print(stats_sueno)

stats_cuest <- cuestionario %>%
  summarise(
    meq_mean    = round(mean(MEQ, na.rm=TRUE), 1),
    meq_sd      = round(sd(MEQ, na.rm=TRUE), 1),
    psqi_mean   = round(mean(Pittsburgh, na.rm=TRUE), 1),
    psqi_sd     = round(sd(Pittsburgh, na.rm=TRUE), 1),
    stai2_mean  = round(mean(STAI2, na.rm=TRUE), 1),
    stai2_sd    = round(sd(STAI2, na.rm=TRUE), 1),
    estres_mean = round(mean(Daily_stress, na.rm=TRUE), 1),
    estres_sd   = round(sd(Daily_stress, na.rm=TRUE), 1)
  )
print(stats_cuest)

# CIRCULAR ANALYSIS

# convert clock times to radians: 0h at top, clockwise
sueno <- sueno %>%
  mutate(
    h_acostarse  = as.numeric(substr(`In Bed Time`, 1, 2)) +
      as.numeric(substr(`In Bed Time`, 4, 5))/60,
    h_despertar  = as.numeric(substr(`Out Bed Time`, 1, 2)) +
      as.numeric(substr(`Out Bed Time`, 4, 5))/60,
    # hours before noon shifted to next day
    h_acostarse  = ifelse(h_acostarse<12, h_acostarse+24, h_acostarse),
    rad_acostarse = pi/2-(h_acostarse/24)*2*pi,
    rad_despertar = pi/2-(h_despertar/24)*2*pi
  )

# circular statistics for bedtime and wake-up
acostarse_ok <- na.omit(sueno$rad_acostarse)
despertar_ok <- na.omit(sueno$rad_despertar)

acostarse_circ <- circular(acostarse_ok, units="radians")
despertar_circ <- circular(despertar_ok, units="radians")

# helper to convert radians back to readable time
rad_a_hora <- function(rad) {
  h <- (pi/2-rad)*24/(2*pi)
  h <- h%%24
  sprintf("%02d:%02d", floor(h), round((h-floor(h))*60))
}

tabla_circ <- tibble(
  Variable          = c("Bedtime", "Wake-up time"),
  Mean_direction    = c(rad_a_hora(mean(acostarse_circ)),
                        rad_a_hora(mean(despertar_circ))),
  Mean_resultant_R  = c(round(rho.circular(acostarse_circ), 3),
                        round(rho.circular(despertar_circ), 3)),
  Circular_variance = c(round(1-rho.circular(acostarse_circ), 3),
                        round(1-rho.circular(despertar_circ), 3)),
  Kappa             = c(round(mle.vonmises(acostarse_circ)$kappa, 2),
                        round(mle.vonmises(despertar_circ)$kappa, 2))
)
print(tabla_circ)

# axis labels: 0h top, 6h right, 12h bottom, 18h left
ejes_pos    <- circular(c(pi/2, 0, -pi/2, pi))
ejes_etiq   <- c("0h", "6h", "12h", "18h")

#  Bandwidth comparison for bedtime
par(mfrow=c(2, 2), mar=c(2, 2, 3, 2))

plot(kern.den.circ(acostarse_circ, bw="rt"),
     plot.type="circle", shrink=2,
     main="RT", points.plot=FALSE,
     xlab="", ylab="", axes=FALSE)
axis.circular(at=ejes_pos, labels=ejes_etiq)

plot(kern.den.circ(acostarse_circ, bw="dpi"),
     plot.type="circle", shrink=2,
     main="DPI", points.plot=FALSE,
     xlab="", ylab="", axes=FALSE)
axis.circular(at=ejes_pos, labels=ejes_etiq)

plot(kern.den.circ(acostarse_circ, bw="AA"),
     plot.type="circle", shrink=2,
     main="STE", points.plot=FALSE,
     xlab="", ylab="", axes=FALSE)
axis.circular(at=ejes_pos, labels=ejes_etiq)

plot(kern.den.circ(acostarse_circ, bw="CV"),
     plot.type="circle", shrink=2,
     main="LCV", points.plot=FALSE,
     xlab="", ylab="", axes=FALSE)
axis.circular(at=ejes_pos, labels=ejes_etiq)

# BANDWIDTH COMPARISON FOR WAKE-UP TIME 

# 1. Set up the 2x2 grid layout and margins
par(mfrow=c(2, 2), mar=c(2, 2, 3, 2))

# 2. Quadrant 1: Rule of Thumb (RT)
plot(kern.den.circ(despertar_circ, bw="rt"),
     plot.type="circle", shrink=2,
     main=" RT", points.plot=FALSE,
     xlab="", ylab="", axes=FALSE)
axis.circular(at=ejes_pos, labels=ejes_etiq)

# 3. Quadrant 2: Direct Plug-In (DPI)
plot(kern.den.circ(despertar_circ, bw="dpi"),
     plot.type="circle", shrink=2,
     main="DPI", points.plot=FALSE,
     xlab="", ylab="", axes=FALSE)
axis.circular(at=ejes_pos, labels=ejes_etiq)

# 4. Quadrant 3: Solve-The-Equation (STE)
plot(kern.den.circ(despertar_circ, bw="AA"),
     plot.type="circle", shrink=2,
     main="STE", points.plot=FALSE,
     xlab="", ylab="", axes=FALSE)
axis.circular(at=ejes_pos, labels=ejes_etiq)

# 5. Quadrant 4: Likelihood Cross-Validation (LCV)
plot(kern.den.circ(despertar_circ, bw="CV"),
     plot.type="circle", shrink=2,
     main=" LCV", points.plot=FALSE,
     xlab="", ylab="", axes=FALSE)
axis.circular(at=ejes_pos, labels=ejes_etiq)

# EUCLIDEAN KDE

# TST: compare four bandwidth selectors
tst <- na.omit(sueno$`Total Sleep Time (TST)`)
h_rt   <- bw.nrd(tst)
h_dpi  <- bw.SJ(tst, method="dpi")
h_lscv <- bw.ucv(tst)
h_bcv  <- suppressWarnings(bw.bcv(tst))

{
  par(mfrow=c(2, 2), mar=c(4, 4, 2, 1))
  plot(density(tst, bw=h_rt),   main="", xlab="TST (min)", ylab="Density")
  legend("topright", legend=paste("RT, h =",   round(h_rt, 2)),   bty="n")
  plot(density(tst, bw=h_dpi),  main="", xlab="TST (min)", ylab="Density")
  legend("topright", legend=paste("DPI, h =",  round(h_dpi, 2)),  bty="n")
  plot(density(tst, bw=h_lscv), main="", xlab="TST (min)", ylab="Density")
  legend("topright", legend=paste("LSCV, h =", round(h_lscv, 2)), bty="n")
  plot(density(tst, bw=h_bcv),  main="", xlab="TST (min)", ylab="Density")
  legend("topright", legend=paste("BCV, h =",  round(h_bcv, 2)),  bty="n")
}

# Daily Stress: build panas scores then compare bandwidths
cuestionario <- cuestionario %>%
  mutate(
    panas_pos = panas_pos_10+panas_pos_14+panas_pos_18+
      panas_pos_22+`panas_pos_9+1`,
    panas_neg = panas_neg_10+panas_neg_14+panas_neg_18+
      panas_neg_22+`panas_neg_9+1`
  )

estres_q  <- na.omit(cuestionario$Daily_stress)
h_rt_e   <- bw.nrd(estres_q)
h_dpi_e  <- bw.SJ(estres_q, method="dpi")
h_lscv_e <- bw.ucv(estres_q)
h_bcv_e  <- bw.bcv(estres_q)

par(mfrow=c(2, 2), mar=c(4, 4, 2, 1))
plot(density(estres_q, bw=h_rt_e),   main="", xlab="Daily Stress", ylab="Density")
legend("topright", legend=paste("RT, h =",   round(h_rt_e, 2)),   bty="n")
plot(density(estres_q, bw=h_dpi_e),  main="", xlab="Daily Stress", ylab="Density")
legend("topright", legend=paste("DPI, h =",  round(h_dpi_e, 2)),  bty="n")
plot(density(estres_q, bw=h_lscv_e), main="", xlab="Daily Stress", ylab="Density")
legend("topright", legend=paste("LSCV, h =", round(h_lscv_e, 2)), bty="n")
plot(density(estres_q, bw=h_bcv_e),  main="", xlab="Daily Stress", ylab="Density")
legend("topright", legend=paste("BCV, h =",  round(h_bcv_e, 2)),  bty="n")

# PANAS Positive: same bandwidth comparison
panas_p   <- na.omit(cuestionario$panas_pos)
h_rt_p   <- bw.nrd(panas_p)
h_dpi_p  <- bw.SJ(panas_p, method="dpi")
h_lscv_p <- bw.ucv(panas_p)
h_bcv_p  <- bw.bcv(panas_p)

par(mfrow=c(2, 2), mar=c(4, 4, 2, 1))
plot(density(panas_p, bw=h_rt_p),   main="", xlab="PANAS Positive", ylab="Density")
legend("topright", legend=paste("RT, h =",   round(h_rt_p, 2)),   bty="n")
plot(density(panas_p, bw=h_dpi_p),  main="", xlab="PANAS Positive", ylab="Density")
legend("topright", legend=paste("DPI, h =",  round(h_dpi_p, 2)),  bty="n")
plot(density(panas_p, bw=h_lscv_p), main="", xlab="PANAS Positive", ylab="Density")
legend("topright", legend=paste("LSCV, h =", round(h_lscv_p, 2)), bty="n")
plot(density(panas_p, bw=h_bcv_p),  main="", xlab="PANAS Positive", ylab="Density")
legend("topright", legend=paste("BCV, h =",  round(h_bcv_p, 2)),  bty="n")

# REGRESSION
# merge sleep and questionnaire, average per user
datos <- inner_join(sueno, cuestionario, by="user_id")

datos <- datos %>%
  group_by(user_id) %>%
  summarise(across(where(is.numeric), mean, na.rm=TRUE)) %>%
  ungroup()

datos <- datos %>%
  mutate(
    panas_pos = panas_pos_10+panas_pos_14+panas_pos_18+
      panas_pos_22+`panas_pos_9+1`,
    panas_neg = panas_neg_10+panas_neg_14+panas_neg_18+
      panas_neg_22+`panas_neg_9+1`
  )

cat("Participants:", nrow(datos), "\n")

# MODEL 1: CIRCULAR-LINEAR 
datos_reg <- datos %>%
  select(rad_acostarse, Daily_stress, panas_pos) %>%
  na.omit()

bedtime_reg <- circular(datos_reg$rad_acostarse, units="radians")
estres      <- datos_reg$Daily_stress
panas       <- datos_reg$panas_pos

# select bandwidth via cross-validation
bw_reg <- bw.reg.circ.lin(bedtime_reg, estres)
cat("Circular-linear bandwidth:", bw_reg, "\n")

# fit NW and LL estimators for both responses
est_nw_estres <- kern.reg.circ.lin(bedtime_reg, estres, bw=bw_reg, method="NW")
est_ll_estres <- kern.reg.circ.lin(bedtime_reg, estres, bw=bw_reg, method="LL")
est_nw_panas  <- kern.reg.circ.lin(bedtime_reg, panas,  bw=bw_reg, method="NW")
est_ll_panas  <- kern.reg.circ.lin(bedtime_reg, panas,  bw=bw_reg, method="LL")

# line plot: NW vs LL
par(mfrow=c(1, 2))

# 1. Align the domain
grid_x <- as.numeric(est_nw_estres$x) - 2*pi

# 2.cut the lines exactly at the min and max data points (removes the unstable tails)
valid_idx <- grid_x >= min(datos_reg$rad_acostarse) & grid_x <= max(datos_reg$rad_acostarse)

# Plot 1: Daily stress (NW vs LL)
plot(as.numeric(datos_reg$rad_acostarse), datos_reg$Daily_stress,
     xlab="Bedtime", ylab="Daily Stress",
     main="",
     xlim=c(min(datos_reg$rad_acostarse)-0.1, max(datos_reg$rad_acostarse)+0.1),
     ylim=c(10, 80), col="gray60", pch=1)

lines(grid_x[valid_idx], est_nw_estres$y[valid_idx], col="black", lwd=2)
lines(grid_x[valid_idx], est_ll_estres$y[valid_idx], col="red", lwd=2)
# Moved legend to top-left to avoid overlaps
legend("topleft", legend=c("NW", "LL"), col=c("black", "red"), lty=1, bty="n")

# Plot 2: PANAS Positive (NW vs LL)
plot(as.numeric(datos_reg$rad_acostarse), datos_reg$panas_pos,
     xlab="Bedtime", ylab="PANAS Positive",
     main="",
     xlim=c(min(datos_reg$rad_acostarse)-0.1, max(datos_reg$rad_acostarse)+0.1),
     ylim=c(70, 170), col="gray60", pch=1)

lines(grid_x[valid_idx], est_nw_panas$y[valid_idx], col="black", lwd=2)
lines(grid_x[valid_idx], est_ll_panas$y[valid_idx], col="red", lwd=2)
# Moved legend to top-left to avoid overlaps
legend("topleft", legend=c("NW", "LL"), col=c("black", "red"), lty=1, bty="n")

# circle plot: NW estimator (Watson-style)
{
  par(mfrow=c(1, 2), mar=c(2, 2, 3, 2))
  plot(est_nw_estres,
       plot.type="circle", points.plot=TRUE,
       main="", line.col=1, lwd=2, points.col="gray50")
  plot(est_nw_panas,
       plot.type="circle", points.plot=TRUE,
       main="", line.col=1, lwd=2, points.col="gray50")
}

#MODEL 2: PARAMETRIC Johnson & Wehrly 
mod_estres <- lm(Daily_stress~sin(rad_acostarse)+cos(rad_acostarse),
                 data=datos_reg)
summary(mod_estres)

mod_panas <- lm(panas_pos~sin(rad_acostarse)+cos(rad_acostarse),
                data=datos_reg)
summary(mod_panas)

# predict over observed range
teta_seq <- seq(min(datos_reg$rad_acostarse),
                max(datos_reg$rad_acostarse),
                length.out=200)

pred_estres <- predict(mod_estres,
                       newdata=data.frame(rad_acostarse=teta_seq))
pred_panas  <- predict(mod_panas,
                       newdata=data.frame(rad_acostarse=teta_seq))

# plot parametric fit over raw data
{
  par(mfrow=c(1, 2), mar=c(4, 4, 2, 2))
  plot(datos_reg$rad_acostarse, datos_reg$Daily_stress,
       pch=1, col="gray40",
       xlab="Bedtime", ylab="Daily Stress", main="")
  lines(teta_seq, pred_estres, col="black", lwd=2)
  
  plot(datos_reg$rad_acostarse, datos_reg$panas_pos,
       pch=1, col="gray40",
       xlab="Bedtime", ylab="PANAS Positive", main="")
  lines(teta_seq, pred_panas, col="black", lwd=2)
}
