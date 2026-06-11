# TFG-Paula Garcia

# circular statistics analysis and sleep and activity data

this repo contains the R code used for a university tfg on circular statistics applied to sleep timing data from the MMASH dataset (multilevel monitoring of activity and sleep in healthy people).

the project covers the main topics of circular statistics: descriptive measures, common circular distributions, kernel density estimation and non parametric regression.

---

## files

- FigureSection2.R :  figures for the circular statistics theory section: circular mean vs linear mean, and density plots for von mises, wrapped normal, wrapped cauchy, wrapped skew normal, mixture of von mises and projected normal distributions

- FigureSection3.R:  circular KDE visualization built with von mises kernel components, showing the effect of bandwidth choice (oversmoothing vs undersmoothing)

- FiguresSection5.R:  full analysis on real data: descriptive stats, circular analysis of bedtime and wake-up times, euclidean KDE with different bandwidth selectors (RT, DPI, LSCV, BCV), and circular-linear regression (Nadaraya-Watson, local linear and johnson & wehrly parametric model)

---

## data

the dataset used is the [MMASH dataset](https://physionet.org/content/mmash/1.0.0/), publicly available on PhysioNet. it contains actigraphy, sleep, saliva, RR intervals and questionnaire data from 22 healthy participants.

you need to download it and update the `ruta_base` path in `FiguresSection5.R` to match your local folder.

---

## requirements

```r
install.packages(c("tidyverse", "data.table", "circular", "NPCirc", "ggplot2"))
```

---

## notes

- bedtime hours are converted to radians with 0h at the top and clockwise orientation, to match a clock face
- hours before noon are shifted to the next day to handle the midnight crossing correctly
- the circular KDE uses a von mises kernel via `kern.den.circ()` from the NPCirc package
