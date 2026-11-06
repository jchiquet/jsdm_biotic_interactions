## Theme, colours and labels shared by the figure scripts. Figures are written
## to figures/; their captions are in the manuscript.
suppressMessages(library(ggplot2))
FIG.DIR <- "figures"
dir.create(FIG.DIR, showWarnings = FALSE)
fig_path <- function(name) file.path(FIG.DIR, paste0(name, ".png"))

theme_set(theme_minimal(base_size = 12))
th <- theme(plot.title = element_text(face = "bold", size = 14),
            plot.subtitle = element_text(size = 10.5, colour = "grey30"),
            strip.text = element_text(face = "bold"),
            panel.grid.minor = element_blank(),
            legend.position = "bottom")

## mechanisms (Okabe-Ito colours)
MECH <- c(lottery = "lottery", gompertz = "Gompertz", ricker = "Ricker")
MECH.COL <- setNames(c("#0072B2", "#D55E00", "#009E73"), MECH)
mech_factor <- function(x) factor(MECH[x], MECH)

METRIC.LAB <- c(AUC_competition          = "AUC, competition",
                AUC_facilitation         = "AUC, facilitation",
                precision_k_competition  = "top-k precision, competition",
                precision_k_facilitation = "top-k precision, facilitation")
metric_factor <- function(x) factor(METRIC.LAB[x], METRIC.LAB)

BASE.LAB <- c(no.covar = "none", cat.2 = "categorical (10)", quadratic = "quadratic",
              poly.3 = "cubic polynomial", spline.2 = "spline (10 df)")

## the two estimators most figures compare
LAB.PLNNET <- "PLNnetwork"
LAB.QUICK  <- "log-GGM"
EST.LAB <- c(plnnetwork = LAB.PLNNET, quick = LAB.QUICK)
EST.COL <- setNames(COLOURS[c("plnnetwork", "quick")], EST.LAB)
EST.SHP <- setNames(c(16, 17), EST.LAB)
LAB.INC <- "increment"

## facet row strips on two lines
row_wrap <- labeller(.rows = label_wrap_gen(16))
