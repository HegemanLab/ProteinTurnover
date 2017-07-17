## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
library(ProteinTurnover)

## ---- examples-----------------------------------------------------------
ex_ra <- relAbForTimes(data=isocounts)
# plot(ex_ra)
# regressionPlot(ex_ra)

ex_di <- pepfit(TimePoint, RelAb, Channel, data=isodata,
            Elements=list(N=12,C=45,H=73,O=15))

ex_in <- pepfit(TimePoint, RelAb, Channel, data=isoincorp,
                type="incorporation", 
                Elements=list(N=12,C=45,H=73,O=15))

