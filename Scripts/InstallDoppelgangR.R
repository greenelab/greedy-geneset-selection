source("http://bioconductor.org/biocLite.R")
biocLite(c("sva", "impute", "mnormt", "digest", "knitr", "ROCR", "pROC"), suppressUpdates=TRUE)
install.packages("doppelgangR-master", repos=NULL)