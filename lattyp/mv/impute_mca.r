options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
# print(args)
# print(args[[1]])

# q()

library("missMDA")
d <- read.table(args[[1]], sep="\t", header=TRUE, na.strings="NA", colClass=c('factor'))
res.impute <- imputeMCA(d, threshold=0.015, maxiter=9000)
write.table(res.impute$completeObs, file=args[[2]], sep="\t", quote=FALSE)
