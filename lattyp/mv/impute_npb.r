options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
kmax <- as.integer(args[[3]])
# print(args)
# print(args[[1]])

# q()

#library("missMDA")
d <- read.table(args[[1]], sep="\t", header=TRUE, na.strings="NA", colClass=c('factor'))

require(NPBayesImpute)
model <- CreateModel(d, NULL, kmax, 0, 0.25, 0.25)


# burn-in
model$Run(10000,0,0)

for (y in 0:99) {
        model$Run(0, 100, 1)
        result <- model$snapshot
        imputed <- GetDataFrame(result$ImputedX,d)
        write.table(imputed, file=paste(args[[2]], ".", y, sep=""), sep="\t", quote=FALSE)
}
# res.impute <- imputeMCA(d, threshold=0.015, maxiter=9000)
# write.table(res.impute$completeObs, file=args[[2]], sep="\t", quote=FALSE)
