debugData <- readRDS("dev/debugDataTraces.rds")
C_sub <- debugData$C_sub
nvarGi.sparse <- debugData$nvarGi.sparse
nvarGi.dense <- debugData$nvarGi.dense
dim.var.corrected <- debugData$dim.var.corrected
dim.nvar <- debugData$dim.nvar

all.equal(as.matrix(nvarGi.sparse), as.matrix(nvarGi.dense)) # TRUE

traces.sparse <- forceSymmetric(kmmr:::getTracesSparse(C_sub@p,
                                                       C_sub@i,
                                                       C_sub@x,
                                                       nvarGi.sparse@p,
                                                       nvarGi.sparse@i,
                                                       nvarGi.sparse@x,
                                                       dim.var.corrected,
                                                       dim.nvar))

traces.dense <- forceSymmetric(kmmr:::getTraces(C_sub@p,
                                                C_sub@i,
                                                C_sub@x,
                                                as.matrix(nvarGi.dense),
                                                dim.var.corrected,
                                                dim.nvar))

all.equal(traces.dense, traces.sparse) # "Mean relative difference: 8.445986"
