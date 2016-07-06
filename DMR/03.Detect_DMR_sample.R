library('BSgenome.NMR.BGI.sample')
library(MEDIPS)
#genome
#genome$PseudoC_23
BM=c(paste('/home/enjoycloud/Documents/MyLocal/NakedMole/Meth2/SAMPLE/BAM/BM_brain',c(1:3)
           ,'.sort.bam',sep=''))
WM=c(paste('/home/enjoycloud/Documents/MyLocal/NakedMole/Meth2/SAMPLE/BAM/WM_brain',c(1:3)
           ,'.sort.bam',sep=''))
BSgenome='BSgenome.NMR.BGI.sample'
uniq=TRUE
extend=300
shift=0
ws=100
chr.select=c(paste('PseudoC',c(0:23),sep='_'))
#sr = MEDIPS.saturation(file = BM[1], BSgenome = BSgenome,
#                       uniq = uniq, extend = extend, shift = shift, window_size = ws,
#                       chr.select = chr.select, nit = 10, nrit = 1, empty_bins = TRUE,
#                       rank = FALSE)
#sr
#MEDIPS.plotSaturation(sr)
#cr = MEDIPS.seqCoverage(file = BM[1], pattern = "CG",
#                        BSgenome = BSgenome, chr.select = chr.select, extend = extend,
#                        shift = shift, uniq = uniq)
#MEDIPS.plotSeqCoverage(seqCoverageObj = cr, type = "pie", cov.level = c(0,1,5,10,50))
#MEDIPS.plotSeqCoverage(seqCoverageObj = cr, type = "hist", t = 200,main = "Sequence pattern coverage, histogram")
#er = MEDIPS.CpGenrich(file = BM[1], BSgenome = BSgenome,chr.select = chr.select, extend = extend, shift = shift,
#                      uniq = uniq)
#er$enrichment.score.relH
BM_MID=vector("list", 3)
for (i in seq(1,3,1)){BM_MID[i]= MEDIPS.createSet(file = BM[i],
                                                  BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
                                                  window_size = ws, chr.select = chr.select)}
WM_MID=vector("list", 3)
for (i in seq(1,3,1)){WM_MID[i]= MEDIPS.createSet(file = WM[i],
                                                  BSgenome = BSgenome, extend = extend, shift = shift, uniq = uniq,
                                                  window_size = ws, chr.select = chr.select)}
CS = MEDIPS.couplingVector(pattern = "CG", refObj = BM_MID[[1]])
#6.2 Coverage, methylation 6.2 Coverage, methylation p proles and diferential 
##Detect  Coverage ####
mr.edgeR = MEDIPS.meth(MSet1 = BM_MID, MSet2 = WM_MID,
                       CSet = CS, p.adj = "fdr",
                       diff.method = "edgeR", prob.method = "poisson", MeDIP = T,
                       CNV = F, type = "rpkm", minRowSum = 1)
## Detect DMR , based on your criteria ####
mr.edgeR.s = MEDIPS.selectSig(results = mr.edgeR, p.value = 0.1,
                              adj = T, ratio = NULL, bg.counts = NULL, CNV = F)
mr.edgeR.s.gain = mr.edgeR.s[which(mr.edgeR.s[, grep("logFC",colnames(mr.edgeR.s))] > 0), ]
mr.edgeR.s.lost = mr.edgeR.s[which(mr.edgeR.s[, grep("logFC",colnames(mr.edgeR.s))] < 0), ]
### Merging####
mr.edgeR.s.gain.m = MEDIPS.mergeFrames(frames = mr.edgeR.s.gain,distance = 1)
mr.edgeR.s.lost.m = MEDIPS.mergeFrames(frames = mr.edgeR.s.lost,distance = 1)
### Summary DMR ####
columns = names(mr.edgeR)[grep("counts", names(mr.edgeR))]
rois_gain = MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.gain.m,columns = columns, summarize = NULL)
rois_gain.s = MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.gain.m,columns = columns, summarize = "avg")
rois_lost = MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.lost.m,columns = columns, summarize = NULL)
rois_lost.s = MEDIPS.selectROIs(results = mr.edgeR, rois = mr.edgeR.s.lost.m,columns = columns, summarize = "avg")
### Miscellaneous###
##Plot##
#MEDIPS.plotCalibrationPlot(CSet = CS, main = "Calibration Plot",
#                           MSet = BM_MID[[1]], plot_chr = 'PseudoC_3', rpkm = TRUE,
#                           xrange = TRUE)
cor.matrix = MEDIPS.correlation(MSets = c(BM_MID,WM_MID),plot=F)
save(BM_MID,WM_MID,mr.edgeR,file='Sample_DMR_Feb_12.Rdata')
