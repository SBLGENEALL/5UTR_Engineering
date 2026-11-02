#!/usr/bin/env Rscript
# Plot evidence-source composition across four nested selection stages.
# Usage: Rscript plot_evidence_source_composition.R --results results/four_stage
args <- commandArgs(trailingOnly=TRUE)
opt <- function(k,d) { i<-match(k,args); if(!is.na(i)&&i<length(args)) args[i+1] else d }
indir <- opt("--results","results/four_stage")
input <- file.path(indir,"candidate_groups_and_features.csv")
if(!file.exists(input)) stop("Missing ",input,"; run analyze_four_stages.R first")
d <- read.csv(input,check.names=FALSE,stringsAsFactors=FALSE)
if(!all(c("analysis_category","sortseq_rank","stage") %in% names(d))) stop("Required columns missing")
if(nrow(d)!=2000L) stop("Expected 2000 library rows, found ",nrow(d))
catcode <- toupper(sub("^\\s*([A-H]|K[1-4])(?:[^A-Za-z0-9].*)?$","\\1",trimws(as.character(d$analysis_category)),perl=TRUE))
# Also accept category codes embedded in labels such as A_publicTE...
evidence <- c("A","C","D","E","K1","K2")
predictor <- c("B","K3")
exploratory <- c("F","G","H","K4")
valid <- c(evidence,predictor,exploratory)
unknown <- unique(d$analysis_category[!catcode %in% valid])
if(length(unknown)) stop("Unrecognized categories: ",paste(unknown,collapse=", "),". Review category mapping before plotting.")
d$evidence_source <- ifelse(catcode %in% evidence,"Experimental evidence",ifelse(catcode %in% predictor,"Predictor-based","Exploratory / control"))
stage_names <- c("Library 2000","Top 200","Tested 61","Hits 14")
groups <- list(rep(TRUE,nrow(d)),!is.na(d$sortseq_rank),!is.na(d$sortseq_rank)&d$sortseq_rank<=61,d$stage=="Hit 14")
stage_n <- vapply(groups,sum,integer(1))
if(!identical(unname(stage_n),c(2000L,200L,61L,14L))) stop("Unexpected stage sizes: ",paste(stage_n,collapse=", "))
source_names <- c("Experimental evidence","Predictor-based","Exploratory / control")
count <- sapply(groups,function(g) as.integer(table(factor(d$evidence_source[g],levels=source_names))))
rownames(count) <- source_names
colnames(count) <- stage_names
pct <- sweep(count,2,colSums(count),"/")*100
long <- do.call(rbind,lapply(seq_along(stage_names),function(j) data.frame(stage=stage_names[j],stage_n=stage_n[j],evidence_source=source_names,n=count[,j],percent=round(pct[,j],2),stringsAsFactors=FALSE)))
write.csv(long,file.path(indir,"evidence_source_composition.csv"),row.names=FALSE)
colors <- c("#1428A0","#00AFC1","#A6B2C5")
out <- file.path(indir,"evidence_source_composition.pdf")
pdf(out,width=11.8,height=6.6,pointsize=12)
par(mar=c(5.4,9,4.7,2.5),xpd=FALSE)
bp <- barplot(pct[,4:1,drop=FALSE],horiz=TRUE,col=colors,border="white",space=.65,
              xlim=c(0,100),names.arg=stage_names[4:1],las=1,cex.names=1.08,
              xlab="Composition (%)",main="Selection-source composition across four stages",
              cex.main=1.2)
# White text on dark blue, dark text on cyan/grey; omit labels for tiny segments.
for(j in 1:4) {
  acc <- 0
  for(i in 1:3) {
    v <- pct[i,5-j]
    if(v>=7) {
      label <- sprintf("%d (%.1f%%)",count[i,5-j],v)
      text(acc+v/2,bp[j],label,col=if(i==1)"white" else "#14213D",font=2,cex=.83)
    }
    acc <- acc+v
  }
}
legend("bottomright",legend=source_names,fill=colors,border=NA,bty="n",cex=.85,inset=c(.01,.02))
mtext("Category reflects original library selection route, not exclusive evidence availability.",side=1,line=4.1,cex=.76,col="#4C5871")
dev.off()
cat("Created: ",out,"\nCSV: ",file.path(indir,"evidence_source_composition.csv"),"\n",sep="")
print(long,row.names=FALSE)
