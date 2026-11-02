#!/usr/bin/env Rscript
# Publication-style data-driven summary; base R only.
args <- commandArgs(trailingOnly=TRUE)
opt <- function(flag, default) { i <- match(flag,args); if(!is.na(i) && i<length(args)) args[i+1] else default }
out <- opt("--results","results")
alpha <- as.numeric(opt("--alpha","0.05"))
ft <- read.csv(file.path(out,"feature_tests.csv"),check.names=FALSE)
ct <- read.csv(file.path(out,"category_enrichment.csv"),check.names=FALSE)
hits <- read.csv(file.path(out,"hit_feature_table.csv"),check.names=FALSE)
if(!nrow(ft)) stop("feature_tests.csv has no rows")
cmp <- if("tested_only"%in%ft$comparison) "tested_only" else "whole_library"
f <- ft[ft$comparison==cmp,,drop=FALSE]
c <- ct[ct$comparison==cmp,,drop=FALSE]
f <- f[order(f$p_holm,f$p),,drop=FALSE]
c <- c[order(c$p_holm,c$p),,drop=FALSE]
sigf <- f[is.finite(f$p_holm)&f$p_holm<alpha,,drop=FALSE]
sigc <- c[is.finite(c$p_holm)&c$p_holm<alpha,,drop=FALSE]
write.csv(sigf,file.path(out,"significant_features.csv"),row.names=FALSE)
write.csv(sigc,file.path(out,"significant_categories.csv"),row.names=FALSE)
fmt <- function(x,d=3) ifelse(is.finite(x),formatC(x,format="f",digits=d),"NA")
dir.create(out,showWarnings=FALSE,recursive=TRUE)
lines <- c(paste0("Comparison: ",cmp),paste0("Features tested: ",nrow(f)),
           paste0("Significant features (Holm-adjusted p < ",alpha,"): ",nrow(sigf)),
           paste0("Significant categories: ",nrow(sigc)),"",
           "SIGNIFICANT FEATURES")
if(nrow(sigf)) for(i in seq_len(nrow(sigf))) {
 z<-sigf[i,]
 lines<-c(lines,sprintf("%s | hit mean %s | other mean %s | difference %s | Holm p %s",
  z$feature,fmt(z$hit_mean),fmt(z$other_mean),fmt(z$mean_difference),fmt(z$p_holm,4)))
} else lines<-c(lines,"None after Holm correction")
lines<-c(lines,"","SIGNIFICANT CATEGORIES")
if(nrow(sigc)) for(i in seq_len(nrow(sigc))) {
 z<-sigc[i,]
 lines<-c(lines,sprintf("%s | hit n %s | other n %s | Holm p %s",
 z$category,z$hit_n,z$other_n,fmt(z$p_holm,4)))
} else lines<-c(lines,"None after Holm correction")
lines<-c(lines,"","Interpretation: higher means are not automatically better or significant.",
 "Whole-library comparison is exploratory; tested-only is preferred when validation data are supplied.",
 "Category and feature p-values are Holm-adjusted within each comparison family.",
 "Associations do not establish causal effects.")
writeLines(lines,file.path(out,"significant_summary.txt"))

navy<-"#162744"; purple<-"#7155B8"; grey<-"#A7B0BF"; ink<-"#263447"; pale<-"#F2EFFA"
pdf(file.path(out,"top_utr_summary.pdf"),width=15,height=9.5,useDingbats=FALSE)
par(mar=c(0,0,0,0),family="sans")
plot.new();plot.window(xlim=c(0,15),ylim=c(0,9.5))
rect(0,8.55,15,9.5,col=navy,border=NA)
text(.55,9.12,"TOP 5'UTR CANDIDATES | DATA-DRIVEN SUMMARY",adj=0,col="white",cex=1.6,font=2)
text(.55,8.75,paste("Comparison:",if(cmp=="tested_only")"validated hits vs tested non-hits" else "validated hits vs remaining library"),adj=0,col="#E6E8F4",cex=1)
panel<-function(x0,y0,x1,y1,title) {
 rect(x0,y0,x1,y1,col="white",border="#D8DCE5",lwd=1.1)
 text(x0+.18,y1-.30,title,adj=0,col=navy,font=2,cex=1.12)
 segments(x0+.18,y1-.48,x1-.18,y1-.48,col="#E6E7ED")
}
panel(.25,6.55,4.45,8.30,"A  Study design")
panel(4.65,6.55,14.75,8.30,"B  Significant findings (Holm-adjusted)")
n_hit<-length(unique(hits$analysis_id))
n_other<-if(nrow(f))f$other_n[1] else NA
text(.55,7.50,paste0("Top UTRs: ",n_hit),adj=0,col=purple,font=2,cex=1.3)
text(.55,7.06,paste0("Comparator: ",n_other," (feature-specific n may vary)"),adj=0,col=ink,cex=.88)
text(4.95,7.58,paste0("Significant features: ",nrow(sigf)," / ",nrow(f)),adj=0,col=purple,font=2,cex=1.2)
text(4.95,7.10,paste0("Significant categories: ",nrow(sigc)," / ",nrow(c)),adj=0,col=ink,cex=1)
text(4.95,6.77,"Higher average alone does not imply statistical significance.",adj=0,col=ink,cex=.87)
panel(.25,2.50,9.15,6.30,"C  Feature comparison (mean; adjusted p-value)")
panel(9.35,2.50,14.75,6.30,"D  Significant categories")
# show up to 8 top features by Holm p, even if none significant
showf<-head(f,8)
headings<-c("Feature","Top mean","Other mean","Delta","Holm p")
xx<-c(.55,4.35,5.72,7.12,8.73)
for(j in seq_along(headings)) text(xx[j],5.50,headings[j],adj=if(j==1)0 else 1,col=navy,font=2,cex=.82)
if(nrow(showf)) for(i in seq_len(nrow(showf))) {
 z<-showf[i,]; y<-5.15-(i-1)*.34
 if(i%%2==1) rect(.43,y-.17,8.98,y+.16,col="#F5F6FA",border=NA)
 if(is.finite(z$p_holm) && z$p_holm<alpha) rect(.43,y-.17,.50,y+.16,col=purple,border=NA)
 vals<-c(substr(z$feature,1,29),fmt(z$hit_mean),fmt(z$other_mean),fmt(z$mean_difference),fmt(z$p_holm,4))
 for(j in seq_along(vals)) text(xx[j],y,vals[j],adj=if(j==1)0 else 1,col=ink,cex=.82)
}
if(nrow(sigc)) {
 showcat<-head(sigc,7)
 for(i in seq_len(nrow(showcat))) {
 z<-showcat[i,];y<-5.44-(i-1)*.46
 text(9.65,y,substr(as.character(z$category),1,24),adj=0,col=purple,font=2,cex=.86)
 text(14.4,y,sprintf("Top %d | p.adj %s",z$hit_n,fmt(z$p_holm,3)),adj=1,col=ink,cex=.82)
 }
} else text(9.65,5.3,"No categories pass adjusted p < 0.05",adj=0,col=ink,cex=.95)
panel(.25,.45,14.75,2.27,"E  Interpretation")
notes<-c(if(nrow(sigf))paste("Significant feature(s):",paste(head(sigf$feature,6),collapse=", ")) else "No feature passed multiple-testing correction.",
         if(nrow(sigc))paste("Significant category(s):",paste(head(sigc$category,5),collapse=", ")) else "No category passed multiple-testing correction.",
         if(cmp=="whole_library")"Exploratory: the remaining library was not individually tested; avoid calling it validated non-hits." else "Tested-only comparison reduces the bias of using untested library members as controls.",
         "Direction, effect size, sample size and batch effects must be considered together.")
for(i in seq_along(notes)) text(.57,1.91-(i-1)*.36,notes[i],adj=0,col=ink,cex=.91)
dev.off()
cat("Created:",file.path(out,"top_utr_summary.pdf"),"\n")
cat(paste(lines,collapse="\n"),"\n")
