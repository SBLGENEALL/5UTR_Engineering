#!/usr/bin/env Rscript
# Cleaner four-stage plots from existing analysis outputs; base R only.
args<-commandArgs(trailingOnly=TRUE)
getopt<-function(k,d){i<-match(k,args);if(!is.na(i)&&i<length(args))args[i+1] else d}
indir<-getopt("--results","results/four_stage")
datafile<-file.path(indir,"candidate_groups_and_features.csv")
if(!file.exists(datafile))stop("Missing ",datafile,"; run analyze_four_stages.R first")
d<-read.csv(datafile,check.names=FALSE,stringsAsFactors=FALSE)
group_names<-c("Library 2000","Top 200","Tested 60","Hits 14")
group_col<-c("#24345C","#1428A0","#0689D8","#00AFC1")
groups<-list(rep(TRUE,nrow(d)),!is.na(d$sortseq_rank),!is.na(d$sortseq_rank)&d$sortseq_rank<=60,d$stage=="Hit 14")
stopifnot(identical(as.integer(vapply(groups,sum,integer(1))),c(2000L,200L,60L,14L)))
out<-file.path(indir,"four_stage_figures_clean.pdf")
pdf(out,width=12,height=7.5,onefile=TRUE,pointsize=11)
op<-par(no.readonly=TRUE)
# Page 1: horizontal category composition with external legend
cats<-sort(unique(d$analysis_category[!is.na(d$analysis_category)]))
tab<-sapply(groups,function(g){t<-table(factor(d$analysis_category[g],levels=cats));if(sum(t)>0) 100*as.numeric(t)/sum(t) else rep(0,length(t))})
rownames(tab)<-cats
# stable category palette, color-blind-friendly qualitative, many distinct hues
catcols<-grDevices::colorRampPalette(c("#24345C","#1428A0","#315FC3","#0689D8","#00AFC1","#50C8D1","#8DA8D8","#7387B2","#4C6B94","#1B819B","#8AC8E8","#B7C8DF"))(length(cats))
layout(matrix(c(1,2),1,2),widths=c(3.7,1.25))
par(mar=c(5,10,5,1),xpd=FALSE)
barplot(tab,horiz=TRUE,beside=FALSE,col=catcols,border="white",space=.5,
 names.arg=group_names,las=1,xlim=c(0,100),xlab="Composition (%)",
 main="Selection category composition",cex.names=1.1,cex.main=1.3)
abline(v=seq(0,100,20),col=adjustcolor("grey75",alpha.f=.45),lty=3)
par(mar=c(4,0,5,1),xpd=NA)
plot.new()
legend("topleft",legend=cats,fill=catcols,border=NA,cex=.75,bty="n",title="Library category",ncol=1,y.intersp=1.15)
# Page 2: category differences as a heatmap; readable even with 12 categories
layout(1)
par(mar=c(9,12,5,2))
z<-tab
image(x=seq_len(ncol(z)),y=seq_len(nrow(z)),z=t(z),col=grDevices::hcl.colors(100,"YlGnBu"),
 axes=FALSE,xlab="",ylab="",main="Category frequency by selection stage")
axis(1,at=seq_len(ncol(z)),labels=group_names,las=2)
axis(2,at=seq_len(nrow(z)),labels=rownames(z),las=2,cex.axis=.75)
for(i in seq_len(nrow(z)))for(j in seq_len(ncol(z))){
 v<-z[i,j]; if(!is.finite(v)) v<-0; mx<-suppressWarnings(max(z[is.finite(z)],na.rm=TRUE)); if(!is.finite(mx)) mx<-0; text(j,i,sprintf("%.1f%%",v),cex=.65,col=if(isTRUE(v>mx*.60))"white" else "black")
}
# Features: consistent four-stage colors; label medians, no claims of significance
feature_names<-setdiff(names(d),c("id","sortseq_rank","stage","analysis_category"))
feature_names<-feature_names[vapply(d[feature_names],is.numeric,logical(1))]
for(f in feature_names){
 x<-suppressWarnings(as.numeric(d[[f]]))
 vals<-lapply(groups,function(g){v<-x[g];v[is.finite(v)]})
 if(any(vapply(vals,length,integer(1))<3))next
 par(mar=c(6.5,5.5,5,2),xpd=FALSE)
 boxplot(vals,names=group_names,col=adjustcolor(group_col,alpha.f=.75),border=group_col,
 outline=FALSE,boxwex=.55,las=1,ylab=f,main=paste("Feature distribution:",f),cex.main=1.25)
 abline(h=pretty(range(unlist(vals)),n=5),col=adjustcolor("grey75",alpha.f=.35),lty=3)
 boxplot(vals,add=TRUE,axes=FALSE,outline=FALSE,boxwex=.55,col=adjustcolor(group_col,alpha.f=.65),border=group_col)
 # use modest subsampling for display; deterministic for reproducibility
 set.seed(101)
 for(j in seq_along(vals)){
   v<-vals[[j]];if(length(v)>150)v<-sample(v,150)
   points(j+runif(length(v),-.16,.16),v,pch=16,cex=.55,col=adjustcolor(group_col[j],alpha.f=.78))
 }
 m<-vapply(vals,median,numeric(1))
 # Draw the median LAST so it remains visible above jittered points.
 for(j in seq_along(m)){
   segments(j-.275,m[j],j+.275,m[j],col="white",lwd=5)
   segments(j-.275,m[j],j+.275,m[j],col="#14213D",lwd=2.5)
   points(j,m[j],pch=21,bg="white",col="#14213D",cex=1.0,lwd=1.6)
 }
 mtext(paste(sprintf("%s median: %.3g",group_names,m),collapse="   |   "),side=1,line=4.2,cex=.72)
}
# Rank map
par(mar=c(6,5,5,2))
r<-d$sortseq_rank
plot(1,1,type="n",xlim=c(1,200),ylim=c(.3,2.2),yaxt="n",xlab="Sort-seq Rank (lower is better)",ylab="",main="Validated hits within the Top 200")
rect(1,.3,60,2.2,col=adjustcolor(group_col[3],alpha.f=.10),border=NA)
abline(v=60,lty=2,col="grey55")
points(r[!is.na(r)&r<=60&d$stage!="Hit 14"],rep(1,sum(!is.na(r)&r<=60&d$stage!="Hit 14")),pch=16,col=group_col[3],cex=.9)
points(r[d$stage=="Hit 14"],rep(1.5,sum(d$stage=="Hit 14")),pch=19,col=group_col[4],cex=1.4)
legend("topright",legend=c("Tested non-hit (46)","Validated hit (14)"),pch=16,col=group_col[3:4],bty="n")
par(op);dev.off()
cat("Created:",out,"\n")
