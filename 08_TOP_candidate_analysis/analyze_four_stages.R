#!/usr/bin/env Rscript
# Four-stage 5'UTR comparison. Base R only. All proprietary data remain on local server.
args <- commandArgs(trailingOnly=TRUE)
opt <- function(k,d=NULL) { i<-match(k,args); if(!is.na(i)&&i<length(args)) args[i+1] else d }
libfile<-opt("--library"); topfile<-opt("--top200","top200.csv"); hitfile<-opt("--hits","hits.txt"); out<-opt("--out","results/four_stage")
if(is.null(libfile)) stop("Provide --library /absolute/path/library.csv")
for(p in c(libfile,topfile,hitfile)) if(!file.exists(p)) stop("Missing file: ",p)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
lib<-read.csv(libfile,check.names=FALSE,stringsAsFactors=FALSE)
top<-read.csv(topfile,check.names=FALSE,stringsAsFactors=FALSE)
normalize<-function(x) toupper(gsub("\\s+","",trimws(as.character(x))))
findcol<-function(d,poss) { z<-which(tolower(gsub("[^a-z0-9]","",names(d)))%in%tolower(gsub("[^a-z0-9]","",poss))); if(length(z)) names(d)[z[1]] else NA_character_ }
if(!("library_index"%in%names(lib))) stop("Library CSV missing library_index")
lib$id<-sprintf("UTR%04d",as.integer(lib$library_index))
if(anyNA(lib$id)||anyDuplicated(lib$id)) stop("Invalid or duplicated library_index")
idcol<-findcol(top,c("UTR ID","UTR_ID","utr_id"))
rankcol<-findcol(top,c("Rank","sortseq_rank"))
if(is.na(idcol)||is.na(rankcol)) stop("Top200 needs Rank and UTR ID columns")
top$id<-normalize(top[[idcol]])
top$sortseq_rank<-suppressWarnings(as.integer(as.character(top[[rankcol]])))
if(nrow(top)!=200L||anyNA(top$sortseq_rank)||anyDuplicated(top$sortseq_rank)||!setequal(top$sortseq_rank,1:200)) stop("Top200 must contain exactly ranks 1..200, once each")
if(anyDuplicated(top$id)) stop("Duplicate UTR IDs in top200")
if(length(setdiff(top$id,lib$id))) stop("Top200 IDs absent from library: ",paste(setdiff(top$id,lib$id),collapse=", "))
hits<-unique(normalize(readLines(hitfile,warn=FALSE)));hits<-hits[nzchar(hits)&!grepl("^#",hits)]
if(length(hits)!=14) stop("Expected 14 hit IDs; got ",length(hits))
tested<-top$id[top$sortseq_rank<=60]
if(length(setdiff(hits,tested))) stop("Hit IDs outside tested rank 1..60: ",paste(setdiff(hits,tested),collapse=", "))
lib$sortseq_rank<-top$sortseq_rank[match(lib$id,top$id)]
lib$stage<-ifelse(lib$id%in%hits,"Hit 14",ifelse(lib$id%in%tested,"Tested non-hit 46",ifelse(!is.na(lib$sortseq_rank),"Top200 untested 140","Other library 1800")))
lib$in_top200<-!is.na(lib$sortseq_rank)
lib$in_tested<-lib$id%in%tested
lib$is_hit<-lib$id%in%hits
seqcol<-findcol(lib,c("utr5_sequence_tss_corrected","sequence","utr_sequence","seq"))
if(!is.na(seqcol)) {
 s<-toupper(gsub("[^ACGTU]","",as.character(lib[[seqcol]])))
 lib$sequence_length<-nchar(s)
 lib$sequence_GC_percent<-ifelse(nchar(s)>0,100*nchar(gsub("[^GC]","",s))/nchar(s),NA_real_)
 lib$sequence_uAUG<-vapply(s,function(z){ z<-chartr("T","U",z); m<-gregexpr("AUG",z,fixed=TRUE)[[1]]; if(m[1]<0) 0L else length(m) },integer(1))
}
catcol<-findcol(lib,c("library_group","category","selection_group","analysis_category"))
if(!is.na(catcol)) lib$analysis_category<-as.character(lib[[catcol]]) else lib$analysis_category<-"Category unavailable"
feature_candidates<-c("sequence_length","sequence_GC_percent","sequence_uAUG","gc_content","GC_percent","length","mfe","MFE","mfe_kcal_mol","rnafold_mfe","uaug_count","uorf_count","heavy_ensemble_score","robust_public_te_rank","model_support_score","multi_omics_utr_rank")
features<-intersect(feature_candidates,names(lib))
# Avoid duplicating derived features with source columns.
if("sequence_GC_percent"%in%features) features<-setdiff(features,c("gc_content","GC_percent"))
if("sequence_length"%in%features) features<-setdiff(features,"length")
if("sequence_uAUG"%in%features) features<-setdiff(features,"uaug_count")
write.csv(lib[,unique(c("id","sortseq_rank","stage","analysis_category",features))],file.path(out,"candidate_groups_and_features.csv"),row.names=FALSE)
groups<-list(Library_2000=rep(TRUE,nrow(lib)),Top_200=lib$in_top200,Tested_60=lib$in_tested,Hit_14=lib$is_hit)
counts<-do.call(rbind,lapply(names(groups),function(g) data.frame(group=g,n=sum(groups[[g]]))))
write.csv(counts,file.path(out,"group_counts.csv"),row.names=FALSE)
cat_table<-do.call(rbind,lapply(names(groups),function(g){t<-table(lib$analysis_category[groups[[g]]]);data.frame(group=g,category=names(t),n=as.integer(t),percent=100*as.numeric(t)/sum(t))}))
write.csv(cat_table,file.path(out,"category_composition.csv"),row.names=FALSE)
comparisons<-list(Top200_vs_remaining1800=list(a=lib$in_top200,b=!lib$in_top200),Hits14_vs_tested_nonhits46=list(a=lib$is_hit,b=lib$in_tested&!lib$is_hit),Hits14_vs_other_top200_186=list(a=lib$is_hit,b=lib$in_top200&!lib$is_hit))
results<-list()
for(nm in names(comparisons)) for(f in features){
 x<-suppressWarnings(as.numeric(lib[[f]]));a<-x[comparisons[[nm]]$a];b<-x[comparisons[[nm]]$b];a<-a[is.finite(a)];b<-b[is.finite(b)]
 if(length(a)<3||length(b)<3) next
 p<-if(length(unique(c(a,b)))>1) tryCatch(wilcox.test(a,b,exact=FALSE)$p.value,error=function(e) NA_real_) else 1
 results[[length(results)+1]]<-data.frame(comparison=nm,feature=f,n_a=length(a),n_b=length(b),median_a=median(a),median_b=median(b),mean_a=mean(a),mean_b=mean(b),median_difference=median(a)-median(b),p=p)
}
stats<-if(length(results))do.call(rbind,results) else data.frame()
if(nrow(stats)){stats$p_holm<-ave(stats$p,stats$comparison,FUN=function(z)p.adjust(z,"holm"));stats$interpretation<-ifelse(stats$p_holm<0.05,"Holm-adjusted p < 0.05","Not significant after correction")}
write.csv(stats,file.path(out,"feature_comparisons.csv"),row.names=FALSE)
cats<-list()
for(nm in names(comparisons))for(cat in sort(unique(lib$analysis_category))){
 a<-comparisons[[nm]]$a;b<-comparisons[[nm]]$b;v<-lib$analysis_category==cat
 mat<-matrix(c(sum(a&v),sum(a&!v),sum(b&v),sum(b&!v)),nrow=2,byrow=TRUE)
 z<-fisher.test(mat)
 cats[[length(cats)+1]]<-data.frame(comparison=nm,category=cat,n_a=sum(a&v),total_a=sum(a),n_b=sum(b&v),total_b=sum(b),odds_ratio=unname(z$estimate),p=z$p.value)
}
ct<-do.call(rbind,cats);ct$p_holm<-ave(ct$p,ct$comparison,FUN=function(z)p.adjust(z,"holm"))
write.csv(ct,file.path(out,"category_comparisons.csv"),row.names=FALSE)
pdf(file.path(out,"four_stage_figures.pdf"),width=11,height=7.5,onefile=TRUE)
par(mar=c(9,4,4,2))
tab<-sapply(groups,function(ix){ t<-table(factor(lib$analysis_category[ix],levels=sort(unique(lib$analysis_category))));100*as.numeric(t)/sum(t)})
rownames(tab)<-sort(unique(lib$analysis_category))
barplot(tab,col=grDevices::hcl.colors(nrow(tab),"Set 3"),ylab="Composition (%)",main="Selection-category composition across four stages",las=2,legend.text=rownames(tab),args.legend=list(x="topright",cex=.6,bty="n"))
for(f in features){
 x<-suppressWarnings(as.numeric(lib[[f]]));if(sum(is.finite(x))<10)next
 vals<-lapply(groups,function(ix)x[ix&is.finite(x)])
 par(mar=c(6,5,4,2))
 boxplot(vals,names=c("Library 2000","Top 200","Tested 60","Hits 14"),col=c("grey85","#a6cee3","#b2df8a","#fdbf6f"),main=f,ylab=f,outline=FALSE)
 stripchart(vals,vertical=TRUE,method="jitter",add=TRUE,pch=16,cex=.35,col=adjustcolor("black",alpha.f=.20))
}
par(mar=c(5,5,4,2))
plot(top$sortseq_rank,rep(0,nrow(top)),type="n",xlim=c(1,200),ylim=c(0,2),yaxt="n",xlab="Sort-seq Rank (1 = highest)",ylab="",main="Top 200 ranking and validated hits")
abline(v=60,lty=2,col="grey50")
points(lib$sortseq_rank[lib$in_tested&!lib$is_hit],rep(1,sum(lib$in_tested&!lib$is_hit)),pch=16,col="#6baed6")
points(lib$sortseq_rank[lib$is_hit],rep(1.4,sum(lib$is_hit)),pch=19,col="#e6550d",cex=1.2)
legend("topright",legend=c("Tested non-hit (46)","Validated hit (14)","Test cutoff Rank 60"),col=c("#6baed6","#e6550d","grey50"),pch=c(16,19,NA),lty=c(NA,NA,2),bty="n")
dev.off()
writeLines(c("Four-stage 5UTR analysis completed.",paste("Top200:",sum(lib$in_top200),"Tested:",sum(lib$in_tested),"Hits:",sum(lib$is_hit)),paste("Features tested:",length(features)),paste("Significant feature tests after Holm:",if(nrow(stats))sum(stats$p_holm<0.05,na.rm=TRUE) else 0),paste("Significant category tests after Holm:",sum(ct$p_holm<0.05,na.rm=TRUE)),"Top200 vs remaining 1800 and hits vs tested nonhits are primary comparisons.","Hits vs other Top200 is descriptive/exploratory: 140 are untested.","No non-significant comparison is interpreted as equivalence.","Top200 was preselected by Sort-seq; results are associative, not causal."),file.path(out,"summary.txt"))
cat("Completed:",normalizePath(out),"\n")
