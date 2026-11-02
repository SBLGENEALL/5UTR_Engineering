#!/usr/bin/env Rscript
# 08 TOP candidate analysis; base R only; no model training
args <- commandArgs(trailingOnly=TRUE)
opt <- function(flag, default=NULL) { i<-match(flag,args); if(!is.na(i) && i<length(args)) args[i+1] else default }
library_file <- opt("--library")
hits_file <- opt("--hits","hits.txt")
validation_file <- opt("--validation")
outdir <- opt("--out","results")
id_mode <- opt("--id-mode","auto")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(hits_file))
read_tab <- function(p) {
  sep <- if(grepl("\\.tsv$|\\.txt$",p,ignore.case=TRUE)) "\t" else ","
  read.table(p,header=TRUE,sep=sep,quote="\"",comment.char="",check.names=FALSE,
             stringsAsFactors=FALSE,fill=TRUE)
}
if(is.null(library_file)) {
  candidates <- c("07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.csv",
                  "07_library_design/tables/selected_2000_50_100bp_cluster_diverse_evidence_balanced_library.tsv")
  candidates <- candidates[file.exists(candidates)]
  if(!length(candidates)) stop("Supply --library /path/to/final_2000_table.csv (or .tsv)")
  library_file <- candidates[1]
}
lib<-read_tab(library_file)
hit<-trimws(readLines(hits_file,warn=FALSE))
hit<-toupper(hit[nzchar(hit) & !grepl("^#",hit)])
hit<-unique(hit)
if(!length(hit)) stop("No hit IDs in hits.txt")
norm<-function(x) toupper(trimws(as.character(x)))
id_cols<-intersect(c("UTR_ID","utr_id","utr_name","variant_id","name","id"),names(lib))
# Library index -> UTRxxxx is NOT assumed unless user explicitly confirms --id-mode index
if(id_mode=="index") {
  if(!("library_index"%in%names(lib))) stop("No library_index")
  lib$analysis_id<-sprintf("UTR%04d",as.integer(lib$library_index))
} else if(length(id_cols)) {
  lib$analysis_id<-norm(lib[[id_cols[1]]])
} else stop("No explicit UTR ID column. Verify mapping, then rerun with --id-mode index.")
if(anyDuplicated(lib$analysis_id)) stop("Duplicate library UTR IDs")
missing<-setdiff(hit,lib$analysis_id)
if(length(missing)) stop(paste("Unmatched hit IDs:",paste(missing,collapse=", ")))
lib$is_hit<-lib$analysis_id%in%hit
# sequence-derived GC and length, only if absent
seqcol<-intersect(c("utr5_sequence_tss_corrected","sequence","utr_sequence","seq"),names(lib))
if(length(seqcol)) {
  seqs<-toupper(gsub("[^ACGTU]","",as.character(lib[[seqcol[1]]])))
  if(!("length"%in%names(lib))) lib$length<-nchar(seqs)
  if(!("gc_content"%in%names(lib))) lib$gc_content<-ifelse(nchar(seqs)>0,(nchar(gsub("[^GC]","",seqs))/nchar(seqs)),NA)
}
if("gc_content"%in%names(lib)) {
  x<-suppressWarnings(as.numeric(lib$gc_content))
  if(all(is.na(x)|x<=1)) x<-x*100
  lib$GC_percent<-x
}
category<-intersect(c("library_group","category","Category","selection_group"),names(lib))
if(!length(category)) warning("Category column not found")
if(length(category)) lib$analysis_category<-as.character(lib[[category[1]]])
write.csv(lib[lib$is_hit,,drop=FALSE],file.path(outdir,"hit_feature_table.csv"),row.names=FALSE)
write.csv(data.frame(UTR_ID=hit),file.path(outdir,"hit_ids_used.csv"),row.names=FALSE)
val<-NULL
if(!is.null(validation_file)) {
  val<-read_tab(validation_file)
  vid<-intersect(c("UTR_ID","utr_id","analysis_id"),names(val))
  if(!length(vid)) stop("Validation file requires UTR_ID")
  val$analysis_id<-norm(val[[vid[1]]])
  if(anyDuplicated(val$analysis_id)) stop("Validation IDs must be unique (one row per candidate)")
  unknown<-setdiff(val$analysis_id,lib$analysis_id)
  if(length(unknown)) stop(paste("Unknown validation IDs:",paste(unknown,collapse=", ")))
  val$is_hit<-val$analysis_id%in%hit
  if(!all(hit%in%val$analysis_id)) stop("Some hits absent from validation data")
  write.csv(val,file.path(outdir,"validation_checked.csv"),row.names=FALSE)
}
# Primary comparison is tested hit vs tested non-hit; secondary is hit vs entire library non-hit.
comparisons<-list(whole_library=lib)
if(!is.null(val)) comparisons$tested_only<-lib[lib$analysis_id%in%val$analysis_id,,drop=FALSE]
feature_names<-intersect(c("GC_percent","length","mfe","MFE","mfe_kcal_mol","rnafold_mfe",
                           "uaug_count","uorf_count","heavy_ensemble_score",
                           "robust_public_te_rank","model_support_score","multi_omics_utr_rank"),names(lib))
# Include additional numeric feature columns without identifiers, ranks only if numeric
extra<-names(lib)[vapply(lib,is.numeric,logical(1))]
extra<-setdiff(extra,c("library_index","is_hit"))
feature_names<-unique(c(feature_names,extra))
features<-list()
for(cmp in names(comparisons)) {
 d<-comparisons[[cmp]]
 for(f in feature_names) {
   x<-suppressWarnings(as.numeric(d[[f]]))
   a<-x[d$is_hit]; b<-x[!d$is_hit]; a<-a[is.finite(a)]; b<-b[is.finite(b)]
   if(length(a)<3 || length(b)<3 || length(unique(c(a,b)))<2) next
   p<-tryCatch(wilcox.test(a,b,exact=FALSE)$p.value,error=function(e) NA_real_)
   features[[length(features)+1]]<-data.frame(comparison=cmp,feature=f,
     hit_n=length(a),other_n=length(b),hit_mean=mean(a),other_mean=mean(b),
     hit_median=median(a),other_median=median(b),mean_difference=mean(a)-mean(b),p=p)
 }
}
ft<-if(length(features)) do.call(rbind,features) else data.frame()
if(nrow(ft)) ft$p_holm<-ave(ft$p,ft$comparison,FUN=function(x)p.adjust(x,"holm"))
write.csv(ft,file.path(outdir,"feature_tests.csv"),row.names=FALSE)
cats<-list()
if(length(category)) for(cmp in names(comparisons)) {
 d<-comparisons[[cmp]]; d<-d[!is.na(d$analysis_category),,drop=FALSE]
 for(cat in sort(unique(d$analysis_category))) {
   a<-sum(d$is_hit & d$analysis_category==cat);b<-sum(!d$is_hit & d$analysis_category==cat)
   c<-sum(d$is_hit & d$analysis_category!=cat);e<-sum(!d$is_hit & d$analysis_category!=cat)
   z<-fisher.test(matrix(c(a,b,c,e),nrow=2))
   cats[[length(cats)+1]]<-data.frame(comparison=cmp,category=cat,hit_n=a,other_n=b,
     hit_fraction=ifelse(a+c>0,a/(a+c),NA),other_fraction=ifelse(b+e>0,b/(b+e),NA),
     odds_ratio=unname(z$estimate),p=z$p.value)
 }
}
ct<-if(length(cats)) do.call(rbind,cats) else data.frame()
if(nrow(ct)) ct$p_holm<-ave(ct$p,ct$comparison,FUN=function(x)p.adjust(x,"holm"))
write.csv(ct,file.path(outdir,"category_enrichment.csv"),row.names=FALSE)
# Expression correlations only if a continuous expression column is supplied.
corrows<-list()
if(!is.null(val)) {
 expcol<-intersect(c("relative_gMFI","relative_gmfi","gMFI_ratio","gmfi_ratio","gMFI","gmfi"),names(val))
 if(length(expcol)) {
   joined<-merge(val,lib[,c("analysis_id",feature_names),drop=FALSE],by="analysis_id",all.x=TRUE,suffixes=c("",".lib"))
   y<-suppressWarnings(as.numeric(joined[[expcol[1]]]))
   for(f in feature_names) {
     if(!f%in%names(joined)) next
     x<-suppressWarnings(as.numeric(joined[[f]])); ok<-is.finite(x)&is.finite(y)
     if(sum(ok)<5 || length(unique(x[ok]))<3 || length(unique(y[ok]))<3) next
     z<-suppressWarnings(cor.test(x[ok],y[ok],method="spearman",exact=FALSE))
     corrows[[length(corrows)+1]]<-data.frame(feature=f,n=sum(ok),rho=unname(z$estimate),p=z$p.value)
   }
 }
}
cr<-if(length(corrows)) do.call(rbind,corrows) else data.frame()
if(nrow(cr)) cr$p_holm<-p.adjust(cr$p,"holm")
write.csv(cr,file.path(outdir,"expression_correlations.csv"),row.names=FALSE)
pdf(file.path(outdir,"feature_boxplots.pdf"),width=8,height=5)
for(f in feature_names) {
 x<-suppressWarnings(as.numeric(lib[[f]]))
 if(sum(is.finite(x))<10) next
 boxplot(x~factor(lib$is_hit,levels=c(FALSE,TRUE),labels=c("Other library","Selected hits")),
         main=paste(f,"(whole library)"),ylab=f,col=c("grey85","skyblue"))
 stripchart(x~factor(lib$is_hit,levels=c(FALSE,TRUE)),vertical=TRUE,method="jitter",
            pch=16,cex=.35,add=TRUE,col=rgb(0,0,0,.3))
}
dev.off()
if(length(category)) {
 pdf(file.path(outdir,"category_composition.pdf"),width=10,height=6)
 tab<-table(lib$analysis_category,lib$is_hit)
 barplot(prop.table(tab,2)*100,beside=FALSE,legend.text=TRUE,
         args.legend=list(x="topright",cex=.65),ylab="Percent",main="Library vs selected hits")
 dev.off()
}
writeLines(c(paste("Library:",normalizePath(library_file)),
             paste("Library rows:",nrow(lib)),paste("Hit IDs:",length(hit)),
             paste("Validated candidates:",if(is.null(val)) "not supplied" else nrow(val)),
             "Primary inference should use tested_only; whole_library is descriptive.",
             "Hits selected by batch top-7: results are exploratory; selection bias and batch effects remain.",
             "Missing MFE/category fields are NOT fabricated.",
             "Check ID mapping: --id-mode index only if UTRxxxx equals library_index."),
           file.path(outdir,"analysis_summary.txt"))
cat("Completed:",normalizePath(outdir),"\n")
