options(stringsAsFactors=FALSE)
TAB <- "/private/tmp/claude-502/-Users-ls1020/cdc231a7-d353-43cf-ac37-fd155cce1deb/scratchpad/tabs3"
dir.create(TAB, showWarnings=FALSE, recursive=TRUE)
sc <- readRDS("results/iter003/combined_scenario_metrics.rds")
sc$ld <- factor(ifelse(is.na(sc$n_ref),"in-sample",paste0("ref",sc$n_ref)),
                levels=c("in-sample","ref500","ref200"))
CORE <- c("fb_xregion","fb_pooled","functional_beatrice","beatrice","finemap","susie","polyfun_est","polyfun_oracle")
d <- sc[sc$method %in% CORE & sc$annotation_type %in% c("binary","continuous"),]
med <- function(df,m,by){a<-aggregate(df[[m]],df[by],function(x) median(x,na.rm=TRUE));names(a)[ncol(a)]<-m;a}
W <- function(x,f) write.csv(x, file.path(TAB,f), row.names=FALSE)
tab <- function(metric, dp=3){
  a <- med(d, metric, c("annotation_type","ld","method"))
  w <- reshape(a, idvar=c("annotation_type","ld"), timevar="method", direction="wide")
  names(w) <- gsub(paste0(metric,"."),"",names(w),fixed=TRUE)
  keep <- intersect(CORE, names(w)); w[keep] <- lapply(w[keep], function(z) round(z,dp))
  w[order(w$annotation_type,w$ld), c("annotation_type","ld",keep)] }
W(tab("total_mass_ratio",2),   "t_mass.csv")
W(tab("hi_pip_reliab",3),      "t_reliab.csv")
W(tab("max_fdr_violation_n20",3), "t_fdr.csv")
W(tab("ap",3),                 "t_auprc.csv")

# vs finemap, paired
key <- c("job_dir","S","phi","region_size")
pw <- reshape(d[,c(key,"method","ap","max_fdr_violation_n20")], idvar=key, timevar="method", direction="wide")
meta <- unique(sc[,c("job_dir","annotation_type","n_ref")]); pw <- merge(pw, meta, by="job_dir")
pw$ld <- factor(ifelse(is.na(pw$n_ref),"in-sample",paste0("ref",pw$n_ref)),
                levels=c("in-sample","ref500","ref200"))
out <- do.call(rbind, lapply(c("fb_xregion","fb_pooled","functional_beatrice","beatrice"), function(j){
  z <- aggregate(cbind(a=pw[[paste0("ap.",j)]] > pw[["ap.finemap"]],
                       f=pw[[paste0("max_fdr_violation_n20.",j)]] < pw[["max_fdr_violation_n20.finemap"]]) ~ ld,
                 pw, function(x) round(100*mean(x,na.rm=TRUE),1))
  data.frame(method=j, ld=z$ld, beats_finemap_AUPRC=z$a, beats_finemap_FDR=z$f) }))
W(out, "t_vs_finemap.csv")

# by S (ref500)
a <- med(d[d$ld=="ref500",], "ap", c("S","method"))
w <- reshape(a, idvar="S", timevar="method", direction="wide"); names(w) <- gsub("ap.","",names(w),fixed=TRUE)
keep <- intersect(CORE,names(w)); w[keep] <- lapply(w[keep], function(z) round(z,3))
W(w[order(as.numeric(w$S)), c("S",keep)], "t_byS_auprc.csv")
a <- med(d[d$ld=="ref500",], "total_mass_ratio", c("S","method"))
w <- reshape(a, idvar="S", timevar="method", direction="wide"); names(w) <- gsub("total_mass_ratio.","",names(w),fixed=TRUE)
keep <- intersect(CORE,names(w)); w[keep] <- lapply(w[keep], function(z) round(z,2))
W(w[order(as.numeric(w$S)), c("S",keep)], "t_byS_mass.csv")

# enrichment slopes
E <- d[!is.na(d$enrichment_fold),]
a <- med(E,"ap",c("annotation_type","ld","enrichment_fold","method"))
sl <- do.call(rbind, lapply(split(a, list(a$method,a$annotation_type,a$ld), drop=TRUE), function(z){
  if (nrow(z)<3) return(NULL)
  data.frame(method=as.character(z$method[1]), annotation=z$annotation_type[1], ld=as.character(z$ld[1]),
             slope=round(unname(coef(lm(ap~enrichment_fold,z))[2]),4))}))
w <- reshape(sl, idvar=c("method","annotation"), timevar="ld", direction="wide")
names(w) <- gsub("slope.","",names(w),fixed=TRUE)
W(w[order(w$annotation,-w$`in-sample`),], "t_enrichment.csv")
cat("tables written\n"); print(list.files(TAB))
