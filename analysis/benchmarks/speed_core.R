args <- commandArgs(TRUE); n <- as.integer(args[1]); tag <- args[2]
options(keep.source=TRUE, keep.source.pkgs=TRUE, Ncpus=1)
data.table::setDTthreads(1)
devtools::load_all(Sys.getenv('COVERAGESIM_REPO', 'repo'), quiet=TRUE)
dest <- file.path('bench', tag); dir.create(dest, recursive=TRUE, showWarnings=FALSE)
set.seed(6618); leader <- sample(50:400,n,TRUE); cds <- sample(100:1000,n,TRUE)*3L; trailer <- sample(100:1500,n,TRUE)
step <- function(name, expr) {
  cat('START',name,format(Sys.time()),'\n'); flush.console()
  Rprof(file.path(dest,paste0(name,'.prof')),interval=.02)
  t <- system.time(value <- force(expr)); Rprof(NULL)
  write.table(data.frame(step=name,t=unname(t['elapsed'])),file.path(dest,'times.tsv'),sep='\t',row.names=FALSE,col.names=!file.exists(file.path(dest,'times.tsv')),append=TRUE,quote=FALSE)
  capture.output(summaryRprof(file.path(dest,paste0(name,'.prof'))),file=file.path(dest,paste0(name,'_summary.txt')))
  cat('END',name,t['elapsed'],'\n'); flush.console(); value
}
g <- step('genome',simGenome(n=n,out_dir=dest,architecture='all',leader_length=leader,cds_length=cds,trailer_length=trailer,cds_exons=3,debug_on=FALSE))
c <- step('counts',simCountTables(ORFik::loadRegion(g['txdb'],'cds'),libtypes=c('RFP','RNA'),conditions=c('WT','Mutant'),replicates=3,plot_PCA=FALSE,print_statistics=FALSE))
r <- step('regions',simCountTablesRegions(c))
step('coverage',simNGScoverage(g,r,out_dir=dest,exp_save_dir=dest,cores=1,ground_truth=TRUE,fragment_mode='simulated_rpf',validate=FALSE))
