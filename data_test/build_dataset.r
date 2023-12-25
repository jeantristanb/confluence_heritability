library(data.table)
Data<-fread('gctb_2.0_tutorial/gwas/sim_1.assoc.linear')
DataBim<-fread('gctb_2.0_tutorial/data/1000G_eur_chr22.bim')
system('~/bin/plink -bfile gctb_2.0_tutorial/data/1000G_eur_chr22 --keep-allele-order --freq -out 1000G_eur_chr22')
DataBim<-DataBim[,c(1,2,4,5,6)]
names(DataBim)<-c('CHR','SNP','BP','A1','A2')
Data2<-merge(Data,DataBim,by=c('CHR','SNP','BP','A1'),all=T)
DataFreq<-fread('1000G_eur_chr22.frq')
DataWithFrq<-merge(Data2,DataFreq,by=c('CHR','SNP','A1','A2'))
DataWithFrq$N<-DataWithFrq$NCHROBS/2
DataWithFrq$SE <- DataWithFrq$BETA/DataWithFrq$STAT
write.table(DataWithFrq,file='sumstat_1.tsv', row.names=F, col.names=T, sep='\t', quote=F)
