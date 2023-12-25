Dirgit=/home/jeantristan/Travail/git
sumstat=$Dirgit/confluence_heritability/data_test/sumstat_1.tsv
output=testnf
hchr=CHR;hbp=BP;ha1=A1;ha2=A2;hrs=SNP;hrs=SNP;hp=P;hbeta=BETA;hse=SE;hfreq=MAF;hn=N
listbin=listbin
listinfo=listinfo

ls $Dirgit//confluence_heritability/data_test/chr22/*.bin > $listbin
ls $Dirgit/confluence_heritability/data_test/chr22/*.info > $listinfo

echo "$Dirgit/confluence_heritability/main.nf  --sumstat $sumstat   --sumstat_head_chr $hchr --sumstat_head_bp $hbp --sumstat_head_a1 $ha1 --sumstat_head_a2 $ha2 --sumstat_head_rs $hrs --sumstat_head_pval $hp --sumstat_head_beta $hbeta --sumstat_head_se $hse --sumstat_head_freq $hfreq    --gctb_ld_bin $listbin --gctb_ld_info $listinfo -resume -profile slurmSingularity --memory_gctb "130GB"   --output_dir $output --output_pat $output --sumstat_head_n $hn"



