# Run GCTB to estimate heritability
## run :
* used by default singularity image see [dockerhub](https://hub.docker.com/repository/docker/jeantristanb/confluence) and [github](https://github.com/jeantristanb/confluence_docker)
* profile :
 * slurmSingularity
 * batch
 * docker
 * dnanexus 

## model 
* option to defined action of script : `--model` :
  * `gctb` : run gctb to estimated ld [default]
  * `gctb_ld` : build ld matrix
## general option 
* general option
 * `output_pat` : output of pattern
 * `output_dir` : directory outpuyt
 * `gctb_bin` [gctb]
 * `memory_gctb` : memory used gctb [20.Gb]
 
## heritability using sumstat
 * `sumstat` : sumstat file can  be gzip
   * `sumstat_head_chr` : chromosome header
   * `sumstat_head_bp` : bp header
   * `sumstat_head_a1` : effect allele header
   * `sumstat_head_a2` : non effect allele header
   * `sumstat_head_rs` : rs header
   * `sumstat_head_pval` : pval header
   * `sumstat_head_se` : standard error header
   * `sumstat_head_beta` : beta header
   * `sumstat_head_freq` : frequences
   * nvalue can be done as argument or header of sumstat:
    * `sumstat_head_n` : one or n header, separated by header
    * `sumstat_n` :  n values [int]
 * `gctb_impute_n` : option to impute n, see `--impute-n` in gctb [default 0]
 * `gctb_exclude_mhc` : option to impute n, see `--exclude-mhc` in gctb [default 0]
 * matrice ld :
  * `gctb_ld_bin` :  file contained list of ld build with gctb, extension bin [default : ""]
  * `gctb_ld_bin_dir` : directory containing bin file [default : "" ]
 * info ld :
  * `gctb_ld_info` :  file contained list of info of `gctb_ld_bin` build with gctb, extension bin  [default : ""]
  * `gctb_ld_info_dir` :   directory containing info file of bin file [default : "" ]
 * model [see more information](https://cnsgenomics.com/software/gctb/#Bayesianalphabet) :
  * `gctb_bayesmod` : model S, R [default S], tested just with S
  * `gctb_hsqinit` : heritability (see option `--hsq` of gctb)
 * `update_rsid` : file contained positions to update, file one position by line with chr bp rs, no header

### ressource 
 * ld using 50kb of 2.3 mb: see [gctb website](https://cnsgenomics.com/software/gctb/#LDmatrices) or [zenodo](https://zenodo.org/records/3375373#.XyFgOS17G8o)


## run gctb to compute heritabilty

### Created file for matrice
```
ls  ../ressource/ukb_50k_bigset_2.8M/*.bin > list_bin
ls  ../ressource/ukb_50k_bigset_2.8M/*.info > list_info
```

### run gctb
```
#ID CHR POS ALT REF AF_ALT beta.ALT SE P Rsq N_case N_control
hchr=CHR;hbp=POS;ha1=ALT;ha2=REF;hrs=ID;hp=P;hbeta=beta.ALT;hse=SE;hn="N_case,N_control";hfreq=AF_ALT
#gctb_listld  
confluence_heritability/main.nf  --sumstat $sumstat   --sumstat_head_chr $hchr --sumstat_head_bp $hbp --sumstat_head_a1 $ha1 --sumstat_head_a2 $ha2 --sumstat_head_rs $hrs --sumstat_head_pval $hp --sumstat_head_beta $hbeta --sumstat_head_se $hse --sumstat_head_n $hn --sumstat_head_freq $hfreq   --gctb_ld_bin list_bin --gctb_ld_info list_info -resume -profile slurmSingularity
```

## build LD ressourxce using plink file as input
* see ressource (here)[https://github.com/jeantristanb/confluence/tree/main/buildld/afr_1000g]
build ld used for heritability, bin info, clean SNPs
* plink option and clean :
 * `bfile` : plink without ext [ "" ]
 * `cpu_plink` : cpu number for plink [ 2 ]
 * `memory_plink` :  memory of plink [10.Gb]
* clean plink file
 * `plink_maf` :  clean minor allele frequence, val 0 - 0.5 [ .01]
 * `plink_geno` : missing for genotype, ratio val 0 - 1 [-1]
 * `plink_hwe` : p-value for hardy wendberg equilibrium by position[-1]
 * `gctb_keepind` : individual to extract from plink and analyse  [""]
 * `plink_shuffle` : number positions randomly choose in plink [-1]
* genome reference :
 * `list_map` : csv file contained file with genetics map, header chro,file

### example run pipeline 
 * using `bfile` from 1000 Genome vcf hg19 in plink
 * extract individual without relatdness using information

###Created file without inbredding
see :  `confluence/buildld/afr_1000g/extract_ind_afr_v2.r`

```
fileind=confluence/buildld/afr_1000g/list_afr_noappar.ind
```

### genomic reference

```
wget -c ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/ASW_omni_recombination_20130507.tar
tar -xf ASW_omni_recombination_20130507.tar
```

```
echo "chro,file" > list_map.csv
ls confluence/buildld/afr_1000g/ASW/*.gz|awk -F "-" '{print $(NF-1)","$0}' >> list_map.csv
```

### run pipeline
 * see [data\_test](./data_test/) folder

```
sumstat=$Dirgit/confluence_heritability/data_test/sumstat_1.tsv
output=testnf
hchr=CHR;hbp=BP;ha1=A1;ha2=A2;hrs=SNP;hrs=SNP;hp=P;hbeta=BETA;hse=SE;hfreq=MAF;hn=N
listbin=listbin
listinfo=listinfo

ls $Dirgit//confluence_heritability/data_test/chr22/*.bin > $listbin
ls $Dirgit/confluence_heritability/data_test/chr22/*.info > $listinfo

echo "$Dirgit/confluence_heritability/main.nf  --sumstat $sumstat   --sumstat_head_chr $hchr --sumstat_head_bp $hbp --sumstat_head_a1 $ha1 --sumstat_head_a2 $ha2 --sumstat_head_rs $hrs --sumstat_head_pval $hp --sumstat_head_beta $hbeta --sumstat_head_se $hse --sumstat_head_freq $hfreq    --gctb_ld_bin $listbin --gctb_ld_info $listinfo -resume -profile slurmSingularity --memory_gctb "130GB"   --output_dir $output --output_pat $output --sumstat_head_n $hn"
```


## ressource :
*  `interpolate_ld` : script modify of TODO

# run on dna-nexus

see manual [here](https://documentation.dnanexus.com/user/running-apps-and-workflows/running-nextflow-pipelines)


## install and prepared folder
```
# install dx command line
#pip3 install dxpy --user
#connection
dx login or dx login --token xxxxx
#create a directory 
dx mkdir Aim2_Polygenicity/jeantristan
#what file are in one directory 
dx ls Aim2_Polygenicity
## upload data_test in polygenicity
#dx upload -r $dirgit/confluence_heritability/data_test --destination ./Aim2_Polygenicity/jeantristan/
```
### build nextflow
* command line take 1 to 2 minutes : 
* when command finish you will received a confirmation emel
* command line return applet project reference to used on the next steps (otherwise you can find in folder where you have built nextflow applet)
```
dx login 
#deleted a previous workflow 
dx rm Aim2_Polygenicity/jeantristan/heritability/confluence_heritability
# @e build with profile dnanexus contained singularity image
dx build --nextflow   --repository https://github.com/jeantristanb/confluence_heritability   --destination Aim2_Polygenicity/jeantristan/heritability --profile  dnanexus
```

### run example 

```
dx login --token XXX
# variable to run pipeline
# general direction where are data using project name (found on dnanexus)
dirdata=dx://project-GYBZ3yj4BzFgB61K7vX9Fq9K:./Aim2_Polygenicity/jeantristan/
# dir ld with bin file and info
dirbin=$dirdata/data_test/chr22/
# sumstat
sumstat=$dirdata/data_test/sumstat_1.tsv
#header of sumstat
hchr=CHR;hbp=BP;ha1=A1;ha2=A2;hrs=SNP;hrs=SNP;hp=P;hbeta=BETA;hse=SE;hfreq=MAF;hn=N
#Output of result
outputdir=$dirdata/out_data_test
# output pattern
output=test

# run applet, using project id and applet.
dx run project-GYBZ3yj4BzFgB61K7vX9Fq9K:applet-Gf5Vp4842p63vz54544qyKP2 -i nextflow_pipeline_params="--sumstat $sumstat   --sumstat_head_chr $hchr --sumstat_head_bp $hbp --sumstat_head_a1 $ha1 --sumstat_head_a2 $ha2 --sumstat_head_rs $hrs --sumstat_head_pval $hp --sumstat_head_beta $hbeta --sumstat_head_se $hse --sumstat_head_freq $hfreq    --gctb_ld_bin_dir $dirbin --gctb_ld_info_dir $dirbin -resume  --memory_gctb 130GB   --output_dir $outputdir --output_pat $output --sumstat_head_n $hn"
```
