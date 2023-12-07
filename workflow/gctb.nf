include {format_sumstat} from '../process/gctb.nf'
include {run_gctb_sumstat} from '../process/gctb.nf'
include {interpolate_ld} from '../process/gctb.nf'
include {errormess} from '../process/utils.nf'
include {extractindfile} from '../process/plink.nf'
include {shuffle} from '../process/plink.nf'
include {clean_plink} from '../process/plink.nf'
include {computefreq_plink} from '../process/plink.nf'
include {shrunkld} from '../process/gctb.nf'
include {sparceld_aftshrunk} from '../process/gctb.nf'

workflow build_ld {
   if(params.bfile==""){
      errormess("bfile must be initialised to build ld using gctb")
  }  
   bfile=params.bfile
   bedfileI=Channel.fromPath("${bfile}.bed",checkIfExists:true).combine(Channel.fromPath("${bfile}.bim",checkIfExists:true)).combine(Channel.fromPath("${bfile}.fam",checkIfExists:true))
   if(params.gctb_keepind!='')extractindfile(Channel.fromPath(params.gctb_keepind,checkIfExists:true), channel.from("${params.output_pat}.ind"))
   else extractindfile(Channel.fromPath("${bfile}.bim",checkIfExists:true), channel.from("${params.output_pat}.ind"))
   //    tuple path(bed), path(bim), path(fam), path(list_ind), val(output)
   clean_plink(bedfileI, channel.fromPath("01"),  channel.fromPath("02"), extractindfile.out, channel.from("${params.output_dir}/cleanplk"), channel.from("${params.output_pat}"))
   if(params.plink_shuffle>0){
     shuffle(clean_plink.out.all, channel.from(params.plink_shuffle), channel.from("${params.output_dir}/cleanplk"), channel.from("${params.output_pat}"))
     allplk=shuffle.out.all
     bimplk=shuffle.out.bim
   }else{
     allplk=clean_plink.out.all
     bimplk=clean_plink.out.bim
   }
   computefreq_plink(allplk, extractindfile.out,channel.from("${params.output_dir}/freq"),channel.from(params.output_pat+"_frq"))
   map_ch=Channel.fromPath(params.list_map, checkIfExists:true).splitCsv(header: true, sep: ',')
        .map{row ->
            def chro= row['chro']
            def file= row['file']
            return [ chro, file]
        }
  interpolate_ld(map_ch.combine(bimplk).combine(Channel.from("${params.output_dir}/interpolate_ld")).combine(channel.from(params.output_pat+'')))
 shrunkld (interpolate_ld.out.combine(allplk).combine(Channel.from("${params.output_dir}/ld/shrunk/")).combine(channel.from(params.output_pat)))
 sparceld_aftshrunk(shrunkld.out.res.combine(Channel.from("${params.output_dir}/ld/sparce/")).combine(channel.from(params.output_pat+'_sparse')))

}


workflow gctb {
  println "other option : "+ params.gctb_otheroption
  sumstat=channel.fromPath(params.sumstat, checkIfExists:true)
  if(params.update_rsid!='')rsidfile=channel.fromPath(params.update_rsid, checkIfExists:true)
  format_sumstat(sumstat.combine(rsidfile).combine(channel.from("${params.output_dir}/sumstat_format")).combine(channel.from(params.output_pat+'.sumstat')))
  listfile_gtc_bin=Channel.fromPath(file(params.gctb_ld_bin, checkIfExists:true).readLines(), checkIfExists:true).collect()
  listfile_gtc_info=Channel.fromPath(file(params.gctb_ld_info, checkIfExists:true).readLines(), checkIfExists:true).collect()
  run_gctb_sumstat(format_sumstat.out, listfile_gtc_bin, listfile_gtc_info, channel.from(params.output_pat))

}
