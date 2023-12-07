#!/usr/bin/env python3
import argparse
from math import sqrt
import sys
from scipy import stats
import gzip
import os
EOL=chr(10)

def GetSep(Sep):
   ListOfSep=["\t"," ",","]
   if len(Sep)>2 :
      Sep=Sep.upper()[:3]
      if Sep=='COM' :
         Sep=','
      elif Sep=='TAB' :
         Sep='\t'
      elif Sep=='WHI' :
         Sep=' '
   if Sep not in ListOfSep :
      return None
   return Sep


def formatrs(chro, bp,a1,a2) :
  a1=a1.upper()
  a2=a2.upper()
  if a1 > a2 :
     AA1=a1
     AA2=a2
  else :
     AA1=a2
     AA2=a1
  newrs=chro+'_'+bp+'_'+AA1+'_'+AA2
  return(newrs)

def spl_gz(z,sep) :
   return z.decode('utf-8').replace('\n','').split(sep)

def spl_nogz(z,sep) :
   return z.replace('\n','').split(sep)



def openf(File) :
  def is_gz_file(filepath):
      with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'
  def checkexists(path_to_file) :
    if exists(path_to_file)==False :
     sys.exist('file '+ path_to_file+' doesn t exist')
  balisegz=False
  if is_gz_file(File) :
    readf=gzip.open(File)
    balisegz=True
  else :
    readf=open(File)
  spl=spl_nogz
  if balisegz:
     spl=spl_gz
  return (readf, balisegz,spl)


def extractinfobim(bimfile):
   readbim=open(bimfile)
   dicrskey={}
   dickeyrs={}
   listkey=set([])
   listchrbp=set([])
   listdup=set([])
   for line in readbim :
      splline=line.replace('\n','').split()
      chro=splline[0]
      rs=splline[1]
      bp=splline[3]
      a1=splline[4]
      a2=splline[5]
      keyposalt=formatrs(chro, bp,a1,a2)
      keychr=chro+':'+bp
      if keychr in listchrbp :
         listdup.add(keychr)
      else :
         listchrbp.add(keychr)
      listkey.add(keyposalt) 
      dicrskey[rs]=keyposalt
      dickeyrs[keyposalt]=rs
   return (listkey, listdup, dicrskey,dickeyrs)

def readsummarystat(sumstat, dickeyrs,listdup, listkey, rs_header,n_header, chr_header, bp_header, beta_header, z_header, a1_header, a2_header, se_header, af_header, p_header,n_value, maf, nlim, keep_genet,sep, exclude_pos) :
    def getposheader( head, headersumstat):
        nheader=None
        if head is not None:
          head=head.lower()
          try :
           nheader=headersumstat.index(head)
          except :
            print(head+" not found in "+" ".join(headersumstat))
            sys.exit(2)
        print('head '+str(head)+' found in column '+str(nheader))
        return nheader
    def gv(infogwas, posheader, othervalue="NA", formatfct=str) :
       if posheader is not None:
          try :
             val=infogwas[posheader]
          except :
            print(infogwas, posheader, othervalue, formatfct)
            sys.exit(2)
          try :
            return formatfct(val)
          except :
            return othervalue 
       return othervalue
    (readsummstat, balisegz,spl)=openf(sumstat)
    headersumstat=[x.lower() for x in spl(readsummstat.readline(), sep)]
    n_headerp=getposheader(n_header, headersumstat)
    se_headerp=getposheader(se_header, headersumstat)
    a2_headerp=getposheader(a2_header, headersumstat)
    a1_headerp=getposheader(a1_header, headersumstat)
    z_headerp=getposheader(z_header, headersumstat)
    beta_headerp=getposheader(beta_header, headersumstat)
    bp_headerp=getposheader(bp_header, headersumstat)
    chr_headerp=getposheader(chr_header, headersumstat)
    rs_headerp=getposheader(rs_header, headersumstat)
    af_headerp=getposheader(af_header, headersumstat)
    p_headerp=getposheader(p_header, headersumstat)
    ncol=len(headersumstat)
    if maf is not None:
      mafupper= 1 - maf
    if n_value is None:
       n_value_default="NA"
    dicres={}
    listnewkey=set([])
    balisen=(nlim  is not None) and (n_headerp is not None)
    for line in readsummstat :
        splline=spl(line, sep)
        if len(splline)!=ncol :
           print(splline)
           print("column number different between header and line\n")
           continue
        key=formatrs(splline[chr_headerp],splline[bp_headerp],splline[a1_headerp],splline[a2_headerp])
        key2=splline[chr_headerp]+':'+splline[bp_headerp]
        if key2 in exclude_pos :
          continue
        balisersbim=(dickeyrs is not None) and (key in dickeyrs)
        if (key2 not in listdup)  and (keep_genet==0 or balisersbim):
           freq=gv(splline,af_headerp)
           balise=True
           if freq!="NA" and (maf is not None):
             freq=float(freq)
             balise=freq>=maf and freq<=mafupper
           if balisen and balise:
             Nval=float(splline[n_headerp])
             balise=balise and Nval > nlim
           if balise :
             rsbim="NA"
             if balisersbim==False:
                keybim='NA'
             else :
                 keybim=dickeyrs[key]
             dicres[key]=[keybim,gv(splline, rs_headerp, key), gv(splline, chr_headerp), gv(splline, bp_headerp), gv(splline, a1_headerp).upper(), gv(splline, a2_headerp).upper(), gv(splline, z_headerp,formatfct=float), gv(splline, beta_headerp,formatfct=float), gv(splline, se_headerp,formatfct=float), gv(splline,af_headerp,formatfct=float), gv(splline, n_headerp, n_value_default,formatfct=float), gv(splline, p_headerp,formatfct=float)]  
             listnewkey.add(key)
    return (dicres, listnewkey)

def ExtractFreqN(bfile, dicsumstat, listkey,dicbim,freq_header,rs_header, n_header, chr_header,bp_header, bin_plk, keep, threads, memory, nvalue, maf) :
   if (n_header or nvalue) and (freq_header is not None):
       return (dicsumstat, listkey)
   if bfile == None :
      return (dicsumstat, listkey)
   plkfreqfil=os.path.basename(bfile)
   out_range="tmp.frq"
   if bfile is None :
     print("no header for n or freq and bfile")
     return (dicsumstat, listkey)
     sys.exit()
   Cmd=bin_plk+" -bfile "+bfile+" --freq --keep-allele-order "
   rangelist=out_range+".bed"
   writebed=open(rangelist, 'w') 
   [writebed.write("\t".join([dicsumstat[key][2],dicsumstat[key][3], dicsumstat[key][3], dicsumstat[key][0]])+'\n') for key in listkey]
   writebed.close()
   if keep is not None:
     Cmd+=" --keep "+keep
   Cmd+=" --threads "+str(threads)
   Cmd+=" --extract range "+rangelist
   Cmd+=" --memory "+str(memory)
   Cmd+=" --out "+plkfreqfil
   os.system(Cmd)
   #data_n=pd.read_csv(plkfreqfil+".frq",delim_whitespace=True)
   readfreq=open(plkfreqfil+".frq")
   #CHR	Chromosome code
   #SNP	Variant identifier
   #A1	Allele 1 (usually minor)
   #A2	Allele 2 (usually major)
  #MAF	Allele 1 frequency
  #NCHROBS	Number of allele observations
   header=readfreq.readline().replace('\n', '').split()
   # CHR               SNP   A1   A2          MAF  NCHROBS
   posrs=1
   posa1=2
   posa2=3
   posaf=4
   posn=5
   if maf is not None:
      mafupper= 1 - maf
   newlistkey=set([])
   for line in readfreq :
      splline=line.replace('\n', '').split()
      key=dicbim[splline[posrs]]
      # 0 SNP1, 1 SNP2, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
      if not freq_header :
        af=float(splline[posaf])
        balise=True
        if maf is not None:
          balise=af>=maf and af<=mafupper
        if splline[posa1]==dicsumstat[key][5] and splline[posa2]==dicsumstat[key][4] :
         af = 1 - af
        if balise :
         dicsumstat[key][9]=af
         newlistkey.add(key)
        else :
         print("key "+key)
         del dicsumstat[key]
      else : 
         newlistkey.add(key)
      if not n_header : 
         dicsumstat[key][10]=int(splline[posn])/2
   return (dicsumstat,newlistkey)

def addz(dicsumstat, listkey ,z_header) :
 # 0 SNP1, 1 SNP2, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
  if z_header :
     return dicsumstat
  for key in listkey :
    if dicsumstat[key][7]!='NA' and dicsumstat[key][8]!='NA' :
      dicsumstat[key][6]=float(dicsumstat[key][7])/float(dicsumstat[key][8])
    else :
      dicsumstat[key][6]='NA'
  return dicsumstat
 
def addbetase(dicsumstat, listkey ,beta_header, se_header, z_header) :
   def fn(x, fct) :
       try :
         fct(x)
       except :
         return 'NA'
   if beta_header and se_header :
       return dicsumstat
   if not z_header :
     print('no z header and beta and se header')
     return dicsumstat
   for key in listkey :
    z, p, n=fn(dicsumstat[key][6], float),fn(dicsumstat[key][9], float),fn(dicsumstat[key][10], int)
    if z!='NA' and p!='NA' and n!='NA':
     b=z/sqrt(2*p*(1-p)*(n+z**2))
     se=1/sqrt(2*p*(1-p)*(n+z**2))
     dicsumstat[key][7]=b
     dicsumstat[key][8]=se
    else :
      del dicsumstat[key]
      print("pos "+key+" : z "+str(z)+" p : "+str(p)+" n "+str(n)+" is na")
   return(dicsumstat)

def compute_usingp(dicsumstat, listkey,used_p, beta_header, se_header, z_header) :
  if used_p==0 :
    return(dicsumstat)
  for key in listkey :
    if dicsumstat[key][7]!='NA':
      if float(dicsumstat[key][7])> 0 :
        sens=1
      else :
        sens = -1 
    elif dicsumstat[key][6]!='NA' :
      if float(dicsumstat[key][6]) > 0 :
        sens = 1
      else :
        sens = -1
    else :
      print("oos "+ key+" z and beta NA")
      continue
    dicsumstat[key][6]=abs(stats.norm.ppf(1-dicsumstat[key][11]/2))*sens
    z, p, n=float(dicsumstat[key][6]),float(dicsumstat[key][9]),int(dicsumstat[key][10])
    b=z/sqrt(2*p*(1-p)*(n+z**2))
    se=1/sqrt(2*p*(1-p)*(n+z**2))
    dicsumstat[key][7]=b
    dicsumstat[key][8]=se
  return dicsumstat

def checkrs(dicsumstat, listkey) :
   for key in listkey :
      #0 SNP1, 1 SNP2, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
      if dicsumstat[key][0]=='NA' or (not dicsumstat[key][0]) or dicsumstat[key][0]=="":
         dicsumstat[key][0]=dicsumstat[key][1]
   return dicsumstat





def parseArguments():
    parser = argparse.ArgumentParser(description='extract rs in gwas file')
    parser.add_argument('--inp_resgwas',type=str,required=True)
    parser.add_argument('--rs',type=str,required=False, help="list rs", default="")
    parser.add_argument('--chro_header',type=str,required=False,help="chro header in inp files", default=None)
    parser.add_argument('--n_header',type=str,required=False,help="chro header in inp files", default=None)
    parser.add_argument('--pos_header',type=str,required=False,help="pos header in inp files", default=None)
    parser.add_argument('--beta_header',type=str,required=False,help="beta header in inp files", default=None)
    parser.add_argument('--z_header',type=str,required=False,help="Z header in inp files",default=None)
    parser.add_argument('--se_header',type=str,required=False,help="se header in inp files",default=None)
    parser.add_argument('--p_header',type=str,required=False,help="p-value header in inp files", default=None)
    parser.add_argument('--rs_header',type=str,required=False,help="beta header in inp files", default=None)
    parser.add_argument('--n_lim',type=int,required=False,help="beta header in inp files",default=None)
    parser.add_argument('--a1_header',type=str,required=False,help="beta header in inp files", default=None)
    parser.add_argument('--keep_genet',type=str,help="beta header in inp files", default=0)
    parser.add_argument('--a2_header',type=str,required=False,help="beta header in inp files", default=None)
    parser.add_argument('--freq_header',type=str,required=False,help="frequencies header in inp files", default=None)
    parser.add_argument('--maf',type=float,default=0.0,help="minor allele frequencies")
    parser.add_argument('--bin_plk',type=str,required=False,help="plink binary", default="plink")
    parser.add_argument('--exclude_pos',type=str,required=False,help="exclude position", default=None)
    parser.add_argument('--bfile',type=str,required=False,help="bfile if need to compute frequency or N", default=None)
    parser.add_argument('--keep',type=str,required=False,help="file of data used for if need to compute frequency or N", default=None)
    parser.add_argument('--used_p',type=int,required=False,help="file of data used for if need to compute frequency or N", default=0)
    parser.add_argument('--threads',type=int,required=False,help="", default=1)
    parser.add_argument('--memory',type=int,required=False,help="", default="1000")
    parser.add_argument('--n',required=False, help="bim file ")
    parser.add_argument('--info_file',required=False, default="")
    parser.add_argument('--out_head',type=str,default="out",help="out header")
    args = parser.parse_args()
    return args


def check_args(args):
   def check1args(arg,head, need=False):
      head=head.lower()
      if head in l_infohead :
         headtmp=l_filehead[l_infohead.index(head)].replace(' ','')
         if len(headtmp)> 0:
            print(headtmp,head)
            return headtmp
      if arg :
        print(arg,head)
        return arg
      if need :    
        print("args "+ head +" not found") 
        exit(2)
      print("None : ",head)
      return None
   infofile=args.info_file
   infohead=infofile.split(",")
   infohead=[x for x in infohead if len(x)>2]
   if len(infohead) >1 :
     l_infohead=[x.split(":")[0].lower() for x in infohead if len(x.split(":"))==2]
     l_filehead=[x.split(":")[1] for x in infohead if len(x.split(":"))==2]
   else : 
     l_infohead=[]
     l_filehead=[] 
   return (check1args(args.rs_header, 'rsID'), check1args(args.n_header, 'N'),check1args(args.chro_header, 'Chro'), check1args(args.pos_header, 'Pos'), check1args(args.beta_header,'Beta'), check1args(args.z_header,'Z'), check1args(args.a1_header,'A1', True), check1args(args.a2_header, 'A2', True), check1args(args.se_header,'Se'), check1args(args.freq_header,'freqA1'), check1args(args.p_header, 'Pval'))
    

def check_args(args):
#chr     pos     snpid   effect_allele   other_allele    eaf     ntotal  beta    se      p_value rsq
#MarkerName	Chromosome	Position	EA	NEA	EAF	Nsample	Ncohort	Effects	beta_0	se_0	beta_1	se_1	beta_2	se_2	chisq_association	ndf_association	P-value_association	chisq_ancestry_het	ndf_ancestry_het	P-value_ancestry_het	chisq_residual_het	ndf_residual_het	P-value_residual_het	lnBF	Comments
   listp=["p_wald","p-value", 'p', 'p_value', 'pvalue_fe','p-value_association']
   listchr=['chr','chro', 'chromosome']
   listpos=['pos_b37', 'ps', 'bp', 'pos', 'position']
   listse=['se', 'stderr', 'std_re']
   listn=['n', 'nvalue', 'n_val', 'n_value', 'n_total_sum', 'ntotal','n_total']
   lista0=['allele0','allele2', 'a2', 'other_allele','nea']
   lista1=['allele1', 'a1', 'effect_allele', 'ea']
   listz=['z']
   listbeta=['beta', 'B', 'effect', "beta_fe",'b']
   listrs=['rs', 'rsid', 'snp', 'markername', 'snpid']
   listaf=['af', 'frq', 'freq1', 'effect_allele_freq', 'eaf',"af_weight", 'freq']
   def check1args(arg,header,l_filehead,l_infohead,infohead,listhead_check, need) :
    headfinal=None
    l_infohead=[x.lower() for x in l_infohead]
    infohead=infohead.lower()
    #print(l_infohead, infohead)
    if arg :
        headfinal=arg
    elif infohead in l_infohead :
      headfinal=l_filehead[l_infohead.index(infohead)]
    else :
      head=[x for x in header if x in listhead_check]
      if len(head)>1:
        print(" ".join(listhead_check )+"\n"+"\t".join(listhead_check ))
        print(head)
        print("more than one header "+infohead )
        if need :
          sys.exit(3)
        else :
          head=[head[1]]
      if len(head) ==0 :
        if need :
          print('not found '+infohead)
          sys.exit(2)
        headfinal=None
      else :
        headfinal=head[0]
    print(infohead+' used head '+ str(headfinal))
    return(headfinal)
   infofile=args.info_file
   infohead=infofile.split(",")
   infohead=[x for x in infohead if len(x)>2]
   if len(infohead) >1 :
     l_infohead=[x.split(":")[0].lower() for x in infohead if len(x.split(":"))==2]
     l_filehead=[x.split(":")[1] for x in infohead if len(x.split(":"))==2]
   else : 
     l_infohead=[]
     l_filehead=[]
   (readhead,balisegz, spl)=openf(args.inp_resgwas)
   line10i=[readhead.readline() for x in range(0,10)]
   readhead.close() 
   if balisegz :
     line10=[x.decode('utf-8') for x in line10i]
   else :
     line10=line10i
   if "sep" not in l_infohead :
      sep=None
      for lsep in [',',' ', '\t', ';'] :
        countsep=[x.count(lsep) for x in line10]
        print(lsep,countsep)
        if countsep[0]>3 and countsep.count(countsep[0]) > 5 :
          if sep ==None :
            print(countsep[0], countsep.count(countsep[0]))
            sep=lsep
          else :
            print('not found separator between space, comma, tab and ;\n exit')
            sys.exit(3)
      if not sep :
         print('not found separator between space, comma, tab and ;\n exit')
         sys.exit(3)
   else :
      sep=GetSep(l_filehead[l_infohead.index('sep')])
   print('sep used '+str(sep))
   header=[x.lower() for x in spl(line10i[0], sep)]
   return (check1args(args.rs_header,header, l_filehead,l_infohead,'rsID', listrs, False), check1args(args.n_header,header, l_filehead,l_infohead, 'N',listn, False), check1args(args.chro_header, header,l_filehead,l_infohead,'Chro',listchr, False), check1args(args.pos_header,header,l_filehead,l_infohead, 'Pos',listpos, False), check1args(args.beta_header,header ,l_filehead,l_infohead,'Beta', listbeta, False), check1args(args.z_header,header,l_filehead,l_infohead,'Z',listz, False), check1args(args.a1_header,header,l_filehead,l_infohead,'A1', lista1,True), check1args(args.a2_header,header,l_filehead,l_infohead,'A2', lista0,True), check1args(args.se_header, header,l_filehead,l_infohead,'Se', listse,False), check1args(args.freq_header, header,l_filehead,l_infohead,'freqA1',listaf, False), check1args(args.p_header, header,l_filehead,l_infohead,'Pval', listp, True),sep)
        


      
      

args = parseArguments()


(rs_header, n_header, chro_header, pos_header, beta_header, z_header, a1_header,a2_header, se_header, freq_header, p_header ,sep)=check_args(args)
## keep information bim bam
print(" reading bfile :begin")
#(listkey, listdup, dicrskey,dickeyrs)
keep_genet=args.keep_genet
if args.bfile :
  (listkey, listdup, dicrskey,dickeyrs)=extractinfobim(args.bfile+".bim")
else :
  keep_genet=False 
  (listkey, listdup, dicrskey,dickeyrs)=(None, [], None, None)

print(" reading bfile : end")
print(" reading sumstat: begin")
# sumstat, dickeyrs,listdup, listkey, rs_header,n_header, chr_header, bp_header, beta_header, z_header, a1_header, a2_header, se_header, af_header, p_header,n_value, maf, nlim, keep_genet,sep
list_exclude=[]
if args.exclude_pos :
  reada=open(args.exclude_pos)
  list_exclude=set([':'.join(a.replace('\n','').split()[0:2]) for a in reada])
  reada.close()

(dicsumstat,listkey)=readsummarystat(args.inp_resgwas, dickeyrs,listdup, listkey, rs_header, n_header, chro_header, pos_header, beta_header, z_header, a1_header, a2_header, se_header, freq_header, p_header,args.n, args.maf, args.n_lim ,keep_genet,sep, list_exclude)

print(len(listkey))
print(" reading sumstat: end")
print(" reading add N and freq: begin")
(dicsumstat,listkey)=ExtractFreqN(args.bfile, dicsumstat, listkey,dicrskey,freq_header,rs_header,n_header,chro_header,pos_header, args.bin_plk, args.keep, args.threads,args.memory, args.n, args.maf)
print(len(listkey))
print(" reading add N and freq: End")
print(" reading add beta se")
dicsumstat=addbetase(dicsumstat, listkey ,beta_header, se_header, z_header)
print(" reading add beta se :end")
dicsumstat=addz(dicsumstat, listkey ,args.z_header)
dicsumstat=compute_usingp(dicsumstat, listkey,args.used_p, beta_header, se_header, z_header)
dicsumstat=checkrs(dicsumstat, listkey)

writeplink=open(args.out_head+".plink",'w')
writelz=open(args.out_head+".lz",'w')
writegcta=open(args.out_head+".gcta",'w')
#0 SNP1, 1 SNP2, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
# 0 SNPBim, 1 SNPSumStat, 2 CHR, 3 BP, 4 A1, 5 A2, 6 Z, 7 BETA, 8 SE, 9 AF, 10 N, P 11 
headplink=["SNP","SNPI","SNPSS","CHR","BP","A1","A2", "Z","BETA","SE", "FRQ","N","P"]
headgcta=["SNPKEY","SNP","SNPSS","chr","bp","A1","A2", "z","b","se", "freq","N","p"]
#["#CHROM","BEGIN","END","MARKER_ID","PVALUE"]
headlz=["#CHROM","BEGIN","END","MARKER_ID","PVALUE"]

writeplink.write("\t".join(headplink)+"\n")
writegcta.write("\t".join(headgcta)+"\n")
writelz.write("\t".join(headlz)+"\n")
nbcol=len(headplink)-1
for key in listkey :
    info=dicsumstat[key]
    if len(info)!=nbcol :
      print(info)
      sys.exit('size of key should be size of output ')
    strdata=[str(x) for x in info]
    if strdata[11] == 'NA' :
      continue
    writegcta.write(key+"\t"+"\t".join(strdata)+"\n")



