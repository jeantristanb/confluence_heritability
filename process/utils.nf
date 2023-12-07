/*return dummy dir for file null*/
def errormess(message,exitn=0){
    if(message=="")return(0)
    println(message)
    System.exit(exitn)
}

def getdummy_dir(){
  filescript=file(workflow.scriptFile)
  projectdir="${filescript.getParent()}"
  dummy_dir="${projectdir}/data/"
  return dummy_dir
}

/*change string in memory for nextflow*/
def strmem(val){
 return val as nextflow.util.MemoryUnit
}

