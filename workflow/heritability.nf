#!/usr/bin/env nextflow

nextflow.enable.dsl=2


/*
 * Authors       :
 *
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of Cancer Group of Sydney Brenner Institute an
 *  2015-2021
 *
 *
 * Description  : Nextflow pipeline to computed heritabilty
 *
 *(C) University of the Witwatersrand, Johannesburg, 2016-2021
 *This is licensed under the XX licence. See the "LICENSE" file for details
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;

include {getdummy_dir} from '../process/utils.nf'
include {strmem} from '../process/lz.nf'

dummy_dir=getdummy_dir()


workflow {



}


