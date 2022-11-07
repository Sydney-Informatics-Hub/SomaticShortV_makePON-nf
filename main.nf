#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2

// =================================================================
// main.nf is the pipeline script for a nextflow pipeline
// Should contain the following sections:
	// Import subworkflows
	// Log info function
	// Help function 
	// Main workflow structure
	// Some tests to check input data, essential arguments

// Examples are included for each section. Remove them and replace
// with project-specific code. For more information on nextflow see:
// https://www.nextflow.io/docs/latest/index.html and the SIH Nextflow
// upskilling repo @ INSERT REPO PATH 
//
// ===================================================================

// Import subworkflows to be run in the workflow
// Each of these is a separate .nf script saved in modules/
// Add as many of these as you need. The example below will
// look for the process called process in modules/moduleName.nf
// Include { process } from './modules/moduleName'
  
  refdir='/scratch/wz54/gs5517/sarek_testing/Reference/v0'

// Could not set it to "" or "./" :-(
base_path="/scratch/wz54/npd561/PIPE-2629_thyroid_carcinoma/nextflow_pipelines/nfcore/sarek/Somatic-ShortV/nextflow/make_PON_and_run_mutect2/final_scripts_production/make_PoN"


// Set PATH pointing to the 'bam' files
params.bams = "/scratch/er01/PIPE-2629-ThyroidCancer/nf_sarek/preprocess_*/Preprocessing/*/Recalibrated/*-{N,T}.recal.bam"
// bam pair channel
bam_pair_ch=Channel.fromFilePairs( params.bams )

intervalList=['a','b','c','d','e','f','g','h','i','j','k','l','m','n']

/*
# This file needs to be created on the fly in bash!! - sample_map_vcf.txt
*/

params.sample_map_vcfs = "$base_path/sample_map_vcf.txt"


//mkdir temp_folder
params.temp_dir_folder="$base_path/temp_folder"

params.outdir="$base_path/results"


  
  


/// Print a header for your pipeline 

log.info """\

      ========================================================================================
      ========================================================================================
         Make a Panel Of Normals (PON) for Somatic Short-variant calling 
      ========================================================================================
      ========================================================================================

 -._    _.--'"`'--._    _.--'"`'--._    _.--'"`'--._    _  
    '-:`.'|`|"':-.  '-:`.'|`|"':-.  '-:`.'|`|"':-.  '.` :    
  '.  '.  | |  | |'.  '.  | |  | |'.  '.  | |  | |'.  '.:    
  : '.  '.| |  | |  '.  '.| |  | |  '.  '.| |  | |  '.  '.  
  '   '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.:_ | :_.' '.  `.  
         `-..,..-'       `-..,..-'       `-..,..-'       `       


             ~~~~ Version: ${params.version} ~~~~
 

 Created by the Sydney Informatics Hub, University of Sydney

 Find documentation and more info @ GITHUB REPO DOT COM

 Cite this pipeline @ INSERT DOI

 Log issues @ GITHUB REPO DOT COM

 All of the default parameters are set in `nextflow.config`

=======================================================================================
Workflow run parameters 
=======================================================================================

input       : ${params.input}
outDir      : ${params.outDir}
workDir     : ${workflow.workDir}

=======================================================================================

 """

/// Help function 
// This is an example of how to set out the help function that 
// will be run if run command is incorrect (if set in workflow) 
// or missing/  

def helpMessage() {
    log.info"""
  Usage:  nextflow run <PATH TO REPO>/myPipeline-nf <args> --outDir

  Required Arguments:
	--outDir	Specify path to output directory

	--input		Specify full path and name of sample
			input file (tab separated).
    """.stripIndent()
}

/// Main workflow structure. Include some input/runtime tests here.
// Make sure to comment what each step does for readability. 



include { run_Mutect2_eachNormalSample_splitGatherApproach; GatherVcfs_step; Create_GenomicsDB_from_normalMutect2Calls_GenomicsDBImport; Combine_normalCallsUsing_CreateSomaticPanelOfNormals} from './make_PON_all_steps.nf'


workflow {

// Show help message if --help is run or if any required params are not 
// provided at runtime

        if ( params.help || params.input == false ){   
        // Invoke the help function above and exit
              helpMessage()
              exit 1

        // consider adding some extra contigencies here.
        // could validate path of all input files in list?
        // could validate indexes for input files exist?
        // could validate indexes for reference exist?
        // confirm with each tool, any requirements for their run?

// if none of the above are a problem, then run the workflow
	} else {
	
  // Define input channels 
  cohort_ch = Channel.fromPath("${params.cohort}")
  outDir_ch = Channel.fromPath("${params.outDir}")

	run_Mutect2_eachNormalSample_splitGatherApproach(bam_pair_ch,intervalList,params.temp_dir_folder,base_path)
	
	GatherVcfs_step(run_Mutect2_eachNormalSample_splitGatherApproach.out.collect(),bam_pair_ch,base_path)

	Create_GenomicsDB_from_normalMutect2Calls_GenomicsDBImport(params.sample_map_vcfs,params.temp_dir_folder,GatherVcfs_step.out[0].collect(),base_path)

	Combine_normalCallsUsing_CreateSomaticPanelOfNormals(Create_GenomicsDB_from_normalMutect2Calls_GenomicsDBImport.out,params.temp_dir_folder,base_path)

}}

workflow.onComplete {
  summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
outDir      : ${params.outDir}

=======================================================================================
  """
  println summary

}
