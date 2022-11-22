#!/usr/bin/env nextflow

/// To use DSL-2 will need to include this
nextflow.enable.dsl=2




params.refFasta="$refdir_path/Homo_sapiens_assembly38.fasta"
params.refFastaIndex="$refdir_path/Homo_sapiens_assembly38.fasta.fai"
params.refFastaDict="$refdir_path/Homo_sapiens_assembly38.fasta.dict"



refFastaFile = file(params.refFasta)
refFastaIndexFile = file(params.refFastaIndex)
refFastaDictFile = file(params.refFastaDict)


/*
(1)
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon.sh

(2)
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon_gathervcfs.sh

(3)
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon_genomicsdbimport.sh

(4)
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_cohort_pon.sh
*/


/*
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/README.md

# PON-Step1 - Mutect2: Each Normal Sample
#To create a PoN, call on each normal sample in this mode
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon.sh
#To be parallelised

#Scatter and create PoN across genomic intervals
#https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon.sh
#

*/


process run_Mutect2_eachNormalSample_splitGatherApproach  {


        tag "Mutect2 $bam_id $splitIntervalNumber"

        publishDir "$params.outdir/", mode:'copy'

        input:
                tuple val(bam_id) , file(bams)

                each splitIntervalNumber

		path temp_dir
		path base_path
	

        output:
                path ("${bam_id}_${splitIntervalNumber}_out.vcf.gz")



        script:

        """

        gatk Mutect2 \
             -R $refdir_path/Homo_sapiens_assembly38.fasta \
             -I ${bams[0]} \
	     --max-mnp-distance 0 \
             -L "$base_path/100M_primary_interval_${splitIntervalNumber}.list" \
             -O ${bam_id}_${splitIntervalNumber}_out.vcf.gz
        """


}

// #PON-Step2 - GatherVcfs
// #https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon_gathervcfs.sh

process GatherVcfs_step {

        tag "GatherVcfs_step $bam_id"
        publishDir "$params.outdir/", mode:'copy'


        input:
		path ('*') 
		tuple val(bam_id) , file(bams)

		path base_path

        output:
                path ("${bam_id}_gathered_vcfs_across_subintervals_sorted.vcf.gz") 
                path ("${bam_id}_gathered_vcfs_across_subintervals_sorted.vcf.gz.tbi")
 
                path ("${bam_id}_gathered_vcfs_across_subintervals.list") 


        script:

        """
        ls $base_path/results/${bam_id}*_out.vcf.gz   >${bam_id}_gathered_vcfs_across_subintervals.list

        # GatherVcfs requires intervals in order, so add chrM using MergeVcfs
        gatk GatherVcfs \
                -I  ${bam_id}_gathered_vcfs_across_subintervals.list \
                -O  ${bam_id}_gathered_vcfs_across_subintervals.vcf.gz

        #Sort
        gatk SortVcf \
                -I ${bam_id}_gathered_vcfs_across_subintervals.vcf.gz \
                -O ${bam_id}_gathered_vcfs_across_subintervals_sorted.vcf.gz
        """

}




// #Consolidate samples with GenomicsDBImport
// #https://github.com/Sydney-Informatics-Hub/Somatic-ShortV/blob/master/gatk4_pon_genomicsdbimport.sh



process Create_GenomicsDB_from_normalMutect2Calls_GenomicsDBImport {

        tag "GenomicsDBImport"
        publishDir "$params.outdir/", mode: 'copy'


        input:

	        path sample_vcfs_input
        	path temp_dir

		path ('*')

		path base_path


        output:
        	path 'pon_db'


        script:

        """


         gatk GenomicsDBImport \
                --sample-name-map $sample_vcfs_input \
                --overwrite-existing-genomicsdb-workspace \
                --genomicsdb-workspace-path pon_db \
                --reader-threads ${task.cpus} \
                --tmp-dir ${temp_dir} \
                --intervals "$base_path/100M_primary_interval.list"


        """


}



process Combine_normalCallsUsing_CreateSomaticPanelOfNormals {

	tag "CreateSomaticPanelOfNormals"

        publishDir "$params.outdir/", mode:'copy'

        input:
	        path pon_db
        	path temp_dir

		path base_path

        output:
       		path 'pon.vcf.gz'
       		path 'pon.vcf.gz.tbi'

        script :

        """

        gatk CreateSomaticPanelOfNormals \
                -R $refdir_path/Homo_sapiens_assembly38.fasta \
                -V gendb://${pon_db} \
                -O pon.vcf.gz


        """


}

