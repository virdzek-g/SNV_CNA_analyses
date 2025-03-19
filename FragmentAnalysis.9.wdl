workflow FragmentAnalysis {
	
    # a string for the sample name
    String sample_name     
    File input_bam
    File input_bam_index
    File annotation
    File plot_list
    String? size_range = "15 500"
    Int? map_quality= 20
    Float cram_disk_multiplier = 4
    Int? additional_disk
    Int disk_pad = select_first([additional_disk, 20])
    String? Docker = "gvirdz/fragmentsize:0.3"
    Float input_size = size(input_bam, "GB")
    
    call Fragment_analysis {
        input:
        	input_bam=input_bam,
            input_bam_index=input_bam_index,
            sample_name=sample_name,
            annotation=annotation,
            plot_list=plot_list,
            map_quality=map_quality,
            size_range=size_range,
            Docker=Docker,
            disk_size = ceil(input_size * cram_disk_multiplier) + disk_pad
       }
       
      output {
      File Fragment_analysis_FragmentFM = Fragment_analysis.FragmentFM
	}
}
task Fragment_analysis {
    File input_bam
    File input_bam_index
    String sample_name
    File annotation
   	File plot_list
	Int? map_quality
    String? size_range
    String? Docker
	Int disk_size
    Int? maxRetries = 1
    String? runtime_cpu = "8"
    String? runtime_memory = "32"
    Float? runtime_disk_multiplier = 1
    Int? runtime_disk = ceil((size(input_bam, 'G')) * runtime_disk_multiplier) +  1
    Int? runtime_preemptible = 2
    
    
    command <<<

    set -euxo pipefail
    
    #mkdir /processing_dir
        
    #bamname=$(basename ${input_bam})
    #bainame=$(basename ${input_bam_index})
    
    #ln -s ${input_bam} /processing_dir/$bamname
    #ln -s ${input_bam_index} /processing_dir/$bainame.bai

	mkdir results
    
    python ../../GenerateFragmentFeatures.py --input_bam ${input_bam} \
    	--sample_name ${sample_name} \
        --annotation ${annotation} \
        --plot_list ${plot_list} \
        --results_dir ./results \
        --map_quality ${map_quality}  \
        --size_range ${size_range} \
        --cpus ${runtime_cpu}
   
    echo "Sample Name: ${sample_name}"
    ls -l ./results/
   	mv ./results/${sample_name}_FragmentFM.tsv .

    
   >>>
   
    output {
      File FragmentFM="${sample_name}_FragmentFM.tsv"
   }
   
    runtime {
            docker: "${Docker}"	 
        	cpu: runtime_cpu
        	memory: '${runtime_memory} GB'
        	preemptible: runtime_preemptible
            disks: "local-disk " + disk_size + " HDD"
            maxRetries: maxRetries
    }

}
    