version 1.0

workflow count_reads_workflow {
  input {
    File bam_file
  }

  call count_reads {
    input: 
      bam_file = bam_file
    }

  output {
    Int reads_0_120 = count_reads.reads_0_120
    Int reads_120_150 = count_reads.reads_120_150
    Int reads_150_plus = count_reads.reads_150_plus
  }
}


task count_reads {
  input {
  File bam_file
  Int? extra_disk
  Float? extra_mem
  Int disk_size = ceil(size(bam_file,"GB") * 2.5 + 50) + select_first([extra_disk,0])
  Float memory = 7.5 + select_first([extra_mem, 0])
 }
  command {
    samtools view -c -F 4 -F 256 -F 2048 ${bam_file} | awk '$10>=0 && $10<=120' > reads_0_120
    samtools view -c -F 4 -F 256 -F 2048 ${bam_file} | awk '$10>120 && $10<=150' > reads_120_150
    samtools view -c -F 4 -F 256 -F 2048 ${bam_file} | awk '$10>150' > reads_150_plus
  }

  output {
    Int reads_0_120 = read_int("reads_0_120")
    Int reads_120_150 = read_int("reads_120_150")
    Int reads_150_plus = read_int("reads_150_plus")
  }

  runtime {
   docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
   cpu: 2
   memory: memory + " GB"
   disks: "local-disk " + disk_size + " HDD"
  }
}