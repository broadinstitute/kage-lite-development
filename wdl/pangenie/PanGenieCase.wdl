version 1.0

struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}

workflow PanGenieCase {
    input {
        File panel_vcf_gz
        File panel_vcf_gz_tbi
        File input_fasta
        File reference_fasta
        Array[String] chromosomes
        String sample_name
        String output_prefix

        String docker
        File? monitoring_script
    }

    call PanGenie {
        input:
            panel_vcf_gz = panel_vcf_gz,
            panel_vcf_gz_tbi = panel_vcf_gz_tbi,
            input_fasta = input_fasta,
            reference_fasta = reference_fasta,
            chromosomes = chromosomes,
            sample_name = sample_name,
            output_prefix = output_prefix,
            docker = docker,
            monitoring_script = monitoring_script
    }

    output {
        File genotyping_vcf_gz = PanGenie.genotyping_vcf_gz
        File genotyping_vcf_gz_tbi = PanGenie.genotyping_vcf_gz_tbi
        File histogram = PanGenie.histogram
        File path_segments_fasta = PanGenie.path_segments_fasta
    }
}

task PanGenie {
    input {
        File panel_vcf_gz
        File panel_vcf_gz_tbi
        File input_fasta
        File reference_fasta
        Array[String] chromosomes
        String sample_name
        String output_prefix

        String docker
        File? monitoring_script
        Int? kmer_length = 31
        String? extra_args

        RuntimeAttributes runtime_attributes = {}
    }

    Int num_chromosomes = length(chromosomes)

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # subset reference and panel VCF to chromosomes
        samtools faidx -r <(echo -e "~{sep="\n" chromosomes}") ~{reference_fasta} > reference.subset.fa
        bcftools view ~{panel_vcf_gz} -r ~{sep="," chromosomes} > panel.subset.vcf

        NPROC=$(nproc)
        NUM_CHROMOSOMES=~{num_chromosomes}
        NUM_THREADS=$(( NPROC < NUM_CHROMOSOMES ? NPROC : NUM_CHROMOSOMES ))

        /pangenie/build/src/PanGenie -g \
            -o ~{output_prefix} \
            -r reference.subset.fa \
            -s ~{sample_name} \
            -t $NUM_THREADS \
            -j $NUM_THREADS \
            -v panel.subset.vcf \
            -i ~{input_fasta} \
            -k ~{kmer_length} \
            ~{extra_args}

        bgzip -c ~{output_prefix}_genotyping.vcf > ~{output_prefix}_genotyping.vcf.gz
        bcftools index -t ~{output_prefix}_genotyping.vcf.gz
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File genotyping_vcf_gz = "~{output_prefix}_genotyping.vcf.gz"
        File genotyping_vcf_gz_tbi = "~{output_prefix}_genotyping.vcf.gz.tbi"
        File histogram = "~{output_prefix}_histogram.histo"
        File path_segments_fasta = "~{output_prefix}_path_segments.fasta"
        File monitoring_log = "monitoring.log"
    }
}