version 1.0

import "../kage/KAGEPanel.wdl" as KAGEPanel

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

workflow LeaveOneOutEvaluation {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File reference_fasta
        File reference_fasta_fai
        String output_prefix
        Array[String] chromosomes
        Array[String] leave_one_out_sample_names

        # TODO we require the alignments to subset by chromosome; change to start from raw reads
        Array[File] leave_one_out_bams
        Array[File] leave_one_out_bam_idxs

        File evaluation_script

        String docker
        String kage_docker
        String pangenie_docker
        File? monitoring_script

        RuntimeAttributes? runtime_attributes
        RuntimeAttributes? medium_runtime_attributes
        RuntimeAttributes? large_runtime_attributes
    }

    call PreprocessPanelVCF {
        # TODO move into KAGEPanel?
        input:
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            chromosomes = chromosomes,
            output_prefix = output_prefix,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = runtime_attributes
    }

    scatter (i in range(length(leave_one_out_sample_names))) {
        String leave_one_out_sample_name = leave_one_out_sample_names[i]
        String leave_one_out_bam = leave_one_out_bams[i]
        String leave_one_out_bam_idx = leave_one_out_bam_idxs[i]
        String leave_one_out_output_prefix = output_prefix + ".LOO-" + leave_one_out_sample_name

        call PreprocessCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_bam = leave_one_out_bam,
                input_bam_idx = leave_one_out_bam_idx,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                output_prefix = leave_one_out_sample_name,
                chromosomes = chromosomes,
                docker = docker,
                monitoring_script = monitoring_script
        }

        call CreateLeaveOneOutPanelVCF {
            input:
                input_vcf_gz = PreprocessPanelVCF.preprocessed_panel_vcf_gz,
                input_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_vcf_gz_tbi,
                output_prefix = leave_one_out_output_prefix,
                leave_one_out_sample_name = leave_one_out_sample_name,
                docker = docker,
                monitoring_script = monitoring_script,
                runtime_attributes = runtime_attributes
        }

        call KAGEPanel.KAGEPanel as KAGELeaveOneOutPanel {
            input:
                input_vcf_gz = CreateLeaveOneOutPanelVCF.leave_one_out_panel_vcf_gz,
                input_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_one_out_panel_vcf_gz_tbi,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                output_prefix = leave_one_out_output_prefix,
                chromosomes = chromosomes,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }

        # KAGE+GLIMPSE case
        call KAGECase {
            input:
                input_fasta = PreprocessCaseReads.preprocessed_fasta,
                panel_index = KAGELeaveOneOutPanel.index,
                panel_kmer_index_only_variants_with_revcomp = KAGELeaveOneOutPanel.kmer_index_only_variants_with_revcomp,
                output_prefix = leave_one_out_sample_name,
                sample_name = leave_one_out_sample_name,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }

        # KAGE evaluation

        # PanGenie case

        # PanGenie evaluation
    }

    output {
    }
}

task PreprocessPanelVCF {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        Array[String] chromosomes
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools view --no-version ~{input_vcf_gz} -r ~{sep="," chromosomes} | \
            bcftools norm --no-version -m+ | \
            bcftools view --no-version --max-alleles 2 |
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.vcf.gz -- -t AF
        bcftools index -t ~{output_prefix}.preprocessed.vcf.gz
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File preprocessed_panel_vcf_gz = "~{output_prefix}.preprocessed.vcf.gz"
        File preprocessed_panel_vcf_gz_tbi = "~{output_prefix}.preprocessed.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }
}

task PreprocessCaseReads {
    input {
        File input_bam
        File input_bam_idx
        File reference_fasta
        File reference_fasta_fai
        Array[String] chromosomes
        String output_prefix

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    String filter_N_regex = "/^>/{N;/^>.*\\n.*N.*/d}"

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # hacky way to get chromosomes into bed file
        grep -P '~{sep="\\t|" chromosomes}\t' ~{reference_fasta_fai} | cut -f 1,2 | sed -e 's/\t/\t1\t/g' > chromosomes.bed

        # subset cram to chromosomes
        samtools view -L chromosomes.bed -C -o ~{output_prefix}.cram --write-index ~{input_bam} -T ~{reference_fasta}

        # filter out read pairs containing N nucleotides
        # TODO move functionality into KAGE code
        samtools fasta --reference ~{reference_fasta} ~{output_prefix}.cram | sed -E '~{filter_N_regex}' > ~{output_prefix}.preprocessed.fasta
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File preprocessed_fasta = "~{output_prefix}.preprocessed.fasta"
        File monitoring_log = "monitoring.log"
    }
}

task CreateLeaveOneOutPanelVCF {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        String output_prefix
        String leave_one_out_sample_name

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        bcftools view --no-version ~{input_vcf_gz} -s ^~{leave_one_out_sample_name} | bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LOO.vcf.gz -- -t AF
        bcftools index -t ~{output_prefix}.preprocessed.LOO.vcf.gz
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File leave_one_out_panel_vcf_gz = "~{output_prefix}.preprocessed.LOO.vcf.gz"
        File leave_one_out_panel_vcf_gz_tbi = "~{output_prefix}.preprocessed.LOO.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }
}

task KAGECase {
    input {
        File input_fasta
        File panel_index
        File panel_kmer_index_only_variants_with_revcomp
        String output_prefix
        String sample_name
        Int average_coverage

        String docker
        File? monitoring_script

        Int? num_threads = 1
        String kmer_mapper_args = "-d True -c 100000000 -t ~{num_threads}"

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        kmer_mapper map \
            ~{kmer_mapper_args} \
            -i ~{panel_kmer_index_only_variants_with_revcomp} \
            -f ~{input_fasta} \
            -o ~{output_prefix}.kmer_counts.npy

        kage genotype \
            -i ~{panel_index} \
            -c ~{output_prefix}.kmer_counts.npy \
            --average-coverage ~{average_coverage} \
            -s ~{sample_name} \
            -I true \
            -o ~{output_prefix}.kage.vcf
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File kmer_counts = "~{output_prefix}.kmer_counts.npy"
        File kage_vcf = "~{output_prefix}.kage.vcf"
        #File glimpse_vcf = "~{output_prefix}.kage.glimpse.vcf"
        File monitoring_log = "monitoring.log"
    }
}