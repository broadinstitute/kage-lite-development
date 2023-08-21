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

workflow KAGECase {
    input {
        File input_cram
        File panel_index
        File panel_kmer_index_only_variants_with_revcomp
        File panel_split_vcf_gz # for GLIMPSE
        File panel_split_vcf_gz_tbi
        File panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        File panel_multi_split_vcf_gz_tbi
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
        Array[String] chromosomes
        Boolean subset_reads = true
        String sample_name
        Float average_coverage

        String docker
        String kage_docker
        File? monitoring_script
    }

    call IndexCaseReads {
        # TODO we require the alignments to subset by chromosome; change to start from raw reads
        input:
            input_cram = input_cram,
            docker = docker,
            monitoring_script = monitoring_script
    }

    if (subset_reads) {
        call PreprocessCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = input_cram,
                input_cram_idx = IndexCaseReads.cram_idx,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                reference_dict = reference_dict,
                output_prefix = sample_name,
                chromosomes = chromosomes,
                docker = docker,
                monitoring_script = monitoring_script
        }
    }

    if (!subset_reads) {
        call PreprocessCaseReadsWithoutSubsetting {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = input_cram,
                input_cram_idx = IndexCaseReads.cram_idx,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                reference_dict = reference_dict,
                output_prefix = sample_name,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }
    }

    call KAGECase {
        input:
            input_fasta = select_first([PreprocessCaseReads.preprocessed_fasta, PreprocessCaseReadsWithoutSubsetting.preprocessed_fasta]),
            panel_index = panel_index,
            panel_kmer_index_only_variants_with_revcomp = panel_kmer_index_only_variants_with_revcomp,
            panel_multi_split_vcf_gz = panel_multi_split_vcf_gz,
            panel_multi_split_vcf_gz_tbi = panel_multi_split_vcf_gz_tbi,
            reference_fasta_fai = reference_fasta_fai,
            output_prefix = sample_name,
            sample_name = sample_name,
            average_coverage = average_coverage,
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    scatter (j in range(length(chromosomes))) {
        call GLIMPSECaseChromosome {
            input:
                kage_vcf_gz = KAGECase.kage_vcf_gz,
                kage_vcf_gz_tbi = KAGECase.kage_vcf_gz_tbi,
                panel_split_vcf_gz = panel_split_vcf_gz,
                panel_split_vcf_gz_tbi = panel_split_vcf_gz_tbi,
                reference_fasta_fai = reference_fasta_fai,
                chromosome = chromosomes[j],
                output_prefix = sample_name,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }
    }

    call GLIMPSECaseGather {
        input:
            chromosome_glimpse_vcf_gzs = GLIMPSECaseChromosome.chromosome_glimpse_vcf_gz,
            chromosome_glimpse_vcf_gz_tbis = GLIMPSECaseChromosome.chromosome_glimpse_vcf_gz_tbi,
            output_prefix = sample_name,
            docker = kage_docker,
            monitoring_script = monitoring_script
    }

    output {
    }
}

# some 1000G CRAM indices have issues due to htsjdk version, so we reindex; see e.g. https://github.com/broadinstitute/gatk/issues/7076
task IndexCaseReads {
    input {
        File input_cram

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    # if we instead simply use File cram_idx = "~{input_cram}.crai" in the output block,
    # Terra tries to localize an index adjacent to the CRAM in PreprocessCaseReads???
    String output_prefix = basename(input_cram, ".cram")

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        samtools index -@ $(nproc) ~{input_cram} ~{output_prefix}.cram.crai
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, true]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File cram_idx = "~{output_prefix}.cram.crai"
    }
}

task PreprocessCaseReads {
    input {
        File input_cram
        File input_cram_idx
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
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

        # filter out read pairs containing N nucleotides
        # TODO move functionality into KAGE code
        samtools view --reference ~{reference_fasta} -@ $(nproc) -L chromosomes.bed -u ~{input_cram} | \
            samtools fasta --reference ~{reference_fasta} -@ $(nproc) | \
            sed -E '~{filter_N_regex}' > ~{output_prefix}.preprocessed.fasta
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File preprocessed_fasta = "~{output_prefix}.preprocessed.fasta"
    }
}

task PreprocessCaseReadsWithoutSubsetting {
    input {
        File input_cram
        File input_cram_idx
        File reference_fasta
        File reference_fasta_fai
        File reference_dict
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

        # filter out read pairs containing N nucleotides
        # TODO move functionality into KAGE code
        samtools fasta --reference ~{reference_fasta} -@ $(nproc) ~{input_cram} | sed -E '~{filter_N_regex}' > ~{output_prefix}.preprocessed.fasta
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 2])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File preprocessed_fasta = "~{output_prefix}.preprocessed.fasta"
    }
}

task KAGECase {
    input {
        File input_fasta
        File panel_index
        File panel_kmer_index_only_variants_with_revcomp
        File panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        File panel_multi_split_vcf_gz_tbi
        File reference_fasta_fai
        String output_prefix
        String sample_name
        Float average_coverage

        String docker
        File? monitoring_script

        String kmer_mapper_args = "-c 100000000"
        Boolean? ignore_helper_model = true
        String? kage_genotype_extra_args


        RuntimeAttributes runtime_attributes = {}
        Int? cpu = 8
    }

    Int cpu_resolved = select_first([runtime_attributes.cpu, cpu])

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        kmer_mapper map \
            ~{kmer_mapper_args} \
            -t ~{cpu_resolved} \
            -i ~{panel_kmer_index_only_variants_with_revcomp} \
            -f ~{input_fasta} \
            -o ~{output_prefix}.kmer_counts.npy

        kage genotype \
            -i ~{panel_index} \
            -c ~{output_prefix}.kmer_counts.npy \
            --average-coverage ~{average_coverage} \
            -s ~{sample_name} \
            -I ~{ignore_helper_model} \
            ~{kage_genotype_extra_args} \
            -o ~{output_prefix}.kage.bi.vcf

        # we need to add split multiallelics to biallelic-only KAGE VCF
        # create single-sample header from LOO panel w/ split multiallelics
        bcftools view --no-version -h -G ~{panel_multi_split_vcf_gz} | \
            sed 's/##fileformat=VCFv4.2/##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods.">/g' | \
            sed 's/INFO$/INFO\tFORMAT\t~{sample_name}/g' > ~{output_prefix}.multi.split.header.txt
        # create single-sample missing genotypes from LOO panel w/ split multiallelics
        bcftools view --no-version -H -G ~{panel_multi_split_vcf_gz} | \
            sed 's/$/\tGT:GL\t.\/.:nan,nan,nan/g' > ~{output_prefix}.multi.split.GT.txt
        # create single-sample VCF w/ split multiallelics
        bgzip -c <(cat ~{output_prefix}.multi.split.header.txt ~{output_prefix}.multi.split.GT.txt) > ~{output_prefix}.multi.split.vcf.gz
        bcftools index -t ~{output_prefix}.multi.split.vcf.gz

        bgzip -c ~{output_prefix}.kage.bi.vcf > ~{output_prefix}.kage.bi.vcf.gz
        bcftools index -t ~{output_prefix}.kage.bi.vcf.gz

        bcftools concat --no-version -a ~{output_prefix}.kage.bi.vcf.gz ~{output_prefix}.multi.split.vcf.gz -Oz -o ~{output_prefix}.kage.vcf.gz
        bcftools index -t ~{output_prefix}.kage.vcf.gz
    }

    runtime {
        docker: docker
        cpu: cpu_resolved
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 100]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File monitoring_log = "monitoring.log"
        File kmer_counts = "~{output_prefix}.kmer_counts.npy"
        File kage_vcf_gz = "~{output_prefix}.kage.vcf.gz"
        File kage_vcf_gz_tbi = "~{output_prefix}.kage.vcf.gz.tbi"
    }
}

task GLIMPSECaseChromosome {
    input {
        File kage_vcf_gz
        File kage_vcf_gz_tbi
        File panel_split_vcf_gz       # for GLIMPSE
        File panel_split_vcf_gz_tbi
        File reference_fasta_fai
        String chromosome
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

        bcftools view --no-version -r ~{chromosome} ~{kage_vcf_gz} | \
            sed -e 's/nan/-1000000.0/g' | sed -e 's/-inf/-1000000.0/g' | sed -e 's/inf/-1000000.0/g' | bgzip > ~{output_prefix}.kage.nonan.~{chromosome}.vcf.gz
        bcftools index -t ~{output_prefix}.kage.nonan.~{chromosome}.vcf.gz

        # TODO update to GLIMPSE2; first figure out why it complains about AC/AN and GT being inconsistent?
        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static
        chmod +x GLIMPSE_phase_static

        CHROMOSOME_LENGTH=$(grep -P "~{chromosome}\t" ~{reference_fasta_fai} | cut -f 2)
        ./GLIMPSE_phase_static \
            -I ~{output_prefix}.kage.nonan.~{chromosome}.vcf.gz \
            -R ~{panel_split_vcf_gz} \
            --input-region ~{chromosome}:1-$CHROMOSOME_LENGTH \
            --output-region ~{chromosome}:1-$CHROMOSOME_LENGTH \
            --input-GL \
            --thread $(nproc) \
            --output ~{output_prefix}.kage.glimpse.~{chromosome}.vcf.gz
        bcftools index -t ~{output_prefix}.kage.glimpse.~{chromosome}.vcf.gz
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
        File monitoring_log = "monitoring.log"
        File chromosome_glimpse_vcf_gz = "~{output_prefix}.kage.glimpse.~{chromosome}.vcf.gz"
        File chromosome_glimpse_vcf_gz_tbi = "~{output_prefix}.kage.glimpse.~{chromosome}.vcf.gz.tbi"
    }
}

task GLIMPSECaseGather {
    input {
        Array[File] chromosome_glimpse_vcf_gzs
        Array[File] chromosome_glimpse_vcf_gz_tbis
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

        bcftools concat --no-version ~{sep=" " chromosome_glimpse_vcf_gzs} -Oz -o ~{output_prefix}.kage.glimpse.vcf.gz
        bcftools index -t ~{output_prefix}.kage.glimpse.vcf.gz
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
        File monitoring_log = "monitoring.log"
        File glimpse_vcf_gz = "~{output_prefix}.kage.glimpse.vcf.gz"
        File glimpse_vcf_gz_tbi = "~{output_prefix}.kage.glimpse.vcf.gz.tbi"
    }
}