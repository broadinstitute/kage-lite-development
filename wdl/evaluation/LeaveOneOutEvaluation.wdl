version 1.0

import "../kage/KAGEPanel.wdl" as KAGEPanel
import "../pangenie/PanGenieCase.wdl" as PanGenieCase

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
        File reference_dict
        File repeat_mask_bed
        File segmental_duplications_bed
        File simple_repeats_bed
        File challenging_medically_relevant_genes_bed
        String output_prefix
        Array[String] chromosomes
        Array[String] leave_one_out_sample_names
        Boolean do_pangenie

        # TODO we require the alignments to subset by chromosome; change to start from raw reads
        Array[File] leave_one_out_crams

        String docker
        String gatk_docker
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
            repeat_mask_bed = repeat_mask_bed,
            segmental_duplications_bed = segmental_duplications_bed,
            simple_repeats_bed = simple_repeats_bed,
            challenging_medically_relevant_genes_bed = challenging_medically_relevant_genes_bed,
            chromosomes = chromosomes,
            output_prefix = output_prefix,
            docker = docker,
            monitoring_script = monitoring_script,
            runtime_attributes = runtime_attributes
    }

    scatter (i in range(length(leave_one_out_sample_names))) {
        String leave_one_out_sample_name = leave_one_out_sample_names[i]
        String leave_one_out_cram = leave_one_out_crams[i]
        String leave_one_out_output_prefix = output_prefix + ".LOO-" + leave_one_out_sample_name

        call IndexCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = leave_one_out_cram,
                docker = docker,
                monitoring_script = monitoring_script
        }

        call PreprocessCaseReads {
            # TODO we require the alignments to subset by chromosome; change to start from raw reads
            input:
                input_cram = leave_one_out_cram,
                input_cram_idx = IndexCaseReads.cram_idx,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                reference_dict = reference_dict,
                output_prefix = leave_one_out_sample_name,
                chromosomes = chromosomes,
                docker = gatk_docker,
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
                input_vcf_gz = CreateLeaveOneOutPanelVCF.leave_one_out_panel_bi_vcf_gz,
                input_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_one_out_panel_bi_vcf_gz_tbi,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                output_prefix = leave_one_out_output_prefix,
                chromosomes = chromosomes,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }

        # KAGE+GLIMPSE case
        call KAGEPlusGLIMPSECase {
            input:
                input_fasta = PreprocessCaseReads.preprocessed_fasta,
                panel_index = KAGELeaveOneOutPanel.index,
                panel_kmer_index_only_variants_with_revcomp = KAGELeaveOneOutPanel.kmer_index_only_variants_with_revcomp,
                panel_multi_split_vcf_gz = CreateLeaveOneOutPanelVCF.leave_one_out_panel_multi_split_vcf_gz,
                panel_multi_split_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_one_out_panel_multi_split_vcf_gz_tbi,
                panel_split_vcf_gz = CreateLeaveOneOutPanelVCF.leave_one_out_panel_split_vcf_gz,
                panel_split_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_one_out_panel_split_vcf_gz_tbi,
                reference_fasta_fai = reference_fasta_fai,
                chromosomes = chromosomes,
                output_prefix = leave_one_out_sample_name,
                sample_name = leave_one_out_sample_name,
                docker = kage_docker,
                monitoring_script = monitoring_script
        }

        # KAGE evaluation
        call CalculateMetrics as CalculateMetricsKAGE {
            input:
                case_vcf_gz = KAGEPlusGLIMPSECase.kage_vcf_gz,
                case_vcf_gz_tbi = KAGEPlusGLIMPSECase.kage_vcf_gz_tbi,
                panel_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                panel_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                label = "KAGE",
                sample_name = leave_one_out_sample_name,
                docker = docker,
                monitoring_script = monitoring_script
        }

        # KAGE+GLIMPSE evaluation
        call CalculateMetrics as CalculateMetricsKAGEPlusGLIMPSE {
            input:
                case_vcf_gz = KAGEPlusGLIMPSECase.glimpse_vcf_gz,
                case_vcf_gz_tbi = KAGEPlusGLIMPSECase.glimpse_vcf_gz_tbi,
                panel_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                panel_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                label = "KAGE+GLIMPSE",
                sample_name = leave_one_out_sample_name,
                docker = docker,
                monitoring_script = monitoring_script
        }

        if (do_pangenie) {
            # PanGenie case
            call PanGenieCase.PanGenie as PanGenieCase {
                input:
                    panel_vcf_gz = CreateLeaveOneOutPanelVCF.leave_one_out_panel_vcf_gz,
                    panel_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_one_out_panel_vcf_gz_tbi,
                    input_fasta = PreprocessCaseReads.preprocessed_fasta,
                    reference_fasta = reference_fasta,
                    chromosomes = chromosomes,
                    sample_name = leave_one_out_sample_name,
                    output_prefix = leave_one_out_sample_name,
                    docker = pangenie_docker,
                    monitoring_script = monitoring_script
            }

            # PanGenie evaluation
            call CalculateMetrics as CalculateMetricsPanGenie {
                input:
                    case_vcf_gz = PanGenieCase.genotyping_vcf_gz,
                    case_vcf_gz_tbi = PanGenieCase.genotyping_vcf_gz_tbi,
                    panel_vcf_gz = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz,
                    panel_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_split_vcf_gz_tbi,
                    label = "PanGenie",
                    sample_name = leave_one_out_sample_name,
                    docker = docker,
                    monitoring_script = monitoring_script
            }
        }
    }

    output {
    }
}

task PreprocessPanelVCF {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File repeat_mask_bed
        File segmental_duplications_bed
        File simple_repeats_bed
        File challenging_medically_relevant_genes_bed
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

        bcftools view --no-version ~{input_vcf_gz} -r ~{sep="," chromosomes} -Ou | \
            bcftools norm --no-version -m+ -Ou | \
            bcftools plugin fill-tags --no-version -Ou -- -t AF,AC,AN | \
            truvari anno svinfo | \
            bcftools annotate --no-version -a ~{repeat_mask_bed} -c CHROM,FROM,TO -m +RM -Ou | \
            bcftools annotate --no-version -a ~{segmental_duplications_bed} -c CHROM,FROM,TO -m +SD -Ou | \
            bcftools annotate --no-version -a ~{simple_repeats_bed} -c CHROM,FROM,TO -m +SR -Ou | \
            bcftools annotate --no-version -a ~{challenging_medically_relevant_genes_bed} -c CHROM,FROM,TO -m +CMRG -Oz -o ~{output_prefix}.preprocessed.vcf.gz
        bcftools index -t ~{output_prefix}.preprocessed.vcf.gz

        bcftools view --no-version --min-alleles 3  ~{output_prefix}.preprocessed.vcf.gz -Ou | \
            bcftools norm --no-version -m- -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.multi.split.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.multi.split.vcf.gz

        bcftools norm --no-version -m- ~{output_prefix}.preprocessed.vcf.gz -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.split.temp.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.split.temp.vcf.gz

        bcftools annotate --no-version -a ~{output_prefix}.preprocessed.multi.split.vcf.gz ~{output_prefix}.preprocessed.split.temp.vcf.gz -m +MULTIALLELIC \
            -Oz -o ~{output_prefix}.preprocessed.split.vcf.gz
        bcftools index -t ~{output_prefix}.preprocessed.split.vcf.gz
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
        File preprocessed_panel_split_vcf_gz = "~{output_prefix}.preprocessed.split.vcf.gz"
        File preprocessed_panel_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.split.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
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
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File cram_idx = "~{output_prefix}.cram.crai"
        File monitoring_log = "monitoring.log"
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

        # subset cram to chromosomes
        gatk PrintReads \
            -L chromosomes.bed \
            -I ~{input_cram} \
            --read-index ~{input_cram_idx} \
            -R ~{reference_fasta} \
            --disable-sequence-dictionary-validation \
            -O ~{output_prefix}.bam

        # filter out read pairs containing N nucleotides
        # TODO move functionality into KAGE code
        samtools fasta ~{output_prefix}.bam | sed -E '~{filter_N_regex}' > ~{output_prefix}.preprocessed.fasta
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

        bcftools view --no-version ~{input_vcf_gz} -s ^~{leave_one_out_sample_name} --trim-alt-alleles -Ou | \
            bcftools view --no-version --min-alleles 2 -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LOO.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LOO.vcf.gz

        bcftools view --no-version --min-alleles 2 --max-alleles 2  ~{output_prefix}.preprocessed.LOO.vcf.gz -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LOO.bi.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LOO.bi.vcf.gz

        bcftools view --no-version --min-alleles 3  ~{output_prefix}.preprocessed.LOO.vcf.gz -Ou | \
            bcftools norm --no-version -m- -Ou | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LOO.multi.split.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LOO.multi.split.vcf.gz

        bcftools norm --no-version -m- ~{output_prefix}.preprocessed.LOO.vcf.gz -Oz | \
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LOO.split.vcf.gz -- -t AF,AC,AN
        bcftools index -t ~{output_prefix}.preprocessed.LOO.split.vcf.gz
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
        File leave_one_out_panel_split_vcf_gz = "~{output_prefix}.preprocessed.LOO.split.vcf.gz"
        File leave_one_out_panel_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.LOO.split.vcf.gz.tbi"
        File leave_one_out_panel_bi_vcf_gz = "~{output_prefix}.preprocessed.LOO.bi.vcf.gz"
        File leave_one_out_panel_bi_vcf_gz_tbi = "~{output_prefix}.preprocessed.LOO.bi.vcf.gz.tbi"
        File leave_one_out_panel_multi_split_vcf_gz = "~{output_prefix}.preprocessed.LOO.multi.split.vcf.gz"
        File leave_one_out_panel_multi_split_vcf_gz_tbi = "~{output_prefix}.preprocessed.LOO.multi.split.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }
}

task KAGEPlusGLIMPSECase {
    input {
        File input_fasta
        File panel_index
        File panel_kmer_index_only_variants_with_revcomp
        File panel_multi_split_vcf_gz # for filling in biallelic-only VCFs produced by KAGE
        File panel_multi_split_vcf_gz_tbi
        File panel_split_vcf_gz       # for GLIMPSE
        File panel_split_vcf_gz_tbi
        File reference_fasta_fai
        Array[String] chromosomes
        String output_prefix
        String sample_name
        Int average_coverage

        String docker
        File? monitoring_script

        String kmer_mapper_args = "-c 100000000"

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        NPROC=$(nproc)

        kmer_mapper map \
            ~{kmer_mapper_args} \
            -t $NPROC \
            -i ~{panel_kmer_index_only_variants_with_revcomp} \
            -f ~{input_fasta} \
            -o ~{output_prefix}.kmer_counts.npy

        kage genotype \
            -i ~{panel_index} \
            -c ~{output_prefix}.kmer_counts.npy \
            --average-coverage ~{average_coverage} \
            -s ~{sample_name} \
            -I true \
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

        bcftools view --no-version ~{output_prefix}.kage.vcf.gz | \
            sed -e 's/nan/-1000000.0/g' | sed -e 's/-inf/-1000000.0/g' | sed -e 's/inf/-1000000.0/g' | bgzip > ~{output_prefix}.kage.nonan.vcf.gz
        bcftools index -t ~{output_prefix}.kage.nonan.vcf.gz

        # TODO update to GLIMPSE2; first figure out why it complains about AC/AN and GT being inconsistent?
        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static
        chmod +x GLIMPSE_phase_static

        # TODO parallelize
        for CHROMOSOME in ~{sep=" " chromosomes}
        do
            CHROMOSOME_LENGTH=$(grep -P "$CHROMOSOME\t" ~{reference_fasta_fai} | cut -f 2)
            ./GLIMPSE_phase_static \
                -I ~{output_prefix}.kage.nonan.vcf.gz \
                -R ~{panel_split_vcf_gz} \
                --input-region $CHROMOSOME:1-$CHROMOSOME_LENGTH \
                --output-region $CHROMOSOME:1-$CHROMOSOME_LENGTH \
                --input-GL \
                --thread $NPROC \
                --output $CHROMOSOME.vcf.gz
            bcftools index -t $CHROMOSOME.vcf.gz
        done

        bcftools concat --no-version ~{sep=".vcf.gz " chromosomes}.vcf.gz -Oz -o ~{output_prefix}.kage.glimpse.vcf.gz
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
        File kmer_counts = "~{output_prefix}.kmer_counts.npy"
        File kage_vcf_gz = "~{output_prefix}.kage.vcf.gz"
        File kage_vcf_gz_tbi = "~{output_prefix}.kage.vcf.gz.tbi"
        File glimpse_vcf_gz = "~{output_prefix}.kage.glimpse.vcf.gz"
        File glimpse_vcf_gz_tbi = "~{output_prefix}.kage.glimpse.vcf.gz.tbi"
        File monitoring_log = "monitoring.log"
    }
}

task CalculateMetrics {
    input {
        File case_vcf_gz        # bi+multi split
        File case_vcf_gz_tbi
        File panel_vcf_gz       # bi+multi split
        File panel_vcf_gz_tbi
        String label
        String sample_name

        String docker
        File? monitoring_script

        RuntimeAttributes runtime_attributes = {}
    }

    command <<<
        set -e

        # Create a zero-size monitoring log file so it exists even if we don't pass a monitoring script
        touch monitoring.log
        if [ -s ~{monitoring_script} ]; then
            bash ~{monitoring_script} > monitoring.log &
        fi

        # split multiallelics in case (may be redundant)
        bcftools norm --no-version -m- ~{case_vcf_gz} -Oz -o case.split.vcf.gz
        bcftools index -t case.split.vcf.gz

        # mark case variants in panel
        bcftools annotate --no-version -a case.split.vcf.gz ~{panel_vcf_gz} -m +CASE -Oz -o panel.annot.vcf.gz
        bcftools index -t panel.annot.vcf.gz

        conda install -y seaborn

        python - --case_vcf_gz case.split.vcf.gz \
                 --panel_vcf_gz panel.annot.vcf.gz \
                 --label ~{label} \
                 --sample_name ~{sample_name} \
                 <<-'EOF'
        import argparse
        import allel
        import numpy as np
        import pandas as pd
        import sklearn.metrics
        import matplotlib
        import matplotlib.pyplot as plt
        import seaborn as sns

        matplotlib.use('Agg')

        def load_callset(callset_file, **kwargs):
            return allel.read_vcf(callset_file, **kwargs)

        def get_gt_vsp(callset):
            return allel.GenotypeDaskArray(callset['calldata/GT'])

        def calculate_metrics_and_plot(case_vcf_gz, panel_vcf_gz, label, sample_name):
            samples = [sample_name]
            panel_callset = load_callset(panel_vcf_gz, samples=samples, fields='*', alt_number=1)
            is_case_V = panel_callset['variants/CASE']
            panel_gt_vp = get_gt_vsp(panel_callset)[is_case_V, 0, :].compute()


            is_snp_v = panel_callset['variants/is_snp'][is_case_V]
            is_multiallelic_v = panel_callset['variants/MULTIALLELIC'][is_case_V]

            altfreq_v = panel_callset['variants/AF'][is_case_V]
            is_altfreq_v = [['[0%, 1%)', (0. <= altfreq_v) & (altfreq_v < 0.01)],
                            ['[1%, 5%)', (0.01 <= altfreq_v) & (altfreq_v < 0.05)],
                            ['[5%, 10%)', (0.05 <= altfreq_v) & (altfreq_v < 0.1)],
                            ['[10%, 50%)', (0.1 <= altfreq_v) & (altfreq_v < 0.50)],
                            ['[50%, 100%]', (0.5 <= altfreq_v) & (altfreq_v <= 1.)]]

            altlen_v = panel_callset['variants/altlen'][is_case_V]
            is_altlen_v = [['(-inf,-500]', altlen_v <= -500],
                           ['(-500,-50]', (-500 < altlen_v) & (altlen_v <= -50)],
                           ['(-50,-1]', (-50 < altlen_v) & (altlen_v <= -1)],
                           ['0 (SNP)', is_snp_v],
                           ['[1,50)', (1 <= altlen_v) & (altlen_v < 50)],
                           ['[50,500)', (50 <= altlen_v) & (altlen_v < 500)],
                           ['[500,5000)', (500 <= altlen_v) & (altlen_v < 5000)],
                           ['[5000,inf)', 10000 <= altlen_v]]
            is_sv_v = (altlen_v <= -50) | (altlen_v >= 50)

            callset = load_callset(case_vcf_gz, samples=samples, fields='*', alt_number=1)
            gt_vp = get_gt_vsp(callset)[:, 0, :].compute()

            is_missing_v = np.any(gt_vp == -1, axis=1) | np.any(panel_gt_vp == -1, axis=1)

            metrics_dicts = []

            num_i = len(is_altfreq_v)
            num_j = len(is_altlen_v)
            num_evals = np.zeros((num_i, num_j))
            precisions = np.zeros((num_i, num_j))
            recalls = np.zeros((num_i, num_j))
            f1s = np.zeros((num_i, num_j))

            for context in ['ALL', 'CMRG', 'US', 'RM', 'SD', 'SR']:

                if context == 'ALL':
                    is_context_v = True
                elif context == 'US':
                    is_context_v = ~(panel_callset[f'variants/RM'][is_case_V] |
                                     panel_callset[f'variants/SD'][is_case_V] |
                                     panel_callset[f'variants/SR'][is_case_V])
                else:
                    is_context_v = panel_callset[f'variants/{context}'][is_case_V]

                for is_multiallelic in [True, False]:
                    is_allelic_v = is_multiallelic_v if is_multiallelic else ~is_multiallelic_v
                    allelic = 'multiallelic' if is_multiallelic else 'biallelic'

                    for i, (filter_name_i, is_v_i) in enumerate(is_altfreq_v):
                        for j, (filter_name_j, is_v_j) in enumerate(is_altlen_v):
                            is_eval_v = ~is_missing_v & is_v_i & is_v_j & is_allelic_v & is_context_v

                            enc_gt_n = np.sum(gt_vp[is_eval_v], axis=1)
                            panel_enc_gt_n = np.sum(panel_gt_vp[is_eval_v], axis=1)

                            num_eval = np.sum(is_eval_v)
                            precision = np.nan if panel_enc_gt_n.size == 0 else sklearn.metrics.precision_score(panel_enc_gt_n, enc_gt_n, average='weighted')
                            recall = np.nan if panel_enc_gt_n.size == 0 else sklearn.metrics.recall_score(panel_enc_gt_n, enc_gt_n, average='weighted')
                            f1 = np.nan if panel_enc_gt_n.size == 0 else sklearn.metrics.f1_score(panel_enc_gt_n, enc_gt_n, average='weighted')

                            num_evals[i][j] = num_eval
                            precisions[i][j] = precision
                            recalls[i][j] = recall
                            f1s[i][j] = f1

                            metrics_dicts.append({
                                'LABEL': label,
                                'SAMPLE_NAME': sample_name,
                                'CONTEXT': context,
                                'MULTIALLELIC': is_multiallelic,
                                'ALTFREQ': filter_name_i,
                                'ALTLEN': filter_name_j,
                                'NUM_EVAL': num_eval,
                                'PRECISION': precision,
                                'RECALL': recall,
                                'F1': f1
                            })

                    if num_evals.sum() == 0:
                        continue

                    non_sv_enc_gt_n = np.sum(gt_vp[~is_missing_v & ~is_sv_v & is_allelic_v & is_context_v], axis=1)
                    non_sv_panel_enc_gt_n = np.sum(panel_gt_vp[~is_missing_v & ~is_sv_v & is_allelic_v & is_context_v], axis=1)
                    non_sv_precision = np.nan if non_sv_panel_enc_gt_n.size == 0 else sklearn.metrics.precision_score(non_sv_panel_enc_gt_n, non_sv_enc_gt_n, average='weighted')
                    non_sv_recall = np.nan if non_sv_panel_enc_gt_n.size == 0 else sklearn.metrics.recall_score(non_sv_panel_enc_gt_n, non_sv_enc_gt_n, average='weighted')
                    non_sv_f1 = np.nan if non_sv_panel_enc_gt_n.size == 0 else sklearn.metrics.f1_score(non_sv_panel_enc_gt_n, non_sv_enc_gt_n, average='weighted')
                    non_sv_count = np.sum(~is_sv_v & is_allelic_v & is_context_v)

                    sv_enc_gt_n = np.sum(gt_vp[~is_missing_v & is_sv_v & is_allelic_v & is_context_v], axis=1)
                    sv_panel_enc_gt_n = np.sum(panel_gt_vp[~is_missing_v & is_sv_v & is_allelic_v & is_context_v], axis=1)
                    sv_precision = sklearn.metrics.precision_score(sv_panel_enc_gt_n, sv_enc_gt_n, average='weighted')
                    sv_recall = sklearn.metrics.recall_score(sv_panel_enc_gt_n, sv_enc_gt_n, average='weighted')
                    sv_f1 = sklearn.metrics.f1_score(sv_panel_enc_gt_n, sv_enc_gt_n, average='weighted')
                    sv_count = np.sum(is_sv_v & is_allelic_v & is_context_v)

                    fig, ax = plt.subplots(4, 1, figsize=(12, 16), sharex=True)

                    ax[0] = sns.heatmap(num_evals, ax=ax[0], linewidths=1, linecolor='k', annot=True,
                                        norm=matplotlib.colors.LogNorm(), cmap='Blues')
                    cbar = ax[0].collections[0].colorbar
                    cbar.ax.set_ylabel('number of variants', rotation=270, labelpad=20)
                    ax[0].set_title(f'{label}\n{sample_name}\ncontext = {context}, {allelic}\n' +
                                    f'non-SV allele count = {non_sv_count}, SV allele count = {sv_count}')

                    ax[1] = sns.heatmap(precisions, ax=ax[1], linewidths=1, linecolor='k', annot=True,
                                        vmin=0.7, cmap='Greens')
                    cbar = ax[1].collections[0].colorbar
                    cbar.ax.set_ylabel('precision\n(weighted GC)', rotation=270, labelpad=30)
                    ax[1].set_title(f'non-SV precision = {non_sv_precision:.4f}\nSV precision = {sv_precision:.4f}')

                    ax[2] = sns.heatmap(recalls, ax=ax[2], linewidths=1, linecolor='k', annot=True,
                                        vmin=0.7, cmap='Greens')
                    cbar = ax[2].collections[0].colorbar
                    cbar.ax.set_ylabel('recall\n(weighted GC)', rotation=270, labelpad=30)
                    ax[2].set_title(f'non-SV recall = {non_sv_recall:.4f}\nSV recall = {sv_recall:.4f}')

                    ax[3] = sns.heatmap(f1s, ax=ax[3], linewidths=1, linecolor='k', annot=True,
                                        vmin=0.7, cmap='Greens')
                    cbar = ax[3].collections[0].colorbar
                    cbar.ax.set_ylabel('F1\n(weighted GC)', rotation=270, labelpad=30)
                    ax[3].set_title(f'non-SV F1 = {non_sv_f1:.4f}\nSV F1 = {sv_f1:.4f}')

                    ax[3].set_xlabel('ALT length - REF length (bp)')
                    ax[3].set_xticklabels([filter_name for filter_name, _ in is_altlen_v], rotation=0)
                    for a in ax:
                        a.set_ylabel('AF')
                        a.set_yticklabels([filter_name for filter_name, _ in is_altfreq_v], rotation=0)

                    plt.savefig(f'{sample_name}.{context}.{allelic}.{label}.metrics.png')
                    plt.show()

            metrics_df = pd.DataFrame.from_dict(metrics_dicts)
            metrics_df.to_csv(f'{sample_name}.{label}.metrics.tsv', sep='\t', index=False)

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--case_vcf_gz',
                                type=str)

            parser.add_argument('--panel_vcf_gz',
                                type=str)

            parser.add_argument('--label',
                                type=str)

            parser.add_argument('--sample_name',
                                type=str)

            args = parser.parse_args()

            calculate_metrics_and_plot(args.case_vcf_gz,
                                       args.panel_vcf_gz,
                                       args.label,
                                       args.sample_name)

        if __name__ == '__main__':
            main()
        EOF
    >>>

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
        File metrics_tsv = "~{sample_name}.~{label}.metrics.tsv"
        Array[File] metrics_plots = glob("~{sample_name}.*.png")
        File monitoring_log = "monitoring.log"
    }
}