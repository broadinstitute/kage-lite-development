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
        File repeat_mask_bed
        File segmental_duplications_bed
        File simple_repeats_bed
        String output_prefix
        Array[String] chromosomes
        Array[String] leave_one_out_sample_names

        # TODO we require the alignments to subset by chromosome; change to start from raw reads
        Array[File] leave_one_out_bams
        Array[File] leave_one_out_bam_idxs

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
            repeat_mask_bed = repeat_mask_bed,
            segmental_duplications_bed = segmental_duplications_bed,
            simple_repeats_bed = simple_repeats_bed,
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
        call KAGEPlusGLIMPSECase {
            input:
                input_fasta = PreprocessCaseReads.preprocessed_fasta,
                panel_index = KAGELeaveOneOutPanel.index,
                panel_kmer_index_only_variants_with_revcomp = KAGELeaveOneOutPanel.kmer_index_only_variants_with_revcomp,
                glimpse_panel_vcf_gz = CreateLeaveOneOutPanelVCF.leave_one_out_panel_vcf_gz,
                glimpse_panel_vcf_gz_tbi = CreateLeaveOneOutPanelVCF.leave_one_out_panel_vcf_gz_tbi,
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
                panel_vcf_gz = PreprocessPanelVCF.preprocessed_panel_vcf_gz,
                panel_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_vcf_gz_tbi,
                label = "KAGE",
                sample_name = leave_one_out_sample_name,
                docker = docker,
                monitoring_script = monitoring_script,
        }

        # KAGE+GLIMPSE evaluation
        call CalculateMetrics as CalculateMetricsKAGEPlusGLIMPSE {
            input:
                case_vcf_gz = KAGEPlusGLIMPSECase.glimpse_vcf_gz,
                case_vcf_gz_tbi = KAGEPlusGLIMPSECase.glimpse_vcf_gz_tbi,
                panel_vcf_gz = PreprocessPanelVCF.preprocessed_panel_vcf_gz,
                panel_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_vcf_gz_tbi,
                label = "KAGE+GLIMPSE",
                sample_name = leave_one_out_sample_name,
                docker = docker,
                monitoring_script = monitoring_script,
        }

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
                panel_vcf_gz = PreprocessPanelVCF.preprocessed_panel_vcf_gz,
                panel_vcf_gz_tbi = PreprocessPanelVCF.preprocessed_panel_vcf_gz_tbi,
                label = "PanGenie",
                sample_name = leave_one_out_sample_name,
                docker = docker,
                monitoring_script = monitoring_script,
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
            bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.subset.vcf.gz -- -t AF

        truvari anno svinfo -o ~{output_prefix}.subset.svinfo.vcf.gz ~{output_prefix}.subset.vcf.gz

        bcftools annotate --no-version -a ~{repeat_mask_bed} -c CHROM,FROM,TO -m +RM -Oz -o ~{output_prefix}.subset.svinfo.RM.vcf.gz ~{output_prefix}.subset.svinfo.vcf.gz
        bcftools annotate --no-version -a ~{segmental_duplications_bed} -c CHROM,FROM,TO -m +SD -Oz -o ~{output_prefix}.subset.svinfo.RM.SD.vcf.gz ~{output_prefix}.subset.svinfo.RM.vcf.gz
        bcftools annotate --no-version -a ~{simple_repeats_bed} -c CHROM,FROM,TO -m +SR -Oz -o ~{output_prefix}.preprocessed.vcf.gz ~{output_prefix}.subset.svinfo.RM.SD.vcf.gz
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

        bcftools view --no-version ~{input_vcf_gz} -s ^~{leave_one_out_sample_name} | bcftools plugin fill-tags --no-version -Oz -o ~{output_prefix}.preprocessed.LOO.vcf.gz -- -t AF,AC,AN
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

task KAGEPlusGLIMPSECase {
    input {
        File input_fasta
        File panel_index
        File panel_kmer_index_only_variants_with_revcomp
        File glimpse_panel_vcf_gz
        File glimpse_panel_vcf_gz_tbi
        File reference_fasta_fai
        Array[String] chromosomes
        String output_prefix
        String sample_name
        Int average_coverage

        String docker
        File? monitoring_script

        Int? kmer_mapper_num_threads = 1
        String kmer_mapper_args = "-d True -c 100000000 -t ~{kmer_mapper_num_threads}"

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

        bgzip -c ~{output_prefix}.kage.vcf > ~{output_prefix}.kage.vcf.gz
        bcftools index -t ~{output_prefix}.kage.vcf.gz

        bcftools view --no-version ~{output_prefix}.kage.vcf | sed -e 's/nan/-0.01/g' | sed -e 's/inf/100000/g' | bgzip > ~{output_prefix}.kage.nonan.vcf.gz
        bcftools index -t ~{output_prefix}.kage.nonan.vcf.gz

        # TODO update to GLIMPSE2; first figure out why it complains about AC/AN and GT being inconsistent?
        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_phase_static
        chmod +x GLIMPSE_phase_static

        for CHROMOSOME in ~{sep=" " chromosomes}
        do
            CHROMOSOME_LENGTH=$(grep -P "$CHROMOSOME\t" ~{reference_fasta_fai} | cut -f 2)
            ./GLIMPSE_phase_static \
                -I ~{output_prefix}.kage.nonan.vcf.gz \
                -R ~{glimpse_panel_vcf_gz} \
                --input-region $CHROMOSOME:1-$CHROMOSOME_LENGTH \
                --output-region $CHROMOSOME:1-$CHROMOSOME_LENGTH \
                --input-GL \
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
        File case_vcf_gz
        File case_vcf_gz_tbi
        File panel_vcf_gz
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

        python - --case_vcf_gz ~{case_vcf_gz} \
                 --panel_vcf_gz ~{panel_vcf_gz} \
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

        matplotlib.use('Agg')

        def load_callset(callset_file, **kwargs):
            return allel.read_vcf(callset_file, **kwargs)

        def get_gt_vsp(callset):
            return allel.GenotypeDaskArray(callset['calldata/GT'])

        def calculate_metrics_and_plot(case_vcf_gz, panel_vcf_gz, label, sample_name):
            samples = [sample_name]
            panel_callset = load_callset(panel_vcf_gz, samples=samples, fields='*', alt_number=1)
            panel_gt_vp = get_gt_vsp(panel_callset)[:, 0, :].compute()

            fig, ax = plt.subplots(3, 1, figsize=(12, 12), sharex=True)

            ax[0].set_title(f'{sample_name}, {label}')
            ax[0].set_ylabel('number of variants')
            ax[1].set_ylabel('precision\n(weighted)')
            ax[1].set_ylim([0., 1.005])
            ax[2].set_ylabel('recall\n(weighted)')
            ax[2].set_ylim([0., 1.005])
            ax[2].set_xlabel('ALT length - REF length (bp)')

            is_snp_v = panel_callset['variants/is_snp']
            altlen_v = panel_callset['variants/altlen']
            is_fv = [['(-inf,-500]', altlen_v <= -500],
                     ['(-500,-100]', (-500 < altlen_v) & (altlen_v <= -100)],
                     ['(-100,-50]', (-100 < altlen_v) & (altlen_v <= -50)],
                     ['(-50,-10]', (-50 < altlen_v) & (altlen_v <= -10)],
                     ['(-10,-5]', (-10 < altlen_v) & (altlen_v <= -5)],
                     ['(-5,-1]', (-5 < altlen_v) & (altlen_v <= -1)],
                     ['0 (SNP)', is_snp_v],
                     ['[1,5)', (1 <= altlen_v) & (altlen_v < 5)],
                     ['[5,10)', (5 <= altlen_v) & (altlen_v < 10)],
                     ['[10,50)', (10 <= altlen_v) & (altlen_v < 50)],
                     ['[50,100)', (50 <= altlen_v) & (altlen_v < 100)],
                     ['[100,500)', (100 <= altlen_v) & (altlen_v < 500)],
                     ['[500,1000)', (500 <= altlen_v) & (altlen_v < 1000)],
                     ['[1000,5000)', (1000 <= altlen_v) & (altlen_v < 5000)],
                     ['[5000,10000)', (5000 <= altlen_v) & (altlen_v < 10000)],
                     ['[10000,inf)', 10000 <= altlen_v]]

            callset = load_callset(case_vcf_gz, samples=samples, fields='*', alt_number=1)
            gt_vp = get_gt_vsp(callset)[:, 0, :].compute()

            is_missing_v = np.any(gt_vp == -1, axis=1) | np.any(panel_gt_vp == -1, axis=1)

            metrics_dicts = []
            for context in ['US', 'RM', 'SD', 'SR']:
                num_evals = []
                precisions = []
                recalls = []
                for f, (filter_name, is_v) in enumerate(is_fv):
                    if context == 'US':
                        is_context_v = ~(panel_callset[f'variants/RM'] | panel_callset[f'variants/SD'] | panel_callset[f'variants/SR'])
                    else:
                        is_context_v = panel_callset[f'variants/{context}']
                    is_eval_v = ~is_missing_v & is_v & is_context_v

                    enc_gt_n = np.sum(gt_vp[is_eval_v], axis=1)
                    panel_enc_gt_n = np.sum(panel_gt_vp[is_eval_v], axis=1)

                    num_eval = np.sum(is_eval_v)
                    precision = sklearn.metrics.precision_score(panel_enc_gt_n, enc_gt_n, average='weighted')
                    recall = sklearn.metrics.recall_score(panel_enc_gt_n, enc_gt_n, average='weighted')

                    num_evals.append(num_eval)
                    precisions.append(precision)
                    recalls.append(recall)

                    metrics_dicts.append({
                        'LABEL': label,
                        'SAMPLE_NAME': sample_name,
                        'CONTEXT': context,
                        'ALTLEN': filter_name,
                        'NUM_EVAL': num_eval,
                        'PRECISION': precision,
                        'RECALL': recall
                    })

                ax[0].semilogy(num_evals, label=f'{context}')
                ax[0].legend(loc='upper right')
                ax[1].plot(precisions, label=f'{label} {context}')
                ax[2].plot(recalls, label=f'{label} {context}')

            metrics_df = pd.DataFrame.from_dict(metrics_dicts)
            metrics_df.to_csv(f'{sample_name}.{label}.metrics.tsv', sep='\t', index=False)

            ax[2].set_xticks(range(len(is_fv)))
            ax[2].set_xticklabels([filter_name for filter_name, _ in is_fv], rotation=270)
            plt.savefig(f'{sample_name}.{label}.metrics.png')


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
        File metrics_plot = "~{sample_name}.~{label}.metrics.png"
        File monitoring_log = "monitoring.log"
    }
}