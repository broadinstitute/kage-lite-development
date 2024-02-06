from dataclasses import dataclass


@dataclass
class GenotypingConfig:
    avg_coverage: int = 15
    use_naive_priors: bool = 15
    n_threads: int = 4
    ignore_helper_model: bool = False
    min_genotype_quality: float = 0
    sample_name_output: str = "DONOR"
    ignore_homo_ref: bool = False
    use_gpu: bool = False


    @classmethod
    def from_command_line_args(cls, args):
        return cls(
            args.average_coverage,
            args.use_naive_priors,
            args.n_threads,
            args.ignore_helper_model,
            args.min_genotype_quality,
            args.sample_name_output,
            args.ignore_homo_ref,
            args.gpu
        )