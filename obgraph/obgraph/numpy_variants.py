import logging

import dill
import numpy as np


class NumpyVariants:
    properties = {"header", "variants"}
    def __init__(self, header, variants):
        self.header = header
        self.variants = variants

    @classmethod
    def from_file(cls, file_name):
        try:
            with open(file_name, 'rb') as f:
                data = dill.load(f)
        except FileNotFoundError:
            with open(file_name + ".npz", 'rb') as f:
                data = dill.load(f)

        logging.info("Loaded from %s" % file_name)
        return data

    def to_file(self, file_name):
        if not file_name.endswith(".pkl"):
            file_name = file_name + ".pkl"

        with open(file_name, 'wb') as f:
            dill.dump(self, f)
        logging.info("Saved to %s" % file_name)

    def to_vcf_with_genotypes(self, file_name, sample_name, genotypes, add_header_lines=None, ignore_homo_ref=False, add_genotype_likelihoods=None):
        logging.info("Writing to file %s" % file_name)
        if ignore_homo_ref:
            logging.info("Will not write variants with genotype 0/0")

        if add_genotype_likelihoods is not None:
            p = add_genotype_likelihoods
            genotype_likelihoods = p * np.log10(np.e)

        with open(file_name, "w") as f:
            # f.write(self._header_lines)
            for header_line in self.header:  # last element is empty
                if header_line.startswith("#CHROM"):
                    if sample_name != "":
                        header_line = header_line.strip() + "\t" + sample_name + "\n"
                    if add_header_lines is not None:
                        for additional_header_line in add_header_lines:
                            f.writelines([additional_header_line + "\n"])
                f.writelines([header_line])

            lines = []
            for i, (variant, genotype) in enumerate(zip(self.variants, genotypes)):
                genotype = genotype.decode("utf-8")
                if ignore_homo_ref and genotype == "0/0":
                    continue

                variant = variant.decode("utf-8").strip()
                if i % 1000000 == 0:
                    logging.info("%d variants written to file." % i)

                line = "%s\t%s\n" % (variant, genotype)
                if add_genotype_likelihoods is not None:
                    """
                    #likelyhoods = add_genotype_likelyhoods[i]
                    phred_likelyhoods = [
                        #float(round(np.maximum(0, np.minimum(255, -prob * np.log10(np.e))), 4))
                        int(np.minimum(2550000000000, -10 * prob * np.log10(np.e)))
                        for prob in likelyhoods
                    ]
                    """
                    gl_str = ",".join(str(round(p, 4)) for p in genotype_likelihoods[i])
                    #print(likelyhoods, phred_likelyhoods, np.exp(likelyhoods))
                    line = "%s:GL\t%s:%s\n" % (variant, genotype, gl_str)
                    #print(line)

                lines.append(line)

            f.writelines(lines)

    @classmethod
    def from_vcf(cls, vcf_file_name):
        header = []
        variants = []
        with open(vcf_file_name, "r") as f:
            for i, line in enumerate(f):
                #line_decoded = line.decode("utf-8")
                if line.startswith("#"):
                    header.append(line)
                else:
                    # add byte line
                    # remove info field to save space
                    l = line.split()
                    l[7] = "."
                    line = "\t".join(l)
                    variants.append(str.encode(line))  # add as byte

                if i % 10000 == 0:
                    logging.info("%d variants processed" % i)

        return cls(header, variants)


