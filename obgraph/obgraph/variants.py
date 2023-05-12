import gzip
import io
import logging
import time

import numpy as np

from .util import encode_chromosome


def get_variant_type(vcf_line):

    l = vcf_line.split()
    if len(l[3]) == len(l[4]):
        return "SNP"
    elif "VT=SNP" in vcf_line:
        return "SNP"
    elif len(l[3]) > 1 and len(l[4]) > 1:
        return "SUBSTITUTION"
    elif len(l[3]) > len(l[4]) and len(l[4]) == 1:
        return "DELETION"
    elif len(l[3]) < len(l[4]) and len(l[3]) == 1:
        return "INSERTION"
    elif "VT=INDEL" in vcf_line:
        if len(l[3]) > len(l[4]):
            return "DELETION"
        else:
            return "INSERTION"
    else:
        return ""
        raise Exception("Unsupported variant type on line %s" % vcf_line)


class TruthRegions:
    def __init__(self, file_name):
        self.regions = []
        f = open(file_name)

        for line in f:
            l = line.split()
            start = int(l[1])
            end = int(l[2])

            self.regions.append((start, end))

    def is_inside_regions(self, position):
        for region in self.regions:
            if position >= region[0] and position < region[1]:
                return True
        return False


class VcfVariant:
    def __init__(self, chromosome, position, ref_sequence="", variant_sequence="", genotype=None, type="", vcf_line=None, vcf_line_number=None):
        self.chromosome = chromosome
        self.position = position
        self.ref_sequence = ref_sequence
        self.variant_sequence = variant_sequence
        self.genotype = genotype
        self.vcf_line = vcf_line
        if self.genotype == "1|0":
            self.genotype = "0|1"

        self.type = type
        self.vcf_line_number = vcf_line_number

        self._filter = "PASS"

    def get_genotype(self):
        return self.genotype.replace("|", "/").replace("1/0", "0/1").replace("./.", "0/0").replace(".", "0/0")

    def copy(self):
        return VcfVariant(self.chromosome, self.position, self.ref_sequence, self.variant_sequence, self.genotype, self.type, self.vcf_line, self.vcf_line_number)

    def set_filter_by_prob(self, prob_correct, criteria_for_pass=0.99):
        if prob_correct > criteria_for_pass:
            self._filter = "PASS"
        else:
            self._filter = "LowQUAL"

    def get_numeric_genotype(self):
        g = self.get_genotype()
        if g == "0/0":
            return 1
        elif g == "1/1":
            return 2
        elif g == "0/1":
            return 3
        else:
            raise Exception("Invalid genotype. Genotype is: %s" % g)


    def set_genotype(self, genotype, is_numeric=False):
        if is_numeric:
            assert genotype in [0, 1, 2, 3]
            if genotype == 1 or genotype == 0:
                self.genotype = "0/0"
            elif genotype == 2:
                self.genotype = "1/1"
            else:
                self.genotype = "0/1"
        else:
            assert genotype in ["0|0", "0/0", "0|1", "0/1", "1/1", "1|1"], "Invalid genotype %s" % genotype
            self.genotype = genotype

    def id(self):
        #return (self.chromosome, self.position)
        return (self.chromosome, self.position, self.ref_sequence.lower(), self.variant_sequence.lower())

    def get_vcf_line(self):
        chromosome = self.chromosome
        if chromosome == 23:
            chromosome = "X"
        elif chromosome == 24:
            chromosome = "Y"
        else:
            chromosome = str(chromosome)

        return "%s\t%d\t.\t%s\t%s\t.\t%s\t.\tGT\t%s\n"  % (chromosome, self.position, self.ref_sequence, self.variant_sequence, self._filter, self.genotype if self.genotype is not None else ".")

    def __str__(self):
        return "%s:%d %s/%s %s %s" % (self.chromosome, self.position, self.ref_sequence, self.variant_sequence, self.genotype, self.type)

    def length(self):
        if self.type == "SNP":
            return 1
        elif self.type == "DELETION":
            return len(self.get_deleted_sequence())
        elif self.type == "INSERTION" or self.type == "SUBSTITUTION":
            return len(self.variant_sequence) - 1
        else:
            raise Exception("Variant type %s not implemented" % self.type)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if self.id() == other.id():
            if self.genotype == other.genotype:
                return True
        else:
            logging.error("ID mismatch: %s != %s" % (self.id(), other.id()))

        return False

    def get_deleted_sequence(self):
        if self.type == "DELETION":
            return self.ref_sequence[1:]

        raise Exception("Variant %s is not a deletion" % self)

    def get_variant_sequence(self):
        if self.type == "DELETION":
            return ""
        elif self.type == "INSERTION" or self.type == "SUBSTITUTION":
            return self.variant_sequence[1:]
        else:
            return self.variant_sequence

    def get_inserted_sequence(self):
        # different from get_variant_sequence for deletions
        if self.type == "DELETION":
            inserted_sequence = self.ref_sequence[1:]
        elif self.type == "INSERTION":
            inserted_sequence = self.variant_sequence[1:]
        elif self.type == "SUBSTITUTION":
            inserted_sequence = self.variant_sequence[1:]
        else:
            raise Exception("Unsupported variant")

        return inserted_sequence

    def get_reference_sequence(self):
        if self.type == "INSERTION":
            return ""
        elif self.type == "INSERTION" or self.type == "SUBSTITUTION":
            return self.ref_sequence[1:]
        else:
            return self.ref_sequence

    def get_variant_allele_frequency(self):
        try:
            info_field = self.vcf_line.split("\t")[7]
            af = float(info_field.split("AF=")[1].split(";")[0])
        except IndexError:
            logging.error("Could not get allele frequency from vcf line. VCF needs to contain an info column with AF= some allele frequency")
            raise
        return af

    def get_reference_allele_frequency(self):
        return 1 - self.get_variant_allele_frequency()

    def get_individuals_and_numeric_haplotypes(self):
        # genotypes are 1, 2, 3 for homo ref, homo alt and hetero
        for individual_id, genotype_string in enumerate(self.vcf_line.split()[9:]):
            genotype_string = genotype_string.split(":")[0].replace("/", "|")
            if genotype_string == "0|0":
                numeric = (1, 1)
            elif genotype_string == "1|1":
                numeric = (2, 2)
            elif genotype_string == "0|1":
                numeric = (1, 2)
            elif genotype_string == "1|0":
                numeric = (2, 1)
            else:
                logging.error("Could not parse genotype string %s" % genotype_string)
                raise Exception("Unknown genotype")

            yield (individual_id, numeric)

    def _numeric_genotype(self, genotype_string):
        numeric = 0
        if genotype_string == "0|0":
            numeric = 1
        elif genotype_string == "1|1":
            numeric = 2
        elif genotype_string == "0|1" or genotype_string == "1|0":
            numeric = 3
        else:
            logging.error("Could not parse genotype string %s" % genotype_string)
            raise Exception("Unknown genotype")

        return numeric

    def _numeric_genotype2(self, genotype_string):
        numeric = 4
        if genotype_string == "0|0":
            numeric = 0
        elif genotype_string == "1|1":
            numeric = 2
        elif genotype_string == "0|1" or genotype_string == "1|0":
            numeric = 1
        elif genotype_string == ".":
            numberic = 4
        else:
            logging.error("Could not parse genotype string %s" % genotype_string)
            raise Exception("Unknown genotype")

        return numeric

    def _numeric_genotype3(self, genotype_string):
        numeric = 3
        if genotype_string == "0|0":
            numeric = 0
        elif genotype_string == "1|1":
            numeric = 2
        elif genotype_string == "0|1" or genotype_string == "1|0":
            numeric = 1
        else:
            logging.error("Could not parse genotype string %s" % genotype_string)
            raise Exception("Unknown genotype")

        return numeric

    def get_individuals_and_numeric_genotypes(self, encoding_version="1"):
        # genotypes are 1, 2, 3 for homo ref, homo alt and hetero
        for individual_id, genotype_string in enumerate(self.vcf_line.split()[9:]):
            genotype_string = genotype_string.split(":")[0].replace("/", "|")
            if encoding_version == "1":
                numeric = self._numeric_genotype(genotype_string)
            elif encoding_version == "2":
                numeric = self._numeric_genotype2(genotype_string)
            elif encoding_version == "3":
                numeric = self._numeric_genotype3(genotype_string)
            else:
                raise Exception("Invalid encoding version. Must be 1 or 2")
            yield (individual_id, numeric)

    @classmethod
    def from_vcf_line(cls, line, vcf_line_number=None, dont_encode_chromosome=False):
        l = line.split()
        chromosome = l[0]
        if not dont_encode_chromosome:
            chromosome = int(encode_chromosome(l[0]))

        position = int(l[1])
        ref_sequence = l[3].lower().replace("n", "a")  # we only want actg
        variant_sequence = l[4].lower().replace("n", "a")

        if len(l) >= 10:
            genotype = l[9].split(":")[0].replace("/", "|")
        else:
            genotype = ""

        #logging.info("Genotype is %s. Line length is %d (%s) " % (genotype, len(l), l))

        return cls(chromosome, position, ref_sequence, variant_sequence, genotype, get_variant_type(line), line, vcf_line_number=vcf_line_number)

    def get_reference_position_before_variant(self):
        # Returns the position of the last base pair before the variant starts (will be end of node before variant)
        # For SNPS, the position is the actual SNP node
        # FOr deletions, the position is one basepair before the deleted sequence
        # For insertion, the position is one baseiapr before the inserted sequence
        # Subtract -1 in the end to make it 0-based
        if self.type == "SNP":
            return self.position - 1 - 1
        else:
            return self.position - 1

    def get_reference_position_after_variant(self):
        # Returns the next reference position after the variant is finished
        start = self.get_reference_position_before_variant()
        if self.type == "SNP":
            return start + 2
        elif self.type == "INSERTION":
            return start + 1
        elif self.type == "DELETION":
            return start + len(self.ref_sequence) - 1 + 1
        elif self.type == "SUBSTITUTION":
            return start + len(self.ref_sequence)
        else:
            raise Exception("Not able to get ref pos after variant of type %s" % self.type)


class VcfVariants:
    def __init__(self, variant_genotypes=[], skip_index=False, header_lines=""):
        self._header_lines = header_lines
        self.variant_genotypes = variant_genotypes
        self._index = {}
        self._position_index = {}
        if not skip_index:
            self.make_index()

    def __getitem__(self, item):
        return self.variant_genotypes[item]

    def get_header(self):
        return self._header_lines

    def n_individuals(self):
        return len(list(self.variant_genotypes[0].get_individuals_and_numeric_genotypes()))

    def n_variants(self):
        assert isinstance(self.variant_genotypes, list), "Only possible when variants is a list"
        return len(self.variant_genotypes)

    def get_variant_by_line_number(self, line_number):
        return self.variant_genotypes[line_number]

    def add_variants(self, variant_list):
        self.variant_genotypes.extend(variant_list)

    def get_chunks(self, chunk_size=5000, add_variants_to_list=None):
        logging.info("Getting variant chunks of size %d" % chunk_size)
        out = []
        prev_time = time.time()
        i = 0
        for variant_number, variant in enumerate(self.variant_genotypes):
            if variant_number % 100000 == 0:
                logging.info("%d variants read" % variant_number)
            out.append(variant)
            if add_variants_to_list is not None:
                add_variants_to_list.append(variant)
            i += 1
            if i >= chunk_size and chunk_size > 0:
                logging.info("Reading chunk of %d variants took %.5f sec" % (chunk_size, time.time()-prev_time))
                prev_time = time.time()
                yield out
                out = []
                i = 0
        logging.info("Reading chunk of %d variants took %.5f sec" % (chunk_size, time.time() - prev_time))
        yield out

    def make_index(self):
        logging.info("Making vcf index")
        for variant in self.variant_genotypes:
            self._index[variant.id()] = variant
        logging.info("Done making vcf index")

    def apply_on_sequence(self, sequence, chromosome, sequence_start, sequence_end):
        assert self._position_index is not None, "Must create position index first"

        variants = self.get_variants_in_region(chromosome, sequence_start+1, sequence_end)

        sequence = list(sequence)
        n_changed = 0
        for variant in variants:
            if variant.get_genotype() == "0/1" or variant.get_genotype() == "1/1":
                if variant.type == "SNP":
                    position = variant.position - sequence_start - 1
                    assert sequence[position].lower() == variant.ref_sequence.lower(), "Variant %s does not match sequence %s. Position %s. Sequence at position: %s" % (variant, sequence, position, sequence[position-1:position+2])
                    sequence[position] = variant.variant_sequence
                    n_changed += 1

        #logging.info("Sequence changed %d times" % n_changed)
        return ''.join(sequence)

    def make_position_index(self):
        logging.info("Finding unique chromosomes")
        unique_chromosomes = np.unique([variant.chromosome for variant in self])

        logging.info("Making index")
        self._position_index = {chromosome: [] for chromosome in unique_chromosomes}
        self._variants_by_chromosome = {chromosome: [] for chromosome in unique_chromosomes}

        for variant in self:
            self._position_index[variant.chromosome].append(variant.position)
            self._variants_by_chromosome[variant.chromosome].append(variant)

        logging.info("Converting to numpy")
        for chrom in unique_chromosomes:
            self._position_index[chrom] = np.array(self._position_index[chrom])
        """
        for chrom in unique_chromosomes:
            logging.info("Making for chromosome %s" % chrom)
            self._position_index[chrom] = np.array([variant.position for variant in self if variant.chromosome == chrom])
            self._variants_by_chromosome[chrom] = [variant for variant in self if variant.chromosome == chrom]
        """

    @staticmethod
    def translate_numeric_genotypes_to_literal(genotypes):
        numeric_to_literal = {0: "0/0", 1: "0/0", 2: "1/1", 3: "0/1"}
        return [numeric_to_literal[g] for g in genotypes]

    def has_variant_left_of_variant(self, variant, window=31):
        other = self.get_variants_in_region(variant.chromosome, variant.position - window, variant.position)
        assert variant.position not in [v.position for v in other], "Variants between %d and %d contain %s position: %s" % (variant.position-window, variant.position, variant, other)
        if len(other) > 0:
            return True
        return False

    def has_variant_right_of_variant(self, variant, window=31):
        other = self.get_variants_in_region(variant.chromosome, variant.position+1, variant.position+window)
        assert variant.position not in [v.position for v in other]
        if len(other) > 0:
            return True
        return False

    def get_variants_in_region(self, chromosome, start, end):
        chromosome = int(chromosome)
        if chromosome not in self._variants_by_chromosome:
            #logging.info("Invalid chromosome %s of type %s. Possible chromosomes are %s" % (chromosome, type(chromosome), self._variants_by_chromosome.keys()))
            return []
        start_index = np.searchsorted(self._position_index[chromosome], start)
        end_index = np.searchsorted(self._position_index[chromosome], end)
        return self._variants_by_chromosome[chromosome][start_index:end_index]

    def has_variant(self, variant_genotype):
        if variant_genotype.id() in self._index:
            return True

        return False

    def has_variant_genotype(self, variant_genotype):
        if self.has_variant(variant_genotype) and self._index[variant_genotype.id()].genotype == variant_genotype.genotype:
            return True

        return False

    def get(self, variant):
        if variant.id() not in self._index:
            return None
        return self._index[variant.id()]

    def __getitem__(self, item):
        if item >= len(self.variant_genotypes):
            logging.error("Variant with id %d does not exist. There are only %d variants" % (item, len(self.variant_genotypes)))

        return self.variant_genotypes[item]

    def __len__(self):
        return len(self.variant_genotypes)

    def __iter__(self):
        return self.variant_genotypes.__iter__()

    def __next__(self):
        return self.variant_genotypes.__next__()

    def to_vcf_file(self, file_name, add_individual_to_header="DONOR", ignore_homo_ref=False, add_header_lines=None, sample_name_output="DONOR"):
        logging.info("Writing to file %s" % file_name)
        with open(file_name, "w") as f:
            #f.write(self._header_lines)
            for header_line in self._header_lines.split("\n")[0:-1]:  # last element is empty
                if header_line.startswith("#CHROM"):
                    if sample_name_output != "":
                        header_line += "\t" + sample_name_output
                    # first write additional header line, the chrom-line should be the last one
                    if add_header_lines is not None:
                        for additional_header_line in add_header_lines:
                            f.writelines([additional_header_line + "\n"])

                f.writelines([header_line + "\n"])


            for i, variant in enumerate(self):
                if i % 1000000 == 0:
                    logging.info("%d variants written to file" % i)

                if ignore_homo_ref and (variant.get_genotype() == "0/0" or variant.get_genotype() == ""):
                    continue

                f.writelines([variant.get_vcf_line()])


    @staticmethod
    def read_header_from_vcf(vcf_file_name):
        f = open(vcf_file_name)
        is_bgzipped = False
        if vcf_file_name.endswith(".gz"):
            logging.info("Assuming gzipped file")
            is_bgzipped = True
            gz = gzip.open(vcf_file_name)
            f = io.BufferedReader(gz, buffer_size=1000 * 1000 * 2)  # 2 GB buffer size?
            logging.info("Made gz file object")

        header_lines = ""
        for line in f:
            if is_bgzipped:
                line = line.decode("utf-8")

            if line.startswith("#"):
                header_lines += line
            else:
                break

        return header_lines

    @classmethod
    def from_vcf(cls, vcf_file_name, skip_index=False, limit_to_n_lines=None, make_generator=False, limit_to_chromosome=None):
        logging.info("Reading variants from file")
        variant_genotypes = []

        if limit_to_chromosome is not None:
            logging.info("Will only read variants from chromsome %s" % limit_to_chromosome)

            limit_to_chromosome = encode_chromosome(limit_to_chromosome)

        # Get header
        header_lines = VcfVariants.read_header_from_vcf(vcf_file_name)

        f = open(vcf_file_name)
        is_bgzipped = False
        if vcf_file_name.endswith(".gz"):
            logging.info("Assuming gzipped file")
            is_bgzipped = True
            gz = gzip.open(vcf_file_name)
            f = io.BufferedReader(gz, buffer_size=1000 * 1000 * 2)  # 2 GB buffer size?
            logging.info("Made gz file object")

        if make_generator:
            assert limit_to_chromosome is None, "Cannot both limit to chromosome and make generator (not implemented)"
            logging.info("Returning variant generator")
            if is_bgzipped:
                f = (line for line in f if not line.decode("utf-8").startswith("#"))
                return cls((VcfVariant.from_vcf_line(line.decode("utf-8"), vcf_line_number=i) for i, line in enumerate(f) if not line.decode("utf-8").startswith("#")), skip_index=skip_index, header_lines=header_lines)
            else:
                f = (line for line in f if not line.startswith("#"))
                return cls((VcfVariant.from_vcf_line(line, vcf_line_number=i) for i, line in enumerate(f)), skip_index=skip_index, header_lines=header_lines)

        n_variants_added = 0
        prev_time = time.time()
        variant_number = -1
        for i, line in enumerate(f):
            if is_bgzipped:
                line = line.decode("utf-8")

            if line.startswith("#"):
                continue

            variant_number += 1

            if i % 5000000 == 0:
                logging.info("Read %d variants from file (time: %.3f). %d variants added" % (i, time.time()-prev_time, n_variants_added))
                prev_time = time.time()

            if limit_to_n_lines is not None and i >= limit_to_n_lines:
                logging.warning("Limited to %d lines" % limit_to_n_lines)
                break


            variant = VcfVariant.from_vcf_line(line, vcf_line_number=variant_number)
            if limit_to_chromosome is not None and variant.chromosome != limit_to_chromosome:
                if len(variant_genotypes) > 0:
                    logging.info("Stopping reading file since limiting to chromosome and now on new chromosome")
                    break
                else:
                    continue

            n_variants_added += 1
            variant_genotypes.append(variant)

        return cls(variant_genotypes, skip_index, header_lines)

    def __str__(self):
        return str(list(self))

    def __repr__(self):
        return self.__str__()

    def copy(self):
        return VcfVariants([v.copy() for v in self])

    def intersect(self, other_variants):
        new_variants = []
        for variant in self:
            if other_variants.has_variant(variant):
                new_variants.append(variant)

        logging.info("Variants before intersection: %d Variants after: %d" % (len(self.variant_genotypes), len(new_variants)))
        return VcfVariants(new_variants, skip_index=True, header_lines=self._header_lines)

    def print_variants_not_in_other(self, other):
        for variant in self:
            if not other.has_variant(variant):
                print("Only in set 1: %s" % variant)

    def print_diff(self, other):
        for variant in self:
            if not other.has_variant(variant):
                print("Only in set 1: %s" % variant)
            elif not other.has_variant_genotype(variant):
                print("Genotype mismatch between %s and %s" % (variant, other.get(variant)))

    def compute_similarity_to_other_variants_unsorted(self, other_variants):
        n_identical = 0
        for variant in self:
            if other_variants.has_variant(variant):
                n_identical += 1
        return n_identical / len(self)

    def compute_similarity_to_other_variants(self, other_variants):
        n_identical = 0
        for variant, other_variant in zip(self, other_variants):
            assert variant.position == other_variant.position
            if variant.genotype == other_variant.genotype:
                n_identical += 1

        return n_identical / len(self)