import dill


# Wrapper for all indexes required by kage genotyping
class IndexBundle:
    index_names = ["VariantToNodes", "NumpyVariants", "NodeCountModelAdvanced", "HelperVariants", "CombinationMatrix", "TrickyVariants", "KmerIndex"]
    def __init__(self, indexes):
        self.indexes = indexes

    @classmethod
    def from_file(cls, file_name, skip=None):
        with open(file_name, 'rb') as f:
            return dill.load(f)

    def to_file(self, file_name, compress=True):
        with open(file_name, 'wb') as f:
            dill.dump(self, f)
