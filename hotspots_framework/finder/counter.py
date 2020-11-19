from collections import defaultdict

import bgreference as bgref


class MutationCounter:

    def __init__(self, genome):
        self.genome = genome

        self.n_mutations_by_cohort = defaultdict(lambda: defaultdict(int))
        self.samples_and_alternates = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))

        self.cohort_to_sample = defaultdict(set)
        self.mutations_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.original_reference = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

        self.reference_mismatch = 0

    def add_mutation(self, chromosome, position, ref, alt, sample, cohort):
        chr_position = f'{chromosome}_{position}'
        # Read substitutions of any length
        if ref != '-' and alt != '-':
            if len(alt) == 1:
                muttype = 'snv'
                # Check reference
                if ref == bgref.refseq(self.genome, chromosome, position, 1):
                    self.mutations_dict['snv'][sample][chr_position].append(alt)
                else:
                    self.reference_mismatch += 1
                    return
            else:
                muttype = 'mnv'
                # TODO reference is not checked
                self.mutations_dict['mnv'][sample][chr_position].append(alt)

        # Read indels of any length
        else:
            # Insertions
            if ref == '-':
                muttype = 'ins'
                self.mutations_dict['ins'][sample][chr_position].append(alt)
            # Deletion
            elif alt == '-':
                muttype = 'del'
                # TODO reference is not checked
                self.mutations_dict['del'][sample][chr_position].append(alt)
                # TODO check that all alternates have the same ref when adding
                self.original_reference['del'][sample][chr_position].add(ref)
            else:
                # TODO is this possible?
                return

        self.samples_and_alternates[cohort][muttype][chr_position][sample].append(alt)
        self.n_mutations_by_cohort[cohort][muttype] += 1
        # Keep track of samples from each group
        self.cohort_to_sample[cohort].add(sample)

    def get_samples_per_cohort(self):
        for cohort, set_of_samples in self.cohort_to_sample.items():
            yield cohort, set_of_samples

    def get_samples_per_mutation(self, mut, muttype, cohort):
        return set(self.samples_and_alternates[cohort][muttype][mut].keys())

    def get_alternates_per_mutation(self, mut, muttype, cohort):
        return self.samples_and_alternates[cohort][muttype][mut]

    def get_mutations(self, sample, muttype):
        return self.mutations_dict[muttype][sample]

    def n_samples_cohort(self, cohort):
        return len(self.cohort_to_sample[cohort])

    def n_mutations_cohort(self, cohort, muttype=None):
        if muttype is None:
            # Compute the total
            n_muts = sum(self.n_mutations_by_cohort[cohort].values())
        else:
            n_muts = self.n_mutations_by_cohort[cohort][muttype]
        return n_muts

    def n_mutations(self):
        for cohort, data in self.n_mutations_by_cohort.items():
            for muttype, nmuts in data.items():
                yield cohort, muttype, nmuts
            yield cohort, 'total', self.n_mutations_cohort(cohort)

    def reference(self, mut, sample):
        return self.original_reference['del'][sample][mut]

