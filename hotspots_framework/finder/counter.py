from collections import defaultdict

import bgreference as bgref


class MutationCounter:

    def __init__(self, genome):
        self.genome = genome

        self.samples_and_alternates = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
        self.mutations_by_sample = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: set())))

        self.original_reference = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

        self.reference_mismatch = 0

    def add_mutation(self, chromosome, position, ref, alt, sample, cohort):
        chr_position = f'{chromosome}_{position}'
        # Read substitutions of any length
        if ref != '-' and alt != '-':
            if len(alt) == 1:
                muttype = 'snv'
                # Check reference
                if ref != bgref.refseq(self.genome, chromosome, position, 1):
                    self.reference_mismatch += 1
                    return
            else:
                muttype = 'mnv'
                # TODO reference is not checked

        # Read indels of any length
        else:
            # Insertions
            if ref == '-':
                muttype = 'ins'
            # Deletion
            elif alt == '-':
                muttype = 'del'
                # TODO reference is not checked
                # TODO check that all alternates have the same ref when adding
                self.original_reference['del'][sample][chr_position].add(ref)
            else:
                # TODO is this possible?
                return

        self.samples_and_alternates[cohort][muttype][chr_position][sample].append(alt)
        self.mutations_by_sample[cohort][sample][muttype].add(chr_position)

    def get_cohorts(self):
        return self.mutations_by_sample.keys()

    def get_samples(self, cohort):
        return self.mutations_by_sample[cohort].keys()

    def discard_mutation(self, cohort, sample, muttype, mut):
        self.samples_and_alternates[cohort][muttype][mut].pop(sample)
        self.mutations_by_sample[cohort][sample][muttype].discard(mut)
        if len(self.mutations_by_sample[cohort][sample][muttype]) == 0:
            self.mutations_by_sample[cohort][sample].pop(muttype)
            if len(self.mutations_by_sample[cohort][sample]) == 0:
                self.mutations_by_sample[cohort].pop(sample)
                if len(self.mutations_by_sample[cohort]) == 0:
                    self.mutations_by_sample[cohort].pop(cohort)

    def get_samples_per_mutation(self, mut, muttype, cohort):
        return set(self.samples_and_alternates[cohort][muttype][mut].keys())

    def get_alternates_per_mutation(self, mut, muttype, cohort):
        return self.samples_and_alternates[cohort][muttype][mut]

    def get_mutations_and_alternates(self, cohort, sample):
        for muttype in ('snv', 'mnv', 'ins', 'del'):
            for mut in self.mutations_by_sample[cohort][sample][muttype]:
                yield mut, muttype, self.samples_and_alternates[cohort][muttype][mut][sample]

    def n_mutations_cohort(self, cohort, muttype):
        n_muts = 0
        for sample in self.mutations_by_sample[cohort]:
            for mut in self.mutations_by_sample[cohort][sample][muttype]:
                n_muts += len(self.samples_and_alternates[cohort][muttype][mut][sample])
        return n_muts

    def reference(self, mut, sample):
        return self.original_reference['del'][sample][mut]
