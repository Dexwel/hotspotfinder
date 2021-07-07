from collections import defaultdict

import bgreference as bgref


class MutationCounter:
    """Class to parse and filter mutations"""

    def __init__(self, genome):

        self.nucleotides = 'ACGT'
        self.genome = genome

        self.samples_and_alternates = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list))))
        self.mutations_by_sample = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: set())))

        self.original_reference = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

        self.reference_mismatch = 0
        self.unknown_nucleotides_in_mutation = 0
        self.unknown_nucleotides_in_context = 0

    def add_mutation(self, chromosome, position, ref, alt, alt_type, sample, cohort):
        """Filter unique mutation"""

        chr_position = f'{chromosome}_{position}'

        # Read substitutions
        if alt_type == 'snp':
            muttype = 'snv'
            # Reference filter
            if ref != bgref.refseq(self.genome, chromosome, position, 1):
                self.reference_mismatch += 1
                return
            # Alternate filter
            if alt not in self.nucleotides:
                self.unknown_nucleotides_in_mutation += 1
                return
            # Context filter
            pentamer_sequence = bgref.refseq(self.genome, chromosome, int(position) - 2, 5)
            if len(pentamer_sequence) != len([i for i in pentamer_sequence if i in self.nucleotides]):
                self.unknown_nucleotides_in_context += 1
                return

        # Read MNVs
        elif alt_type == 'mnp':
            muttype = 'mnv'
            # Reference filter
            if ref != bgref.refseq(self.genome, chromosome, position, len(ref)):
                self.reference_mismatch += 1
                return
            # Alternate filter
            if len(alt) != len([i for i in alt if i in self.nucleotides]):
                self.unknown_nucleotides_in_mutation += 1
                return
            # Context filter
            pentamer_sequence = bgref.refseq(self.genome, chromosome, int(position) - 2, 5)
            if len(pentamer_sequence) != len([i for i in pentamer_sequence if i in self.nucleotides]):
                self.unknown_nucleotides_in_context += 1
                return

        # Read indels
        elif alt_type == 'indel':

            # Simple insertions (e.g, G>GA, G>GAAA)
            if ref == '-':
                muttype = 'ins'
                # Alternate filter
                if len(alt) != len([i for i in alt if i in self.nucleotides]):
                    self.unknown_nucleotides_in_mutation += 1
                    return
                # Context filter
                pentamer_sequence = bgref.refseq(self.genome, chromosome, int(position) - 2, 5)
                if len(pentamer_sequence) != len([i for i in pentamer_sequence if i in self.nucleotides]):
                    self.unknown_nucleotides_in_context += 1
                    return

            # Simple deletions (e.g, GT>G, GTG>G)
            elif alt == '-':
                muttype = 'del'
                # Reference filter
                if ref != bgref.refseq(self.genome, chromosome, position, len(ref)):
                    self.reference_mismatch += 1
                    return
                # Alternate filter
                if len(ref) != len([i for i in ref if i in self.nucleotides]):
                    self.unknown_nucleotides_in_mutation += 1
                    return
                # Context filter
                pentamer_sequence = bgref.refseq(self.genome, chromosome, int(position) - 2, 5)
                if len(pentamer_sequence) != len([i for i in pentamer_sequence if i in self.nucleotides]):
                    self.unknown_nucleotides_in_context += 1
                    return
                # todo remove original_reference
                self.original_reference['del'][sample][chr_position].add(ref)
                alt = ref
            else:
                # TODO implement parsing of complex indels (GTG>GA, GTG>GAAA)?
                return
        else:
            return

        self.samples_and_alternates[cohort][muttype][chr_position][sample].append(alt)
        self.mutations_by_sample[cohort][sample][muttype].add(chr_position)

    def get_cohorts(self):
        return self.mutations_by_sample.keys()

    def get_samples(self, cohort):
        return self.mutations_by_sample[cohort].keys()

    def discard_mutation(self, cohort, sample, muttype, mut):
        if sample in self.samples_and_alternates[cohort][muttype][mut].keys():
            self.samples_and_alternates[cohort][muttype][mut].pop(sample)  # Error here
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
