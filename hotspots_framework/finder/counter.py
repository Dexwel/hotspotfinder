from collections import defaultdict

import bgreference as bgref


class MutationCounter:

    def __init__(self, genome):
        self.genome = genome

        self.cohort_to_sample = defaultdict(set)
        self.cohort_total_mutations = defaultdict(lambda: defaultdict(int))
        self.mutations_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        self.hotspots_samples = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
        self.original_reference = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))

        self.reference_mismatch = 0

    def add_mutation(self, chromosome, position, ref, alt, sample, cohort):
        chr_position = f'{chromosome}_{position}'
        if ref != alt:
            # Read substitutions of any length
            if ref != '-' and alt != '-':
                if len(alt) == 1:
                    # Check reference
                    if ref == bgref.refseq(self.genome, chromosome, position, 1):
                        self.mutations_dict['snv'][sample][chr_position].append(alt)
                        self.cohort_total_mutations['snv'][cohort] += 1
                        self.hotspots_samples[cohort]['snv'][chr_position].add(sample)
                    else:
                        self.reference_mismatch += 1
                        return
                else:
                    # TODO reference is not checked
                    self.mutations_dict['mnv'][sample][chr_position].append(alt)
                    self.cohort_total_mutations['mnv'][cohort] += 1
                    self.hotspots_samples[cohort]['mnv'][chr_position].add(sample)

            # Read indels of any length
            else:
                # Insertions
                if ref == '-':
                    self.mutations_dict['ins'][sample][chr_position].append(alt)
                    self.cohort_total_mutations['ins'][cohort] += 1
                    self.hotspots_samples[cohort]['ins'][chr_position].add(sample)
                # Deletion
                elif alt == '-':
                    # TODO reference is not checked
                    self.mutations_dict['del'][sample][chr_position].append(alt)
                    self.cohort_total_mutations['del'][cohort] += 1
                    self.hotspots_samples[cohort]['del'][chr_position].add(sample)
                    self.original_reference['del'][sample][chr_position].add(ref)

        # Keep track of samples from each group
        self.cohort_to_sample[cohort].add(sample)
        # Keep track of how many mutations each group has
        self.cohort_total_mutations['total'][cohort] += 1

    def get_samples_per_cohort(self):
        for cohort, set_of_samples in self.cohort_to_sample.items():
            yield cohort, set_of_samples

    def get_samples_per_mutation(self, mut, muttype, cohort):
        return self.hotspots_samples[cohort][muttype][mut]

    def get_alternates_per_mutation(self, mut, muttype, sample):
        return self.mutations_dict[muttype][sample][mut]

    def get_mutations(self, sample, muttype):
        for chr_pos, alts in self.mutations_dict[muttype][sample].items():
            yield chr_pos, alts

    def n_samples_cohort(self, cohort):
        return len(self.cohort_to_sample[cohort])

    def n_mutations_cohort(self, cohort, type_=None):
        return self.cohort_total_mutations[type_ or 'total'][cohort]

    def n_mutations(self):
        for muttype, data in self.cohort_total_mutations.items():
            for cohort, nmuts in data.items():
                yield cohort, muttype, nmuts

    def reference(self, mut, sample):
        return self.original_reference['del'][sample][mut]

