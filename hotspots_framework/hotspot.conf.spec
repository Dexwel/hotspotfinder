genome = string(default='hg38')
cores = integer(default=None)
output_format = option('tsv', 'tsv.gz', default='tsv.gz')

mappable_regions = string()
blacklisted_regions = string()
population_variants = string()

genomic_elements = string(default='all')


[finder]
# TODO can we add a minimum number? 3
samples_cutoff = integer(default=3)
mutations_cutoff = integer(default=3)
# TODO remove "remove"
remove_unknown_reference_nucleotides = boolean(default=True)
remove_nonannotated_hotspots = boolean(default=True)
groupby = option('GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE', 'PLATFORM', default=None)
split_alternates = boolean(default=False)
annotate = boolean(default=True)


