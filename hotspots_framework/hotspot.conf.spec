genome = string(default='hg38')
cores = integer(default=None)
output_format = option('tsv', 'tsv.gz' default='tsv.gz')

mappable_regions = string()
blacklisted_regions = string()
population_variants = string()

genomic_elements = string(default='all')
remove_nonannotated_hotspots = boolean(default=True)

[finder]
mutations_cutoff = integer(default=3)
# TODO remove "remove"
remove_unknown_reference_nucleotides = boolean(default=True)
split_alternates = boolean(default=False)
groupby = option('none', 'GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE', 'PLATFORM', default='none')


