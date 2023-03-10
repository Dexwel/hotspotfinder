genome = string(default='hg38')
cores = integer(default=None)
output_format = option('tsv', 'tsv.gz', default='tsv.gz')

mappable_regions = string()
blacklisted_regions = string()
population_variants = string()
repeats = string()
ig_tr_regions = string()
genomic_elements = string()

[finder]
samples_cutoff = integer(default=2)
mutations_cutoff = integer(default=2)
remove_nonannotated_hotspots = boolean(default=True)
groupby = option('GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE', 'PLATFORM', 'TYPE', default=None)
split_alternates = boolean(default=False)
annotate = boolean(default=True)


