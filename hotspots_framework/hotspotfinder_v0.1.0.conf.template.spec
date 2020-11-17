[genome]
build = string(default='hg38')

[hotspot_mutations]
cutoff = integer(default=3)

[reference_nucleotides]
remove_unknowns = boolean(default=True)

[group]
groupby = option('none', 'GROUP', 'GROUP_BY', 'COHORT', 'CANCER_TYPE', 'PLATFORM', default='none')

[alternates]
split = boolean(default=False)

[mappability]
mappable_regions = string()
blacklisted_regions = string()

[polymorphisms]
population_variants = string()

[genomic_regions]
genomic_elements = string(default='all')
remove_nonannotated_hotspots = boolean(default=True)

[settings]
output_format = option('tsv', default='tsv')
gzip = boolean(default=False)
cores = integer(default=None)