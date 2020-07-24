[genome]
build = string(default='hg38')

[hotspot_mutations]
cutoff = integer(default=3)

[reference_nucleotides]
remove_unknowns = boolean(default=True)

[group]
groupby = option('none', 'GROUP', 'GROUP_BY', 'COHORT' and 'CANCER_TYPE' default='none')

[alternates]
split = boolean(default=False)

[mappability]
mappable_regions = string(default='bgdata')
blacklisted_regions = string(default='bgdata')

[polymorphisms]
population_variants = string(default='bgdata')

[genomic_regions]
genomic_elements = string(default='all')
remove_unnanotated_hotspots = boolean(default=False)

[settings]
output_format = option('tsv', 'vcf', default='tsv')
gzip = boolean(default=False)
cores = integer(default=None)