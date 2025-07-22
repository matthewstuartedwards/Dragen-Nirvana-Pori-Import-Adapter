def determineZygosity(genotype):
    """
    Determine the zygosity from a genytype string.
    The genotype is expected to be in the format '0/0', '0/1', '1/0', or '1/1'.
    Returns 'hom' for homozygous and 'het' for heterozygous.
    If the genotype is not recognized, it returns an empty string.
    """
    if genotype == '0/0':
        return 'hom'
    elif genotype == '0/1' or genotype == '1/0':
        return 'het'
    elif genotype == '1/1':
        return 'hom'
    else:
        return ''