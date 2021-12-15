def sort_bgc(product):
    """Sort BGC by its type. Uses AntiSMASH annotations
    (see https://docs.antismash.secondarymetabolites.org/glossary/#cluster-types)"""
    
    # TODO: according with current (2021-05) antiSMASH rules:
    # prodigiosin and PpyS-KS -> PKS
    # CDPS -> NRPS
    pks1_products = {'t1pks', 'T1PKS'}
    pksother_products = {'transatpks', 't2pks', 't3pks', 'otherks', 'hglks', 
                     'transAT-PKS', 'transAT-PKS-like', 'T2PKS', 'T3PKS', 
                     'PKS-like', 'hglE-KS'}
    nrps_products = {'nrps', 'NRPS', 'NRPS-like', 'thioamide-NRP', 'NAPAA'}
    ripps_products = {'lantipeptide', 'thiopeptide', 'bacteriocin', 'linaridin', 
                     'cyanobactin', 'glycocin', 'LAP', 'lassopeptide', 
                     'sactipeptide', 'bottromycin', 'head_to_tail', 'microcin', 
                     'microviridin', 'proteusin', 'lanthipeptide', 'lipolanthine', 
                     'RaS-RiPP', 'fungal-RiPP', 'TfuA-related', 'guanidinotides', 
                     'RiPP-like', 'lanthipeptide-class-i', 'lanthipeptide-class-ii', 
                     'lanthipeptide-class-iii', 'lanthipeptide-class-iv',
                     'lanthipeptide-class-v', 'ranthipeptide', 'redox-cofactor'
                     'thioamitides' 'epipeptide', 'cyclic-lactone-autoinducer',
                     'spliceotide', 'RRE-containing'}
    saccharide_products = {'amglyccycl', 'oligosaccharide', 'cf_saccharide', 
                     'saccharide'}
    others_products = {'acyl_amino_acids', 'arylpolyene', 'aminocoumarin', 
                       'ectoine', 'butyrolactone', 'nucleoside', 'melanin', 
                       'phosphoglycolipid', 'phenazine', 'phosphonate', 'other', 
                       'cf_putative', 'resorcinol', 'indole', 'ladderane', 
                       'PUFA', 'furan', 'hserlactone', 'fused', 'cf_fatty_acid', 
                       'siderophore', 'blactam', 'fatty_acid', 'PpyS-KS', 'CDPS', 
                       'betalactone', 'PBDE', 'tropodithietic-acid', 'NAGGN', 
                       'halogenated', 'pyrrolidine'}
    
    # PKS_Type I
    if product in pks1_products:
        return("PKSI")
    # PKS Other Types
    elif product in pksother_products:
        return("PKSother")
    # NRPs
    elif product in nrps_products:
        return("NRPS")
    # RiPPs
    elif product in ripps_products:
        return("RiPPs")
    # Saccharides
    elif product in saccharide_products:
        return("Saccharides")
    # Terpenes
    elif product == 'terpene':
        return("Terpene")
    # PKS/NRP hybrids
    elif len(product.split(".")) > 1:
        #print("  Possible hybrid: (" + cluster + "): " + product)
        # cf_fatty_acid category contains a trailing empty space
        subtypes = set(s.strip() for s in product.split("."))
        if len(subtypes - (pks1_products | pksother_products | nrps_products)) == 0:
            if len(subtypes - nrps_products) == 0:
                return("NRPS")
            elif len(subtypes - (pks1_products | pksother_products)) == 0:
                return("PKSother") # pks hybrids
            else:
                return("PKS-NRP_Hybrids")
        elif len(subtypes - ripps_products) == 0:
            return("RiPPs")
        elif len(subtypes - saccharide_products) == 0:
            return("Saccharide")
        else:
            return("Others") # other hybrid
    # Others
    elif product in others_products:
        return("Others")
    # ??
    elif product == "":
        # No product annotation. Perhaps not analyzed by antiSMASH
        return("Others")
    else:
        print("  Warning: unknown product '{}'".format(product))
        return("Others")

def get_composite_bgc_similarities(bgcs_1, bgcs_2, sim_matrix):
    num_pairs = 0
    sum_sim = 0.00
    min_sim = (1.00, -1, -1)
    max_sim = (0.00, -1, -1)
    for bgc_1 in bgcs_1:
        for bgc_2 in bgcs_2:
            sim = 0.00 if (bgc_1 == bgc_2) else (sim_matrix[bgc_1][bgc_2] if ((bgc_1 in sim_matrix) and (bgc_2 in sim_matrix[bgc_1])) else sim_matrix[bgc_2][bgc_1])
            sum_sim += sim
            if (sim < min_sim[0]):
                min_sim = (sim, bgc_1, bgc_2) if (bgc_2 > bgc_1) else (sim, bgc_2, bgc_1)            
            if (sim > max_sim[0]):
                max_sim = (sim, bgc_1, bgc_2) if (bgc_2 > bgc_1) else (sim, bgc_2, bgc_1)
            num_pairs += 1
    return ((sum_sim / num_pairs), min_sim, max_sim) # let it throw infinite division error by itself


def BGC_dic_gen(filtered_matrix):
    """Generates the: { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } part of the BGCs variable."""
    
    # It would be tempting to do `bgc_dict = defaultdict(list)` but then the 
    # dictionary keeps the defaultdict definition. Later on in the distance
    # calculation function, we do a 
    #  ```try: unshared_occurrences = BGCs[A][unshared_domain]
    #  except KeyError: unshared_occurrences = BGCs[B][unshared_domain]```
    # Problem is the `try` would always be successful because of defaultdict and
    # an empty list is added for the unshared_domain!
    bgc_dict = {}
    for row in filtered_matrix:
        header = row[-1] + ":" + row[3] + ":" + row[4] # add domain positions
        try: #Should be faster than performing `if key in dictionary.keys()`
            bgc_dict[row[5]].append(header)
        except KeyError: # First occurrence of this domain in the cluster
            bgc_dict[row[5]] = [header]
            
    return bgc_dict
