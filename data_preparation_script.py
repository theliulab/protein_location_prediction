import pandas as pd
import numpy as np
import re


def get_localization_marker_information(data, uniprot):
    gene_list = []
    protein_list = []
    subcellular_location_list = []

    for i in data.index:
        protein_name = re.sub("^[^\\|]*\\|([^\\|]+)\\|.*$", "\\1", data['Protein1'][i])

        # skip if protein not in uniprot list
        if (uniprot['Entry'].eq(protein_name).any()) is False:
            continue

        if len(uniprot.loc[uniprot['Entry'] == protein_name, 'Subcellular location [CC]']) == 0:
            continue

        location_uniprot = uniprot.loc[uniprot['Entry'] == protein_name, 'Subcellular location [CC]'].values[0]

        if pd.isnull(location_uniprot):
            continue

        occ = location_uniprot.count("{")
        if occ > 1:
            continue

        location = re.search("SUBCELLULAR LOCATION: (.+?)(\\.|{)", location_uniprot).group(1)

        # only take mitochondrion proteins as lm
        if 'itoch' in location and 'soform' not in location and 'Mitochondrion' != location.strip():
            # gene_a, protein_a, crosslinks_ab, subcellular_location, transmembrane
            gene_list.append(data.iloc[i]['gene_a'])
            protein_list.append(protein_name)
            subcellular_location_list.append(location)

    new_data = pd.DataFrame({'gene': gene_list, 'protein': protein_list,
                             'subcellular_location': subcellular_location_list})

    return new_data


# function assigning location markers if proteins are in one location and dont have transmembrane region
def get_transmembrane_information(data, uniprot):
    gene_list = []
    protein_list = []
    transmembrane_list = []

    for i in data.index:
        protein_name = re.sub("^[^\\|]*\\|([^\\|]+)\\|.*$", "\\1", data['Protein1'][i])

        # skip if protein not in uniprot list
        if (uniprot['Entry'].eq(protein_name).any()) is False:
            continue

        if len(uniprot.loc[uniprot.Entry == protein_name, 'Transmembrane']) == 0:
            continue

        transmem_uniprot = uniprot.loc[uniprot.Entry == protein_name, 'Transmembrane'].values[0]

        if pd.isnull(transmem_uniprot):
            continue
        else:
            transmem = re.findall("TRANSMEM (.+?)\\;", transmem_uniprot)

            # gene_a, protein_a, crosslinks_ab, subcellular_location, transmembrane
            gene_list.append(data.iloc[i]['gene_a'])
            protein_list.append(protein_name)
            transmembrane_list.append(transmem)

    new_data = pd.DataFrame({'gene': gene_list, 'protein': protein_list,
                             'transmembrane': transmembrane_list})

    return new_data


# function splitting proteins+crosslinks into multiple rows if they have transmembrane regions
def combine_lm_transmem_and_proteins(data, lm_data, transmem_data):
    # if gene in this list continue
    gene_helper_list = []

    gene_list = []
    protein_list = []
    crosslinks_list = []
    subcellular_location_list = []
    transmembrane_list = []

    genes = ['gene_a', 'gene_b']
    proteins = ['Protein1', 'Protein2']

    for i in data.index:
        for q in range(len(genes)):
            gene = data.iloc[i][genes[q]]

            protein = re.sub("^[^\\|]*\\|([^\\|]+)\\|.*$", "\\1", data[proteins[q]][i])

            # get all crosslinks as list
            if (genes[q] == genes[0]):
                crosslinks = data.iloc[i]['crosslinks_a']
            elif (genes[q] == genes[1]):
                crosslinks = data.iloc[i]['crosslinks_b']

            crosslinks_split = crosslinks.split('#')

            if gene in gene_helper_list:
                continue

            gene_helper_list.append(gene)

            # get all crosslinks for this gene by adding crosslinks where its gene_b
           # gene_b_rows = data.loc[data['gene_b'] == gene]

            # if crosslinked gene is not gene_a/gene_b and not in data table, add it
            # gather all its crosslinks for prediction
            # if gene is gene_a, it has all its possible crosslinks. not the case if gene is gene_b
            # for n in list(range(gene_b_rows.shape[0])):
            #     xlink = str(gene_b_rows.iloc[n]['crosslinks_b'])
            #     xlink_split = xlink.split('#')
            #
            #     for o in xlink_split:
            #         sp = o.split('-')
            #
            #         if len(sp) < 4:
            #             print('h')
            #
            #
            #         if gene == sp[0] and o not in crosslinks_split:
            #             crosslinks_split.append(o)
            #         elif gene == sp[2]:
            #             reversed_link = sp[2] + '-' + sp[3] + '-' + sp[0] + '-' + sp[1]
            #             if reversed_link not in crosslinks_split:
            #                 crosslinks_split.append(reversed_link)

            # # if its already in gene list
            # if gene in gene_list:
            #     xlinks_found = crosslinks_list[gene_list.index(gene)]
            #     for n in xlinks_found:
            #         sp = n.split('-')
            #
            #         if gene == sp[0] and n not in crosslinks_split:
            #             crosslinks_split.append(n)
            #         elif gene == sp[2]:
            #             reversed_link = sp[2] + '-' + sp[3] + '-' + sp[0] + '-' + sp[1]
            #             if reversed_link not in crosslinks_split:
            #                 crosslinks_split.append(reversed_link)


            crosslinks_inter = []
            # get inter-links
            for m in crosslinks_split:
                xl = m.split('-')
                if gene == xl[0] and xl[0] != xl[2] and m not in crosslinks_inter:
                    crosslinks_inter.append('-'.join(xl))
                elif gene == xl[2]:
                    rev_xlink = xl[2] + '-' + xl[3] + '-' + xl[0] + '-' + xl[1]
                    if rev_xlink not in crosslinks_inter:
                        crosslinks_inter.append(rev_xlink)

            if lm_data['gene'].eq(gene).any():
                # exclude brackets if they exist in string
                location = lm_data.loc[lm_data.gene == gene, 'subcellular_location'].values[0]
                if location.count("{") > 0:
                    # get anything until bracket
                    re1 = re.compile("(.*?)\s*\(")
                    location = re1.match(location)
                elif 'ytop' in location or 'ucle' in location or 'soform' in location:
                    location = np.nan
            else:
                location = np.nan

            # check if protein has transmembrane region
            # if so, add rows for each crosslinks before, between and after tm
            if transmem_data['gene'].eq(gene).any():
                transmem_regions = transmem_data.loc[transmem_data.gene == gene, 'transmembrane'].values[0]

                crosslinks_before_tm = []
                crosslinks_after_tm = []
                crosslinks_in_tm = []
                if len(transmem_regions) == 1:
                    transmem_all = transmem_regions[0].split('..')
                    transmem_start = int(transmem_all[0])
                    transmem_end = int(transmem_all[1])
                    for k in crosslinks_inter:
                        link = k.split('-')

                        # also before and after will be the same if there are multiple transmembrane regions
                        if int(link[1]) < transmem_start:
                            crosslinks_before_tm.append(k)
                        elif int(link[1]) > transmem_end:
                            crosslinks_after_tm.append(k)
                        elif transmem_end >= int(link[1]) >= transmem_start:
                            crosslinks_in_tm.append(k)
                    # concatenate list and add them
                    gene_list.extend((gene, gene, gene))
                    protein_list.extend((protein, protein, protein))
                    subcellular_location_list.extend((location, location, location))
                    crosslinks_list.extend(
                        ('#'.join(crosslinks_before_tm), '#'.join(crosslinks_in_tm), '#'.join(crosslinks_after_tm)))
                    transmembrane_list.extend((np.nan, transmem_regions[0], np.nan))

                elif len(transmem_regions) == 2:
                    for j in range(len(transmem_regions)):
                        transmem_all = transmem_regions[j].split('..')
                        transmem_start = int(transmem_all[0])
                        transmem_end = int(transmem_all[1])

                        crosslinks_before_tm = []
                        crosslinks_after_tm = []
                        crosslinks_in_tm = []
                        # get all crosslinks
                        # if first tm region, take links before and in tm
                        if transmem_regions[j] == transmem_regions[0]:
                            for k in crosslinks_inter:
                                link = k.split('-')
                                if int(link[1]) < transmem_start:
                                    crosslinks_before_tm.append(k)
                                elif transmem_end >= int(link[1]) >= transmem_start:
                                    crosslinks_in_tm.append(k)

                            # concatenate list and add them
                            gene_list.extend((gene, gene))
                            protein_list.extend((protein, protein))
                            subcellular_location_list.extend((location, location))
                            crosslinks_list.extend(
                                ('#'.join(crosslinks_before_tm), '#'.join(crosslinks_in_tm)))
                            transmembrane_list.extend((np.nan, transmem_regions[j]))
                        elif transmem_regions[j] == transmem_regions[len(transmem_regions) - 1]:
                            transmem_before = transmem_regions[j - 1].split('..')
                            transmem_before_end = int(transmem_before[1])

                            for k in crosslinks_inter:
                                link = k.split('-')
                                # also before and after will be the same if there are multiple transmembrane region
                                if transmem_start > int(link[1]) > transmem_before_end:
                                    crosslinks_before_tm.append(k)
                                elif int(link[1]) > transmem_end:
                                    crosslinks_after_tm.append(k)
                                elif transmem_end >= int(link[1]) >= transmem_start:
                                    crosslinks_in_tm.append(k)

                            gene_list.extend((gene, gene, gene))
                            protein_list.extend((protein, protein, protein))
                            subcellular_location_list.extend((location, location, location))
                            crosslinks_list.extend(('#'.join(crosslinks_before_tm), '#'.join(crosslinks_in_tm),
                                                    '#'.join(crosslinks_after_tm)))
                            transmembrane_list.extend((np.nan, transmem_regions[j], np.nan))

                else:
                    for j in range(len(transmem_regions)):
                        transmem_all = transmem_regions[j].split('..')
                        transmem_start = int(transmem_all[0])
                        transmem_end = int(transmem_all[1])

                        crosslinks_before_tm = []
                        crosslinks_after_tm = []
                        crosslinks_in_tm = []
                        # get all crosslinks
                        # if first tm region, take links before and in tm
                        if transmem_regions[j] == transmem_regions[0]:
                            for k in crosslinks_inter:
                                link = k.split('-')
                                if int(link[1]) < transmem_start:
                                    crosslinks_before_tm.append(k)
                                elif transmem_end >= int(link[1]) >= transmem_start:
                                    crosslinks_in_tm.append(k)

                            # concatenate list and add them
                            gene_list.extend((gene, gene))
                            protein_list.extend((protein, protein))
                            subcellular_location_list.extend((location, location))
                            crosslinks_list.extend(
                                ('#'.join(crosslinks_before_tm), '#'.join(crosslinks_in_tm)))
                            transmembrane_list.extend((np.nan, transmem_regions[j]))
                        # if neither first nor last tm region, take links before, in and after tm
                        elif transmem_regions[j] != transmem_regions[0] and transmem_regions[j] != transmem_regions[
                            len(transmem_regions) - 1]:
                            transmem_before = transmem_regions[j - 1].split('..')
                            transmem_before_end = int(transmem_before[1])

                            transmem_after = transmem_regions[j + 1].split('..')
                            transmem_after_start = int(transmem_after[1])

                            for k in crosslinks_inter:
                                link = k.split('-')
                                # also before and after will be the same if there are multiple transmembrane regions
                                if transmem_start > int(link[1]) > transmem_before_end:
                                    crosslinks_before_tm.append(k)
                                elif transmem_end >= int(link[1]) >= transmem_start:
                                    crosslinks_in_tm.append(k)
                            # concatenate list and add them
                            gene_list.extend((gene, gene))
                            protein_list.extend((protein, protein))
                            subcellular_location_list.extend((location, location))
                            crosslinks_list.extend(('#'.join(crosslinks_before_tm), '#'.join(crosslinks_in_tm)))
                            transmembrane_list.extend((np.nan, transmem_regions[j]))
                        # if last tm, take links in and after tm region
                        elif transmem_regions[j] == transmem_regions[len(transmem_regions) - 1]:
                            transmem_before = transmem_regions[j - 1].split('..')
                            transmem_before_end = int(transmem_before[1])

                            for k in crosslinks_inter:
                                link = k.split('-')
                                # also before and after will be the same if there are multiple transmembrane region
                                if transmem_start > int(link[1]) > transmem_before_end:
                                    crosslinks_before_tm.append(k)
                                elif int(link[1]) > transmem_end:
                                    crosslinks_after_tm.append(k)
                                elif transmem_end >= int(link[1]) >= transmem_start:
                                    crosslinks_in_tm.append(k)

                            gene_list.extend((gene, gene, gene))
                            protein_list.extend((protein, protein, protein))
                            subcellular_location_list.extend((location, location, location))
                            crosslinks_list.extend(('#'.join(crosslinks_before_tm), '#'.join(crosslinks_in_tm),
                                                    '#'.join(crosslinks_after_tm)))
                            transmembrane_list.extend((np.nan, transmem_regions[j], np.nan))

            else:
                # if there are non just add all crosslinks
                gene_list.append(gene)
                protein_list.append(protein)
                crosslinks_list.append('#'.join(crosslinks_inter))
                subcellular_location_list.append(location)
                transmembrane_list.append(np.nan)

    new_data = pd.DataFrame({'gene': gene_list, 'protein': protein_list,
                             'subcellular_location': subcellular_location_list,
                             'crosslinks': crosslinks_list,
                             'transmembrane': transmembrane_list})

    return new_data


def change_crosslinks_for_MTproteins(data):
    gene_a = []
    gene_b = []
    protein_a = []
    protein_b = []
    crosslinks_ab = []

    for i in data.index:
        genea = data.iloc[i]['gene_a']
        geneb = data.iloc[i]['gene_b']

        crosslinks = []
        xlinks_ab = (data.iloc[i]['crosslinks_ab']).split('#')
        xlinks_ba = (data.iloc[i]['crosslinks_ba']).split('#')

        # take the highest number of crosslinks for this gene
        if len(xlinks_ab) > len(xlinks_ba):
            xlinks = xlinks_ab
        elif len(xlinks_ba) > len(xlinks_ab):
            xlinks = xlinks_ba
        else:
            xlinks = xlinks_ab

        if 'MT-' in genea and 'MT-' in geneb:
            genea_split = genea.split('-')
            geneb_split = geneb.split('-')

            for link in xlinks:
                ls = link.split('-')
                if len(ls) > 4:
                    if ls[0] == 'MT' and ls[3] == 'MT':
                        cl = ls[0] + ls[1] + '-' + ls[2] + '-' + ls[3] + ls[4] + '-' + ls[5]
                        crosslinks.append(cl)
                    elif ls[0] == 'MT':
                        cl = ls[0] + ls[1] + '-' + ls[2] + '-' + ls[3] + '-' + ls[4]
                        crosslinks.append(cl)
                    elif ls[3] == 'MT':
                        cl = ls[0] + '-' + ls[1] + '-' + ls[2] + ls[3] + '-' + ls[4]
                        crosslinks.append(cl)
                else:
                    crosslinks.append(link)

            gene_a.append(genea_split[0] + genea_split[1])
            gene_b.append(geneb_split[0] + geneb_split[1])
            protein_a.append(data.iloc[i]['protein_a'])
            protein_b.append(data.iloc[i]['protein_b'])
            crosslinks_ab.append('#'.join(crosslinks))
        elif 'MT-' in genea:
            genea_split = genea.split('-')
            for link in xlinks:
                ls = link.split('-')
                if len(ls) > 4:
                    if ls[0] == 'MT' and ls[3] == 'MT':
                        cl = ls[0] + ls[1] + '-' + ls[2] + '-' + ls[3] + ls[4] + '-' + ls[5]
                        crosslinks.append(cl)
                    elif ls[0] == 'MT':
                        cl = ls[0] + ls[1] + '-' + ls[2] + '-' + ls[3] + '-' + ls[4]
                        crosslinks.append(cl)
                    elif ls[3] == 'MT':
                        cl = ls[0] + '-' + ls[1] + '-' + ls[2] + ls[3] + '-' + ls[4]
                        crosslinks.append(cl)
                else:
                    crosslinks.append(link)

            gene_a.append(genea_split[0] + genea_split[1])
            gene_b.append(data.iloc[i]['gene_b'])
            protein_a.append(data.iloc[i]['protein_a'])
            protein_b.append(data.iloc[i]['protein_b'])
            crosslinks_ab.append('#'.join(crosslinks))

        elif 'MT-' in geneb:
            geneb_split = geneb.split('-')
            for link in xlinks:
                ls = link.split('-')
                if len(ls) > 4:
                    if ls[0] == 'MT' and ls[3] == 'MT':
                        cl = ls[0] + ls[1] + '-' + ls[2] + '-' + ls[3] + ls[4] + '-' + ls[5]
                        crosslinks.append(cl)
                    elif ls[0] == 'MT':
                        cl = ls[0] + ls[1] + '-' + ls[2] + '-' + ls[3] + '-' + ls[4]
                        crosslinks.append(cl)
                    elif ls[3] == 'MT':
                        cl = ls[0] + '-' + ls[1] + '-' + ls[2] + ls[3] + '-' + ls[4]
                        crosslinks.append(cl)
                else:
                    crosslinks.append(link)

            gene_a.append(data.iloc[i]['gene_a'])
            gene_b.append(geneb_split[0] + geneb_split[1])
            protein_a.append(data.iloc[i]['protein_a'])
            protein_b.append(data.iloc[i]['protein_b'])
            crosslinks_ab.append('#'.join(crosslinks))
        else:
            gene_a.append(genea)
            gene_b.append(geneb)
            protein_a.append(data.iloc[i]['protein_a'])
            protein_b.append(data.iloc[i]['protein_b'])
            crosslinks_ab.append(data.iloc[i]['crosslinks_ab'])

    new_data = pd.DataFrame({'gene_a': gene_a, 'gene_b': gene_b,
                             'protein_a': protein_a,
                             'protein_b': protein_b,
                             'crosslinks_ab': crosslinks_ab})

    return new_data


if __name__ == '__main__':
    # adjust your data paths
    data = pd.read_csv(
        'lysine_xlinks_xlilo.csv',sep=';',dtype='str')
    uniprot = pd.read_csv('../protein_location_prediction_local/uniprotkb_AND_reviewed_true_AND_model_o_2023_07_18.tsv',
                          sep='\t', header=0)

    lm_data = get_localization_marker_information(data, uniprot)
    lm_data = lm_data.drop_duplicates(subset=['gene'], keep='first')

    transmem_data = get_transmembrane_information(data, uniprot)
    transmem_data = transmem_data.drop_duplicates(subset=['gene'], keep='first')

    #updated_data = change_crosslinks_for_MTproteins(data)

    # combine and extend transmembrane proteins to multiple row
    combined_data = combine_lm_transmem_and_proteins(data, lm_data, transmem_data)

    combined_data.to_csv('combined_data.csv', index=False)
    print('done')
