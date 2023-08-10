import pandas as pd
import numpy as np
import re


def get_localization_marker_information(data, uniprot):
    gene_list = []
    protein_list = []
    subcellular_location_list = []

    for i in data.index:
        protein_name = re.sub("^[^\\|]*\\|([^\\|]+)\\|.*$", "\\1", data['protein_a'][i])

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
        protein_name = re.sub("^[^\\|]*\\|([^\\|]+)\\|.*$", "\\1", data['protein_a'][i])

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

    for i in data.index:
        gene = data.iloc[i]['gene_a']
        if gene == 'MT-ATP8' or gene == 'MT-ND1':
            continue
        protein = re.sub("^[^\\|]*\\|([^\\|]+)\\|.*$", "\\1", data['protein_a'][i])

        if gene in gene_helper_list:
            continue

        gene_helper_list.append(gene)
        # get all crosslinks as list
        crosslinks = data.iloc[i]['crosslinks_ab']
        crosslinks_split = crosslinks.split('#')

        # get all crosslinks for this gene by adding crosslinks where its gene_b
        gene_b_rows = data.loc[data['gene_b'] == gene]
        for n in list(range(gene_b_rows.shape[0])):
            xlink = gene_b_rows.iloc[n]['crosslinks_ab']
            xlink_split = xlink.split('#')

            for o in xlink_split:
                sp = o.split('-')
                if gene == sp[0] and o not in crosslinks_split:
                    crosslinks_split.append(o)
                elif gene == sp[2]:
                    reversed_link = sp[2] + '-' + sp[3] + '-' + sp[0] + '-' + sp[1]
                    if reversed_link not in crosslinks_split:
                        crosslinks_split.append(reversed_link)


        crosslinks_inter = []
        # get inter-links
        for m in crosslinks_split:
            xl = m.split('-')
            if xl[0] != xl[2]:
                crosslinks_inter.append('-'.join(xl))

        # TODO maybe only add lm if its more than just Mitochondrion
        if lm_data['gene'].eq(gene).any():
            # TODO only add if location is mitochondrion and not just mitochondrion
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
                    # TODO check for the before element in list
                    # also before and after will be the same if there are multiple transmembrane regions
                    if int(link[1]) < transmem_start:
                        crosslinks_before_tm.append(k)
                    # TODO check for the next element in list
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
                transmembrane_list.extend((np.nan,transmem_regions[0],np.nan))

            else:
                for j in range(len(transmem_regions)):
                    transmem_all = transmem_regions[j].split('..')
                    transmem_start = int(transmem_all[0])
                    transmem_end = int(transmem_all[1])

                    crosslinks_before_tm = []
                    crosslinks_after_tm = []
                    crosslinks_in_tm = []
                    # get all crosslinks

                    if transmem_regions[j] == transmem_regions[0]:
                        for k in crosslinks_inter:
                            link = k.split('-')
                            if link[1] == 'ATP6' or link[1] == 'ND1':
                                continue
                            if int(link[1]) < transmem_start:
                                crosslinks_before_tm.append(k)
                            elif transmem_end >= int(link[1]) >= transmem_start:
                                crosslinks_in_tm.append(k)
                    elif transmem_regions[j] != transmem_regions[0] and transmem_regions[j] != transmem_regions[
                        len(transmem_regions)-1]:
                        transmem_before = transmem_regions[j - 1].split('..')
                        transmem_before_end = int(transmem_before[1])

                        transmem_after = transmem_regions[j + 1].split('..')
                        transmem_after_start = int(transmem_after[1])

                        for k in crosslinks_inter:
                            link = k.split('-')
                            if link[1] == 'ATP6' or link[1] == 'ND1':
                                continue
                            # TODO check for the before element in list
                            # also before and after will be the same if there are multiple transmembrane regions
                            if transmem_start > int(link[1]) > transmem_before_end:
                                crosslinks_before_tm.append(k)
                            # TODO check for the next element in list
                            elif (int(link[1]) > transmem_end) and (int(link[1]) < transmem_after_start):
                                crosslinks_after_tm.append(k)
                            elif transmem_end >= int(link[1]) >= transmem_start:
                                crosslinks_in_tm.append(k)
                    elif transmem_regions[j] == transmem_regions[len(transmem_regions)-1]:
                        transmem_before = transmem_regions[j - 1].split('..')
                        transmem_before_end = int(transmem_before[1])

                        for k in crosslinks_inter:
                            link = k.split('-')
                            if link[1] == 'ATP6' or link[1] == 'ND1':
                                continue
                            # TODO check for the before element in list
                            # also before and after will be the same if there are multiple transmembrane regions
                            if transmem_start > int(link[1]) > transmem_before_end:
                                crosslinks_before_tm.append(k)
                            # TODO check for the next element in list
                            elif int(link[1]) > transmem_end:
                                crosslinks_after_tm.append(k)
                            elif transmem_end >= int(link[1]) >= transmem_start:
                                crosslinks_in_tm.append(k)
                    # concatenate list and add them
                    gene_list.extend((gene, gene,gene))
                    protein_list.extend((protein, protein,protein))
                    subcellular_location_list.extend((location,location,location))
                    crosslinks_list.extend(
                        ('#'.join(crosslinks_before_tm), '#'.join(crosslinks_in_tm), '#'.join(crosslinks_after_tm)))
                    transmembrane_list.extend((np.nan,transmem_regions[j],np.nan))
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


if __name__ == '__main__':
    data = pd.read_csv('../protein_location_prediction_local/SS_unique_lys_crosslink_targetonly_2FDR_separate_all_exported_XlinkCyNET.csv', sep=' ')
    uniprot = pd.read_csv('../protein_location_prediction_local/uniprotkb_AND_reviewed_true_AND_model_o_2023_07_18.tsv', sep='\t', header=0)

    lm_data = get_localization_marker_information(data, uniprot)
    lm_data = lm_data.drop_duplicates(subset=['gene'], keep='first')

    transmem_data = get_transmembrane_information(data, uniprot)
    transmem_data = transmem_data.drop_duplicates(subset=['gene'], keep='first')

    # combine and extend transmembrane proteins to multiple row
    combined_data = combine_lm_transmem_and_proteins(data, lm_data, transmem_data)

    combined_data.to_csv('combined_data.csv',index=False)
    print('done')
