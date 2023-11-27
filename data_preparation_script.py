"""
@author: Julia Ruta
This script extracts topology information and subcellular location
for mitochondrion proteins.
Proteins either need to have a unique subcellular location or in
case of mitochondrion membrane proteins unique locations for
topological domains in order to be qualified as a localization marker.
"""


import sys
import pandas as pd
import numpy as np
import re
import collections

def get_localization_marker_information(data, uniprot):
    """
    method extracting subcellular location from uniprot for
    mitochdrion proteins
    :param data:
    :param uniprot:
    :return: new_data
    """
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

        location = re.search("SUBCELLULAR LOCATION: (.+?)(\\.|{)", location_uniprot).group(1)
        # only take mitochondrion proteins as lm
        if 'itoch' in location and 'soform' not in location and 'Mitochondrion' != location.strip():
            occ = location_uniprot.count("{")
            if occ > 1:
                # only allowed to have multiple locations if its a membrane protein
                if 'membrane' in location:
                    gene_list.append(data.iloc[i]['gene_a'])
                    protein_list.append(protein_name)
                    subcellular_location_list.append(location)
            # gene_a, protein_a, crosslinks_ab, subcellular_location, transmembrane
            else:
                gene_list.append(data.iloc[i]['gene_a'])
                protein_list.append(protein_name)
                subcellular_location_list.append(location)

    new_data = pd.DataFrame({'gene': gene_list, 'protein': protein_list,
                             'subcellular_location': subcellular_location_list})

    return new_data


# function assigning location markers if proteins are in one location and dont have transmembrane region
def get_transmembrane_information(data, uniprot):
    """
    method extracting transmembrane information for proteins and in
    case of a membrane proteins splitting it into rows depending on
    the number of transmembrane domains
    :param data:
    :param uniprot:
    :return: new_data
    """
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
    """
    method combining subcellular location and transmembrane domains
    for proteins; a unique subcellular location signals a localization
    marker
    :param data:
    :param lm_data:
    :param transmem_data:
    :return: new date
    """
    # if gene in this list continue
    gene_helper_list = []

    gene_list = []
    protein_list = []
    crosslinks_list = []
    subcellular_location_list = []
    transmembrane_list = []
    is_lm_list = []

    genes = ['gene_a', 'gene_b']
    proteins = ['Protein1', 'Protein2']

    for i in data.index:
        for q in range(len(genes)):
            is_lm = ""
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

            # if crosslinked gene is not gene_a/gene_b and not in data table, add it
            # gather all its crosslinks for prediction
            # if gene is gene_a, it has all its possible crosslinks. not the case if gene is gene_b

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
                is_lm = "true"
            else:
                location = np.nan
                is_lm = "false"

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
                    is_lm_list.extend((is_lm,is_lm,is_lm))

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
                            is_lm_list.extend((is_lm,is_lm))
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
                            is_lm_list.extend((is_lm,is_lm,is_lm))

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
                            is_lm_list.extend((is_lm,is_lm))
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
                            is_lm_list.extend((is_lm,is_lm))
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
                            is_lm_list.extend((is_lm,is_lm,is_lm))

            else:
                # if there are non just add all crosslinks
                gene_list.append(gene)
                protein_list.append(protein)
                crosslinks_list.append('#'.join(crosslinks_inter))
                subcellular_location_list.append(location)
                transmembrane_list.append(np.nan)
                is_lm_list.append(is_lm)

    new_data = pd.DataFrame({'gene': gene_list, 'protein': protein_list,
                             'subcellular_location': subcellular_location_list,
                             'crosslinks': crosslinks_list,
                             'transmembrane': transmembrane_list,
                             'is_localization_marker': is_lm_list})

    return new_data

def add_topology_information(data,uniprot):
    """
    add topology information for mitochondrion membrane proteins
    :param data:
    :param uniprot:
    :return: new_data
    """
    protein_helper_list = []

    gene_list = []
    protein_list = []
    crosslinks_list = []
    subcellular_location_list = []
    topology_list = []
    transmembrane_list = []
    is_lm_list = []

    data['transmembrane'] = data['transmembrane'].fillna("")
    data['subcellular_location'] = data['subcellular_location'].fillna("")

    for i in data.index:
        is_lm = ""
        protein = data.iloc[i]['protein']
        if protein in protein_helper_list:
            continue
        else:
            protein_helper_list.append(protein)

        sub = data.loc[data['protein'] == protein]

        # if no mito protein, continue
        mask = [bool(re.search(p, s)) for p, s in zip('itoch',sub['subcellular_location'])]

        # only add topology for mito proteins with unique subcellular location
        # if no transmembrane regions, dont consider topological domains
        if sum(mask) == 0 or len(sub.index) < 2:
            gene_list.extend(sub['gene'].tolist())
            protein_list.extend(sub['protein'].tolist())
            crosslinks_list.extend(sub['crosslinks'].tolist())
            subcellular_location_list.extend(sub['subcellular_location'].tolist())
            topology_list.extend([""] * len(sub.index))
            transmembrane_list.extend(sub['transmembrane'].tolist())
            is_lm_list.extend(sub['is_localization_marker'].tolist())
            continue

        if (uniprot['Entry'].eq(protein).any()) is False:
            gene_list.extend(sub['gene'].tolist())
            protein_list.extend(sub['protein'].tolist())
            crosslinks_list.extend(sub['crosslinks'].tolist())
            subcellular_location_list.extend(sub['subcellular_location'].tolist())
            topology_list.extend([""]*len(sub.index))
            transmembrane_list.extend(sub['transmembrane'].tolist())
            is_lm_list.extend(sub['is_localization_marker'].tolist())
            continue

        if len(uniprot.loc[uniprot['Entry'] == protein, 'Topological domain']) == 0:
            gene_list.extend(sub['gene'].tolist())
            protein_list.extend(sub['protein'].tolist())
            crosslinks_list.extend(sub['crosslinks'].tolist())
            subcellular_location_list.extend(sub['subcellular_location'].tolist())
            topology_list.extend([""] * len(sub.index))
            transmembrane_list.extend(sub['transmembrane'].tolist())
            is_lm_list.extend(sub['is_localization_marker'].tolist())
            continue

        topology_uniprot = uniprot.loc[uniprot['Entry'] == protein, 'Topological domain'].values[0]

        if pd.isnull(topology_uniprot) or topology_uniprot == "":
            gene_list.extend(sub['gene'].tolist())
            protein_list.extend(sub['protein'].tolist())
            crosslinks_list.extend(sub['crosslinks'].tolist())
            subcellular_location_list.extend(sub['subcellular_location'].tolist())
            topology_list.extend([""] * len(sub.index))
            transmembrane_list.extend(sub['transmembrane'].tolist())
            is_lm_list.extend(sub['is_localization_marker'].tolist())
            continue

        topologies = re.findall("TOPO_DOM (.+?); \/evidence", topology_uniprot)

        # rewrite topologies into data frame
        topo_start = []
        topo_end = []
        topo_locations = []
        topo_is_lm = []
        for j in topologies:
            topo_split = j.split(";")
            if ".." in topo_split[0]:
                topo_start.append(int(topo_split[0].split("..")[0]))
                topo_end.append(int(topo_split[0].split("..")[1]))
                topo_loc = re.findall('"(.*?)"', topo_split[1])[0]
                topo_locations.append(topo_loc)
                if topo_loc.count('note') < 2:
                    topo_is_lm.append("true")
                else:
                    topo_is_lm.append("false")
            else:
                topo_start.append(int(topo_split[0]))
                topo_end.append(sys.maxsize)
                ttopo_loc = re.findall('"(.*?)"', topo_split[1])[0]
                topo_locations.append(topo_loc)
                if topo_loc.count('note') < 2:
                    topo_is_lm.append("true")
                else:
                    topo_is_lm.append("false")

        topology_info = pd.DataFrame({'start': topo_start, 'end': topo_end, 'topo': topo_locations,'is_lm': topo_is_lm})

        # loop over subset of protein and add topology based on transmembrane regions
        #test = range(len(sub.index))
        for k in range(len(sub.index)):
            # dont add topology if it has a transmembrane region in this row
            if sub.iloc[k]['transmembrane'] != "":
                gene_list.append(sub.iloc[k]['gene'])
                protein_list.append(sub.iloc[k]['protein'])
                crosslinks_list.append(sub.iloc[k]['crosslinks'])
                subcellular_location_list.append(sub.iloc[k]['subcellular_location'])
                topology_list.append("")
                transmembrane_list.append(sub.iloc[k]['transmembrane'])
                is_lm_list.append(sub.iloc[k]['is_localization_marker'])
                continue

            if k == len(sub.index):
                gene_list.append(sub.iloc[k-1]['gene'])
                protein_list.append(sub.iloc[k-1]['protein'])
                crosslinks_list.append(sub.iloc[k-1]['crosslinks'])
                subcellular_location_list.append(sub.iloc[k-1]['subcellular_location'])
                topology_list.append(topology_info.iloc[len(topology_info.index)]['topo'])
                transmembrane_list.append(sub.iloc[k-1]['transmembrane'])
                is_lm_list.append(topology_info.iloc[len(topology_info.index)]['is_lm'])
                continue

            # consider cases for first and last row again
            elif k == 0:
                # check if any end of topology region is smaller than following transmembrane start
                transmembrane_start = int((sub.iloc[k+1]['transmembrane']).split('..')[0])
                #res = list(filter(lambda i: i < transmembrane_start, list(topology_info["end"])))[0]
                res = topology_info.index[topology_info["end"] < transmembrane_start].tolist()
                if not res:
                    gene_list.append(sub.iloc[k]['gene'])
                    protein_list.append(sub.iloc[k]['protein'])
                    crosslinks_list.append(sub.iloc[k]['crosslinks'])
                    subcellular_location_list.append(sub.iloc[k]['subcellular_location'])
                    topology_list.append("")
                    transmembrane_list.append(sub.iloc[k]['transmembrane'])
                    is_lm_list.append(sub.iloc[k]['is_localization_marker'])
                else:
                    gene_list.append(sub.iloc[k]['gene'])
                    protein_list.append(sub.iloc[k]['protein'])
                    crosslinks_list.append(sub.iloc[k]['crosslinks'])
                    subcellular_location_list.append(sub.iloc[k]['subcellular_location'])
                    topology_list.append(topology_info.iloc[res[len(res)-1]]['topo'])
                    transmembrane_list.append(sub.iloc[k]['transmembrane'])
                    is_lm_list.append(topology_info.iloc[res[len(res)-1]]['is_lm'])

            elif k != 0 and k != (len(sub.index)-1):
                transmembrane_end = int((sub.iloc[(k-1)]['transmembrane']).split('..')[0])
                transmembrane_start = int((sub.iloc[(k+1)]['transmembrane']).split('..')[1])

                # check if any start of topology region is larger than tm end before
                #res_start = list(filter(lambda i: i > transmembrane_end, list(topology_info["start"])))[0]
                res_start = topology_info.index[topology_info["start"] > transmembrane_end].tolist()
                # check if any end of topology region is smaller than following tm start
                #res_end = list(filter(lambda i: i < transmembrane_start, list(topology_info["end"])))[0]
                res_end = topology_info.index[topology_info["end"] > transmembrane_start].tolist()

                # intersection of residues
                set_start = set(res_start)
                set_end = set(res_end)
                idx = set(res_start)^set(res_end)

                #print(idx)
                if not idx:
                    gene_list.append(sub.iloc[k]['gene'])
                    protein_list.append(sub.iloc[k]['protein'])
                    crosslinks_list.append(sub.iloc[k]['crosslinks'])
                    subcellular_location_list.append(sub.iloc[k]['subcellular_location'])
                    topology_list.append("")
                    transmembrane_list.append(sub.iloc[k]['transmembrane'])
                    is_lm_list.append(sub.iloc[k]['is_localization_marker'])
                    continue
                #idx = collections.Counter(res_start) & collections.Counter(res_end)
                #res = list(idx.elements())
                res = list(idx)
                #print(res)
                # if res_end == res_start then its the same row
                #if int(res_start) == int(res_end):

                gene_list.append(sub.iloc[k]['gene'])
                protein_list.append(sub.iloc[k]['protein'])
                crosslinks_list.append(sub.iloc[k]['crosslinks'])
                subcellular_location_list.append(sub.iloc[k]['subcellular_location'])
                topology_list.append(topology_info.iloc[res[len(res)-1]]['topo'])
                transmembrane_list.append(sub.iloc[k]['transmembrane'])
                is_lm_list.append(topology_info.iloc[res[len(res) - 1]]['is_lm'])

                # TODO what if its not equal

            elif k == len(sub.index)-1:
                # if its last only take tm region before and add another row
                transmembrane_end = int((sub.iloc[k - 1]['transmembrane']).split('..')[0])
                #res = list(filter(lambda i: i > transmembrane_end, list(topology_info["start"])))[0]
                res = topology_info.index[topology_info["start"] > transmembrane_end].tolist()
                if not res:
                    gene_list.append(sub.iloc[k]['gene'])
                    protein_list.append(sub.iloc[k]['protein'])
                    crosslinks_list.append(sub.iloc[k]['crosslinks'])
                    subcellular_location_list.append(sub.iloc[k]['subcellular_location'])
                    topology_list.append("")
                    transmembrane_list.append(sub.iloc[k]['transmembrane'])
                    is_lm_list.append(sub.iloc[k]['is_localization_marker'])
                else:
                    gene_list.append(sub.iloc[k]['gene'])
                    protein_list.append(sub.iloc[k]['protein'])
                    crosslinks_list.append(sub.iloc[k]['crosslinks'])
                    subcellular_location_list.append(sub.iloc[k]['subcellular_location'])
                    topology_list.append(topology_info.iloc[res[len(res)-1]]['topo'])
                    transmembrane_list.append(sub.iloc[k]['transmembrane'])
                    is_lm_list.append(topology_info.iloc[res[len(res) - 1]]['is_lm'])

    new_data = pd.DataFrame({'gene': gene_list, 'protein': protein_list, 'crosslinks': crosslinks_list,
                             'topology': topology_list, 'subcellular_location': subcellular_location_list,
                             'transmembrane': transmembrane_list,'is_localization_marker': is_lm_list})
    return new_data


if __name__ == '__main__':
    # adjust your data paths
    data = pd.read_csv(
        'lysine_xlinks_xlilo.csv',sep=';',dtype='str')
    uniprot = pd.read_csv('../protein_location_prediction_local/uniprotkb_AND_reviewed_true_AND_model_o_2023_07_18.tsv',
                          sep='\t', header=0)

    # gets subcellular location information
    lm_data = get_localization_marker_information(data, uniprot)
    lm_data = lm_data.drop_duplicates(subset=['gene'], keep='first')
    print('lm done')

    transmem_data = get_transmembrane_information(data, uniprot)
    transmem_data = transmem_data.drop_duplicates(subset=['gene'], keep='first')
    print('transmem done')

    #updated_data = change_crosslinks_for_MTproteins(data)

    # combine and extend transmembrane proteins to multiple row
    combined_data = combine_lm_transmem_and_proteins(data, lm_data, transmem_data)
    print('combined data done')

    # a protein can either have a unique subcellular location, or if it has multiple topological domains
    # it must have a unique location for this domain
    combined_data_with_topology = add_topology_information(combined_data,uniprot)

    combined_data.to_csv('combined_data.csv', index=False)
    combined_data_with_topology.to_csv('combined_data_with_topology.csv',index=False)
    print('csv saved')
