# Prediction of protein locations by XL-MS data
## Python script - One Epoch

### Script execution
- so far done in IDE

### Data preparation script

#### Input
- crosslinking (xlink) information in the following format:
  | gene_a| gene_b | Protein1 | Protein2 | n_score_a | n_score_b | LinkPos1 | LinkPos2 | crosslink_score | crosslink_ab | crosslinks_a | crosslinks_b |
- reviewed uniprot informatio downloaded as tsv with the following columns:
  | Entry | Protein names | Gene Names | Subcellular location [CC] | Topological domain | Transmembrane |

#### Method description
- get localization markers by checking which proteins in xlink data have
  only one subcellular location in uniprot data
- get transmembrane region from uniprot data for proteins in xlink data
- combine proteins in xlink data, localizations for localization marker and transmembrane information to one table
  - for each protein_a in data, crosslinks will be accumulated by taking all (inter-protein) interactions where it is either protein_a or protein_b (this will be done for protein_b in data accordingly)
  - if protein is not first part of crosslink (PB-1-PROTEIN-5), the crosslink will be reversed
    (PROTEIN-5-PB-1)
  - if this protein has transmembrane regions, it will get multiple rows splitting crosslinks and
    locations in each row according to residue numbers of the transmembrane region 
      | Protein | Crosslink | Transmembrane region |
      | :--- | :---: | ---: |
      | Protein1 | P1-3-P5-6 | |
      | Protein1 |           | 5..9 |
      | Protein1 | P1-11-P8-5 | |

- topology information for proteins with transmembrane regions will be added from uniprot data and ordered according to its transmembrane regions
- the output is a combined table, which will is saved and can later on be changed manually
  |gene | protein | crosslinks | topology | subcellular_location | transmembrane |

  
### Prediction script
#### Input
- combined table (output for the data preparation script)

#### Method description
- if transmembrane regions were changed, crosslinks, subcellular locations and topology will be reorderd for each protein
- a localization marker column is added, indicating proteins with an already known and thus starting location/topology as True, unknown proteins or transmembrane regions as False
- the prediction is done by using crosslinked residues of each row of a protein to identify the crosslinked residues of another protein, and assigning topology as well as subcellular location to it
  residues to a new data frame
- afterwards, all predicted residues of a protein will be aggregated and ordered via their residue numbers
- the output table in the following format is saved
  | predicted_gene | predicting_gene | predicted_gene_residue | predicting_gene_residue | predicted_subcellular_location | predicted_topology | predicting_crosslinks | predicted_by_transmembrane | transmembrane_regions | predicting_gene_is_lm |
