# Prediction of protein location by cross-linking mass spectrometry data
The scripts in this repository utilize cross-linking mass spectrometry (XL-MS) data to predict the location and topology of proteins by direct interactions and defined localization marker.

### Summary
- script execution done in IDE, e.g. PyCharm or Spyder
- prediction is split into two steps
  1. data preparation script: gets localization and topology information from uniprot for cross-linked proteins in XL-MS data and assigns localization marker; resulting script can be edited manually
  2. prediction script: predicts the localization and topology of proteins by pre-defined localization marker

### Data preparation script
#### Input
- XL-MS information, the following columns are needed:
  | gene_a | gene_b | Protein1 | Protein2 | n_score_a | n_score_b | LinkPos1 | LinkPos2 | crosslink_score | crosslink_ab | crosslinks_a | crosslinks_b |
  | :--- | :---: | :---: | :---: | :---: |  :---: | :---: | :---: | :---: | ---: |
- crosslinks_a and crosslinks_b columns must contain all crosslinks found in XL-MS data for Protein1 and Protein2
- reviewed uniprot information downloaded as tsv with the following columns:
  | Entry | Protein names | Gene Names | Subcellular location [CC] | Topological domain | Transmembrane |
  | :--- | :---: | ---: | :--- | :---: | ---: |

#### Method description
- define localization marker by checking which proteins in XL-MS data have either only one subcellular location or, in case of membrane proteins, a defined topology and localization for non-transmembrane regions in uniprot data
- accumulate transmembrane regions for membrane proteins
- combine protein information in XL-MS data, subcellular locations and membrane topology
  - for each protein_a, crosslinks will be accumulated by taking all (inter-protein) interactions where it is either protein_a or protein_b (this will be done for protein_b in data accordingly)
  - if protein is not first part of crosslink (PB-1-PROTEIN-5), the crosslink will be reversed (PROTEIN-5-PB-1)
  - if this protein has a transmembrane region, it will get multiple rows splitting crosslinks and
    locations/topology in each row according to residue numbers of the transmembrane region 
      | Protein | Crosslink | Transmembrane region |
      | :--- | :---: | ---: |
      | Protein1 | P1-3-P5-6 | |
      | Protein1 |           | 5..9 |
      | Protein1 | P1-11-P8-5 | |

- topology information for proteins with transmembrane regions will be added from uniprot data and ordered according to their transmembrane regions
- a localization marker column is added, indicating proteins with an already known and thus starting location/topology as True, unknown proteins or transmembrane regions as False
- these localization marker cancidates can be checked for their consistency and crosslinks between them
- the output is a combined table, which will be saved and can be edited manually
  |gene | protein | crosslinks | topology | subcellular_location | transmembrane | is_localization_marker |
  | :--- | :---: | :---: | :--- | :---: | :---: | ---: |
  
### Prediction script
#### Input
- combined table (output of the data preparation script)

#### Method description
- if transmembrane regions were changed, crosslink, subcellular location and topology information will be reorderd and updated for each protein
- the prediction is done by using the crosslinked residue of each row of a protein to identify the residue of a linked protein, and assigning topology as well as subcellular location 
- afterwards, all predicted residues of a protein will be aggregated and ordered via their residue numbers
- the following output table will be saved
  | predicted_gene | predicting_gene | predicted_gene_residue | predicting_gene_residue | predicted_subcellular_location | predicted_topology | predicting_crosslinks | predicted_by_transmembrane | transmembrane_regions | predicting_gene_is_lm |
  | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | ---: |


### Interpretation of results
In the output table, the "predicted_subcellular_location" and "predicted_topology" columns show the sub-compartment localization or the topology information of the protein.
If there is no transmembrane information for a protein in the "Transmembrane region" column, hence a soluble protein, the predicted information of "predicted_subcellular_location" and "predicted_topology" should be consistent. Otherwise, the sub-compartment localization of this protein is ambiguous.
If there is transmembrane information in the "Transmembrane region" column, hence a membrane protein, the predicted information needs to be interpreted according to the transmembrane region, and for ease of interpretation, the predicted information for each protein is ordered via their residue numbers.

### Example
- the example folder contains exemplary input and output tables 
