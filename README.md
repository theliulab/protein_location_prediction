# Prediction of protein locations by XL-MS data
## Python script - One Epoch

### Data preparation script

- takes xlinkx_XLinkCyNet and downloaded uniprot data as input
- get localization markers by checking which proteins in xlink data have
  only one location in uniprot data
- gets transmembrane region for proteins in xlink data
- if proteins include additional seperators (eg '-') in names or crosslinks,
  the seperator will be removed for a decreased difficulty in parsing
- combining of proteins in xlink data, localizations for localization marker and
  transmembrane information to one table
  -  for each protein_a and protein_b in data, crosslinks will be gathered by taking links
    in the specific row as well as all crosslinks its part of the opposite column (eg specific row
    is protein_a, also take links of column protein_b)
  - if protein is not first part of crosslink (PB-1-PROTEIN-5), the crosslink will be inversed
    (PROTEIN-5-PB-1)
  - all inter links for this protein will be accumulated
  - if this protein has transmembrane regions, it will get multiple rows splitting crosslinks and
    locations in each row according to residue numbers of the transmembrane region 
      | Protein | Crosslink | Transmembrane region |
      | :--- | :---: | ---: |
      | Protein1 | P1-3-P5-6 | |
      | Protein1 |           | 5..9 |
      | Protein1 | P1-11-P8-5 | |‚
    
- saving of combined table for user modification

  
### Prediction script
- reading combined table previously prepared and with possible user modification
- if transmembrane regions were changed, crosslinks will be reorderd for each proteins
- prediction will be done by taking crosslinks for each row of a protein and adding the linked protein and their
  residues to a new data frame
- later on the prediction result will be combined by taking each predicted gene once, the genes it is predicted by
  as well as according residue numbers and crosslinks
