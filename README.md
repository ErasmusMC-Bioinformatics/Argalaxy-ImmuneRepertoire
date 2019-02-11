# ARGalaxy Immune Repertoire
This is the GitHub repository for the ARGalaxy Immune repertoire pipeline.  
The Galaxy tool version can be found [here](https://toolshed.g2.bx.psu.edu/repository/browse_repositories_i_own?sort=name&operation=view_or_manage_repository&id=2e457d63170a4b1c).  
The docker version can be found [here](https://github.com/ErasmusMC-Bioinformatics/ARGalaxy-docker).

## Overview

In execution order:

#### imgt_loader or igblast

###### imgt_loader (Recommended)
Start the analysis with [IMGT HighV Quest](https://www.imgt.org/HighV-QUEST/) archives.  
An IMGT archive file holds [multiple tabular files](http://www.imgt.org/IMGT_vquest/share/textes/imgtvquest.html#output3), this script extracts the specific columns relevant to the analysis from several of these files.

`Rscript imgt_loader.r 1_Summary.txt 3_Nt-sequences.txt 5_AA-sequences.txt 6_Junction.txt 4_IMGT-gapped-AA-sequences.txt /path/to/output.txt`


###### igblast
Start the analysis with FASTA files that are aligned with [igblast](https://www.ncbi.nlm.nih.gov/igblast/).  
Note that this method will provide less information than the IMGT archive.

`sh igblast.sh /path/to/sequences.fasta species locus /path/to/output.txt`

#### experimental_design
This script will merge multiple result files (from the last step) into a single file with an additional ID and Replicate column to differentiate the individual samples during the analysis and to allow for analysis across samples.

`Rscript experimental_design.r /path/to/input_1 id_1 [/path/to/input_2 id_2] [/path/to/input_n id_n] /path/to/output`

#### report_clonality
The R script that creates the analysis result.

`sh r_wrapper.sh /path/to/experimental_design/output.txt /path/to/output_dir/output.html /path/to/output_dir "clonaltype" "species" "locus" "filter_productive" "clonality_method"`  

###### parameters
Clonaltype:  
- none
- Top.V.Gene,CDR3.Seq
- Top.V.Gene,CDR3.Seq.DNA
- Top.V.Gene,Top.J.Gene,CDR3.Seq
- Top.V.Gene,Top.J.Gene,CDR3.Seq.DNA
- Top.V.Gene,Top.D.Gene,Top.J.Gene,CDR3.Seq.DNA

Species:
- Homo sapiens functional
- Homo sapiens
- Homo sapiens non-functional
- Bos taurus
- Bos taurus functional
- Bos taurus non-functional
- Camelus dromedarius
- Camelus dromedarius functional
- Camelus dromedarius non-functional
- Canis lupus familiaris
- Canis lupus familiaris functional
- Canis lupus familiaris non-functional
- Danio rerio
- Danio rerio functional
- Danio rerio non-functional
- Macaca mulatta
- Macaca mulatta functional
- Macaca mulatta non-functional
- Mus musculus
- Mus musculus functional
- Mus musculus non-functional
- Mus spretus
- Mus spretus functional
- Mus spretus non-functional
- Oncorhynchus mykiss
- Oncorhynchus mykiss functional
- Oncorhynchus mykiss non-functional
- Ornithorhynchus anatinus
- Ornithorhynchus anatinus functional
- Ornithorhynchus anatinus non-functional
- Oryctolagus cuniculus
- Oryctolagus cuniculus functional
- Oryctolagus cuniculus non-functional
- Rattus norvegicus
- Rattus norvegicus functional
- Rattus norvegicus non-functional
- Sus scrofa
- Sus scrofa functional
- Sus scrofa non-functional

Locus:
- TRA
- TRD
- TRG
- TRB
- IGH
- IGI
- IGK
- IGL

Filter productive:
- yes
- no

Clonality Method:
- none
- old
- boyd

## complete.sh
This script will run all of the above for you, it will detect if you are using FASTA files or IMGT archives and use the appropriate tools.

`sh complete.sh /path/to/input_1 id_1 [/path/to/input_n id_n] /path/to/out_dir/out.html /path/to/out_dir clonaltype species locus filter_productive clonality_method`  
See "report_clonality" for the parameter options.

## Dependencies
- Linux
- R
  - gridExtra
  - ggplot2
  - plyr
  - data.table
  - reshape2
  - lymphclon

#### optional
- Circos
- IgBlast
- igblastwrp
