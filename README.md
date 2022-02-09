

## OMAMO: orthology-based model organism selection


![workflow diagram](logo-omamo.jpg)



OMAMO is a tool that suggests the best model organism to study a biological process based on orthologous relationship between a species and human. 

The user can consider several species as potential model organisms and the algorithm will rank them and report the output for a given biological process (searched as a GO term or a GO ID) is produced in the dataframe format.


### Dependencies
Following Python packages are needed: numpy, matplotlib, pickle and pandas. Besides, you need to install [pyOMA](https://pypi.org/project/pyoma).


## Pipeline

Firstly, download the OMA dataset:

```
wget  https://omabrowser.org/All/OmaServer.h5  -O data/OmaServer.h5  #caution: 94GB
```

Secondly, using the file `data/oma-species.txt` find the five-letter UniProt code for species of interest. For example, consider three species _Dicdyostelium discodeium_ , _Neurospora crassa_ and _Schizosaccharomyces pombe_. Their UniProt codes are `DICDI`, `NEUCR` and `SCHPO`, respectively.

Install omamo from the git checkout:

```bash
pip install <path_to_omamo.git>
```

Once the package is installed, you should be able to run `omamo` as a command. With `omamo -h` see the available options:
```text
usage: omamo [-h] --db DB [--query QUERY] [--ic IC] [--h5-out H5_OUT] [--tsv-out TSV_OUT] --models MODELS [MODELS ...]

Run omamo for a set of model organisms

optional arguments:
  -h, --help            show this help message and exit
  --db DB               Path to the HDF5 database
  --query QUERY         Name of the Query species, defaults to HUMAN
  --ic IC               Path to the information content file (tsv format)
  --h5-out H5_OUT       Path to the HDF5 output file. If omitted, not stored in this format
  --tsv-out TSV_OUT     Path to the TSV output file. If omitted, not stored in this format
  --models MODELS [MODELS ...]
                        List of model species, or a path to a txt file with the model species
```

In order to create the omamo data for _Dicdyostelium discodeium_, _Neurospora crassa_ and _Schizosaccharomyces pombe_, 
we would run omamo with the following parameters:

```
omamo --db OmaServer.h5 --query HUMAN --tsv-out omamo_output_df.csv
```
Finally, the output data frame is ready as a TSV file `omamo_output_df.csv`. For example, for the GO ID of `GO0000012`, OMAMO provides the following ranking for potential model organisms: 


```
head -n 1 omamo_output_df.csv > ranked_organisms.csv
awk '$2 == 12'  omamo_output_df.csv >> ranked_organisms.csv
cat ranked_organisms.csv


row    GO_ID   Species Human_Genes     Species_Genes   No.of_OGs  Average_func.similarity±st.dev      Score

44	12	SCHPO	A0A024R6L5,A8K3W1	Q9USG9,Q9USR0	2	0.8179 ± 0.0986	1.64
55	12	NEUCR	A0A024RE06	DNLI4_NEUCR	1	0.292 ± 0.0	0.29
56	12	DICDI	A0A024RE06	Q54CR9	1	0.2848 ± 0.0	0.28

```



## OMAMO Website

You can also visit the [OMAMO website](https://omabrowser.org/omamo), where you can browse biological processes to study in 50 unicellular species.





## Change log

Version 0.0.1
- Initial release


## Citation

Alina Nicheperovich, Adrian M Altenhoff, Christophe Dessimoz, Sina Majidian, "OMAMO: orthology-based model organism selection", submitted to Bioinformatics journal, [preprint](https://www.biorxiv.org/content/10.1101/2021.11.04.467067v1).



## License

OMAMO is a free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

OMAMO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with OMAMO. If not, see http://www.gnu.org/licenses/.




