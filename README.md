

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
omamo --db OmaServer.h5 --query HUMAN --tsv-out omamo_output_df.csv --models  DICDI NEUCR SCHPO
```

You might face an error about `OSError: ``OmaServer.h5.idx`` does not exist` and `pyoma.browser.db.DBConsistencyError: Suffix index for protein sequences is not available` which you can ignore them. 

Finally, the output data frame is ready as a TSV file `omamo_output_df.csv`. For example, for the GO ID of `GO0000472`, "endonucleolytic cleavage to generate mature 5'-end of SSU-rRNA", OMAMO provides the following ranking for potential model organisms: 


```
head -n 1 omamo_output_df.csv > ranked_organisms.csv
awk '$1 == 472'  omamo_output_df.csv >> ranked_organisms.csv
cat ranked_organisms.csv


GOnr	Species	QuerySpeciesGenes	ModelSpeciesGenes NrOrthologs	FuncSim_Mean	FuncSim_Std	Score
472	DICDI	NOP9;TBL3;ABT1	  Q551Y5;Q7KWS8;esf2	          3  	0.9095	0.1567	2.7286
472	NEUCR	NOP9;TBL3	         nop9;pod-5	          2  	1.0000	0.0000	2.0000
472	SCHPO	NOP9;TBL3	         nop9;utp13	          2  	1.0000	0.0000	2.0000
```



## OMAMO Website

You can also visit the [OMAMO website](https://omabrowser.org/omamo), where you can browse biological processes to study in 50 unicellular species.





## Change log

#### Version 0.2.1
- store ic values in hdf5 database

#### Version 0.2.0
- Overhaul and creating pip package

#### Version 0.0.1
- Initial release


## Citation

Alina Nicheperovich, Adrian M Altenhoff, Christophe Dessimoz, Sina Majidian, "OMAMO: orthology-based model organism selection", submitted to Bioinformatics journal, [preprint](https://www.biorxiv.org/content/10.1101/2021.11.04.467067v1).



## License

OMAMO is a free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

OMAMO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with OMAMO. If not, see http://www.gnu.org/licenses/.




