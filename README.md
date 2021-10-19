

## OMAMO: orthology-based model organism selection


![workflow diagram](logo-omamo.jpg)



OMAMO is a tool that suggests the best model organism to study a biological process based on orthologous relationship between a species and human. 

The user can consider several species as potential model organisms and the algorithm will rank them and report the output for a given biological process (searched as a GO term or a GO ID) is produced in the dataframe format.

## Pipeline

Firstly, download the OMA dataset:

```
wget  https://omabrowser.org/All/OmaServer.h5  -O data/OmaServer.h5  #caution: 94GB
```

Secondly, using the file `data/oma-species.txt` find the five-letter UniProt code for species of interest. For example, consider three species _Dicdyostelium discodeium_ , _Neurospora crassa_ and _Schizosaccharomyces pombe_. Their UniProt codes are `DICDI`, `NEUCR` and `SCHPO`, respectively. 

Then, run the code `omamo_base.py` for each species code (`DICDI`, `NEUCR` and `SCHPO`):

```
species="DICDI"
mkdir output; cd output

python3 ../omamo_base.py ../data/OmaServer.h5 ../data/go_positive_annotations.tsv ${species}
```



Once the code finished running, the outputs include `${species}2.csv` files which should be combined to create a final dataframe using the code `omamo_dataframe.py`:

```
python3 omamo_dataframe.py output
```

where `output` is the name of the directory where the user wishes to save the output. 

Finally, the output data frame is ready as a CSV file `omamo_output_df.csv`.



## OMAMO Website

You can also visit the [OMAMO website](https://omabrowser.org/omamo), where you can browse biological processes to study in 50 unicellular species.





## Change log

Version 0.0.1
- Initial release


## Citation

Alina Nicheperovich, Adrian M Altenhoff, Christophe Dessimoz, Sina Majidian, "OMAMO: orthology-based model organism selection", submitted to Bioinformatics journal.



## License

OMAMO is a free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

OMAMO is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with OMAMO. If not, see http://www.gnu.org/licenses/.




