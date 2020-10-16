Repository for software described in
# Evolution of linkage and genome expansion of protocells: the origin of chromosomes

András Szilágyi (1,2,3); Viktor Péter Kovács (1); Eörs Szathmáry (1,2,3); Mauro Santos (1,4)

1. Institute of Evolution, Centre for Ecological Research, 8237 Tihany, Hungary

2. Department of Plant Systematics, Ecology and Theoretical Biology, Eötvös Loránd University, 1117 Budapest, Hungary

3. Center for the Conceptual Foundations of Science, Parmenides Foundation, 82049 Pullach/Munich, Germany

4. Grup de Genòmica, Bioinformàtica i Biologia Evolutiva (GGBE), Departament de Genètica i de Microbiologia, Universitat Autonòma de Barcelona, Bellaterra, Barcelona, Spain


PLOS Genetics (citation data and URL will be inserted here)

Note: The code available from this repository is provided without any warranty or liability.



### Usage:

There are two compile time `#define`-able code options `FEATURE_NO_MUT` for no mutation and `FEATURE_NO_ASSORT` for a decreased assortment load, giving a total of 4 valid combinations. It's possible to choose between various predefined replication activities stored in the `replaTbl` table.

The command line has the following parameters in this fixed order: `D SD mu randomSeed`

If the command line is empty, the default values are the following:
```
    D = 3
    SD = 30
    mu = 0.001
    randomSeed = 12249
```
