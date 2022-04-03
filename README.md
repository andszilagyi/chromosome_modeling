Repository for software described in
# Evolution of linkage and genome expansion of protocells: the origin of chromosomes

András Szilágyi<sup>1,2,3</sup>; Viktor Péter Kovács<sup>1</sup>; Eörs Szathmáry<sup>1,2,3</sup>; Mauro Santos<sup>1,4</sup>

<sup>1</sup>Institute of Evolution, Centre for Ecological Research, 8237 Tihany, Hungary

<sup>2</sup>Department of Plant Systematics, Ecology and Theoretical Biology, Eötvös Loránd University, 1117 Budapest, Hungary

<sup>3</sup>Center for the Conceptual Foundations of Science, Parmenides Foundation, 82049 Pullach/Munich, Germany

<sup>4</sup>Grup de Genòmica, Bioinformàtica i Biologia Evolutiva (GGBE), Departament de Genètica i de Microbiologia, Universitat Autonòma de Barcelona, Bellaterra, Barcelona, Spain

Szilágyi A, Kovács VP, Szathmáry E, Santos M (2020) Evolution of linkage and genome expansion in protocells: The origin of chromosomes. *PLOS Genetics* **16**(10): e1009155. https://doi.org/10.1371/journal.pgen.1009155

Note: The code available from this repository is provided without any warranty or liability.



### Usage:

There are two compile time `#define`-able code options `FEATURE_NO_MUT` for no mutation and `FEATURE_NO_ASSORT` for a decreased assortment load, giving a total of 4 valid combinations. It's possible to choose between various predefined replication activities stored in the `replaTbl` table.

To build the program under linux or any compatible system using the GCC toolchain, issue the command `make` in the source directory.

The command line has the following parameters in this fixed order: `D SD mu randomSeed`
(For the meaning of these parameters, please see the model description in the article.)

If the command line is empty, the default values are the following:
```
    D = 3
    SD = 30
    mu = 0.001
    randomSeed = 12249
```
