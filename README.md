# filtFst

The C program **filtFst** calculates *F<sub>ST</sub>*, a measure of population differentiation due to genetic structure, on datasets containing genetic data for two populations formatted in a similar style to the output of Richard Hudson's simulation software, [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html).

**filtFst** was originally written by [Jeff Wall](http://profiles.ucsf.edu/jeff.wall) (University of California San Francisco), and later modified by August Woerner (University of Arizona) and [Murray Cox](http://massey.genomicus.com) (Massey University).

**filtFst** implements the variant of *F<sub>ST</sub>* described in:

Hudson, R.R., M. Slatkin and W.P. Maddison. 1992. [Estimation of levels of gene flow from DNA sequence data.](http://www.genetics.org/content/132/2/583.abstract) *Genetics* **132**:583-589.

Specifically, equation 3:

&lang;*F<sub>ST</sub>*&rang; = 1 &minus; (*H<sub>w</sub>* / *H<sub>b</sub>*)

where *H<sub>w</sub>* is the mean pairwise diversity within populations and *H<sub>b</sub>* is the mean pairwise diversity between populations. The implementation here includes a correction factor for unequal sample sizes, as described in the supplementary material of:

Plagnol, V. and J.D. Wall. 2006. [Possible ancestral structure in human populations.](http://doi.org/10.1371/journal.pgen.0020105) *PLoS Genetics* **2**:e105.

**filtFst** can accommodate:

+ Simulated data in a slight modification of Richard Hudson's [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html) format; and
+ Real data in the same format (with missing data encoded as '?').

#### COMPILATION

To compile **filtFst**, install any standard C compiler. For [gcc](https://gcc.gnu.org), type:

```gcc filtFst.c -o filtFst -O3 -lm -Wall```

#### USAGE

The basic usage of the program is:

```filtFst n1 n2```

where

```n1``` is the sample size of population 1 (i.e., the number of sequences)  
```n2``` is the sample size of population 2 (i.e., the number of sequences)


The data file itself is read from stdin (see **EXAMPLES** section below).

#### INPUT FORMAT

The input format is similar to the output format of Richard Hudson's [ms](http://home.uchicago.edu/rhudson1/source/mksamples.html) program, described in:

Hudson, R.R. 2002. [Generating samples under a Wright-Fisher neutral model of genetic variation.](http://doi.org/10.1093/bioinformatics/18.2.337) *Bioinformatics* **18**:337-338.

A similar format is used by Gary Chen's genome-scale Markovian coalescent simulator, [MaCS](https://github.com/gchen98/macs), described in:

Chen, G.K., P. Marjoram and J.D. Wall. 2009. [Fast and flexible simulation of DNA sequence data.](http://doi.org/10.1101/gr.083634.108) *Genome Research* **19**:136-142.

The first line lists the number of segregating sites. The second line lists the positions of each single nucleotide polymorphism.  (These are not used in the *F<sub>ST</sub>* calculation, but one value per segregating site must be present in the data file). Each following line gives a single haplotype with two alleles represented by '0' and '1' and missing data represented by '?'.  Only two character states (excluding the missing data character) are allowed per site and polymorphisms must be point mutations (*i.e.*, the program does not accommodate indels). The first set of ```n1``` sequences should come from population one, while the last set of ```n2``` sequences should come from population two.


```
42
0.176460 0.176470 0.178610 0.179200 0.183260 0.184560 ...
000000000000000100000000000000000000000000
110000100010100101010100011100001010011111
??0000000000000100000000011000000110000000
??0110000000000110100000011000000110000000
000000000000000000000000011000000010000000
000110000000010100100001011000000110000000
??????????????????????????????????????????
??????????????????????????????????????????
000000000000000100000000001000000000000000
000000000010101100010000011000000110000010
.
.
.
```

The *F<sub>ST</sub>* value is returned to stdout.

#### EXAMPLES

The example files ```example_data1.dat``` and ```example_data2.dat``` represent real genomic data for the human loci 1pMB4 and 4qMB105 in French and Mandenka, respectively.

To calculate *F<sub>ST</sub>* between French and Mandenka at 1pMB4:

```
cat example_data1.dat | filtFst 32 34
$ 0.086995
```

To calculate *F<sub>ST</sub>* between French and Mandenka at 4qMB105:

```
cat example_data2.dat | filtFst 32 34
$ 0.139620
```
