## Introduction

pyDIscore is an implementation version of the published Directionality Index score calculation script by Dixon *et al.* [[1]](#1).
Unlike the original written in Perl (<http://renlab.sdsc.edu/yanxiao/download/hic-pip/domaincall_software/perl_scripts/DI_from_matrix.pl>), 
it's written in Python and works compatible with the covNorm (<https://github.com/kaistcbfg/covNormRpkg>) output and numpy array.


## Installlation

pyDIscore works with basic Python packages. Most of them are installed by default when installing Python (such as anaconda).
The code was tested on Python 3.8.10 and 2.7.15.

* numpy
* argparse
* pickle
* gzip
* sys

```bash
git clone https://github.com/kaistcbfg/pyDIscore.git
```

## Usage

```bash
$ python pyDIscore -h
usage: pyDIscore.py [-h] --input-file INPUT_FILE [--input-format INPUT_FORMAT]
                    --chrname CHRNAME --fai-file FAI_FILE
                    [--resolution RESOLUTION] [--window-size WINDOW_SIZE]
                    [--double-count-flag DOUBLE_COUNT_FLAG]
                    [--chrname-number-flag CHRNAME_NUMBER_FLAG]
                    [--fullbin-output-flag FULLBIN_OUTPUT_FLAG]
                    [--output-file OUTPUT_FILE]

python DI score calc.

optional arguments:
  -h, --help            show this help message and exit
  --input-file INPUT_FILE
                        input *.gz file (required)
    --input-format INPUT_FORMAT
                        default cov, format: cov (covNorm) or pkl (numpy
                        pickled array)
  --chrname CHRNAME     target chromosome (required)
  --fai-file FAI_FILE   FAI file for chr size (required)
  --resolution RESOLUTION
                        bin resolution (default 40kb)
  --window-size WINDOW_SIZE
                        window size (default 2Mb)
  --double-count-flag DOUBLE_COUNT_FLAG
                        default False, if True, apply /2 to matrix
  --chrname-number-flag CHRNAME_NUMBER_FLAG
                        default False, if True chr1 -> 1 (X:23, Y:24, M:25
  --fullbin-output-flag FULLBIN_OUTPUT_FLAG
                        default True, if False, print from starter bin
  --output-file OUTPUT_FILE
                        output DI score bedgraph file. if None, print
```

Input format preparation:
* covNorm (cov)
```R
write.table(final_df, file=gzfile("covnorm.gz"), row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
```

* numpy pickled array (pkl): Check [Output](https://github.com/kaistcbfg/covNormRpkg#Output) section of covNorm github for 'contact_map' generation.


```python
import pickle
import gzip

pickle.dump(contact_map, gzip.open('np_pkl_arr.gz', 'wb'))
```


## Citation

<a id="1">[1]</a>  Dixon JR, Selvaraj S, Yue F, Kim A, Li Y, Shen Y, Hu M, Liu JS, Ren B. Topological domains in mammalian genomes identified by analysis of chromatin interactions. *Nature*, **485**, 376–380 (2012). <https://doi.org/10.1038/nature11082>
