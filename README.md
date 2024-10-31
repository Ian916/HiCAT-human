[![DOI](https://zenodo.org/badge/748088016.svg)](https://zenodo.org/doi/10.5281/zenodo.10570850)
# HiCAT-human

We proposed a modified version of our previous HOR annotation tool [HiCAT](https://github.com/xjtu-omics/HiCAT) for automatically annotating centromere HOR patterns from both HiFi reads and assemblies of multiple human samples.

## Dependencies

Python 3.9.13

Development environment: Linux

Development tool: Pycharm

| Packages           | Version |
| ------------------ | ------- |
| biopython          | 1.79    |
| joblib             | 1.1.0   |
| lastz              | 1.04.22 |
| matplotlib         | 3.5.1   |
| numpy              | 1.22.3  |
| pandas             | 1.4.0   |
| python-edlib       | 1.3.9   |
| python-levenshtein | 0.12.2  |
| scikit-learn       | 1.0.2   |
| seqtk              | 1.2     |
| setuptools         | 61.2.0  |

[StringDecomposer](https://github.com/ablab/stringdecomposer)  version 1.1.2. (included in HiCAT-human source code)

## Quick start

`HiCAT-human` is a tool to automatically annotate centromere HOR patterns from both reads and assemblies of multiple human samples.

#### Installation

```
#install
git clone https://github.com/xjtu-omics/HiCAT-human.git
conda install -y --file requirements.txt
cd ./stringdecomposer && make
```

#### Overview

HiCAT-human consists of 4 modules:

- `reads` used for reads HOR annotation.
- `reads_aggregate` used for aggregating reads annotation results.
- `assembly` used for assembly HOR annotation.
- `assembly_match` used for matching assembly annotation results to reads annotation results.

#### Input data

* Reads HOR annotation: whole genome HiFi reads.
* Assembly HOR annotation: haplotype-resolved human genome.

For detail usage, read the docs on the `HiCAT-human` [wiki](https://github.com/Ian916/HiCAT-human/wiki).

## Contact

If you have any questions, please feel free to contact: [gaoxian15002970749@163.com](mailto:gaoxian15002970749@163.com), [xfyang@xjtu.edu.cn](mailto:xfyang@xjtu.edu.cn), [kaiye@xjtu.edu.cn](mailto:kaiye@xjtu.edu.cn)

## Reference

Please cite the following paper when you use HiCAT-human in your work

Gao S, Zhang Y, Bush SJ, Wang B, Yang X, Ye K. Centromere landscapes resolved from hundreds of human genomes. Genomics, Proteomics & Bioinformatics. 2024 Oct 18:qzae071. https://doi.org/10.1093/gpbjnl/qzae071