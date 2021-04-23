# BPPart/BPMax
This repository includes the implementation of BPPart and BPMax, as well as the last stable version of piRNA tool from the Algorithmic Biology Laboratory (with the permission of Dr. Hamidreza Chitsaz)

To see the details about these tools and what they compute, please refer to below manuscripts. Please cite them if you use part of this repository:

1. Ebrahimpour-Boroojeny, A., Rajopadhye, S., & Chitsaz, H. (2019). BPPart and BPMax: RNA-RNA interaction partition function and structure prediction for the base pair counting model. arXiv preprint arXiv:1904.01235.
2. Chitsaz, H., Salari, R., Sahinalp, S. C., & Backofen, R. (2009). A partition function algorithm for interacting nucleic acid strands. Bioinformatics, 25(12), i365-i373.

# Installation
To install the tools after cloning or downloading the repository, simply use the command:

```
  make
```

# Input format
The input format to these tools have beem made the same to make it easier to replace them with one another in a pipeline. The input has to be a fasta file in which each entry contains the the pair of sequence we want to run these tools on. Those sequence should be separated with an `&`. `test_input.fa` is a sample input with two pair of sequences. To make the format of the input clear, here we show the content of the file:

```
  >example1
  ACCGCCGTCTTCGAGGAAAG&CCCGGCTGCTAGCTAGGAGAAATCGCGCATTT
  >example2
  CGCGCTGGATAAATATAGGACCAGGAAT&GCTCGGATAGAGCTAGGAGAAATCGCGCCGCTAGA
```

# Running the tools
To run the tools with the default hyper-parameters use these commands (replace the `test_input.fa` with your desired input file):

```
  ./bppart test_input.fa
```
```
  ./bpmax test_input.fa
```
```
  ./pirna test_input.fa
```

# Tunning the hyper-parameters

To get a list of the hyper-parameters of the each of the tools, simple use `-h`. As an example for bppart:

```
  ./bppart -h
```

As an example, to change the weights of `AU` and `GU` pairs to `1.0` and `2.5`, respectively, use this command (keep in mind that the weight of `CG` pairs are considered to be `3`):
```
  ./bppart -A=1 -G=2.5 test_input.fa
```

As another example, to run piRNA in 25 degrees Celsius, run this:
```
  ./pirna -t=25 -T=25 test_input.fa
```

# Accumulating the results
The results of bppart and bpmax will be automatically accumulated in a tab-separated file, in which ear row contains the results for one entry (one pair of sequences) of the input file.

piRNA generates separate files for each of the entries in the input file. To make it easier to process the data, we have prepared a script to accumulate these results in single tab-separated file. To do so, use this command (replace the `test_input.fa` with your desired input file):
```
   python src_ext_res.py test_input.fa
```

# Precomputed scores
The precomputed scores on the data from RISE database are available in the `pre_computed` folder. `table_x.csv` files have the scores of piRNA at temperature `x`. BPPart and BPMax scores are also available (they are not temperature-dependent). 

To compute the correlations that are presented in the paper and generate the correlation plot of the paper, you can run the command below. Note that `len_data_human_1_101.fa` is the file that has the length information of the RNA pairs in order to normalize the scores.

```
   python correlation.py pre_computed/len_data_human_1_101.fa
```


