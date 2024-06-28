# AMPDesign

Antimicrobial resistance poses a significant threat to global health. AMPs offer a promising alternative to traditional antibiotics, but their design through experimental methods is time-consuming and resource-intensive. Here we proposed AMPDesign, a tool that addresses this challenge by providing a computational approach to AMP design, significantly accelerating the discovery process.

## Overview
![AMPDesign Overview](img.png)

AMPDesign is an innovative computational tool for generating antimicrobial peptide (AMP) sequences. It combines Generative Adversarial Network (GAN) architecture with Reinforcement Learning (RL) and incorporates our previously developed [AMPActiPred](https://onlinelibrary.wiley.com/doi/10.1002/pro.5006) model to produce highly effective antimicrobial peptides.

### Key Features

- First AMP generation toolbox capable of designing AMPs against specific bacterial species
- Generates AMPs targeting 10 distinct bacterial species
- Integrates GAN architecture and Reinforcement Learning
- Incorporates AMPActiPred model for high antimicrobial efficacy

## Quick tutorial

```bash
git clone https://github.com/DamarisDeng/AMPDesign.git
cd AMPDesign
```
### Requirements

We suggest you run the platform under Python 3.7+ with following libs: 

- TensorFlow == 1.13.1
- Numpy >= 1.21
- Scipy >= 1.7.3
- NLTK == 3.7
- deep-forest == 0.1.7
- biopython == 1.81

Or just type `pip install -r requirements.txt` in your terminal.

The basic usage is

```bash
python AMPDesign.py -s <microbial-type> -a <model-loc> -d <training-data> -o <output-data>

# example
python AMPDesign.py -s Ec -a prediction/models/Ec/ -d data/Ec.txt -o save/output.txt
```

### Arguments:

- `-s`: This is used to specify the type of microbial your peptides belong, you can choose from 10 different species, including E. coli (Ec), M. luteus (Ml), K. pneumoniae (Kp), P. aeruginosa (Pa), A. baumannii (Ab), E. feacalis (Ef), S. typhimurium (St), S. aureus (Sa), S. epidermidis (Se), B. subtilis (Bs).
- `-a`: Path to the prediction model. Models can be downloaded from [AMPActiPred website](https://awi.cuhk.edu.cn/~AMPActiPred/download.php).
- `-o`: Output file path.
- `-d`: Training sequence file path.
  - The training data should be a .txt file, with each line containing a peptide sequence.

```
# sample training data
AAGMGFFGAR
AAHCLAIGRR
KKAFAAAAAFAAWAAFAKKKK
AIHKLAHKLTKKTLRAVKKLAN
GWGDTFGKVLKNFAKVAGVKAAK
```

Hope you enjoy using AMPDesign!
