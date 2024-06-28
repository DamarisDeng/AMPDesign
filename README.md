# AMPDesign

Antimicrobial resistance poses a threat to human well-being, especially with the emergence of multidrug-resistant bacteria. Antimicrobial peptides (AMPs), as one of the most promising alternatives to traditional antibiotics, have attracted significant attention in recent years. However, designing AMPs using traditional experimental methods is time-consuming and labor-intensive. Therefore, in this study, we propose an out-of-the-box AMP design toolbox, named AMPDesign. 

AMPDesign is a novel computational tool designed for the generation of antimicrobial peptide (AMP) sequences. It integrates the Generative Adversarial Network (GAN) architecture and Reinforcement Learning (RL), incorporates our previously developed activity predictor [AMPActiPred](https://onlinelibrary.wiley.com/doi/10.1002/pro.5006) model to produce peptides with high antimicrobial efficacy. Notably, AMPDesign is the first AMP generation toolbox capable of designing AMPs against specific bacterial species.

The platform's versatility is demonstrated through its capacity to generate AMPs targeting 10 distinct bacterial species. Researchers can access the comprehensive training datasets and the  prediction model via the [AMPActiPred website](https://awi.cuhk.edu.cn/~AMPActiPred/download.php), facilitating further exploration and application in the field of antimicrobial peptide design.

![img.png](img.png)

## Requirement

**查询requirement**

We suggest you run the platform under Python 3.7+ with following libs: 

- TensorFlow 
- Numpy
- Scipy
- NLTK
- CUDA 7.5+ (Suggested for GPU speed up, not compulsory)
- deep-forest
- Biopython
- 还有什么？

Or just type `pip install -r requirements.txt` in your terminal.

## Quick tutorial

```bash
git clone https://github.com/DamarisDeng/AMPDesign.git
cd AMPDesign
```

The basic usage is

```bash
python AMPDesign.py -s <microbial-type> -a model/Ec/ -d train.fasta -o output.fasta
```

Arguments:

- `-s`: This is used to specify the type of microbial your peptides belong, you can choose from `Ec`, and `Sa`.
- `-a`: Specifiy the input route.
- `-o`: Specify the output filename
- `-d`: the training sequence