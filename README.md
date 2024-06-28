# AMPDesign

AMPDesign is a tool for users to generate antimicrobial peptide (AMP) sequences. Users can specify ...

## Requirement

查询requirement

We suggest you run the platform under Python 3.7+ with following libs: **查询这个**

- TensorFlow 
- Numpy
- Scipy
- NLTK
- CUDA 7.5+ (Suggested for GPU speed up, not compulsory)
- deep-forest
- Biopython

Or just type `pip install -r requirements.txt` in your terminal.

## Quick tutorial

```bash
git clone https://github.com/DamarisDeng/AMPDesign.git
cd AMPDesign
```

The basic usage is

```bash
python3 AMPDesign.py -s <microbial-type> -a model/Ec/  -o output.fasta -d train.fasta
```

Arguments:

- `-s`: This is used to specify the type of microbial your peptides belong, you can choose from `Ec`, and `Sa`.
- `-a`: Specifiy the input route.
- `-o`: Specify the output filename
- `-d`: the training sequence