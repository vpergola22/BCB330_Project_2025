# 🧬 Predicting T-Cell Receptors from Tumor Mutation Data

### 📘 Overview
This project explores the feasibility of predicting **T-cell receptor (TCR)** sequences from **tumour mutation data**, aiming to understand how somatic mutations influence immune recognition and neoantigen presentation.  
Completed as part of **BCB330: Computational Biology** at the University of Toronto.

Using publicly available datasets, I built a pipeline to preprocess and integrate **TCR** and **tumour mutation data**, apply logistic regression models, and identify relationships between recurrent mutations and public TCRs.

---

### 🎯 Objectives
- Curate and preprocess tumour mutation and TCR sequencing data  
- Encode features for model training (one-hot representation)  
- Apply logistic regression to predict public TCR occurrence from mutation bins  
- Evaluate model performance with accuracy and F1 scores  
- Visualize mutation–TCR associations  

---

### 🧠 Methods
**Tools & Packages:**  
`R` (`immunarch`, `maftools`), `Python` (`scikit-learn`, `pandas`, `matplotlib`), `nf-core` (for reproducible data preprocessing pipelines)  

**Pipeline Summary:**  
1. Mutation calling via Sarek + DeepVariant
- Filtered using DeepVariant filter criteria → shown to outperform other variant calling tools
- Analyzed using maftools
3. TCRs derived from single-cell 10x data
- Filtered out barcodes with multiple α or β reads or low confidence
b. Analyzed using immunarch
- Built binary logistic regression models in scikit-learn using public TCRs 

![Pipeline Workflow](docs/workflow.png)

---

### 🧩 Results

The results of this project highlight several key challenges in linking tumour mutation data with T-cell receptor (TCR) sequences:

- **Limited WES availability:** Whole-exome sequencing (WES) data were rare and difficult to access, making it the most restrictive data type encountered.  
- **Scarcity of matched datasets:** Matched WES and TCR sequencing datasets were especially uncommon, presenting a major barrier to building a sufficiently large cohort for analysis.  
- **Incomplete TCR annotation:** The D region was not annotated in most single-cell datasets, reducing clonotype resolution and limiting analyses that rely on full TCR segment information.  
- **Dataset size limitations:** The small number of matched and annotated samples prevented robust model training and evaluation. A larger dataset will be necessary to accurately assess predictive performance.  

All preprocessing and quality control steps were performed using **nf-core** pipelines to ensure reproducibility and standardized data handling across datasets.


---

### 🖥️ Presentation
📑 [**BCB330 Final Presentation (PDF)**](./docs/BCB330_Final_Presentation.pdf)

---

### ⚙️ Repo Structure
├── Data/ # Processed and sample input datasets, also includes literature review of potential datasets
├── Pipeline/             # End-to-end data processing and analysis workflow  
│   ├── fetchngs/         # Retrieves and organizes raw sequencing data from public repositories  
│   ├── ids/              # Handles sample and patient ID matching across datasets  
│   ├── immunarch/        # Processes and analyzes T-cell receptor (TCR) sequencing data  
│   ├── maftools/         # Performs mutation data analysis and visualization using MAF files  
│   ├── regression/       # Contains scripts for logistic regression and predictive modeling  
│   ├── sarek/            # nf-core Sarek pipeline for variant calling and genomic preprocessing  
│   ├── vcf/              # Handles variant call format (VCF) data extraction and standardization  
│   └── vdj/              # Processes VDJ recombination data and extracts clonotype information
├── docs/ # Project presentation and related images
├── nfcore_scripts/ # Analysis and preprocessing scripts
└── README.md # Project overview

---

### 👩‍🔬 Author
**Victoria Pergola**  
Undergraduate Researcher – Campbell Lab  
Interested in computational immunology, ML for cancer research, and drug discovery.


