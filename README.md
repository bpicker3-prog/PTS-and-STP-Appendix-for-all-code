# PTS and STP Appendix for AI-Based Protein Design Benchmarking

This repository contains the computational scripts used for the thesis:

**A Critical Evaluation of Current AI-Based Protein Design Methodologies**

**Brett Picker**
Master of Biochemistry
Arizona State University

The purpose of this repository is to provide the complete computational appendix and reproducibility code associated with the thesis. The study evaluates current protein-to-sequence (PTS) and sequence-to-protein (STP) tools independently and as parts of a combined computational protein-design pipeline.

The repository is intended to supplement the written Methods and Appendix. The thesis describes the experimental reasoning, model comparisons, validation metrics, results, and limitations, while this repository provides the scripts used to perform the computational work.

## Overview of the Computational Pipeline

The study is divided into three major computational stages:

**Protein-to-Sequence (PTS)**

Experimentally determined protein structures are used as starting references for sequence generation.

Models evaluated:

* ProteinMPNN
* LigandMPNN
* SolubleMPNN

**Sequence-to-Protein (STP)**

Protein sequences are folded into predicted tertiary structures and compared with their corresponding wild-type references.

Models and conditions evaluated:

* OmegaFold
* ESMFold
* ColabFold with MSA
* ColabFold without MSA
* Boltz-2 with MSA
* Boltz-2 without MSA

**Combined PTS-to-STP Pipeline**

Sequences generated through each MPNN are independently passed through the STP conditions. The resulting predicted structures are compared against their corresponding wild-type structural references.

The general workflow is:

`PDB reference -> PTS sequence generation -> FASTA -> STP structure prediction -> structural comparison -> metric collection`

## Repository Organization

The code is organized according to the computational order used in the thesis.

### Independent PTS Analysis

The independent PTS scripts generate sequences through ProteinMPNN, LigandMPNN, and SolubleMPNN.

The analysis includes:

* model-generated sequences
* model confidence and available internal scores
* wild-type sequence similarity
* individual amino acid retention
* property-based residue-class retention
* computational time
* secondary-structure comparisons

For the reported secondary-structure analysis, proteins were separated into:

* helix-rich
* sheet-rich
* mixed

The original development scripts may contain additional exploratory buckets, including protein-length and loop-rich categories. These exploratory categories were not included in the final comparative results reported in the thesis.

### Amino Acid and Residue-Class Analysis

Generated sequences were compared with their corresponding wild-type sequences to determine which amino acids were retained or replaced.

For the property-based analysis, amino acids were grouped as:

* Hydrophobic: A, V, L, I, M
* Hydrophilic: S, T, N, Q, C, G
* Aromatic: F, W, Y
* Charged: K, R, H, D, E
* Special: P

These categories are analytical classifications used within this study and are not intended to represent the only possible biochemical classification of amino acids.

The original transition matrices and heatmap-generation scripts are included so that the summarized thesis figures can be traced back to the underlying sequence data.

### Independent STP Analysis

The independent STP scripts evaluate the folding tools using wild-type sequences derived from experimentally determined PDB structures.

The tested conditions include OmegaFold, ESMFold, ColabFold with and without MSA, and Boltz-2 with and without MSA.

This stage served as an initial benchmark of the folding workflow. During analysis, a reference-structure limitation involving biological assemblies was identified. For this reason, the independent STP results should not be treated as the primary standardized comparison of absolute structural performance.

The complete PTS-to-STP pipeline was subsequently modified to retain and standardize the appropriate biological structural references.

### Complete PTS-to-STP Pipeline

The complete pipeline generates sequences through the three MPNNs and passes each generated sequence independently through the six STP conditions.

Biological assemblies are retrieved and standardized before sequence generation. Generated sequences are then processed through:

1. OmegaFold
2. ESMFold
3. ColabFold without MSA
4. ColabFold with MSA
5. Boltz-2 without MSA
6. Boltz-2 with MSA

The resulting structures are compared against the corresponding wild-type reference using a common structural comparison procedure within the complete pipeline.

This combined analysis represents the primary standardized comparison used to determine how the behavior of an STP model changes according to the PTS model that generated its input sequence.

## Structural and Sequence Metrics

Several complementary measurements are used because no individual metric describes the complete protein-design process.

### Wild-Type Sequence Similarity

Generated MPNN sequences are compared with their corresponding wild-type sequences to determine the proportion of positions retaining the original amino acid identity.

### Amino Acid Retention

Individual amino acid transitions are measured to determine whether different MPNNs preferentially retain or replace particular residues.

### Property-Based Residue Retention

Residue substitutions are reduced into the five property-based classes defined above to determine whether sequence changes preserve broader residue characteristics even when exact amino acid identity changes.

### Backbone RMSD

Predicted structures are structurally compared against their corresponding wild-type references.

Within the complete PTS-to-STP worker, Cα positions of the protein backbone are used for structural superposition and RMSD calculation. Lower RMSD represents greater structural agreement with the reference structure.

RMSD is used as a measurement of structural recovery and should not be interpreted as direct evidence of thermodynamic stability or biological function.

### Model Confidence

Available model-confidence information is collected from the structural prediction outputs. Confidence represents the prediction model's assessment of its own output and is interpreted alongside independent structural comparison rather than as a replacement for RMSD.

### Computational Time

Runtime is collected to evaluate the computational cost of each method and to determine how the tools scale when applied to large sequence-generation and structure-prediction datasets.

## Computational Environment

The scripts were developed for execution on the **Arizona State University Sol supercomputer** using SLURM job scheduling.

The workflow uses Bash as the primary scripting environment, with Python programs executed from within Bash where required.

The scripts use multiple software environments because the evaluated models require different dependencies. Examples include environments for:

* LigandMPNN / ProteinMPNN / SolubleMPNN
* OmegaFold
* ESMFold
* ColabFold
* Boltz
* PyMOL

GPU jobs were primarily configured for NVIDIA A100 resources where available.

Exact environment paths, SLURM resource requests, model arguments, and command-line settings are retained in the scripts to document the computational conditions used during the study.

## Reproducibility

The scripts in this repository are provided as a record of the computational procedures used to generate the thesis data.

The code is intentionally retained close to the form in which it was executed rather than rewritten into a generalized software package after completion of the experiments. As a result, some scripts contain environment paths, SLURM settings, directory structures, or variables specific to the ASU Sol environment.

Researchers reproducing the workflow on another computing system will likely need to modify:

* software environment paths
* SLURM partitions and resource requests
* local directory paths
* model installation locations
* database/cache locations
* GPU configuration

These changes should not alter the general computational workflow.

## Important Methodological Notes

The independent STP benchmark and the final combined pipeline should be distinguished when interpreting the code.

A biological-assembly reference issue was identified during the independent STP analysis. The later combined pipeline was modified so that biological assemblies were handled before MPNN sequence generation and subsequent structural prediction.

The repository preserves the scripts associated with both stages so that the development and final implementation of the computational workflow remain transparent.

Some development scripts also contain exploratory structural or protein-length categories that were not used in the final thesis analysis. These portions are retained because the repository serves as a record of the executed computational workflow rather than a rewritten representation of only the final figures.

## Interpretation of the Repository

The code provided here reproduces computational sequence generation, structural prediction, metric collection, and analysis. It does not experimentally validate the biological activity of the generated proteins.

A predicted structure with favorable RMSD or model confidence does not independently demonstrate:

* protein expression
* thermodynamic stability
* solubility
* ligand affinity
* catalytic activity
* biological function

Experimental structural and biochemical validation remain necessary for those conclusions.

## Thesis Relationship

The repository corresponds to the computational Methods, Results, and Appendix of the thesis.

The general organization is:

* Independent PTS analysis -> MPNN sequence behavior
* Independent STP analysis -> folding benchmark
* Combined PTS-to-STP analysis -> complete computational design workflow
* Analysis scripts -> residue transitions, structural metrics, confidence, and computational performance

The thesis should be consulted for the complete experimental rationale, interpretation of the results, limitations, and discussion of appropriate model use cases.

## Citation

If using this repository or its benchmarking framework in academic work, please cite the associated thesis.

A complete thesis citation will be added following final publication through Arizona State University.

## Author

**Brett Picker**
Arizona State University
Master of Biochemistry

## Repository Status

This repository accompanies an academic thesis and may be updated as the thesis is finalized. The version associated with the final submitted thesis should be used when reproducing the reported results.
