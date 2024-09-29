### README for qPCR Primer Design Pipeline (Version 2.0)

This bash script is tailored to automate the process of qPCR primer design for genomic research. It incorporates a series of bioinformatics tools to process genomic sequences, design primers, and assess their specificity. This guide provides detailed instructions on setting up and running this script, along with prerequisites and input requirements.

#### Version 2.0 Updates
- **Improved Pipeline Processes**: The pipeline's internal workflow has been optimized to enhance efficiency.
- **Software Optimization**: Updated usage of several bioinformatics tools to enhance performance.
- **Speed Enhancement**: Primer design speed has increased by 1.5 times due to optimizations in software use and process streamlining.

---

#### Prerequisites
Ensure the following tools are installed and accessible in your environment:
- `BBMap`
- `MPPrimer`
- `e-PCR`
- `Primer3`
- `MFEprimer`
- `BLAT`
- `Python 2.7.18`
- Perl

Set the `PBS_O_WORKDIR` environmental variable to your working directory containing all required files.

#### Expected Input
- **CDS File**: Should be named as `SpeciesName.cds.fa`, e.g., `Bombyx_mori.cds.fa`
- **Genome File**: Should be named as `SpeciesName.fa`, e.g., `Bombyx_mori.fa`
- Files should be placed according to the structured directory setup specified in the script.

#### Pipeline Steps

1. **Gene Name Shortening**: Simplifies gene names to facilitate processing.
2. **Fragmentation**: Splits cDNA sequences into manageable fragments.
3. **Primer3 Input Preparation**: Sets up files for primer design.
4. **Primer Designing**: Designs primers using Primer3.
5. **Primer Formatting**: Formats primers for e-PCR analysis.
6. **Specificity Checks**: Assesses primer specificity via e-PCR.
7. **MFEprimer Analysis**: Indexes and checks primers with MFEprimer for specificity.
8. **Primer File Organization**: Organizes primers by gene into individual files.
9. **Primer Selection**: Filters and selects the best primers based on stringent criteria.
10. **Data Compilation**: Aggregates primer data for database insertion or further analysis.
11. **Best Primers Collection**: Collects the best gene-specific primers, ensuring high specificity.
12. **All Candidate Primers Collection**: Compiles all potential primers that meet initial criteria.

#### Output Files
- Detailed primer lists for each gene with complete parameters.
- Formatted files ready for database use or analytical processes.

#### Usage
Run the script on a UNIX-like system with the necessary computational resources:

```bash
qsub qpcr.sh


