# core-unit-bioinformatics/workflow-nxf-somatic-variantcomp: Citations

Please cite the tools and datasets below when you use results of this pipeline. The
pipeline step(s) that use each tool are given in brackets; the full step-by-tool overview
is in the [README](README.md#pipeline-summary).

## Framework

- [nf-core](https://pubmed.ncbi.nlm.nih.gov/32055031/) (pipeline template, modules)

  > Ewels PA, Peltzer A, Fillinger S, Patel H, Alneberg J, Wilm A, Garcia MU, Di Tommaso P, Nahnsen S. The nf-core framework for community-curated bioinformatics pipelines. Nat Biotechnol. 2020 Mar;38(3):276-278. doi: 10.1038/s41587-020-0439-x. PubMed PMID: 32055031.

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/) (workflow engine)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

- [nf-schema](https://github.com/nextflow-io/nf-schema) (parameter and samplesheet validation)

  > nf-schema, Nextflow plugin, Seqera and nf-core contributors. https://github.com/nextflow-io/nf-schema

- [nf-test](https://www.nf-test.com) (pipeline tests; development only)

  > Forer L, Schönherr S. Improving the reliability, quality, and maintainability of bioinformatics pipelines with nf-test. GigaScience. 2025;14:giaf130. doi: 10.1093/gigascience/giaf130.

## Pipeline tools

- [BCFtools, SAMtools and HTSlib](https://pubmed.ncbi.nlm.nih.gov/33590861/) (sample selection, PASS split, contig renaming, cleaning, reheader, normalization, concatenation, filtering, exact/position matching, read evidence, VEP field extraction with `+split-vep`; alignment indexing; reference index)

  > Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. GigaScience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PubMed PMID: 33590861.

- [BEDTools](https://pubmed.ncbi.nlm.nih.gov/20110278/) (target BED sorting, N-way intersection, merging, region restriction and exclusion; callable regions; repeat and stratification annotation)

  > Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 2010 Mar 15;26(6):841-842. doi: 10.1093/bioinformatics/btq033. PubMed PMID: 20110278.

- [mosdepth](https://pubmed.ncbi.nlm.nih.gov/29096012/) (coverage and callable bases of the alignments)

  > Pedersen BS, Quinlan AR. Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics. 2018 Mar 1;34(5):867-868. doi: 10.1093/bioinformatics/btx699. PubMed PMID: 29096012.

- [RTG Tools (vcfeval)](https://github.com/RealTimeGenomics/rtg-tools) (haplotype-aware matching of truth and query calls)

  > Cleary JG, Braithwaite R, Gaastra K, Hilbush BS, Inglis S, Irvine SA, Jackson A, Littin R, Rathod M, Ware D, Zook JM, Trigg L, De La Vega FM. Comparing variant call files for performance benchmarking of next-generation sequencing variant calling pipelines. bioRxiv. 2015. doi: 10.1101/023754.

- [Ensembl VEP](https://pubmed.ncbi.nlm.nih.gov/27268795/) (consequence, gene, impact and variant class of every compared variant)

  > McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GRS, Thormann A, Flicek P, Cunningham F. The Ensembl Variant Effect Predictor. Genome Biol. 2016 Jun 6;17(1):122. doi: 10.1186/s13059-016-0974-4. PubMed PMID: 27268795.

- [MSIsensor-pro](https://pubmed.ncbi.nlm.nih.gov/32171661/) (`scan`: microsatellite/repeat regions of the reference)

  > Jia P, Yang X, Guo L, Liu B, Lin J, Liang H, Sun J, Zhang C, Ye K. MSIsensor-pro: Fast, Accurate, and Matched-normal-sample-free Detection of Microsatellite Instability. Genomics Proteomics Bioinformatics. 2020 Feb;18(1):65-71. doi: 10.1016/j.gpb.2020.02.001. PubMed PMID: 32171661.

- [GNU Awk](https://www.gnu.org/software/gawk/) (genome BED from the reference index; selection of intervals covered by all BEDs)

  > Free Software Foundation. GNU Awk. https://www.gnu.org/software/gawk/

- [datavzrd](https://github.com/datavzrd/datavzrd) (table browser of the report)

  > Wiegand F, Lähnemann D, Mölder F, Uzuner H, Prinz A, Schramm A, Köster J. datavzrd: Rapid programming- and maintenance-free interactive visualization and communication of tabular data. PLOS ONE. 2025;20(7):e0323079. doi: 10.1371/journal.pone.0323079.

- [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/) (report)

  > Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016 Oct 1;32(19):3047-8. doi: 10.1093/bioinformatics/btw354. PubMed PMID: 27312411; PubMed Central PMCID: PMC5039924.

## Python libraries used by the pipeline scripts (`bin/`)

- [pysam](https://github.com/pysam-developers/pysam) (VCF, alignment and FASTA reading in the reference checks; VCF reading and writing in all VCF scripts; wraps HTSlib, see above)

  > pysam developers. pysam: a Python module for reading and manipulating SAM/BAM/VCF/BCF files. https://github.com/pysam-developers/pysam

- [pandas](https://doi.org/10.25080/Majora-92bf1922-00a) (metric aggregation, report tables)

  > McKinney W. Data structures for statistical computing in Python. Proceedings of the 9th Python in Science Conference. 2010:56-61. doi: 10.25080/Majora-92bf1922-00a.

- [NumPy](https://pubmed.ncbi.nlm.nih.gov/32939066/) (report statistics)

  > Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. Nature. 2020 Sep;585(7825):357-362. doi: 10.1038/s41586-020-2649-2. PubMed PMID: 32939066.

- [SciPy](https://pubmed.ncbi.nlm.nih.gov/32015543/) (Fisher's exact test and Benjamini–Hochberg correction in the enrichment statistics)

  > Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods. 2020 Mar;17(3):261-272. doi: 10.1038/s41592-019-0686-2. PubMed PMID: 32015543.

- [Matplotlib](https://doi.org/10.1109/MCSE.2007.55) (plots)

  > Hunter JD. Matplotlib: a 2D graphics environment. Computing in Science & Engineering. 2007;9(3):90-95. doi: 10.1109/MCSE.2007.55.

- [UpSetPlot](https://github.com/jnothman/UpSetPlot) (UpSet plots), implementing the UpSet technique:

  > Lex A, Gehlenborg N, Strobelt H, Vuillemot R, Pfister H. UpSet: visualization of intersecting sets. IEEE Trans Vis Comput Graph. 2014 Dec;20(12):1983-1992. doi: 10.1109/TVCG.2014.2346248. PubMed PMID: 26356912.

## Statistical methods

- Benjamini-Hochberg correction of the enrichment p-values

  > Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57(1):289-300. doi: 10.1111/j.2517-6161.1995.tb02031.x.

## Variant callers of the tested inputs

These callers produced the VCFs the pipeline was tested with or is designed for (see
[`docs/usage.md`](docs/usage.md#tested-inputs)); cite them for the calls you compare.

- [Illumina DRAGEN](https://help.dragen.illumina.com) (target caller of the real use case, v4.4.4)

  > Behera S, Catreux S, Rossi M, et al. Comprehensive genome analysis and variant detection at scale using DRAGEN. Nat Biotechnol. 2025;43(7):1177-1191. doi: 10.1038/s41587-024-02382-1.

- [GATK MuTect2](https://gatk.broadinstitute.org) (SEQC2 test VCFs)

  > Cibulskis K, Lawrence MS, Carter SL, et al. Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples. Nat Biotechnol. 2013 Mar;31(3):213-219. doi: 10.1038/nbt.2514. PubMed PMID: 23396013.

  > Benjamin D, Sato T, Cibulskis K, Getz G, Stewart C, Lichtenstein L. Calling somatic SNVs and indels with Mutect2. bioRxiv. 2019. doi: 10.1101/861054.

## Datasets

- SEQC2 somatic reference samples HCC1395 / HCC1395BL (test data: WES and WGS VCFs and
  alignments, high-confidence somatic call set v1.2.1)

  > Fang LT, Zhu B, Zhao Y, et al. Establishing community reference samples, data and call sets for benchmarking cancer mutation detection using whole-genome sequencing. Nat Biotechnol. 2021 Sep;39(9):1151-1160. doi: 10.1038/s41587-021-00993-6.

  > Xiao W, Ren L, Chen Z, et al. Toward best practice in cancer mutation detection with whole-genome and whole-exome sequencing. Nat Biotechnol. 2021 Sep;39(9):1141-1150. doi: 10.1038/s41587-021-00994-5.

- GIAB genomic stratifications v3.5 (annotation of the compared variants with the `--giab_<name>` files: homopolymers, tandem repeats, segmental duplications, low mappability, all difficult regions; also in the SEQC2 test data)

  > Dwarshuis N, Kalra D, McDaniel J, et al. The GIAB genomic stratifications resource for human reference genomes. Nat Commun. 2024 Oct 19;15(1):9029. doi: 10.1038/s41467-024-53260-y. PubMed PMID: 39424793.

- GENCODE release 46 (test gene annotation for VEP, `--vep_gtf`)

  > Frankish A, Carbonell-Sala S, Diekhans M, et al. GENCODE: reference annotation for the human and mouse genomes in 2023. Nucleic Acids Res. 2023 Jan 6;51(D1):D942-D949. doi: 10.1093/nar/gkac1072. PubMed PMID: 36420896.

## Software packaging/containerisation tools

- [Anaconda](https://anaconda.com)

  > Anaconda Software Distribution. Computer software. Vers. 2-2.4.0. Anaconda, Nov. 2016. Web.

- [Bioconda](https://pubmed.ncbi.nlm.nih.gov/29967506/)

  > Grüning B, Dale R, Sjödin A, Chapman BA, Rowe J, Tomkins-Tinch CH, Valieris R, Köster J; Bioconda Team. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods. 2018 Jul;15(7):475-476. doi: 10.1038/s41592-018-0046-7. PubMed PMID: 29967506.

- [BioContainers](https://pubmed.ncbi.nlm.nih.gov/28379341/)

  > da Veiga Leprevost F, Grüning B, Aflitos SA, Röst HL, Uszkoreit J, Barsnes H, Vaudel M, Moreno P, Gatto L, Weber J, Bai M, Jimenez RC, Sachsenberg T, Pfeuffer J, Alvarez RV, Griss J, Nesvizhskii AI, Perez-Riverol Y. BioContainers: an open-source and community-driven framework for software standardization. Bioinformatics. 2017 Aug 15;33(16):2580-2582. doi: 10.1093/bioinformatics/btx192. PubMed PMID: 28379341; PubMed Central PMCID: PMC5870671.

- [Seqera Containers](https://seqera.io/containers/) (container of the report step)

  > Seqera. Seqera Containers. https://seqera.io/containers/

- [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241)

  > Merkel, D. (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal, 2014(239), 2. doi: 10.5555/2600239.2600241.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

## Overview figure

- The metro map in the README uses the style of the nf-core metro-map components
  (CC0; James A. Fellows Yates, Maxime U. Garcia, Louis Le Nézet), see
  https://nf-co.re/docs/community/brand/workflow-schematics.
