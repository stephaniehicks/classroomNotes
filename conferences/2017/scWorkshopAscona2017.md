# Notes on talks for the 2017 Ascona Workshop on Statistical Challenges in Single-Cell Biology

#### Organizing committee

* Niko Beerenwinkel
* Peter Bühlmann
* Wolfgang Huber

Notes and slides for the [2017 Ascona Workshop](https://www.bsse.ethz.ch/cbg/cbg-news/ascona-2017.html) 
at Monte Verita, Ascona, Switzerland from Apr 30 - May 5, 2017. Here is a pdf of the 
[program](https://www.ethz.ch/content/dam/ethz/special-interest/bsse/cbg/events/ascona-2017/2017-04-28%20program.pdf) 
and [abstract book](https://www.ethz.ch/content/dam/ethz/special-interest/bsse/cbg/events/ascona-2017/abstract_book.pdf). 



## Monday, May 1

- [John Marioni](http://www.ebi.ac.uk/research/marioni) from [@emblebi](https://twitter.com/emblebi) and [@CRUKresearch](https://twitter.com/crukresearch), Identifying differentially variable genes from single-cell RNA-sequencing data and applications in immunity and aging 
	- Highlighting new computational methods in scRNA-seq
		- Normalization: [Lun et al. (2016)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7) implemented in [scran R/Bioconductor pkg](http://bioconductor.org/packages/release/bioc/html/scran.html), [Bacher et al. (2017)](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4263.html?WT.feed_name=subjects_systems-biology) implemented in [SCnorm](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4263.html?WT.feed_name=subjects_physical-sciences), [Vallejos et al (2017)](https://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4292.html)
		- Confounding factors: Buettner et al. (2015), Finak et al. (2015), Hicks et al. bioRxiv
		- Spatial expression: Achim et al. (2015), Satija et al. (2015)
		- Clustering and pseudotime: [Trapnell et al. (2014)](http://www.nature.com/nbt/journal/v32/n4/full/nbt.2859.html), [Haghverdi et al. (2016)](http://www.nature.com/nmeth/journal/v13/n10/full/nmeth.3971.html), [Kiselev et al. (2017)](https://www.nature.com/nmeth/journal/v14/n5/full/nmeth.4236.html)
	- Interested in identifying Highly Variable Genes (HVGs)
		- HVGs = genes that have more variance in a population than what you expected by chance. Idea being they might hold some substructure in the population vs e.g. housekeeping genes. 
	- Discussed why spike-ins are useful and idea behind them (should be expressed at the same level in every cell and can be used to quantify technical noise)
	- Discussed [BASiCs](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004333) (Bayesian Analysis of Single-Cell Sequencing data): 
		- jointly model endogenous and spike-ins together and inferring normalization and technical noise
		- Can find differentially expressed **and** differentially variable genes
		- uses adaptive MH with Gibbs sampling
	- Demonstrated applications of BASiCs to find genes that are different in variability (but have the same means)
	- Challenges:  
		- Putting in the mean-var relationship into the prior distribution
		- If no spike-ins, put a correlated prior on the mean 
		- Applications to other single-cell data types (e.g. CyTOF)
		- Spatial modeling
		- Data integration and genetics (Using "mult-view" factor analysis models)
	- Nick Navin asked about using UMIs vs ERCC spike-ins in terms of normalization. John describes spike-ins are not perfect in practice. Discussed the idea of using a "control cell" in droplet methods. Also another reason why BASiCs has been adapted for not having spike-ins. 

- [Peter Kharchenko](http://pklab.med.harvard.edu), From one to millions of cells: computational approaches for single-cell analysis
	- Discussing rapid growth in cell numbers, data volume, more multi-patient samples
	- Discussing variability in scRNA-seq data
		- e.g. differences between cells (of the same type), overdispersion, measurement failures, quality of cells
		- can be problems for PCA (non-gaussian, noisy cells from outliers, can makes true heterogeneity)
		- Biological and technical variability (cell-specific noise models)
	- Discussed how to model noise of an individual cell (dropout models)
	- Used known gene sets to get better power to detect patterns of variability 
	- PAGODA: weighted PCA - weights observations that the probability that it's a technical dropout
		- pick up correlated gene sets from de novo gene sets and annotated gene sets
		- "Aspects" = clusters of PCs
	- PAGODA2 - optimized for large sparse measurements
	- Lanczos algorithm. Must be able to break down computations to matrix multiplication. 
	- Challenges: 
		- controlling for sequencing depth. correlates well with cell size. You may kill biological signal if you simply regress it out.
		- sparse measurements. how to process on sparse measurements: nearest neighbors. 
		- visualization on sparse matrices. 
		- making use of reference data. (e.g. organism, tissue-specific data) 
		- distances. most groups use correlation based distance. The problem is there is no universal distance estimates. Depends on question (e.g. if you are looking at cell cycle genes, you should define distance based on cell cycle genes). 

- [Ambrose Carr](http://ambrosejcarr.com), [@ambrosejcarr](https://twitter.com/ambrosejcarr), Bayesian Inference for Single-cell ClUstering and ImpuTing (BISCUIT)
	- Profile tumor immune cells in breast cancer. Used inDrop platform. 
	- cell-type specific capture rates, dropouts
	- Discussed [BISCUIT](https://genomicscomputbiol.org/ojs/index.php/GCB/article/view/46): Bayesian model used to iteratively (1) impute dropouts based on expression or co-expression with other genes and (2) normalize data across library size. 
		- global normalization fails b/c cells have different sizes. Dropouts are not resolved and may introduce spurious correlations. 
	- For each gene, algorithm models clusters of cells using bayesian mixture model (mixtures across mixtures of cell types) 
	- Discussed why downsampling is not appropriate (because you throw out too much data e.g. 90%)

- [Stephanie Hicks](http://www.stephaniehicks.com), [@stephaniehicks](https://twitter.com/stephaniehicks), Missing data and technical variability in single-cell RNA-sequencing experiments

- [Yuanhua Huang](http://homepages.inf.ed.ac.uk/s1333321/), Transcriptome-wide splicing quantification in single cells
	- Two general approaches to quantify isoforms: (1) directly measure junction reads and count, (2) probabilistic method using all reads and infers counts
	- In bulk RNA-seq, if you don't have enough reads, you can increase amplification, but in scRNA-seq the starting number of molecules is very small (dropouts). 
	- [BRIE](https://brie-rna.sourceforge.io) - Bayesian regression for isoform quantification in scRNA-seq (python)
		- uses sequence features as prior information of which exon to include in exon-skipping events
		- This prior information increases 

- [Jean (Zhijin) Wu](http://www.stat.brown.edu/zwu/), The two phases in gene expression regulation
	- Proposes two phases of "off" and "on" genes and the volume is dialed up and down based on cell-wide global effects
	- Method to assess difference in proportion of detected genes and differential expression
		- Empirical Bayes estimation of hyper parameters
		- Method includes a latent indicator Z for whether the gene was expressed
		- Conditional on latent status: expression is modeled as a zero-inflated poisson OR a poisson with parameter as a log normal. 
		- Interested in classification of phases

- [Felix Kurt](https://www.bsse.ethz.ch), Single-cell analysis on microfluidic platforms


## Tuesday, May 2

- [Kobi Benenson](https://www.bsse.ethz.ch/synbio), Synthetic gene circuits for in situ cell classification
	- Interested in finding molecular markers that can discriminate between tumor and healthy cells with a goal to target specific cells using molecular markers 
		- leads to increased potency of tumor cells
	- What types of classifiers? e.g. PCA, SVM, boolean/logic classifiers
		- Uses logic circuits as a platform for biological computing (models molecular interactions in cells)
	- What markers? **microRNAs** (miRNAs) are most informative b/c they have high discrimination power [Golub group, Nature 2005](http://www.nature.com/nature/journal/v435/n7043/full/nature03702.html)
	- proof of concept experiment: Built a logic gene circuit to identify HeLa cells using 6 discriminating miRNAs (identified from literature)
		- The hardest part is to optimize logic gene circuit to give large dynamic range (ratio of HeLa to non-HeLa cells)
	- Applications to [classification of cells in breast cancer](http://www.cell.com/cell-systems/abstract/S2405-4712(17)30003-0) 
	- At single-cell level, there is noise in the input to the logic circuits which propagates to even more noise in output of circuits

- [Bernd Bodenmiller](http://www.bodenmillerlab.org), Highly multiplexed analysis of the tumor ecosystem by mass cytometry
	- Interested in the spatial analysis of cancer tissues
	- Uses **mass cytometry** to identify cell types and functional states based on multiplexed marker expression 
		- Pros: can identify viability, size, DNA content, cell cycle; measures 52 parameters simultaneously and 1000 cells per second
		- Cons: no live cells, relies on antibodies being available
	- Used t-SNE and [PhenoGraph](https://www.c2b2.columbia.edu/danapeerlab/html/phenograph.html) to identify clusters of cell communities
	- Used diffusion maps to describe macrophages as a continuum of cell phenotypes
	- Developed [Imaging Mass Cytometry (IMC)](http://www.nature.com/nmeth/journal/v11/n4/full/nmeth.2869.html) to study cell-to-cell interactions
		- Combines mass cytometry with a laser to give you spatial information in a tissue (to identify cell type and state); can image up to > 50 markers
		- [miCAT](http://biorxiv.org/content/early/2017/02/17/109207) = computational tool to process and analyze multiplexed image cytometry data
	
- [Christof Seiler](http://christofseiler.github.io), CytoGLMM: Bayesian hierarchical linear modeling for **flow cytometry** data
	- Most people related CyTOF data to phenotypic outcomes by clustering cells into subpopulations and then related cluster features to outcomes 
		- e.g. FlowMap-FR, Citrus, flowType-RchyOptimyx, FloReMi, and COMPASS
		- Bad b/c it assumes discrete set of cell subpopulations
	- Introduced R package CytoGLMM (using Stan) to skip clustering step and use GLMM where the response is outcome and explantory variables are 40 protein expression profiles. 
		- Uses GLMM with logistic regression to account for donor-specific variability (e.g. logistic regression: X = 37 transformed protein counts; Y = stimulus condition)

- [Eirini Arvaniti](http://www.imsb.ethz.ch/research/claassen/people/eiriniarvaniti.html), Sensitive detection of rare disease-associated cell subsets via representation learning
	- Interested in detecting subsets of cells associated with phenotypes
		- Identify cell clusters (unsupervised), Compute features (cluster abundances per sample), associate features with phenotype (supervised) 
		- Developed [CellCnn](https://www.ncbi.nlm.nih.gov/pubmed/?term=cellcnn) (detects phenotype-associated cell subsets (gating) using neural networks) which is an improvement on [Citrus](https://www.ncbi.nlm.nih.gov/pubmed/24979804) (data-driven approach for the identification of stratifying subpopulations in multidimensional cytometry datasets). 
			- Good for rare-subpopulations, fast, scales well to larger cohorts
			- Next: will apply to scRNA-seq data and will apply to larger cohorts

- [Raphael Gottardo](http://rglab.org), [@raphg](https://twitter.com/raphg), Identifying novel correlates of protection via single-cell analysis of antigen-specific T-cells
	- Investigate T cell diversity and associate with clinical phenotypes
	- fFow cytometry is still most common assay to study T cell antigens (Intracellular cytokine staining (ICS) assay)
	- Suggests building a gold standard scRNA-seq data set, similar to [FlowCap](http://www.jimmunol.org/content/186/1_Supplement/65.2)
	- Data: flow sort T cells to identify antigen-specific T cells and use **single cell qPCR** or **scRNA-seq** data
	- [MAST](http://bioconductor.org/packages/release/bioc/html/MAST.html): model-based analysis of single cell qPCR, NanoString and RNA-Seq
		- Differential expression, gene set enrichment analysis, CDR: technical and biological [Padovan-Merhar et al. (2015)](https://www.ncbi.nlm.nih.gov/pubmed/25866248) variation
	- Challenges:
		- Integrate multiple data sources (possibly at the single-cell level)
			- surface markers, cytokin secretion and gene expression
		- Increase the dimensionality of single-cell assays 	 

- [Davis McCarthy](https://sites.google.com/site/davismcc/), [@davisjmcc](https://twitter.com/davisjmcc), Mapping genetic effects on interactions with single-cell states 
	- Interested in modeling cell to cell heterogeneity when relating human variation (allele status) to genetic variation 
	- DNA variants in cis with a gene
		`Y = covars + SNP + g + e where g ~ N(0, K \sigma_g), e ~ N(0, | \sigma_e)`
	- Use scRNA-seq for eQTL analyses (e.g. cell differentiation)
	- How to overcome batch in different donors: mix cells from four donors, grow cells up, and use genomic relatedness score to identify donor cell with high accuracy from scRNA-seq reads
	- Deriving phenotypes from scRNA-seq data
		- Questions: should regress out cell-level variables (# detected genes, lib complexity)? normalization? 
	- Look at both mean and variance in eQTLs
		

## Wednesday, May 3

- [Prisca Liberali](http://www.fmi.ch/research/groupleader/?group=135), Mapping genetic interactions during intestinal organoid self- organization
	- Interested in how single cells create organs or tissue
	- Stem cells have an intrinsic self-organization potential
	- Use organoids as models for tissue organization and complex disease
	- [Light Sheet Microscope](https://en.wikipedia.org/wiki/Light_sheet_fluorescence_microscopy) = keeps organoids alive for 5 days and tracks single cells
	- Challenges: 
		- How to follow TF in unsyncronized cells? Use [Wanderlust](https://www.ncbi.nlm.nih.gov/pubmed/24766814) and [Wishbone](http://www.nature.com/nbt/journal/v34/n6/full/nbt.3569.html) to look at trajectory and progression of organioid growth. 
		- How to call phenotypes? Use [Phenograph](https://www.c2b2.columbia.edu/danapeerlab/html/phenograph.html) 
		- How to screen for organoid formation and infer genetic interactions during self-organization? Hierarchical Interaction Score (HIS) to look at enrichment of classes

- [Takashi Hiragi](https://www.embl.de/research/units/dev_biology/hiiragi/), Symmetry breaking and self-organisation in mouse development
	- Use single-cell biology to study early mammalian development 

- [Andre Rendeiro](http://andre-rendeiro.com), [@afrendeiro](https://twitter.com/afrendeiro), Pooled CRISPR screening with single-cell transcriptome readout
	- [CROP-seq](https://www.nature.com/nmeth/journal/v14/n3/fig_tab/nmeth.4177_F1.html) 
		- pooled CRISPR screen = pool guide gRNAs (high throughput, simple readout only)
		- arrayed CRISPR screen = each perturbation happens individually (complex read out, relatively low throughput
	- scRNA-seq microfludics approach. single cells are individual perturbed (genetic knockout) and then phenotype is observed (transcriptomes)
		- How to identify which gRNA perturbed the cell? Integrate a CROPseq-Guide-Puro sgRNA twice (once before gene and once at poly A tail for sequencing
		- More scalable than Perturb-seq, CRISP-seq
	- Validated using a synthetic mixture of cells (controls and knock-out cells)
	- Now moving from drop-seq to 10X genomics platform. 
	- Challenges: 
		- Is variability from perturbation, is it meaningful or can you model this? 

- [Lars Velten](https://www.embl.de/research/units/genome_biology/steinmetz/members/index.php?s_personId=CP-60016232), Quantifying developmental plasticity by integrated single-cell RNA-seq and ex vivo culture
	- Development potential (the potential to change cell states given a perturburation), Developmental fate ()
	- [STEMNET](https://git.embl.de/velten/STEMNET): 
		- Pseudotime with > 2 endpoints
		- Uses a multinomial classifier with elastic net regularization computes class probabilities that HSPCs belonging to different mature lineages. (sparse genes, projects data on simplex)

## Thursday, May 4

- [Valérie Taly](http://recherche.parisdescartes.fr/UMRS1147/Translational-Research-Microfluidics/Group), Droplet-based microfluidic for single cell analysis
	- Droplet based microfluidics.
		- Cells remain viable for days (14) in droplets in perfluorocarbon carrier fluid. Surfactant ensures droplet stability (aka droplets don't combine)
		- Cells can incubate and grow in microdroplets (reproduce every 1.5hr), but need larger droplets. Can identify alive or dead cells and can count alive cells in each droplet. 
	- Can use microfluidics for cancer research (personalized medicine) (e.g. assessing drug resistance at single cell level)
	- Can include one gene per droplet. Multiplexed allows now for up to 7 mutations. 
		- That way you can see which molecule is wild type vs mutant (KRAS mutations are resistant to treatment).  
		- Found patients with KRAS mutations in colorectal cancer
	- Can use droplet-based circulating tumor (Ct) DNA (CtDNA) (early marker for cancer progression) measurement for treatment follow up and progression free survival. 
		- measured 9 plasma samples/patients.  
	- Can use CtDNA measurement for detection of recurrence for early colorectal cancer 
		- Measured plasma before and after surgey, and every 2-5months for three years.  
		- Could predict cancer recurrence prior to 

- [Nicola Aceto](https://biomedizin.unibas.ch/nc/about-us/people/profil/profile/person/aceto/), Single cell analysis of circulating tumor cell clusters
	- Circulating tumor cells (CTCs). One in a billion cells is a CTC (huge technical challenge). Most die in circulation, but a few will survive and form a metastasis. 
	- Technologies for CTC capture. Worked with engineers to design chips. Microfludics devices are required to isolate CTCs from blood samples. 
		- Sample: 10mL whole blood = 50 billion RBCs, 50 million WBCs, 0-100 CTCs 
		- [First gen microfludics chip](https://www.nature.com/nature/journal/v450/n7173/full/nature06385.html)
		- [Second gen herringbone chip](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2972993/). Higher capture rate for CTCs. But hard to capture CTCs because certain cells don't expresses epithelial markers to detect CTCs. 
		- [Third gen i-chip](http://stm.sciencemag.org/content/5/179/179ra47). magnetic bead-tagged CTCs or magnetic bead-tagged WBCs. Tag the WBCs and push into waste channel. CTCs remain unlabeled and live which can be analyzed. 
	- Why do we care about CTC-clusters? 
		- Dogma: Wildly assumed single CTCs enter blood stream and metastasize. Turns out CTCs like to cluster in blood stream, especially with WBCs. 
		- CTC-clusters in breast cancer and prostate cancer. [Aceto et al. Cell 2014](https://www.ncbi.nlm.nih.gov/pubmed/25171411). Patients with clusters had lower overall survival vs patients with single CTCs
		- Mouse models to ask (1) what is the metastatic potential of CTCs compared to single cells. (2) Where to CTC-clusters come from? 
			- Found CTC-clusters are mixed from different intial clusters (colored by red and green). Oligoclonal, not monoclonal. 
			- 2.6% of all CTC events are CTC clusters. 91% of CTC-clusters are multicolored. 53% of all lung foci are multicolored. CTC-clusters are ~50-fold more metastasis-component than single CTCs. 
	- single cell analysis of CTC-clusters (to figure out why they were so much more metastatic vs single cells)
		- identified genes relevant for CTC-clsuters vs CTC single cells. e.g. Plakoglobin helps tether CTCs within clusters, thereby enhancing metastatic spread. 
		- explored how CTC-clusters manage to move through very small human capillaries  
		- [Designing the Cluster-CHIP](http://www.nature.com/nmeth/journal/v12/n7/fig_tab/nmeth.3404_F1.html). Allows the passage for all blood cells, but if you are a cluster of more than 2 cells, then you get trapped at the end. 
		- Measure mutations in patients + CTCs allows you to apply a mutation-based screen to test cancer drug screens. 
	- Heterogeneity of CTC clusters
		- Open questions
			- Are cells within ctc-clusters harboring different mutational profile? 
			- How these mutations impact on gene expression
			- How does CTC profile change during disease progression and therapy resistance
			- What is the identity and role of immune cells attached to CTC clusters
		- Start with blood samples from patient (blue) or mouse models (green). put CTC clusters into tube, put CTC single cells in tube. put WBCs in tubes. 
		- Use single cell Parallel Exome and Transcriptome Squencing (E&T-Seq)
			- from each cell we have a mutational and RNA profile. 

- [Jan Korbel](https://www.embl.de/research/units/genome_biology/korbel/), Single cell-based detection of diverse classes of genomic DNA rearrangements
	- Interested in the formation and selection of structural variation (SVs). SNPs (0.1% of genome) vs heritable SVs (0.3-0.5% of genome)
	- Data: Integrated map of SVs in 2504 human genomes [Sudmant et al. Nature 2015](http://www.nature.com/nature/journal/v526/n7571/full/nature15394.html), [Auton et al. Nature 2015](https://www.nature.com/nature/journal/v526/n7571/full/nature15393.html#group-1)
	- Old technology: use paired-end mapping to discover SVs [Korbel et al. Science 2007](https://www.ncbi.nlm.nih.gov/pubmed/17901297)
		- But method is bad for presence of repeat SV breakpoints, chimeras can complicate SV detection in single cells, method requires high-coverage (>10-fold) to make accurate calls.  
		- Also, SVs associated with repetitive DNA. Inversions are hard to identify and hence mostly overlooked in analyses of human genomes. [Antonacci et al. Hum Mol Gen 2009](https://www.ncbi.nlm.nih.gov/pubmed/19383631). 
	- [Strand-seq](https://www.ncbi.nlm.nih.gov/pubmed/23665005). sequencing [strand-directional DNA reads](https://en.wikipedia.org/wiki/Single-cell_DNA_template_strand_sequencing)
		- complimentary signals can be extracted from data - haplotype, strand, read depth (this is key to call SV in single cells) 
		

- [Tobias Marschall](https://bioinf.mpi-inf.mpg.de/homepage/index.php?&account=marschal), Leveraging single-cell technology for haplotyping
	- Why is this important? to resolve compound heterozygotes, it's important to phase
	- Bulk: to resolve haplotypes from reads, most common solution is the [minimum error correction](https://en.wikipedia.org/wiki/Error_detection_and_correction) 
		- Fixed-parameter tractable (FPT) algorithms. He et al. (2010), WhatsHap [Patterson et al. 2014](https://www.ncbi.nlm.nih.gov/pubmed/25658651), HapCol [Pirola et al. 2015](https://academic.oup.com/bioinformatics/article/32/11/1610/1742594/H-ap-C-ol-accurate-and-memory-efficient-haplotype)
		- Use switch-error rates to measure distance between true haplotype and predicted haplotype. 
			- Could use Hamming distance, but most people use switch distance (distances in switch space). 
			- Wwitch error rates: pacbio (0.13%), 10x (0.025%), illumina 
		- Use Strand-seq as minimum error correction (MEC) 
			- All known FPT algorithms fail 
		- Integrating Strand-seq and PacBio reads
			- As Strand-seq cells increases, the number of SNVs in largest blocks increases
			- Assess hamming distance using ground truth of a 17-family member pedigree [Eberle et al. Genome Res 2016](http://genome.cshlp.org/content/early/2016/11/25/gr.210500.116)
	- Get accurate and dense chromosome-length haplotyping at low cost for single cells. 
	- immediate anchoring of SVs into haplotypes. 
	- WhatsHap: production-quality tool suite
	

- [Nick Navin](http://faculty.mdanderson.org/Nicholas_Navin/Default.asp?SNID=0), [@nicholas_navin](https://twitter.com/nicholas_navin?lang=en), Single Cells, Big Data
	- Models of tumor evolutions:
		- linear, branching, neutral or punctuated evolution [Davis, Gao & Navin 2016](https://www.ncbi.nlm.nih.gov/pubmed/27526321)
	- Methods to resolve intra-tumor heterogeneity: single cell sequencing
		- deep sequencing, multi-region sequencing, single cell sequencing
		- advantages: Can resolve diverse genomes in complex populations of cells, get genomic info on rare subpopulations, can be used to trace lineages, can discover and define new cell types. 
		- limitations: sample size N, technical errors
		- applications: primary tumors, CTCs & metastatsis, therapy resistance 
		- translational applications of SCS in oncology: non-invasive monitoring, prognostic diversity indexes, improving targeted therapy, early detection, scarce clinical samples
	- Development of scDNA-seq
		- Three steps: (1) single cell isolation (2) genome or transcriptome amplification (3) NGS
		- DNA (6 picograms, 2 molecules), RNA (0.1pigograms mRNA 200K-500kK molecules)
		- single nucleus sequencing (SNS) are less sticky than full cells
			- extensions: MALBAC [Ni et al. 2013](), NUC-SEQ [Wang et al. 2014]() - all developed to increase coverage
				- cost now down to $1-2 / cell
			- SNES: Can now sequencing 25,000 exons of the most frequencing mutated genes in cancer 
	- Computational methods for analyzing scDNA-seq
		- Paradigm shift from single genomes from one patient to sequencing hundreds to thousands of genomes from each patients (or regions from a tumor)
			~4TB targeted panels, ~21 TB (exomes), ~600 TB (genomes)
		- technical errors in single cell DNAseq data
			- allelic dropout (now 10-30%) - only one allele is amplified, false positives (1e5 to 1e6) - homozygous mutation to het mutation after sequencing (can filter across cells to get error rates low), no coverage sites (~5-10%) - false negatives, non-uniform coverage
		- Data: observed data genotype matrix (contains FPs, FNs, missing data) - all need to be taken into account in the data analysis
		- Methods: 
			- Monovar (call DNA variants using genotype likelihood model), 
			- SiFit
			- VariaableBinning - using variable bins reduces noise
			- Integer Copy Number Estimation and Ginko
	- metastatic dissemination in colorectal cancer
		- Data suggests it follows a late-dissemination model in which one or more clones seed the liver metastatsis after acquiring most driver mutations in the tumor
	- punctuated copy number evolution in triple negative breast cancer
	- directions/challenges:
		- connect genotypes and phenotypes: live cell imaging & sc seq
		- spatially resolved single sequencing in tissue sections
		- multi-omics measurements from single cells
		- variant calling can take a long time. 



- [Katharina Jahn](https://www.bsse.ethz.ch/cbg/group/people/person-detail.html?persid=208479), Reconstructing tumor mutation histories from single-cell sequencing data

- [Geoffrey Schiebinger](http://www.stat.berkeley.edu/~geoff/), Towards a statistical theory of ontogenetics

- [James Gagnon](https://www.mcb.harvard.edu/directory/james-gagnon/), [@James_Gagnon](https://twitter.com/james_gagnon), Whole organism lineage tracing with genome editing
	- Lineage tracing with genome editing
	- [GESTALT](https://www.ncbi.nlm.nih.gov/pubmed/27229144): genome editing of synthetic target arrays for lineage tracing 
	- Challenges:
		- spatial context = want to read out barcode via imaging that maintains the spatial context in the tissue
		- complexity. each crispr can be mutated 
		- control when the editing is occurring (early editing only); would like a late round of editing. Control the expression of transgeneic Cas9 allows 
		- combine with single cell measurements
			- 13k cells clustered based on gene expression profiles. with Bushra Raj. 
			- only capture barcodes for 10-20% of cells
			- methods for visualizing, interpret and quantity these relationships? 

## Friday, May 5

- [Oliver Stegle](http://www.ebi.ac.uk/research/stegle), Latent variable models for decomposing single-cell expression variation
	- How to model sources of wanted and unwanted variation? 
	- approaches: 
		- (1) use expression residuals 
		- (2) estimate a random effect covariance with low rank structure (assumes genes are independent given the covariance model)
	- Use known annotated cell cycle gene sets to fit cell-to-cell covariances, regress out their effects and then do downstream analyses
		-  then employ latent variable modeling to reconstruct a cell cycle factor (X)
		- tested using ~300 ES cells collected at different stages of the cell cycle
	- [f-scLVM](http://biorxiv.org/content/biorxiv/early/2016/11/15/087775.full.pdf)
		- Use sparse factor analysis model for scRNA-Seq. Use prior gene sets, explicit modeling of unannotated factors, accounting for noise
		- Modeling factor annotations 
		- ARD (relevance prior on the factors). How important is the factor a priori? L2 structure on the factor asking what's the average variance explained by the factor given the genes it's regulating
		- hurdle noise model [Finak et al. 2015](http://link.springer.com/article/10.1186/s13059-015-0844-5)
		- compared to pagoda (one gene set a time)
	- Applications: 
		- To test, use ~300 ES cells collected at different stages of the cell cycle. factorial scLVM automatically discovers the presence of cell cycle
		- modeling sparse factors jointly, you get independent factors (e.g. checkpoint, TP53, etc)
	- Scalability inference
		- factorized variational Bayes (retaining the coupling of the sparsity indicator and factor weights) (Lazaro-Gredilla et al. 2011)
	- How to disentangle the causes of single-cell transcriptome diversity? 
		- Look at other molecular layers (e.g. DNAm) to expression
	- [DeepCpG](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1189-z): Predicting DNAm states from sequencing and methylation context
		- using DNA-seq to learn motifs that are associated or not with methylation heterogeneity at the single-cell level
		- Need long regions (1Kb), but only 2-3 convolution layers in deep learning
	- Challenges: 
		- How to view data from different sources from same cells. Can again use factor models


- [Laurent Modolo](https://scholar.google.fr/citations?user=LlXhId0AAAAJ&hl=en), Sparse Gamma/Poisson PCA to unravel the genomic diversity of single-cell expression data
	- ERCC Spike-ins correlated with CDR
	- count matrix factorization (CMF). X_ij = over dispersed, zero inflated count data. replace norm by approximate likelihood based approaches. 
	- X_ij with poisson rate. factorize E[X] = Gamma = UV^T. 
	- sparsity on the loadings. 
	- Zero inflation on gamma-poisson factor model (poisson-dirac micture)
	- suitable for any count data, especially for NGS. Accounts for over-dispersion (gamma-poisson model), zero-inflation (poisson-diract mixture) and sparsity in V (gamma-dirac mixture). 
	- framework of variational inferences
	- Efficient implementation in C++, incorporated in a R package CMF
	- clustering of the observations accoring to the matrix U (compares Poisson-NMF, PCA, CMF)

- [Magnus Rattray](https://www.research.manchester.ac.uk/portal/Magnus.Rattray.html), [@magnusrattray](https://twitter.com/magnusrattray), Inferring early bifurcation events from single-cell RNA-Seq data
	- Gaussian process recognition
		- flexible non-parametric models. probability distributions over functions. 
		- covariance function (not matrix) has parameters tuning these properties
		- [GPflow package](https://arxiv.org/pdf/1610.08733.pdf) use TensorFlow autodiff (b/c computing derivatives is time-consuming)
	- speeding up inference: GPLVM pseudo time example
		- BEAM (uses splines) more biased towards global branch time
	- allows inference of uncertainty in branching times
	- can identify genes branching earlier than global branching	 
	- challenges: joint inference of branching and pseudo time
	