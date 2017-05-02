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
		- Normalization: [Lun et al. (2016)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0947-7) implemented in [scran R/Bioconductor pkg](http://bioconductor.org/packages/release/bioc/html/scran.html), [Bacher et al. (2017)](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4263.html?WT.feed_name=subjects_systems-biology) implemented in [SCnorm](https://github.com/rhondabacher/SCnorm), Vallejos et al (2017)
		- Confounding factors: Buettner et al. (2015), Finak et al. (2015), Hicks et al. bioRxiv
		- Spatial expression: Achim et al. (2015), Satija et al. (2015)
		- Clustering and pseudotime: Trapnell et al. (2014), Haghvaerdi et al. (2016), Kiselev et al. (2017)
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
		

- [Peter Kharchenko](http://pklab.med.harvard.edu), From one to millions of cells: computational approaches for single--cell analysis
	- Discussing rapid growth in cell numbers, data volume, more multi-patient samples
	- Discussing variability in scRNA-seq data
		- e.g. differences between cells (of the same type), overdispersion, measurement failures, quality of cells, 
		- can be problems for PCA (non-gaussian, noisy cells from outliers, can makes true heterogeneity)
		- Biological and technical variablity (cell-specific noise models)
	- Discussed how to model noise of an individual cell (dropout models)
	- Used known gene sets to get better power to detect patterns of variability 
	- PAGODA: weight PCA - weights observations that the probability that it's a technical dropout
		- pick up correlated gene sets from de novo gene sets and annotated gene sets
		- "Aspects" = clusters of PCs. 
	- PAGODA2 - optimized for large sparse measurements
	- Lanczos algorithm. Must be able to break down computations to matrix multiplication. 
	- Challenges: 
		- controlling for sequencing depth. correlates well with cell size. You may kill biological signal if you simply regress it out.
		- sparse measurements. how to process on sparse measurements: nearest neighbors. 
		- visualization on sparse matricies. 
		- making use of reference data. (e.g. organism, tissue-specific data) 
		- distances. uses mostly correlation based distance. The problem is there is no universal distance estimates. Depends on question e.g. if you are looking at cell cycle genes, you should define distance based on cell cycle genes). 

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
	- What markers? **microRNAs** are most informative b/c they have high discrimination power (Golub group, Nature 2005)
	- proof of concept experiment: Built a logic gene circuit to identify HeLa cells using 6 discriminatng miRNAs (identified from literature)
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
	- 
	
- [Christof Seiler](), CytoGLMM: Bayesian hierarchical linear modeling for flow cytometry data

- [Eirini Arvaniti](), Sensitive detection of rare disease-associated cell subsets via representation learning

- [Raphael Gottardo](), Identifying novel correlates of protection via single-cell analysis of antigen-specific T-cells

- [Davis McCarthy](), Mapping genetic effects on interactions with single-cell states 
		

## Wednesday, May 3

- [Prisca Liberali](), Mapping genetic interactions during intestinal organoid self- organization

- [Takashi Hiragi](), Symmetry breaking and self-organisation in mouse development

- [Andre Rendeiro](), Pooled CRISPR screening with single-cell transcriptome readout

- [Lars Velten](), Quantifying developmental plasticity by integrated single-cell RNA-seq and ex vivo culture


## Thursday, May 4

- [Valérie Taly](), 

- [Nicola Aceto](), 

- [Jan Korbel](), 

- [Tobias Marschall](), 

- [Nick Navin](), 

- [Katharina Jahn](), 

- [Geoffrey Schiebinger](), 

- [James Gagnon](), 

## Friday, May 5

- [Oliver Stegle](), Latent variable models for decomposing single-cell expression variation

- [Laurent Modolo](), Sparse Gamma/Poisson PCA to unravel the genomic diversity of single-cell expression data

- [Magnus Rattray](), Inferring early bifurcation events from single-cell RNA-Seq data

