# Notes on talks for the 2017 Ascona Workshop on Statistical Challenges in Single-Cell Biology

#### Organizing committee

* Niko Beerenwinkel
* Peter Bühlmann
* Wolfgang Huber

Notes and slides for the [2017 Ascona Workshop](https://www.bsse.ethz.ch/cbg/cbg-news/ascona-2017.html) 
at Monte Verita, Ascona, Switzerland from Apr 30 - May 5, 2017. Here is a pdf of the 
[program](https://www.ethz.ch/content/dam/ethz/special-interest/bsse/cbg/events/ascona-2017/2017-04-28%20program.pdf). 


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
		

- [Peter Kharchenko](http://pklab.med.harvard.edu)
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

- [Ambrose Carr](http://ambrosejcarr.com), [@ambrosejcarr](https://twitter.com/ambrosejcarr)
	- Profile tumor immune cells in breast cancer. Used inDrop platform. 
	- cell-type specific capture rates, dropouts
	- Discussed [BISCUIT](https://genomicscomputbiol.org/ojs/index.php/GCB/article/view/46) which is used to iteratively (1) impute dropouts based on expression or co-expression with other genes and (2) normalize data across library size. 
		- global normalization fails b/c cells have different sizes. Dropouts are not resolved and may introduce spurious correlations. 
	- For each gene, algorithm models clusters of cells using bayesian mixture model (mixtures across mixtures of cell types) 
	- Discussed why downsampling is not appropriate (because you throw out too much data e.g. 90%)

- [Stephanie Hicks](http://www.stephaniehicks.com), [@stephaniehicks](https://twitter.com/stephaniehicks), Missing data and technical variability in single-cell RNA-sequencing experiments

- [Yuanhua Huang]()

- [Jean (Zhijin) Wu]() 

- [Felix Kurt]()



## Tuesday, May 2

## Wednesday, May 3

## Thursday, May 4

## Friday, May 5

