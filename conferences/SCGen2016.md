# Notes on talks for the 2016 SCG Conference

Notes and slides for the [2016 Single Cell Genomics (SCG) Conference](https://coursesandconferences.wellcomegenomecampus.org/events/item.aspx?e=596) 
at the Wellcome Genome Campus, Hinxton, Cambridge, UK from Sept 14-16, 2016. Here is a pdf of the 
[SCG Program](http://conf.hinxton.wellcome.ac.uk/advancedcourses/SSG16draftprogrammeFinal.pdf). 
Follow the twitter hashtag [#SCGen16](https://twitter.com/search?f=tweets&vertical=default&q=%23SCGen16&src=typd). 
I missed some of the talks, so pull requests are welcome to fill in gaps! or tweet me
[@stephaniehicks](https://twitter.com/stephaniehicks). 


## Wednesday, Sept 14 

### Keynote Lecture

- [Arnold Kriegstein](https://bms.ucsf.edu/directory/faculty/arnold-kriegstein-md-phd), Genomic insights into human cortical development, lissencephaly, and Zika microcephaly
	- Very interesting story on using scRNA-seq to study the zika virus and congenital microcephaly. Both cause malformations in the brain (but Zika also destroys already developed tissue). The AXL Zika receptor is expressed in radial glia cells (precursors for neuronal expansion in brain cortex) and can be used to predict where damage occurs in the brain. Showed blocking AXL prevents Zika entry into progenitor cells in developing brain. Pregnancy-safe compounds now being tested. 

### Session 1: Neuroscience & Tissue Development

**Chair**: [Rickard Sandberg](http://sandberg.cmb.ki.se) from [@karolinskainst](https://twitter.com/karolinskainst)

- [Sten Linnarsson](http://linnarssonlab.org), [@slinnarsson](https://twitter.com/slinnarsson), Towards a census of mouse brain cell types
	- Example of [identifying neurons that control goosebumps](http://linnarssonlab.org/publications/2016/08/29/sympathetic/). On brain development, "differentiation speed affects how many intermediate types scRNAseq captures". Useful links to data and software available on lab website including [BackSPIN](https://github.com/linnarsson-lab/BackSPIN): biclustering algorithm based on sorting points into neighborhoods (SPIN), implemented in MATLAB and Python and described in [Zeisel et al. (2015)](http://science.sciencemag.org/content/347/6226/1138). 

- [Barbara Treutlein](http://www.treutleinlab.org), Reconstructing human organogenesis using single-cell RNA-seq
	- [Treutlein et al. 2016](http://www.nature.com/nature/journal/v534/n7607/full/nature18323.html) - reprogramming from fibroblast to neuron using scRNA-seq	

- [Maria Kasper](http://ki.se/en/people/markas), Plasticity and heterogeneity of skin cells in health and tissue repair
	- There are 25 (!!) different cell types in a single hair follicle

- [Naomi Habib](http://zlab.mit.edu/team.html), Single nucleus RNA-Seq reveals dynamics of adult neurogenesis
	- Developed [sNuc-Seq and Div-Seq](http://science.sciencemag.org/content/353/6302/925.full) to investigate dynamic transcriptome of rare adult newborn neurons; introduced bi-SNE (biclustering on Stochastic Neighbor Embedding)
		- sNuc-Seq (single nuclei RNA-seq): skips dissociation step and can be done on aging tissue, and fresh/frozen tissue, allows you to obtain neurons. 	
		- Div-Seq (combines sNuc-Seq with pulse labeling of proliferating cells): labels all dividing cells in dynamic processes, no need for marker genes.  

- [Gray Camp](http://www.eva.mpg.de/genetics/staff.html), Human cerebral organoids recapitulate gene expression programs of fetal neocortex development
	- [used scRNA-seq to deconstruct fetal human neocortex & cerebral organoid development](http://www.pnas.org/content/112/51/15672.full). Used [SCDE](http://hms-dbmi.github.io/scde/) identify DE genes between chimp and human (candidates for human-specific function).

- [Marta Rodriguez Orejuela](https://www.mdc-berlin.de/10179514/en/research/research_teams/systems_biology_of_gene_regulatory_elements/team), Characterization of adult neurogenic niches at single cell resolution



## Thursday, Sept 15

### Session 2: Chromatin Structure and Organization 

**Chair**: [Ido Amit](https://www.weizmann.ac.il/immunology/AmitLab/front) from [@WeizmannScience](https://twitter.com/WeizmannScience)

- [Amos Tanay](http://compgenomics.weizmann.ac.il/tanay/), Single cell dynamics of clonal memory
	- [MARS-Seq](http://science.sciencemag.org/content/343/6172/776.abstract) - Based on CEL-Seq technology, but automates processing of cells into 384-well plates and incorporates index sorting up to 10 markers (especially useful for immune cells)

- [Will Greenleaf](http://greenleaf.stanford.edu/index.html), [@WJGreenleaf](https://twitter.com/wjgreenleaf), ATAC-ing regulatory variation in single cells
	- Starts out with both [Waddington's](https://en.wikipedia.org/wiki/C._H._Waddington) Classical Epigenetic Landscape images - [concept of an epigenetic landscape as a visual metaphor](http://www.cell.com/cell/pdf/S0092-8674(07)00186-9.pdf) for the cell (represented by the ball), which can take specific trajectories, leading to different outcomes or cell fates
	- discussed [ATAC-seq](http://www.nature.com/nmeth/journal/v10/n12/full/nmeth.2688.html) (captures open chromatin sites e.g. nucleosomes, TF binding footprints, chromatin states) and [scATAC-Seq](http://www.nature.com/nature/journal/v523/n7561/fig_tab/nature14590_F1.html) using Fluidigm platform

- [Stephen Clark](http://www.babraham.ac.uk/our-research/lymphocyte/geoffrey-butcher/members/65/stephen-clark), Chromatin accessibility, DNA methylation and gene expression from the same single-cell
	- Discussed how chromatin accessibility within genes correlates with expression. Used [scM&T-Seq](http://www.nature.com/nmeth/journal/v13/n3/fig_tab/nmeth.3728_SF1.html) (combines G&T-seq, Smart-Seq2 + scBS-seq) to sequence the methylome and transcriptome in single-cells. 
	
- [Ana Pombo](https://pombolab.wordpress.com), [@apombo1](https://twitter.com/apombo1), Genome Architecture Mapping, new approach to map chromatin contacts
	- Genome Architecture Mapping (GAM): approach to measuring 3-D chromatin topology in a nucleus (spatial information). Maps chromatin contacts using random cryosectioning (slices through a nucleus), extract DNA from sections to sequence, to quantify the frequency of locus co-segregation (currently 30-40 kb resolution) and calculate individual nuclear profiles (NPs). 
	
- [Peter Fraser](http://www.babraham.ac.uk/our-research/nuclear-dynamics/peter-fraser), [@Peter_Fraser1](https://twitter.com/peter_fraser1), Chromosome dynamics revealed by single cell HiC
	- Really cool example of using [single-cell Hi-C](http://www.nature.com/nature/journal/v502/n7469/full/nature12593.html) to investigate chromosomal conformations through stages of the cell cycle. Suggests that TADs may more reflect replication dynamics than transcription regulation. 

- [Amanda Ackerman](http://www.chop.edu/doctors/ackermann-amanda#.V9F6JmX_RGI), Single-cell ATAC-seq identifies epigenetic differences in human pancreatic islet cell subtypes from normal and diabetic donors
	- studies differences between [alpha and beta islet cells](http://www.molmetab.com/article/S2212-8778(16)00003-X/fulltext) using sc and bulk RNA-seq, ChIP-seq, BS-seq and ATAC-seq

- [Jan-Philipp Mallm](https://malone.bioquant.uni-heidelberg.de/people/mallm/index-mallm.html), Dissecting Deregulated Enhancer Activity in Primary Leukemia Cells


### Session 3: Immunology and Cancer

**Chair**: [John Marioni](http://www.ebi.ac.uk/research/marioni) from [@emblebi](https://twitter.com/emblebi) and [@CRUKresearch](https://twitter.com/crukresearch)

- [Sarah Teichmann](http://www.teichlab.org), Understanding Cellular Heterogeneity
	- Discusses [sensitivity, specificity and accuracy of different scRNA-seq protocols](http://biorxiv.org/content/early/2016/09/08/073692) using ERCC spike-ins (work with [@vallens](https://twitter.com/vallens)). Found endogenous genes are more **efficiently** captured than ERCC spike-ins (counterintuitive, but pleasantly surprising). Freeze-thaw cycles decrease RNA content about 20% with each cycle. Protocols are overall very **accurate** (Pearson correlation of expected vs observed), but there are differences in **sensitivity** (especially for detecting lowly expressed genes ~1-10 molecules per cell) - you get a benefit of sequencing up to 1 million reads. 
	- Uses [Gaussian process latent variable model](http://www.cell.com/cell-reports/fulltext/S2211-1247(15)01538-7) for dim reduction
	- Software include [GPfates](https://github.com/Teichlab/GPfates) (models transcriptional cell fates as mixtures of Gaussian Processes) and [TraCeR](http://www.nature.com/nmeth/journal/v13/n4/full/nmeth.3800.html) (reconstructs T-cell clonal relationships from scRNA-seq data	
	- Advertises for [Single Cell Omics Keystone Symposia](http://keystonesymposia.org/17e3) in Stockholm May 26-30, 2017

- [Ido Amit](https://www.weizmann.ac.il/immunology/AmitLab/front), Immunology in the age of single cell genomics
	- Combined CRISPR with [MARS-seq2](https://compgenomics.weizmann.ac.il/tanay/?page_id=672) (~10,000 cells/day, better barcoding, SNR, reproducibility and cost) for pooled index screening - based on gRNA libraries bearing RNA barcodes

- [Timm Schroeder](https://www.bsse.ethz.ch/department/people/detail-person.html?persid=193443), Long-term single cell quantification: New tools for old questions
	- Argued that snapshots of single cells is not sufficient. Need time dimension to really understand what's happening biologically. Critical of pseudotime approaches (not as good as real time in some situations e.g. identifying cycles). Reaffirms the need for good methods software to analyze single cell data
	- [The Tracking Tool]((http://www.nature.com/nbt/journal/v34/n7/full/nbt.3626.html)) - for continuous single cell behavior quantification 

- [Fabian Theis](http://fabian.theis.name), Diffusion pseudotime identifies lineage choice and graded transitions in myeloid progenitors
	- Extracts pseudotemporal ordering and branch points from diffusion maps (examples from early blood development)
	- [Diffusion Pseudotime (DPT)](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3971.html?WT.feed_name=subjects_genetics) software - to estimate order of differentiating cells according to DPT (measures transitions between cells using diffusion-like random walks)

- [Bart Deplancke](http://deplanckelab.epfl.ch), [@BartDeplancke](https://twitter.com/bartdeplancke), Single-cell RNA-seq-based identification and characterisation of somatic stem cells in adipose tissue & beyond
	- [Automated Single-cell Analysis Pipeline (ASAP)](http://asap.epfl.ch)

- [Matan Hofree](https://www.researchgate.net/profile/Matan_Hofree), Unbiased whole tissue analysis of the single cell transcriptional landscape of colon cancer

- [Amir Giladi](http://www.weizmann.ac.il/lifesci/idcards/AmirGiladi0464.html), Transcriptional heterogeneity and lineage commitment in hematopoietic progenitors



## Friday, Sept 16

### Sesssion 4: Transcriptomics

**Chair**: [Sten Linnarsson](http://linnarssonlab.org), [@slinnarsson](https://twitter.com/slinnarsson)

- [Rickard Sandberg](http://sandberg.cmb.ki.se), Single-cell gene expression analyses of allelic transcription and regulation
	- Discussed chromosome-wide evidence of de novo paternal X chromosome inactivation in mouse embryos using scRNA-seq

- [Allon Klein](http://klein.hms.harvard.edu), Population balance reconstruction of differentiation hierarchies in developing and adult tissues by single cell droplet RNA-Seq
	- Discussed [inDrop](http://www.cell.com/cell/abstract/S0092-8674(15)00500-0) (core facility now at 15-20K cells/hr, 80% of cells barcoded, >2K cell input, 7 cents/cell)
	- SPRING - software for visualizing and interacting with single cell data using a force-directed graph approach; publication coming soon. 

- [Arnau Sebé-Pedrós](https://compgenomics.weizmann.ac.il/tanay/?page_id=12), [@ArnauSebe](https://twitter.com/arnausebe), Early metazoan cell type evolution by single cell RNA-seq analysis

- [Omid Faridani](https://www.researchgate.net/profile/Omid_Faridani2), Sequencing Small-RNA transcriptome of individual cells

- [Alexander van Oudenaarden](http://www.hubrecht.eu/onderzoekers/van-oudenaarden-group/), Revealing novel cell types, cell-cell interactions, and cell lineages by single-cell sequencing
	- Used [single-cell strand-specific 5hmC sequencing to explore cell-to-cell variability and cell lineage reconstruction](http://www.nature.com/nbt/journal/v34/n8/full/nbt.3598.html)

- [Stephanie Hicks](http://www.stephaniehicks.com), [@stephaniehicks](https://twitter.com/stephaniehicks), Towards progress in batch effects and biases single-cell RNA-Seq data [[Slides](https://speakerdeck.com/stephaniehicks/towards-progress-in-batch-effects-and-biases-in-single-cell-rna-seq-data)]
	- Pre-print of our [paper on bioRxiv](http://biorxiv.org/content/early/2015/12/27/025528)

- [Marc Wadsworth](https://www.researchgate.net/profile/Marc_Wadsworth), Seq-Well: A Portable Single-Cell RNA-Seq Platform for Low-Input Clinical Samples
	-  Seq-Well: the portable microarray well technology for scRNA-seq (combines microwell and bead-based approaches in a nanowell) - 10K cells/array

- [Eshita Sharma](https://scholar.google.com/citations?user=xljyFDkAAAAJ&hl=en), Single cell preservation for RNAseq

- [Marc Lynch](https://www.fluidigm.com/about/aboutfluidigm), Single-cell transcriptomics and functional analysis of single- cells

- [Manuel Garber](http://garberlab.umassmed.edu), Dissection of T1 Diabetes progression using Single cell RNA sequencing of a RAT model
	- [End sequencing analysis toolkit (ESAT)](https://github.com/garber-lab/ESAT) 


### Sesssion 5: Imaging and Modeling

**Chair**: [Alexander van Oudenaarden](http://www.hubrecht.eu/onderzoekers/van-oudenaarden-group/) from [@_Hubrecht](https://twitter.com/_hubrecht?lang=en)

- [Long Cai](http://singlecell.caltech.edu/cailab/), In situ transcription profiling in tissues by seqFISH
	- [SeqFISH](http://www.nature.com/nmeth/journal/v11/n4/full/nmeth.2892.html) ? multiplexed spatial transcriptomics
	- Clustered 15,000 imaged cells with 125 genes to identify many cell types (mapped back onto slice images)

- [Heather Lee](http://www.babraham.ac.uk/our-research/epigenetics/olivia-casanueva/members/198/heather-lee), Dynamic and heterogeneous DNA methylation in pluripotent cells
	- uses [scBS-seq](http://go.nature.com/2ccGomv) to measure methylation

- [Steffen Rulands](http://www.rulands.net), [@srulands](https://twitter.com/srulands), Dynamic and heterogeneous DNA methylation in pluripotent cells

- [Jeffrey Moffitt](https://scholar.google.com/citations?user=U7eic7AAAAAJ&hl=en), High-throughput, spatially resolved, single-cell transcriptomics with MERFISH
	- [MERFISH](http://www.pnas.org/content/early/2016/09/07/1612826113): an image-based single cell transcriptomics method; now increased throughput for imaging many cells 

- [Jan Philipp Junker](https://scholar.google.com/citations?user=0tt8A_4AAAAJ), Massively parallel clonal analysis using CRISPR/Cas9 induced genetic scars

- [John Marioni](http://www.ebi.ac.uk/research/marioni), Dissecting cell fate choice using single-cell genomics
	- [BASiCs](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004333) (Bayesian Analysis of Single-Cell Sequencing data): finds differentially expressed and differentially variable genes
	- Discussed how bulk RNA-seq normalization tools (e.g. DESeq, TMM) do not work well for scRNA-seq; points to new normalization method in [Bioconductor called scran](http://bioconductor.org/packages/release/bioc/html/scran.html) which calculates cell-specific normalization factors by pooling information across cells to avoid problems due to sparsity in scRNA-seq data
	- Discussed problems with spike-ins in practice; points to [Synthetic Spike-ins Controls (sequins)](http://www.nature.com/nmeth/journal/v13/n9/full/nmeth.3958.html)

- [David van Dijk](https://sciencedavid.wordpress.com), MAGIC: A Diffusion based data imputation method reveals progressions and gene-gene interactions in breast cancer cells undergoing EMT
	- MAGIC: uses random walks/diffusion to generate low-dimension global manifold (data structure) based on local cell similarities; validates method by showing it recovers structure even after inducing a lot of sparsity 

- [Rom Shenhav](http://shalevlab.weizmann.ac.il/group-members/), Single-cell spatial reconstruction reveals global division of labor in the mammalian liver

- [Petra Schwalie](https://scholar.google.com/citations?user=EMSKH8cAAAAJ&hl=en), Accurate identification of somatic stem cells using single-cell RNA-sequencing
	- Used LASSO logistic regression to to identify stem cells in adult mammalian tissue across 14 scRNA-seq adult somatic data sets (train (2/3), tested (1/3))
	


