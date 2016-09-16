# Notes on talks for the 2016 SCG Conference

Notes and slides for the [2016 Single Cell Genomics (SCG) Conference](https://coursesandconferences.wellcomegenomecampus.org/events/item.aspx?e=596) 
at the Wellcome Genome Campus, Hinxton, Cambridge, UK from Sept 14-16, 2016. Here is a pdf of the 
[SCG Program](http://conf.hinxton.wellcome.ac.uk/advancedcourses/SSG16draftprogrammeFinal.pdf). 
Follow the twitter hashtag [#SCGen16](https://twitter.com/search?f=tweets&vertical=default&q=%23SCGen16&src=typd)

Pull requests welcome! or tweet me
[@stephaniehicks](https://twitter.com/stephaniehicks). 


## Wednesday, Sept 14 

### Keynote Lecture

- [Arnold Kriegstein](https://bms.ucsf.edu/directory/faculty/arnold-kriegstein-md-phd), Genomic insights into human cortical development, lissencephaly, and Zika microcephaly

### Session 1: Neuroscience & Tissue Development

**Chair**: [Rickard Sandberg](http://sandberg.cmb.ki.se) from [@karolinskainst](https://twitter.com/karolinskainst)

- [Sten Linnarsson](http://linnarssonlab.org), [@slinnarsson](https://twitter.com/slinnarsson), Towards a census of mouse brain cell types
	- Useful links to data and software available on lab website including [BackSPIN](https://github.com/linnarsson-lab/BackSPIN): biclustering algorithm based on sorting points into neighborhoods (SPIN), implemented in MATLAB and Python and described in [Zeisel et al. (2015)](http://science.sciencemag.org/content/347/6226/1138)

- [Barbara Treutlein](http://www.treutleinlab.org), Reconstructing human organogenesis using single-cell RNA- seq

- [Maria Kasper](http://ki.se/en/people/markas), Plasticity and heterogeneity of skin cells in health and tissue repair

- [Naomi Habib](http://zlab.mit.edu/team.html), Single nucleus RNA-Seq reveals dynamics of adult neurogenesis
	- Developed [sNuc-Seq and Div-Seq](http://science.sciencemag.org/content/353/6302/925.full) to investigate dynamic transcriptome of rare adult newborn neurons; introduced bi-SNE (biclustering on Stochastic Neighbor Embedding)
		- sNuc-Seq (single nuclei RNA-seq): skips dissociation step and can be done on aging tissue, and fresh/frozen tissue, allows you to obtain neurons. 	
		- Div-Seq (combines sNuc-Seq with pulse labeling of proliferating cells): labels all dividing cells in dynamic processes, no need for marker genes.  

- [Gray Camp](http://www.eva.mpg.de/genetics/staff.html), Human cerebral organoids recapitulate gene expression programs of fetal neocortex development

- [Marta Rodriguez Orejuela](https://www.mdc-berlin.de/10179514/en/research/research_teams/systems_biology_of_gene_regulatory_elements/team), Characterization of adult neurogenic niches at single cell resolution



## Thursday, Sept 15

### Session 2: Chromatin Structure and Organization 

**Chair**: [Ido Amit](https://www.weizmann.ac.il/immunology/AmitLab/front) from [@WeizmannScience](https://twitter.com/WeizmannScience)

- [Amos Tanay](http://compgenomics.weizmann.ac.il/tanay/), Single cell dynamics of clonal memory
	- [MARS-Seq](http://science.sciencemag.org/content/343/6172/776.abstract) - Based on CEL-Seq technology, but automates processing of cells into 384-well plates and incorporates index sorting up to 10 markers (especially useful for immune cells)

- [Will Greenleaf](http://greenleaf.stanford.edu/index.html), [@WJGreenleaf](https://twitter.com/wjgreenleaf), ATAC-ing regulatory variation in single cells
	- [scATAC-Seq](http://www.nature.com/nature/journal/v523/n7561/fig_tab/nature14590_F1.html) using Fluidigm platform

- [Stephen Clark](http://www.babraham.ac.uk/our-research/lymphocyte/geoffrey-butcher/members/65/stephen-clark), Chromatin accessibility, DNA methylation and gene expression from the same single-cell
	- Discussed how chromatin accessibility within genes correlates with expression. Used [scM&T-Seq](http://www.nature.com/nmeth/journal/v13/n3/fig_tab/nmeth.3728_SF1.html) to sequence the methylome and transcriptome in single-cells. 
	
- [Ana Pombo](https://pombolab.wordpress.com), [@apombo1](https://twitter.com/apombo1), Genome Architecture Mapping, new approach to map chromatin contacts
	- Genome Architecture Mapping (GAM): approach to measuring 3-D chromatin topology in a nucleus (spatial information). Maps chromatin contacts using random cryosectioning (slices through a nucleus), extract DNA from sections to sequence, to quantify the frequency of locus co-segregation and calculate individual nuclear profiles (NPs). 
	
- [Peter Fraser](http://www.babraham.ac.uk/our-research/nuclear-dynamics/peter-fraser), [@Peter_Fraser1](https://twitter.com/peter_fraser1), Chromosome dynamics revealed by single cell HiC
	- Really cool example of using [single-cell Hi-C](http://www.nature.com/nature/journal/v502/n7469/full/nature12593.html) to investigate chromosomal conformations through stages of the cell cycle. Suggests that TADs may more reflect replication dynamics than transcription regulation. 

- [Amanda Ackerman](http://www.chop.edu/doctors/ackermann-amanda#.V9F6JmX_RGI), Single-cell ATAC-seq identifies epigenetic differences in human pancreatic islet cell subtypes from normal and diabetic donors

- [Jan-Philipp Mallm](https://malone.bioquant.uni-heidelberg.de/people/mallm/index-mallm.html), Dissecting Deregulated Enhancer Activity in Primary Leukemia Cells



### Session 3: Immunology and Cancer

**Chair**: [John Marioni](http://www.ebi.ac.uk/research/marioni) from [@emblebi](https://twitter.com/emblebi) and [@CRUKresearch](https://twitter.com/crukresearch)

- [Sarah Teichmann](http://www.teichlab.org), Understanding Cellular Heterogeneity
	- Discusses [sensitivity, specificity and accuracy of different scRNA-seq protocols](http://biorxiv.org/content/early/2016/09/08/073692) using ERCC spike-ins (work with [@vallens](https://twitter.com/vallens)). Found endogenous genes are more **efficiently** captured than ERCC spike-ins (counterintuitive, but pleasantly surprising). Freeze-thaw cycles decrease RNA content about 20% with each cycle. Protocols are overall very **accurate** (Pearson correlation of expected vs observed), but there are differences in **sensitivity** (especially for detecting lowly expressed genes ~1-10 molecules per cell) - you get a benefit of sequencing up to 1 million reads. 
	- [GPfates](https://github.com/Teichlab/GPfates) - Model transcriptional cell fates as mixtures of Gaussian Processes

- [Ido Amit](https://www.weizmann.ac.il/immunology/AmitLab/front), Immunology in the age of single cell genomics
	- Combined CRISPR with MARS-seq2 for pooled index screening - based on gRNA libraries bearing RNA barcodes

- [Timm Schroeder](https://www.bsse.ethz.ch/department/people/detail-person.html?persid=193443), Long-term single cell quantification: New tools for old questions

- [Amir Giladi](http://www.weizmann.ac.il/lifesci/idcards/AmirGiladi0464.html), Transcriptional heterogeneity and lineage commitment in hematopoietic progenitors

- [Bart Deplancke](http://deplanckelab.epfl.ch), [@BartDeplancke](https://twitter.com/bartdeplancke), Single-cell RNA-seq-based identification and characterisation of somatic stem cells in adipose tissue & beyond
	- [Automated Single-cell Analysis Pipeline (ASAP)](http://asap.epfl.ch)

- [Matan Hofree](https://www.researchgate.net/profile/Matan_Hofree), Unbiased whole tissue analysis of the single cell transcriptional landscape of colon cancer

- [Fabian Theis](http://fabian.theis.name), Diffusion pseudotime identifies lineage choice and graded transitions in myeloid progenitors



## Friday, Sept 16

### Sesssion 4: Transcriptomics

**Chair**: [Sten Linnarsson](http://linnarssonlab.org), [@slinnarsson](https://twitter.com/slinnarsson)

- [Rickard Sandberg](http://sandberg.cmb.ki.se), Single-cell gene expression analyses of allelic transcription and regulation

- [Allon Klein](http://klein.hms.harvard.edu), Population balance reconstruction of differentiation hierarchies in developing and adult tissues by single cell droplet RNA-Seq

- [Arnau Sebé-Pedrós](https://compgenomics.weizmann.ac.il/tanay/?page_id=12), [@ArnauSebe](https://twitter.com/arnausebe), Early metazoan cell type evolution by single cell RNA-seq analysis

- [Omid Faridani](https://www.researchgate.net/profile/Omid_Faridani2), Sequencing Small-RNA transcriptome of individual cells

- [Alexander van Oudenaarden](http://www.hubrecht.eu/onderzoekers/van-oudenaarden-group/), Revealing novel cell types, cell-cell interactions, and cell lineages by single-cell sequencing

- [Stephanie Hicks](http://www.stephaniehicks.com), [@stephaniehicks](https://twitter.com/stephaniehicks), [Towards progress in batch effects and biases single-cell RNA-Seq data](https://speakerdeck.com/stephaniehicks/towards-progress-in-batch-effects-and-biases-in-single-cell-rna-seq-data)

- [Marc Wadsworth](https://www.researchgate.net/profile/Marc_Wadsworth), Seq-Well: A Portable Single-Cell RNA-Seq Platform for Low-Input Clinical Samples
	-  Seq-Well: the portable microarray well technology for scRNA-seq

- [Eshita Sharma](https://scholar.google.com/citations?user=xljyFDkAAAAJ&hl=en), Single cell preservation for RNAseq

- [Marc Unger](https://www.fluidigm.com/about/aboutfluidigm), Single-cell transcriptomics and functional analysis of single- cells

- [Manuel Garber](http://garberlab.umassmed.edu), Dissection of T1 Diabetes progression using Single cell RNA sequencing of a RAT model
	- [End sequencing analysis toolkit (ESAT)](https://github.com/garber-lab/ESAT) 


### Sesssion 5: Imaging and Modeling

**Chair**: [Alexander van Oudenaarden](http://www.hubrecht.eu/onderzoekers/van-oudenaarden-group/) from [@_Hubrecht](https://twitter.com/_hubrecht?lang=en)

- [Long Cai](http://singlecell.caltech.edu/cailab/), In situ transcription profiling in tissues by seqFISH

- [Heather Lee](http://www.babraham.ac.uk/our-research/epigenetics/olivia-casanueva/members/198/heather-lee), Dynamic and heterogeneous DNA methylation in pluripotent cells

- [Steffen Rulands](http://www.rulands.net), [@srulands](https://twitter.com/srulands), Dynamic and heterogeneous DNA methylation in pluripotent cells

- [Jeffrey Moffitt](https://scholar.google.com/citations?user=U7eic7AAAAAJ&hl=en), High-throughput, spatially resolved, single-cell transcriptomics with MERFISH

- [Jan Philipp Junker](https://scholar.google.com/citations?user=0tt8A_4AAAAJ), Massively parallel clonal analysis using CRISPR/Cas9 induced genetic scars

- [John Marioni](http://www.ebi.ac.uk/research/marioni), Dissecting cell fate choice using single-cell genomics

- [David van Dijk](https://sciencedavid.wordpress.com), MAGIC: A Diffusion based data imputation method reveals progressions and gene-gene interactions in breast cancer cells undergoing EMT

- [Rom Shenhav](http://shalevlab.weizmann.ac.il/group-members/), Single-cell spatial reconstruction reveals global division of labor in the mammalian liver

- [Petra Schwalie](https://scholar.google.com/citations?user=EMSKH8cAAAAJ&hl=en), Accurate identification of somatic stem cells using single-cell RNA-sequencing




