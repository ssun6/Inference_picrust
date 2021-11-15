## Manuscript

**Title:**                
_Inference-based accuracy of metagenome prediction tools varies across sample types and functional categories_

**Authors:**               
Shan Sun, Roshonda B Jones, Anthony A Fodor

## Abstract
Background
Despite recent decreases in the cost of sequencing, shotgun metagenome sequencing remains more expensive compared with 16S rRNA amplicon sequencing. Methods have been developed to predict the functional profiles of microbial communities based on their taxonomic composition. In this study, we evaluated the performance of three commonly used metagenome prediction tools (PICRUSt, PICRUSt2, and Tax4Fun) by comparing the significance of the differential abundance of predicted functional gene profiles to those from shotgun metagenome sequencing across different environments.

Results
We selected 7 datasets of human, non-human animal, and environmental (soil) samples that have publicly available 16S rRNA and shotgun metagenome sequences. As we would expect based on previous literature, strong Spearman correlations were observed between predicted gene compositions and gene relative abundance measured with shotgun metagenome sequencing. However, these strong correlations were preserved even when the abundance of genes were permuted across samples. This suggests that simple correlation coefficient is a highly unreliable measure for the performance of metagenome prediction tools. As an alternative, we compared the performance of genes predicted with PICRUSt, PICRUSt2, and Tax4Fun to sequenced metagenome genes in inference models associated with metadata within each dataset. With this approach, we found reasonable performance for human datasets, with the metagenome prediction tools performing better for inference on genes related to “housekeeping” functions. However, their performance degraded sharply outside of human datasets when used for inference.

Conclusion
We conclude that the utility of PICRUSt, PICRUSt2, and Tax4Fun for inference with the default database is likely limited outside of human samples and that development of tools for gene prediction specific to different non-human and environmental samples is warranted.
