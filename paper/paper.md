---
title: 'zalpha: an R package for the identification of regions of the genome under selection'
tags:
- R
- genomics
- genetics
- population genetics
- evolution
- selection
date: "8 August 2020"
output:
  pdf_document: default
  html_document:
    df_print: paged
  word_document: default
authors:
- name: Clare Horscroft
  orcid: 0000-0001-5679-5912
  affiliation: 1, 3
- name: Reuben J Pengelly
  orcid: 0000-0001-7022-645X
  affiliation: 1, 3
- name: Timothy J Sluckin
  orcid: 0000-0002-9163-0061
  affiliation: 2
- name: Andrew Collins
  orcid: 0000-0001-7108-0771
  affiliation: 1, 3
bibliography: paper.bib
affiliations:
- name: Genetic Epidemiology and Bioinformatics, Faculty of Medicine, University of Southampton
  index: 1
- name: Mathematical Sciences, University of Southampton
  index: 2
- name: Institute for Life Sciences, University of Southampton
  index: 3
---

# Summary

Detecting evidence of selection and evolution in population genomes is crucial to understanding the history and the selective pressures experienced by a population. While there are many statistics for identifying regions of the genome under selection, there is a need for software to enable reproducible, standardised results. The statistics implemented in the `zalpha` R package use the relationships and correlations in genetic variation to find patterns that could be indicative of a selective sweep. 

The methods contained within this R package are a development of the statistics published by @Jacobs:2016. This package allows users to run a range of selection statistics on genetic data, which previously were not made publicly available in software. The software is designed to be flexible to allow users to efficiently combine statistics and is open source.

The package also allows users to utilise a linkage disequilibrium (LD) profile, taking into account expected relationships between alleles, ultimately increasing the power of the statistics. This is important as LD varies immensely along the genome, with recombination the biggest contributor to LD fluctuations [@Jeffreys:2001].

# Statement of Need

The purpose of the `zalpha` package is to:

* Allow users to accurately apply the $Z_{\alpha}$ statistic to find candidate regions of the genome for a selective sweep
* Refine $Z_{\alpha}$ results by adjusting for expected correlations between genetic variants
* Further characterise sweeps as ongoing or near fixation using the $Z_{\beta}$ statistic
* Generate results that are reproducible
* Be user-friendly and accessible by using R
	
# Software and Methodology

The `zalpha` package examines correlations between single nucleotide polymorphisms (SNPs) along a chromosome. If SNPs are highly correlated in a region of a chromosome in relation to the rest of the genome, this could indicate the presence of a selective sweep [@Vitti:2013].

Correlation, in the context of genetics, is the ability to predict the value of one SNP, given the value of another. An example is given in \autoref{fig:Figure1}A. The metric used by these statistics to measure correlation is $r^2$ [@Cutter:2019].

![(A) This figure shows a population of seven chromosomes with 10 biallelic SNPs highlighted in orange and yellow to indicate the two alleles of the SNP. The cluster of five SNPs on the left are perfectly correlated. The second cluster on the right are poorly correlated: SNPs cannot be used to predict each other. (B) This figure shows a chromosome with SNPs highlighted. A statistic is being calculated for the target locus in green. The window is defined, and all SNPs falling within the window to the left and right of the target locus are assigned to sets L and R respectively. The correlations between pairs of SNPs are evaluated and represented by the black circles for SNPs within sets L and R, and as purple squares for pairs between the sets. \label{fig:Figure1}](Figure1.png)

When a selective sweep occurs, the locus under selection becomes more frequent in the population, as individuals possessing the beneficial allele are more likely to survive and reproduce. When this happens, variants nearby the selected locus will also sweep, a phenomenon known as “hitchhiking” [@Maynard:1974]. This creates a region of the genome that is highly correlated. Eventually recombination will erode away these correlations. 

`zalpha` allows the user to apply a range of statistics to genetic data. \autoref{fig:Figure1}B shows a target locus with a window, the size of which is set by the user, centred on the locus. Any SNPs either side that fall within the window to the left and right of the target locus are contained within sets L and R respectively. The statistic $Z_{\alpha}$, after which the package is named, is defined as:
\begin{equation}\label{eq:Zalpha}
{Z_{\alpha}=\frac{{|L| \choose 2}^{-1}\sum_{i,j \in L}r^2_{i,j} + {|R| \choose 2}^{-1}\sum_{i,j \in L}r^2_{i,j}}{2}}
\end{equation}
|L| and |R| are the number of SNPs in each set, and $r^2_{i,j}$ is the correlation between two SNPs i and j. \autoref{fig:Figure1}B shows these $r^2$ values as black circles. 

The other base statistic supplied in the `zalpha` package is $Z_{\beta}$, as defined as follows:
\begin{equation}\label{eq:Zbeta}
{Z_{\beta}=\frac{\sum_{i \in L,j \in R}r^2_{i,j}}{|L||R|}}
\end{equation}
In \autoref{fig:Figure1}B the $r^2$ values for $Z_{\beta}$ are represented as purple squares.

Typically, a user will want to find the maximum $Z_{\alpha}$ statistic in a region of a chromosome, and compare this to other regions, to find possible evidence of selection for that region.

The package is designed to be as user-friendly as possible and is reflected in the flexibility of the input requirements. The basic statistics only require three elements:

* vector of physical locations of each SNP,
* a window size, and 
* a matrix of SNP values where the rows are SNPs and the columns are haplotypes. This matrix could be binary, where the 0s represent ancestral alleles and the 1s derived, or it could be nucleotides (i.e. As, Cs, Gs, and Ts), or any other biallelic labelling system.

One of the benefits of this package is the ability to calculate multiple statistics simultaneously. During a selective sweep, the correlations between alleles near to the selected locus increase. This means both $Z_{\alpha}$ and $Z_{\beta}$ should be higher than in other areas of the genome not experiencing selective pressure. Towards the end of a selective sweep however, the correlations between the sets of alleles on the left and the right of the target locus are expected to diminish [@Kim:2004]. This suggests at the end of a sweep, $Z_{\alpha}$ should remain high, but $Z_{\beta}$ will reduce. Thus, it is advantageous to calculate and combine the different statistics to ascertain the strength and stage of sweeps. $Z_{\alpha}$/$Z_{\beta}$ is a simple way to achieve this.

Recombination is a process that has the effect of breaking down the relationship between alleles. However, it is known that recombination does not occur uniformly across the genome. It is therefore imperative to consider recombination when calculating statistics based on LD measures. This package allows the user to supply a population LD profile, providing information on the expected relationships between alleles given the genetic distances between them. Supplying these data increases the power of the statistics and creates more opportunities for combinations and comparisons between statistics. Users can specify whatever units they wish for genetic distance (for example centimorgans (cM)), derived from an appropriate data source. The software contains a function for creating an LD profile from the data. Ideally, an LD profile would be created from a neutral data source without selection, for example from a simulation with relevant population parameters. However, this is not always possible, so creating an LD profile from the same data being analysed is sufficient.

There are many statistics included in the package for adjusting for expected $r^2$ using the LDprofile and genetic distances between SNPs.  It is recommended the user runs all the statistics using the `Zalpha_all()` function and then chooses the ones they are interested in, perhaps even creating their own. For example, $Z_{\alpha}$/${Z_{\alpha}^{E[r^2]}}$ performs well as a simple way to adjust for expected $r^2$. If it is known that the $r^2$ values for each genetic distance are normally distributed, ${Z_{\alpha}^{Zscore}}$ is appropriate, otherwise ${Z_{\alpha}^{BetaCDF}}$ may be useful. For more details of how they are derived see the paper by @Jacobs:2016. This paper also shows how the different statistics perform under a range of demographic scenarios.

The output of the functions is in list format. The SNP positions and the values of the statistic(s) are stored in vectors of equal length in the list. Users can then identify outlying SNPs in their data that are candidate regions for selection.

There are a few other R packages that can be used for selection scans, although none utilise the $Z_{\alpha}$ statistics described here. These include `PopGenome` [@Pfeifer:2014], which calcualtes Kelly's $Z_{nS}$ among other methods, and `rehh` [@Gautier:2017], which implements extended haplotype homozygosity (EHH) and related statistics.

# Conclusion

This new package allows researchers to calculate the $Z_{\alpha}$ suite of selection statistics efficiently using the free, open source R platform. These statistics had previously not been publicly available in software. The package's flexibility allows the user to adjust the statistics for the expected $r^2$ value via an LD profile in a variety of ways, and enables the adjustment of the base statistics to create new and novel methods.

# Acknowledgements

The authors acknowledge the use of the IRIDIS High Performance Computing Facility in the completion of this work. The authors would like to acknowledge Guy Jacobs whose work on these statistics was the foundation for this package.

# References
