# RankComp
Differential expression genes (DEGs) obtaining methods by R language.

RankCompV2

reference: Cai, H., et al. (2018). "Identifying differentially expressed genes from cross-site integrated data based on relative expression orderings." Int J Biol Sci 14(8): 892-900.

Discription: RankComp is based on gene pairs' rank expression order principle so it can find differential expression genes(DEGs) in a small batch datasets.

The functions here may cause much error, just show me in the issues so that I can consummate it.

## The ideas of RankCompV2 are following:

Firstly: for each sample in a dataset (no matter the normal group or disease group), uses the expression  matrix to calculate the genes' expression rank relationship.

Secondly: remaining the same rank expression order gene pairs of all samples which are involve in the same group.

Thirdly: Using the co-occurrence gene pairs of disease group and normal group as background stable REOs.

Fourth: for given gene G1, calculate four parameters: 
1. a is the number of gene pairs which G1's expression is higher than another genes(Gi) in control group's gene pairs which contain G1(pairs G1-Gi).
2. b is similar to a while a is G1's expression higher than another genes(Gj) and b is G1's expression is minor than another gene.
3. x is the number of gene pairs which G1's expression is higher than another gene Gi in normal group while G1's expression is minor than Gi in disease group.
4. y as x, which represents the number of Gene pairs which G1's expression is minor than Gj and higher in disease group.
