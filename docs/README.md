## Maximum Likelihood Estimation for Multi-Site Metabolomics Research 

Studying the statistical relationships between metabolites in a population may help researchers:
- Understand metabolic pathways
- Identify drug interactions
- Conceive new therapeutic approaches
- More!

Metabolomic research may take place in a consortium setting, where members of the consortium share **aggregate-level information** such as:
- Sample size
- Spearman rank correlation matrices

Members may be unwilling to share **subject-level data** due to privacy, data governance or other concerns. 

Therefore, a statistical methodology is proposed which is able to correctly handle conditions typical of consortium metabolomic studies. 

One problem in multi-site studies with large numbers of variables is that not every site is able to take measurements on all variables. Different compounds are best studied using different types of laboratory equipment. Smaller laboratories might have mass spectrometers or NMR systems but perhaps not both. So *different sites do not necessarily detect all the same metabolites*. 

Fundamental insight: 
- If I know something about the relationship between variables $A$ and $B$ based on dataset 1
- And I know something about the relationship between variables $B$ and $C$ based on dataset 2
- For free, I get to know something about the relationship between $A$ and $C$, even though $A$ and $C$ never appear together in either dataset

The following is an illustration of the proposed methodology: 



![Alt text](./img/diagram.svg )

To use this project, run `entry.m`. Explanatory notes are provided as in-line comments. 