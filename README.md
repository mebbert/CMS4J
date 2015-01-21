# CMS4J
A user-friendly version of the composite of multiple signals (CMS) implemented in Java

### SetOperator (SO)
SetOperator will perform intersects, unions, and complements on multi- or
single-sample VCFs while considering sample genotypes. 
Some important features for SetOperator are:

* Handles multi-sample VCFs
* Genotype aware set operations (e.g., only intersect on hets)
* Powerful set operation syntax
* Operation stringing (predefine operations whose result should be passed
  directly to the next operation
* Auto 'chr' handling (i.e., we don't care if some VCFs have 'chr' preprended to
  the #CHROM value and some don't.)
* INDEL fuzzy matching (sometimes INDELS don't align the same even though they
  are the same). We'll tell you about these matches.


