# CMS4J
A user-friendly version of the composite of multiple signals (CMS) implemented in Java

## Hayden will tell you how AWESOME CMS4J is here:
It rocks...
https://help.github.com/articles/github-flavored-markdown/

Output for a functioning program will look like this:

```
======================================================Round 1======================================================================
Starting iHS Analysis
Starting XPEHH Analysis
XPEHH =	SNP [pos=10873592, a0=G, a1=A, rs_number=rs1395681]	-0.7177079065858938	0.9999999999999987

XPEHH =	SNP [pos=10881456, a0=T, a1=C, rs_number=rs28578587]	-0.8013082086861588	-1.0000000000000013

Starting iHH Analysis

======================================================Round 2======================================================================
					    Last 5 lines of every stat

Starting iHS Analysis

iHS =	SNP [pos=14781959, a0=A, a1=G, rs_number=rs446549]	-1.0223001972313361	-1.6838548658693318

iHS =	SNP [pos=14785028, a0=G, a1=A, rs_number=rs431133]	0.5864846531018965	1.7157701287434333

iHS =	SNP [pos=14818410, a0=T, a1=C, rs_number=rs2821870]	-0.6195014088871154	-0.832675276273427

iHS =	SNP [pos=14880003, a0=A, a1=C, rs_number=rs13047108]	-0.36601228485709697	-0.2970113795525913

iHS =	SNP [pos=14931069, a0=C, a1=G, rs_number=rs117062102]	0.44781812890226097	1.4227451282360557

iHS =	SNP [pos=14936462, a0=T, a1=C, rs_number=kgp7953532]	-0.36929558509399946	-0.30394952889423865

Starting XPEHH Analysis

XPEHH =	SNP [pos=14818410, a0=T, a1=C, rs_number=rs2821870]	-1.06020528606459	0.3318146025182604

XPEHH =	SNP [pos=14880003, a0=A, a1=C, rs_number=rs13047108]	-0.7222414566778952	1.9208737366915924

XPEHH =	SNP [pos=14882965, a0=G, a1=A, rs_number=rs144165180]	-1.0105401853071685	0.5653330300714032

XPEHH =	SNP [pos=14931069, a0=C, a1=G, rs_number=rs117062102]	-0.7168941729868864	1.946015924267208

XPEHH =	SNP [pos=14936462, a0=T, a1=C, rs_number=kgp7953532]	-1.027941472371616	0.4835145875288617

Starting iHH Analysis

iHH =	SNP [pos=14785028, a0=G, a1=A, rs_number=rs431133]	1327110.938517952	0.36275875689675086

iHH =	SNP [pos=14818410, a0=T, a1=C, rs_number=rs2821870]	844465.1976076569	-0.37077135068839495

iHH =	SNP [pos=14880003, a0=A, a1=C, rs_number=rs13047108]	1355246.1671102513	0.4055189752244159

iHH =	SNP [pos=14931069, a0=C, a1=G, rs_number=rs117062102]	983494.2243573167	-0.1594735637094967

iHH =	SNP [pos=14936462, a0=T, a1=C, rs_number=kgp7953532]	1071699.7875737464	-0.025417814758463306


======================================================Round 3======================================================================
                                            Last 5 lines of every stat

Starting iHS Analysis

iHS =	SNP [pos=15482605, a0=A, a1=G, rs_number=rs1556276]	-1.8037272191521176	-0.7553394058038247

iHS =	SNP [pos=15491560, a0=G, a1=A, rs_number=rs2403729]	-0.9169210590929179	0.023705225279158363

iHS =	SNP [pos=15493361, a0=G, a1=A, rs_number=rs2155971]	-1.534343807158523	-0.5186905074668807

iHS =	SNP [pos=15494781, a0=T, a1=C, rs_number=rs11701216]	-1.4606836789320814	-0.45398129591171676

iHS =	SNP [pos=15494955, a0=C, a1=T, rs_number=rs2032283]	-0.8606144869173765	0.07316962316416369

iHS =	SNP [pos=15495000, a0=G, a1=A, rs_number=rs2032284]	-2.2655568537089543	-1.1610491465351045


Starting XPEHH Analysis

XPEHH =	SNP [pos=15493361, a0=G, a1=A, rs_number=rs2155971]	0.3857381544408967	1.6719614701450984

XPEHH =	SNP [pos=15494621, a0=G, a1=A, rs_number=rs2263412]	0.39102015701683474	1.6817246194659956

XPEHH =	SNP [pos=15494781, a0=T, a1=C, rs_number=rs11701216]	0.38966805022895307	1.6792254119631398

XPEHH =	SNP [pos=15494955, a0=C, a1=T, rs_number=rs2032283]	0.38972083611853886	1.6793229803597194

XPEHH =	SNP [pos=15495000, a0=G, a1=A, rs_number=rs2032284]	0.2174002298787614	1.3608089725380261


Starting iHH Analysis

iHH =	SNP [pos=15491560, a0=G, a1=A, rs_number=rs2403729]	267644.2342478619	-0.777109149058452

iHH =	SNP [pos=15493361, a0=G, a1=A, rs_number=rs2155971]	578952.2145706662	-0.39462384332500694

iHH =	SNP [pos=15494781, a0=T, a1=C, rs_number=rs11701216]	580826.3380396814	-0.39232122105594325

iHH =	SNP [pos=15494955, a0=C, a1=T, rs_number=rs2032283]	258279.17220160188	-0.7886154349437742

iHH =	SNP [pos=15495000, a0=G, a1=A, rs_number=rs2032284]	1329801.3930938174	0.5278991749070961

```

You can do a 4th round on a much larger Window size but if these are working and that round runs
all the way through I would expect it to be correct as well

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


