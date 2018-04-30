## The code in this repository is sufficient to completely recreate the figures from Schaffer, Stettler, Kato, Choi, Axel, & Abbott (2018). Odor Perception on the Two Sides of the Brain: Consistency Despite Randomness. *Neuron*, https://doi.org/10.1016/j.neuron.2018.04.004.

##### makeFig1.m generates Figure 1. Requires that:

	1.  partsForFig1odorPairComparison.mat has already been generated using calculateFig1parts.m.
	2.  partsForSumFig.mat has been generated using calculateSumFigParts.m and calculateSumFigPart2.m.
	3.  part1ForImportanceFig_pInv.mat and part2ForImportanceFig_pInv.mat have been generated using calculateImportanceFigParts_pInv.m.


##### makeFigure2.m generates Figure 2. Requires that:

	1.  partsForSumFig.mat has been generated using calculateSumFigParts.m and calculateSumFigPart2.m.
	2.  flyPop4thOrderResp…mat has been generated using calculate4thOrderPopResp_Fly_antiHebb.m
	3.  pop4thOrderResp…mat had been generated using calculate4thOrderPopResp_v2.m.
	4.  pop4thOrderRespSmall_…mat had been generated using calculate4thOrderPopResp_v2small.m.