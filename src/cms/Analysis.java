package cms;

import java.util.List;

import tools.WindowStats;

public class Analysis {
	
	public Analysis() {}
	
	public void runCmsAnalysis(List<WindowStats> all_ws) {
		
		//read in all the SNPs from input somehow
		//		-each stat list needs its own snp list that corresponds with it (map???)
		//pull out a positions where all of them correspond with each other
		//for each snp
		//		for each test
		//			-find uniform prior probability; prior=1/tot_num_snps in window
		//			-find the probability that the snp will have the score at that position if it is SELECTED (help)
		//			-find the probability that the snp will have the score at that position if it is NEURTAL (help)
		//			-calculate the posterior probability of that snp (see literature for eqn)
		//			-do the ¹ function to incorporate all the stats into one CMS statistic
		//				just multiply the 3 probabilities together
		//			-NOTE: I don't need ALL 5 stats to do the CMS analysis
		
		
	}

}
