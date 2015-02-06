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
		//			-find the probability that a selected snp will have the score at that position if it is NEUTRAL
		//			-find the probability that a neutral snp will have the score at that position if it is NEURTAL (as well)
		//			-calculate the posterior probability of that snp (see literature for eqn)
		//			-do the ¹ function to incorporate all the statis into one CMS statistic
		
		
	}

}
