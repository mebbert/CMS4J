package cms;

import io.SimulationParser;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

//import org.apache.commons.math3.distribution.*;






import errors.FileParsingException;
import log.Log;
import tools.SNP;
import tools.SimDist;
import tools.WindowStats;

public class Analysis {
	
	private final int NUM_TESTS = 5;
	
	private Log log;
	
//	private NormalDistribution nd;
//	private GammaDistribution gd;
//	private UniformRealDistribution urd;
	
	public Analysis(Log log) {
		
//		nd = new NormalDistribution(0, 1); 			//defined
//		gd = new GammaDistribution(1, 1);			//made up; plot(density(qgamma(pnorm(var), shape=.5, rate=5)))
//		urd = new UniformRealDistribution(-1, 1);	//defined
		
		this.log = log;
	}
	
	/*
	 * Notes for me:
	 * iHS; |iHS| > 2 is significant. To find probability I take the abs(score) of iHS and find nd.cumulativeProbability(score)
	 * iHH is the same as iHS. Need to find nd.cumulativeProbability(abs(iHH))
	 * XPEHH is raw, just find nd.cumulativeProbability(XPEHH)
	 * Fst uses gamma (still figuring out the model parameters). Find gd.cumulativeProbability(fst)
	 * DAF uses urd. urd.cumulativeProbability(daf)
	 */
	public void runCmsAnalysis(List<WindowStats> all_ws) throws FileParsingException {
		
		System.out.println("Starting Analysis");
		
		//read in all the SNPs from input somehow
		//		-each stat list needs its own snp list that corresponds with it (map???)
		//pull out a positions where all of them correspond with each other
		//for each snp
		//		for each test
		//			-find uniform prior probability; prior=1/tot_num_snps in window
		//			-find the probability that the snp will have the score at that position if it is SELECTED (help)
		//			-find the probability that the snp will have the score at that position if it is NEURTAL (help)
		//			-calculate the posterior probability of that snp (see literature for eqn)
		//			-do the pi function to incorporate all the stats into one CMS statistic
		//				just multiply the 3 probabilities together
		//			-NOTE: I don't need ALL 5 stats to do the CMS analysis
		
		System.out.println("Importing Sim Data");
		
		SimulationParser sp = new SimulationParser(log);
		SimDist[] neutral_sim = sp.getNeutralSimulations();
		SimDist[] select_sim = sp.getSelectedSimulations();
		
		for(int i = 0; i < all_ws.size(); i++) {
			WindowStats cur_ws = all_ws.get(i);
			Map<SNP, Double> cms_scores = getCmsScores(cur_ws, neutral_sim, select_sim);//pass simulations too
		}
		
		
	}
	
	/*
	 * Simulations and scores need to be stored into arrays where:
	 * 		[0] = iHS data
	 * 		[1] = iHH data
	 * 		[2] = Fst data
	 * 		[3] = DAF data
	 * 		[4] = XPEHH data
	 */
	private Map<SNP, Double> getCmsScores(WindowStats ws, SimDist[] neut_sim, SimDist[] sel_sim) {//get SimDist too
		
		List<SNP> all_snps = ws.getAllSNPs();
		Double[] tst_scores = new Double[NUM_TESTS];
		
		for(int i = 0; i < all_snps.size(); i++) {
			
			SNP cur_snp = all_snps.get(i);
			
			tst_scores[0] = ws.getIhsScore(cur_snp);
			tst_scores[1] = ws.getIhhScore(cur_snp);
			tst_scores[2] = ws.getFstScore(cur_snp);
			tst_scores[3] = ws.getDafScore(cur_snp);
			tst_scores[4] = ws.getXpehhScore(cur_snp);
			
			Double[] score_probs = new Double[NUM_TESTS];
			
			for(int j = 0; j < NUM_TESTS; j++) {
				
//				Double neut_prob = neut_sim[i].getProb(tst_scores[i]);
//				Double sel_prob = sel_sim[i].getProb(tst_scores[i]);
				
				 
//				 score_probs[i] = cms eqn...
			}
			
//			ws_cms_scores.add(run summation-thing)
			
		}
		
		
		
		return null;
	}
	
/*=======================Work in progress: Theoretical Analysis=================
	private Map<Integer, Double> getScores(WindowStats ws) {
		
		Map<Integer, Double> scores = new TreeMap<Integer, Double>();//key = position; value = cms score
		
//		List<SNP> ihs_snps = ws.getIHSsnps();
//		List<SNP> xpehh_snps = ws.getXPEHHsnps();
//		List<SNP> ihh_snps = ws.getIHHsnps();
//		List<SNP> daf_snps = ws.getDAFsnps();
//		List<SNP> fst_snps = ws.getFSTsnps();
		
		double prior_prob = 1 / (double) ws.getTotNumSNPs();
		
		int cur_pos = ws.getNextPosition(ws.getStPos());
		while(cur_pos > 0) {
			
			List<Double> all_post_prob = new ArrayList<Double>();//to add all the post_prob values
			
			int ihs_indx = getScoreIndex(cur_pos, ws.getIHSsnps());//returns -1 if there is no index
//			double iHS_post_prob = calculatePosteriorProbability(ihs_indx, ws.getIHSstats(), nd);//return -1 if there is no probability
			
			int xpehh_indx = getScoreIndex(cur_pos, ws.getXPEHHsnps());
//			double xpehh_post_prob = calculatePosteriorProbability(xpehh_indx, ws.getXPEHHstats(), nd);
			
			int ihh_indx = getScoreIndex(cur_pos, ws.getIHHsnps());
//			double ihh_post_prob = calculatePosteriorProbability(ihh_indx, ws.getXPEHHstats(), nd);
			
			int daf_indx = getScoreIndex(cur_pos, ws.getDAFsnps());
//			double daf_post_prob = calculatePosteriorProbability(daf_indx, ws.getDAFstats(), urd);
			
			int fst_indx = getScoreIndex(cur_pos, ws.getFSTsnps());
//			double fst_post_prob = calculatePosteriorProbability(fst_indx, ws.getFSTstats(), gd);
			
			//TODO: run CMS on all post_prob
			
			cur_pos = ws.getNextPosition(cur_pos);//will return -1 if there is no more next position
		}
		
		return scores;
	}
	
	private double calculatePosteriorProbability(int indx, 
													List<Double> scores 
//													RealDistribution d
													) {
		
		if(indx == -1)
			return -1;
		
		double score = scores.get(indx);
		
//		return d.cumulativeProbability(score);
		return -1;
	}
	
	private int getScoreIndex(int pos, List<SNP> snps) {
		
		for(int i = 0; i < snps.size(); i++) {
			SNP s = snps.get(i);
			if(s.getPosition() == pos)
				return i;
		}
		
		return -1;
	}
*///============================================================================
	
	
	
	
	
	
	

}
