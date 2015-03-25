package cms;

import io.SimulationParser;
import genTest.HaplotypeTests;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
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
	
	private Map<SNP, Double> cms_scores_broad;
	private Map<SNP, Double> cms_scores_mean;
	
	private Log log;
	
//	private NormalDistribution nd;
//	private GammaDistribution gd;
//	private UniformRealDistribution urd;
	
	public Analysis(Log log) {
		
//		nd = new NormalDistribution(0, 1); 			//defined
//		gd = new GammaDistribution(1, 1);			//made up; plot(density(qgamma(pnorm(var), shape=.5, rate=5)))
//		urd = new UniformRealDistribution(-1, 1);	//defined
		
		cms_scores_broad = new HashMap<SNP, Double>();
		cms_scores_mean = new HashMap<SNP, Double>();
		
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
			calcCmsScores(cur_ws, neutral_sim, select_sim);//pass simulations too
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
	private void calcCmsScores(WindowStats ws, SimDist[] neut_sim, SimDist[] sel_sim) {//get SimDist too
		
		Map<SNP, Double> win_scores_broad = new HashMap<SNP, Double>();
		Map<SNP, Double> win_scores_mean = new HashMap<SNP, Double>();
		
		double prior_prob = 1 / (double) ws.getTotNumSNPs();
		
		List<SNP> all_snps = ws.getAllSNPs();
		Double[] tst_scores = new Double[NUM_TESTS];
		
		for(int i = 0; i < all_snps.size(); i++) {
			
			SNP cur_snp = all_snps.get(i);
			
			tst_scores[0] = ws.getIhsScore(cur_snp);
			tst_scores[1] = ws.getIhhScore(cur_snp);
			tst_scores[2] = ws.getFstScore(cur_snp);
			tst_scores[3] = ws.getDafScore(cur_snp);
			tst_scores[4] = ws.getXpehhScore(cur_snp);
			
			Double broad_cms = calcBroadCMS(cur_snp, tst_scores, neut_sim, sel_sim, prior_prob);
			if(!broad_cms.equals(Double.NaN))
				win_scores_broad.put(cur_snp, broad_cms);
			
			Double byu_cms = calcByuCMS(cur_snp, tst_scores, neut_sim, sel_sim, prior_prob);
			if(!byu_cms.equals(Double.NaN))
				win_scores_mean.put(cur_snp, byu_cms);
			
			//========Testing Mean vs Product (same filter)========
//			Double broad_cms_mean = calcBroadCMSMean(cur_snp, tst_scores, neut_sim, sel_sim, prior_prob);
//			if(!broad_cms_mean.equals(Double.NaN))
//				win_scores_mean.put(cur_snp, broad_cms_mean);
			//=====================================================
		}
		
		win_scores_broad = normalizeData(win_scores_broad);
		for(SNP key : win_scores_broad.keySet())
			cms_scores_broad.put(key, win_scores_broad.get(key));

		win_scores_mean = normalizeData(win_scores_mean);
		for(SNP key : win_scores_mean.keySet())
			cms_scores_mean.put(key, win_scores_mean.get(key));
	}
	
	private Double calcByuCMS(SNP cur_snp, 
								Double[] tst_scores, 
								SimDist[] neut_sim, 
								SimDist[] sel_sim, 
								double prior_prob) {
		
		Double[] score_probs = new Double[NUM_TESTS];
		
		for(int j = 0; j < NUM_TESTS; j++) {
			
			if(!tst_scores[j].equals(Double.NaN)) {
			
				boolean two_sided = false;
				if(j == 0 || j== 1 || j == 3)
					two_sided = true;
				
				Double neut_prob = neut_sim[j].getProb(tst_scores[j], two_sided);
				Double sel_prob = sel_sim[j].getProb(tst_scores[j], two_sided);
				
				/*
				 * To do the CMSgw 
				 * 		I create SimulationParsers from different files
				 * 		Run a similar for loop structure
				 * 		instead of the below steps I do: double cms_bf (bayes factor) = sel_prob / neut_prob
				 * 		ask the question: Do I do this in parallel with the CMSlocal analysis
				 */
				
				double cms_nom = sel_prob * prior_prob;
				double cms_denom = ((sel_prob*prior_prob) + (neut_prob*(1-prior_prob)));
				 
				score_probs[j] = cms_nom / cms_denom;
			}
		}
		
		Double final_score_mean = meanOfScores(score_probs);
		
		return final_score_mean;
	}
	
	private Double calcBroadCMS(SNP cur_snp, 
									Double[] tst_scores, 
									SimDist[] neut_sim, 
									SimDist[] sel_sim, 
									double prior_prob) {
		
		Double[] score_probs = new Double[NUM_TESTS];
		
		boolean complete_data = true;
		for(int j = 0; j < NUM_TESTS; j++) {
			
			if(tst_scores[j].equals(Double.NaN)) {
				complete_data = false;
				break;
			}
			else {
			
				boolean two_sided = false;
				if(j == 0 || j== 1 || j == 3)
					two_sided = true;
				
				Double neut_prob = neut_sim[j].getProb(tst_scores[j], two_sided);
				Double sel_prob = sel_sim[j].getProb(tst_scores[j], two_sided);
				
				/*
				 * To do the CMSgw 
				 * 		I create SimulationParsers from different files
				 * 		Run a similar for loop structure
				 * 		instead of the below steps I do: double cms_bf (bayes factor) = sel_prob / neut_prob
				 * 		ask the question: Do I do this in parallel with the CMSlocal analysis
				 */
				
				double cms_nom = sel_prob * prior_prob;
				double cms_denom = ((sel_prob*prior_prob) + (neut_prob*(1-prior_prob)));
				 
				score_probs[j] = cms_nom / cms_denom;
			}
		}
		
		Double final_score = Double.NaN;
		if(complete_data)
			final_score = productOfScores(score_probs);
		
		return final_score;
	}
	
	private Double calcBroadCMSMean(SNP cur_snp, 
								Double[] tst_scores, 
								SimDist[] neut_sim, 
								SimDist[] sel_sim, 
								double prior_prob) {

		Double[] score_probs = new Double[NUM_TESTS];

		boolean complete_data = true;
		for(int j = 0; j < NUM_TESTS; j++) {

			if(tst_scores[j].equals(Double.NaN)) {
				complete_data = false;
				break;
			}
			else {

				boolean two_sided = false;
				if(j == 0 || j== 1 || j == 3)
					two_sided = true;

				Double neut_prob = neut_sim[j].getProb(tst_scores[j], two_sided);
				Double sel_prob = sel_sim[j].getProb(tst_scores[j], two_sided);

				/*
				 * To do the CMSgw 
				 * 		I create SimulationParsers from different files
				 * 		Run a similar for loop structure
				 * 		instead of the below steps I do: double cms_bf (bayes factor) = sel_prob / neut_prob
				 * 		ask the question: Do I do this in parallel with the CMSlocal analysis
				 */

				double cms_nom = sel_prob * prior_prob;
				double cms_denom = ((sel_prob*prior_prob) + (neut_prob*(1-prior_prob)));

				score_probs[j] = cms_nom / cms_denom;
			}
		}

		Double final_score = Double.NaN;
		if(complete_data)
			final_score = meanOfScores(score_probs);

		return final_score;
	}
	
	
	
	private Map<SNP, Double> normalizeData(Map<SNP, Double> unstd_cms) {
		
		List<SNP> all_keys = new LinkedList<SNP>();
		for(SNP s : unstd_cms.keySet()) 
			all_keys.add(s);
		Collections.sort(all_keys);
		
		List<Double> all_values = new LinkedList<Double>();
		for(int i = 0; i < all_keys.size(); i++) 
			all_values.add(unstd_cms.get(all_keys.get(i)));
		
		all_values = HaplotypeTests.normalizeData(all_values);
		
		Map<SNP, Double> std_cms = new HashMap<SNP, Double>();
		for(int i = 0; i < all_keys.size(); i++)
			std_cms.put(all_keys.get(i), all_values.get(i));
		
		return std_cms;
	}
	
	private Double productOfScores(Double[] score_probs) {
		
		Double score = 1.0;
		
		for(int i = 0; i < NUM_TESTS; i++) {
			if(score_probs[i] != null) 
				score = score*score_probs[i];
		}
		
		return score;
	}
	
	private Double meanOfScores(Double[] score_probs) {
		
		int tot_tests = NUM_TESTS;
		Double score = 0.0;
		for(int i = 0; i < NUM_TESTS; i++) {
			if(score_probs[i] != null)
				score += score_probs[i];
			else
				tot_tests--;
		}
		
		return score / tot_tests;
	}

	@Override
	public String toString() {
		
		List<SNP> all_keys = new LinkedList<SNP>();
		for(SNP s : cms_scores_mean.keySet()) 
			all_keys.add(s);
		Collections.sort(all_keys);
		
		StringBuilder sb = new StringBuilder();
		for(int i = 0; i < all_keys.size(); i++)  {
			SNP key = all_keys.get(i);
			sb.append("\n" + key + "\tCMS Score=\t" + cms_scores_broad.get(key) + "\tMean CMS Score=\t" + cms_scores_mean.get(key));
		}
		
		return sb.toString();
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
