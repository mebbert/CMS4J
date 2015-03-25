package tools;

import java.util.ArrayList;
import java.util.List;

public class SimDist {
	
	private final int BIN_NUM = 60;
	
	private int up_bndry;
	private int low_bndry;
	private double total_prob;
	
	private List<Double> sim_vals;
	
	public SimDist(int low_bndry, int up_bndry) {
		
		this.up_bndry = up_bndry;
		this.low_bndry = low_bndry;
		this.total_prob = 0.0;
		
		sim_vals = new ArrayList<Double>();
	}

	public void addSimValue(double val) {
		
		sim_vals.add(val);
		total_prob += val;
	}
	
	public Double getProb(Double score) {
		
		//TODO: Handle the 2-sided case **************************
		
		if(sim_vals.size() != BIN_NUM) {
			//throw new FileParsingException("Bin numbers don't coincide, error in reading simulated data"
			System.out.println("ERROR WITH BIN NUMBER");
			System.exit(0);
		}
		
		int s_indx = getScoreIndex(score, up_bndry, low_bndry);
		
		return calcProbAtBin(s_indx);
	}
	
	public Double getProb(Double score, boolean two_sided) {
		
		if(!two_sided)
			return getProb(score);
		
		if(sim_vals.size() != BIN_NUM) {
			//throw new FileParsingException("Bin numbers don't coincide, error in reading simulated data"
			System.out.println("ERROR WITH BIN NUMBER");
			System.exit(0);
		}
		
		score = Math.abs(score);
		
		int up_indx = getScoreIndex(score, up_bndry, low_bndry);
		int low_indx = getScoreIndex((-1 * score), up_bndry, low_bndry);
		
		return calcTwoSidedProbAtBin(up_indx, low_indx);
	}
	
	public List<Double> getSimVals() {
		return sim_vals;
	}
	
	private Double calcProbAtBin(int indx) {
		
		Double prob = 0.0;
		
		for(int i = 0; i <= indx; i++)
			prob += sim_vals.get(i);
		
		return prob;
	}
	
	private Double calcTwoSidedProbAtBin(int up_indx, int low_indx)  {
		
		Double prob = 0.0;
		
		for(int i = low_indx; i <= up_indx; i++)
			prob += sim_vals.get(i);
		
		return prob;
	}
	
//	private Double calcTotalProb() {
//		
//		Double prob = 0.0;
//		
//		for(int i = 0; i < sim_vals.size(); i++)
//			prob += sim_vals.get(i);
//		
//		return prob;
//	}
	
	private int getScoreIndex(Double score, int up, int dwn) {
		
		double rng = (double) up - (double) dwn;
		double bin_size = rng / (double) BIN_NUM;
		
//		if(score == 0.513382123798503) {
//			System.out.println("range=" + rng + "\tbin size=" + bin_size);
//		}
		
		for(int i = 0; i < BIN_NUM; i++) {
			if((dwn + bin_size*i) >= score)
				return i;
		}
		
		return BIN_NUM-1;
	}
	
}
