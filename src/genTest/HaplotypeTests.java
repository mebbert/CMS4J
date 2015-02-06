package genTest;

import java.util.ArrayList;
import java.util.List;

import tools.ExtendedHaplotype;
import tools.GeneticMap;
import tools.Individual;
import tools.SNP;

public abstract class HaplotypeTests {
	
	abstract void runStat();
//	abstract List<SNP> getUnusedSNPs();
	abstract List<SNP> getSNPs();
	abstract List<Double> getStats();
	
	protected void setHaplotypeGroups(ExtendedHaplotype anc_eh, 
				ExtendedHaplotype der_eh, 
				Individual[] individuals,
				int snp_index, 
				SNP anc_snp, 
				SNP win_snp) {
		
		//When win_snps's a0 = ancestral type
		if(win_snp.getAllele0().equals(anc_snp.getAllele0())) {
			for(int i = 0; i < individuals.length; i++) {
				
				//adding the index of the individual with the corresponding strand (1)
				int st1_allele = individuals[i].getStrand1Allele(snp_index);
				if(st1_allele == 0)
					anc_eh.add(i, 1);
				else
					der_eh.add(i, 1);
				
				//adding the index of the individual with the corresponding strand (2)
				int st2_allele = individuals[i].getStrand2Allele(snp_index);
				if(st2_allele == 0)
					anc_eh.add(i, 2);
				else
					der_eh.add(i, 2);
			}
		}
		//When win_snps's a1 = ancestral type
		else if(win_snp.getAllele1().equals(anc_snp.getAllele0())) {
			for(int i = 0; i < individuals.length; i++) {
				
				int st1_allele = individuals[i].getStrand1Allele(snp_index);
				if(st1_allele == 0)
					der_eh.add(i, 1);
				else
					anc_eh.add(i, 1);
				
				int st2_allele = individuals[i].getStrand2Allele(snp_index);
				if(st2_allele == 0)
					der_eh.add(i, 2);
				else
					anc_eh.add(i, 2);
			}
		}
	}
	
	protected ExtendedHaplotype setHaplotypeGroup(Individual[] all_indv) {
		
		ExtendedHaplotype eh = new ExtendedHaplotype();
		
		for(int i = 0; i < all_indv.length; i++) {
			eh.add(i, 1);
			eh.add(i, 2);
		}
		
		return eh;
	}
 	
	protected double integrateEhhValues(double[] ehh_vals, int[] ehh_pos, SNP core_snp, GeneticMap gm) {
		
		int core_snp_pos = core_snp.getPosition();
		int core_index = -1;
		
		for(int i = 0; i < ehh_pos.length; i++) {
			if(ehh_pos[i] == core_snp_pos)
				core_index = i;
		}
		
		double ihh_lft = 0.0;
		for(int i = core_index - 1; i >= 0; i--) {
			
			double rate = gm.getRecombRate(ehh_pos[i+1], ehh_pos[i]);
			
			int dist = ehh_pos[i+1] - ehh_pos[i];
			
			double wt_dist = (double) dist * rate;
			
			ihh_lft += wt_dist * ehh_vals[i];
		}
		
		double ihh_rt = 0.0;
		for(int i = core_index + 1; i < ehh_pos.length; i++) {
			
			double rate = gm.getRecombRate(ehh_pos[i], ehh_pos[i-1]);
			
			int dist = ehh_pos[i] - ehh_pos[i-1];
			
			double wt_dist = (double) dist * rate;
			
			ihh_rt += wt_dist * ehh_vals[i];
			
		}
		
		return ihh_lft + ihh_rt;
	}
	
	protected boolean checkValidSnpComparison(SNP core_snp, SNP anc_snp) {
		
		if(anc_snp == null || core_snp == null)
			return false;
		else if(anc_snp.getAllele0().equals(core_snp.getAllele0()))
			return true;
		else if(anc_snp.getAllele0().equals(core_snp.getAllele1())) 
			return true;
		else
			return false;
	}

	protected SNP getAncestralSNP(SNP s, List<SNP> anc_types) {
		
		for(SNP as : anc_types) {
			if(s.getPosition() == as.getPosition())
				return as;
		}
		
		return null;
	}
	
	protected Individual[] combineIndvArrays(Individual[] a, Individual[] b) {
		Individual[] tot = new Individual[a.length + b.length];
		
		for(int i = 0; i < a.length; i++)
			tot[i] = a[i];
		
		for(int j = 0; j < b.length; j++) 
			tot[j + a.length] = b[j];
		
		return tot;
	}
	
	protected List<Double> standardizeData(List<Double> all_unstd_data) {
		//Standardizes the data by calculating the z-score at every position
		//	follows this equation: z = (x-mean[x])/std[x]; 
		//		z = standardized data point
		//		x = unstandardized data point
		//		mean[x] = mean of all data in set
		//		std[x] = standard deviation of all data in set
		
		List<Double> all_data = new ArrayList<Double>();
		
		double win_expect = findMean(all_unstd_data);
		double win_stddv = findStandardDeviation(all_unstd_data, win_expect);
		
		for(Double unstd_data : all_unstd_data) 
			all_data.add(calcZScore(unstd_data, win_expect, win_stddv));
		
		return all_data;
	}
	
	private double calcZScore(double unstd_data, double expect, double stddv) {
		
		return (unstd_data - expect) / stddv;
	}
	
	private double findStandardDeviation(List<Double> all_unstd_data, double mean) {
		
		double[] all_deviations = new double[all_unstd_data.size()];
		
		for(int i = 0; i < all_unstd_data.size(); i++) {
			
			double dev = all_unstd_data.get(i) - mean; 
			dev = Math.pow(dev, 2);//or dev^2
			all_deviations[i] = dev;
		}
		
		double dev_sum = 0.0;
		for(int i = 0; i < all_deviations.length; i++)
			dev_sum += all_deviations[i];
		
		double dev_mean = dev_sum / all_deviations.length;
		
		return Math.sqrt(dev_mean);
	}
	
	private double findMean(List<Double> all_unstd_data) {
		
		double sum = 0.0;
		
		for(Double unstd_data : all_unstd_data)
			sum += unstd_data;
		
		return sum / all_unstd_data.size();
	}
}


