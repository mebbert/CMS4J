package genTest;

import java.util.ArrayList;
import java.util.List;

import log.Log;
import tools.Individual;
import tools.SNP;
import tools.Window;

public class Fst extends HaplotypeTests {
	
	static int NUM_OF_POPULATIONS = 3;
	
	//Same for all populations
	private Window win;
	
	//All three populations intersected with each other
	private Individual[] tp_indv;
	private Individual[] xp_indv;
	private Individual[] op_indv;
	
	//Fst statistic information
	private List<SNP> unused_snps;
	private List<SNP> all_Fst_snps;
	private List<Double> all_Fst;
	
	private Log log;
	
	
	public Fst(Log log,
				Window txin_win, 
				Individual[] tp_inx_indv,
				Individual[] xp_int_indv,
				Individual[] op_inx_indv) {
		
		this.log = log;
		
		this.win = txin_win;
		
		this.tp_indv = tp_inx_indv;
		this.xp_indv = xp_int_indv;
		this.op_indv = op_inx_indv;
		
		unused_snps = new ArrayList<SNP>();
		all_Fst_snps = new ArrayList<SNP>();
		all_Fst = new ArrayList<Double>();
	}
	
	@Override
	public void runStat() {
		
		System.out.println("Starting Fst Analysis");
		log.addLine("Starting Fst Analysis");
		
		List<SNP> win_snps = win.getSNPs();
		for(int i = 0; i < win_snps.size(); i++) {
			
			SNP core_snp = win_snps.get(i);
			int index = win.getSnpIndex(core_snp);
			
			double tp_size = (double) tp_indv.length*2;
			int tp_instance = getInstanceOfAllele(tp_indv, index);
			double tp_freq = (double) tp_instance / tp_size;
			
			double xp_size = (double) xp_indv.length*2;
			int xp_instance = getInstanceOfAllele(xp_indv, index);
			double xp_freq = (double) xp_instance / xp_size;
			
			double op_size = (double) op_indv.length*2;
			int op_instance = getInstanceOfAllele(op_indv, index);
			double op_freq = (double) op_instance /op_size;
			
			double avg_freq = getAverageFrequencyOfAllele(tp_freq, 
															xp_freq, 
															op_freq,
															tp_size,
															xp_size,
															op_size);
			
			double sample_var = calcSampleVariance(tp_freq, 
													xp_freq, 
													op_freq,
													avg_freq,
													tp_size,
													xp_size,
													op_size);
			
			//this is Weirs analysis, I need to do Cockerhams...
			double fst = sample_var / (avg_freq * (1 - avg_freq));
			
			all_Fst_snps.add(core_snp);
			all_Fst.add(fst);
		}
		
//		printStats();
//		logRunStats();
	}
	
	@Override
	public List<SNP> getSNPs() {
		return all_Fst_snps;
	}
	
	@Override
	public List<Double> getStats() {
		return all_Fst;
	}
	
	@Override
	public void printStats() {
//		===============Default Printout===================
		System.out.println("\nShowing Fst Data");
		for(int i = 0; i < all_Fst.size(); i++) {
			System.out.print("Fst =\t");
			System.out.print(all_Fst_snps.get(i) + "\t");
			System.out.println(all_Fst.get(i));	
		}
		
//		===============R Printout==========================
//		StringBuilder fst_sb = new StringBuilder();
//		StringBuilder pos_sb = new StringBuilder();
//		
//		System.out.println("\nShowing R output: Fst");
//		for(int i = 0; i < all_Fst.size(); i++) {
//			
//			fst_sb.append(all_Fst.get(i) + ",");
//			pos_sb.append(all_Fst_snps.get(i).getPosition() + ",");
//		}
//		System.out.println("Fst =\t" + fst_sb.toString());
//		System.out.println("Pos =\t" + pos_sb.toString());
	}

	@Override
	public void logRunStats() {
		
		log.addLine("Out of " + win.getSNPs().size() + " SNPs, " 
				+ all_Fst.size() + " were successful and " + unused_snps.size() 
				+ " SNPs were unsuccessful");
	}
	
	private double calcSampleVariance(double f1, 
										double f2, 
										double f3,
										double f_avg,
										double s1,
										double s2,
										double s3) {
		
		double mean_size = (s1 + s2 + s3) / (double) NUM_OF_POPULATIONS;
		
//		double val1 = s1 * (Math.pow((f1 - f_avg), 2));//equal to s1*(f1 - f_avg)^2
//		double val2 = s2 * (Math.pow((f2 - f_avg), 2));//equal to s2*(f2 - f_avg)^2
//		double val3 = s3 * (Math.pow((f3 - f_avg), 2));//equal to s3*(f3 - f_avg)^2
//		
//		double sample_var = (val1 + val2 + val3) / ((NUM_OF_POPULATIONS - 1) * mean_size);
		
		double val1 = (s1 * (Math.pow((f1 - f_avg), 2))) / ((NUM_OF_POPULATIONS - 1) * mean_size);
		double val2 = (s2 * (Math.pow((f2 - f_avg), 2))) / ((NUM_OF_POPULATIONS - 1) * mean_size);
		double val3 = (s3 * (Math.pow((f3 - f_avg), 2))) / ((NUM_OF_POPULATIONS - 1) * mean_size);
		
		double sample_var = val1 + val2 + val3;
		
		return sample_var;
	}
	
	private double getAverageFrequencyOfAllele(double f1, 
												double f2, 
												double f3,
												double s1,
												double s2,
												double s3) {
		
		double mean_size = (s1 + s2 + s3) / (double) NUM_OF_POPULATIONS;
		
//		double val1 = f1 * s1;
//		double val2 = f2 * s2;
//		double val3 = f3 * s3;
//		
//		double avg_freq = (val1 + val2 + val3) / (NUM_OF_POPULATIONS * mean_size);
		
		double val1 = (f1 * s1) / (NUM_OF_POPULATIONS * mean_size);
		double val2 = (f2 * s2) / (NUM_OF_POPULATIONS * mean_size);
		double val3 = (f3 * s3) / (NUM_OF_POPULATIONS * mean_size);
		
		double avg_freq = val1 + val2 + val3;
		
		return avg_freq;
	}
	
	//returns the instance of a0 (it doesn't matter which allele because Fst is a measure of freq, not instance)
	private int getInstanceOfAllele(Individual[] indv, int index) {
		
		int instance = 0;
		for(int i = 0; i < indv.length; i++) {
			if(indv[i].getStrand1Allele(index) == 0)
				instance++;
			if(indv[i].getStrand2Allele(index) == 0)
				instance++;
		}
		
		return instance;
	}

	

}
