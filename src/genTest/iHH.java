package genTest;

import java.util.ArrayList;
import java.util.List;

import log.Log;
import tools.ExtendedHaplotype;
import tools.GeneticMap;
import tools.Individual;
import tools.SNP;
import tools.Window;

public class iHH extends HaplotypeTests {
	
	private Window win;
	private Individual[] individuals;
	private GeneticMap gm;
	private ExtendedHaplotype anc_eh;
	private ExtendedHaplotype der_eh;
	
	private List<SNP> anc_types;
	private List<Double> all_unstd_iHH;
	private List<Window> all_win;
	
	//iHH statistic information
	private List<SNP> unused_snps;
	private List<SNP> all_iHH_snp;
	private List<Double> all_std_iHH;
	
	private Log log;
	
	/**
	 * For setting up the environment to run the ÆiHH statistic
	 * 
	 * @param log			universal log for progress and error output		
	 * @param win			current Window within the target population (tp)
	 * @param individuals	all Individuals of the target population
	 * @param anc_types		all Ancestral types in the form of SNPs; ancestral type is a0
	 * @param all_win		all Windows in the tested region, usually the chr
	 * @param gm			Genetic Map for the tested region, usually the chr
	 */
	public iHH(Log log, 
			Window win, 
			Individual[] individuals, 
			List<SNP> anc_types, 
			List<Window> all_win, 
			GeneticMap gm) {
		
		this.log = log;
		
		this.win = win;
		this.individuals = individuals;
		this.gm = gm;
		
		this.anc_types = anc_types;
		this.all_win = all_win;
		
		anc_eh = new ExtendedHaplotype();
		der_eh = new ExtendedHaplotype();
		
		unused_snps = new ArrayList<SNP>();
		all_iHH_snp = new ArrayList<SNP>();
		all_unstd_iHH = new ArrayList<Double>();
		all_std_iHH = new ArrayList<Double>();
	}

	/**
	 * Runs the ÆiHH statistic using the environment setup by the constructor. The
	 * below series of evens spans multiple private methods but outlines the 
	 * logic for calculating ÆiHH. This is done for every SNP in the Window (core_snp)
	 * 		-Step 1: Create extended haplotypes for all Individuals (2 per Individual)
	 * 		-Step 2: Separate the pool of haplotypes based upon whether or not they have the Ancestral allele at the core position
	 * 		-Step 3: Calculated EHH values for both groups until reaching a significantly insignificant EHH value (EHH = 0.05)
	 * 		-Step 4: Integrate EHH values from core to ends for both Ancestral and Derived groups; weight value based upon Genetic Map
	 * 		-Step 5: Calculate ÆiHH from previously calculated iHH values
	 * 		-Step 6: Repeat and save these unstandard ÆiHH values for all SNPs in Window
	 * 		-Step 7: Standardize the all the ÆiHH values within the Window
	 * 
	 * Note that many of these functions are extended from HaplotypeTests and 
	 * can't be found in this class.
	 */
	@Override
	public void runStat() {
		
		System.out.println("Starting iHH Analysis");
		log.addLine("Starting iHH Analysis");
		
		//Starting iHH Analysis
		int st_index = win.getStIndex();
		for(int i = 0; i < win.getSNPs().size(); i++) {
			
			Double unstd_iHH = getUnstandardizedIHH(win.getSNPs().get(i), (st_index + i));
			
			//saving the successful unstandardized ÆiHH
			if(unstd_iHH != null)
				all_unstd_iHH.add(unstd_iHH);
		}
		
		//calculating and saving all standardized ÆiHH values
		all_std_iHH = standardizeData(all_unstd_iHH);
		
//		printStats();
//		logRunStats();
	}
	
	@Override
	public List<SNP> getSNPs() {
		return all_iHH_snp;
	}
	
	@Override
	public List<Double> getStats() {
		return all_std_iHH;
	}
	
	@Override
	public void printStats() {
//		===============Default Printout===================
//		System.out.println("\nShowing ÆiHH Data");
//		for(int i = 0; i < all_std_iHH.size(); i++) {
//			System.out.print("iHH =\t");
//			System.out.print(all_iHH_snp.get(i) + "\t");
//			System.out.print(all_unstd_iHH.get(i) + "\t");
//			System.out.println(all_std_iHH.get(i));	
//		}
		
//		===============R Printout==========================
		StringBuilder ihh_sb = new StringBuilder();
		StringBuilder pos_sb = new StringBuilder();
		
		System.out.println("\nShowing R output: iHH");
		for(int i = 0; i < all_std_iHH.size(); i++) {
			
			ihh_sb.append(all_std_iHH.get(i) + ",");
			pos_sb.append(all_iHH_snp.get(i).getPosition() + ",");
		}
		System.out.println("iHH =\t" + ihh_sb.toString());
		System.out.println("Pos =\t" + pos_sb.toString());
		
	}

	@Override
	public void logRunStats() {
		
		log.addLine("Out of " + win.getSNPs().size() + " SNPs, " 
				+ all_std_iHH.size() + " were successful and " + unused_snps.size() 
				+ " SNPs were unsuccessful");
	}
	
	private Double getUnstandardizedIHH(SNP core_snp, int snp_index) {
		
		double unstd_iHH = 0.0;
		
		SNP anc_snp = getAncestralSNP(core_snp, anc_types);
		
		//TODO: check if there are any reverse compliment data being thrown out here
		if(checkValidSnpComparison(core_snp, anc_snp)) {
				
			//Initial Grouping (according to ancestral or derived type)
			setHaplotypeGroups(anc_eh, der_eh, individuals, snp_index, anc_snp, core_snp);
			
			if(anc_eh.size() <= 1 || der_eh.size() <= 1) {
				//No variance and thus no EHH pattern can be found
				unused_snps.add(core_snp);
				return null;
			}
			
			//Starting EHH Analysis
			EHH anc_ehh = new EHH(win, individuals, core_snp, anc_eh, all_win);
			EHH der_ehh = new EHH(win, individuals, core_snp, der_eh, all_win);
			
			//Running Ancestral Analysis
			anc_ehh.calcSignificantEhhValues();
			double[] ehh_values_anc = anc_ehh.getEhhValues();
			int[] ehh_pos_anc = anc_ehh.getEhhPositions();
			
			//Running Derived Analysis
			der_ehh.calcSignificantEhhValues();
			double[] ehh_values_der = der_ehh.getEhhValues();
			int[] ehh_pos_der = der_ehh.getEhhPositions();
			
			//find the area under the curve created by the EHH data
			double anc_ihh = integrateEhhValues(ehh_values_anc, ehh_pos_anc, core_snp, gm);
			double der_ihh = integrateEhhValues(ehh_values_der, ehh_pos_der, core_snp, gm);
			
			//main ÆiHH function; unstandardized
			unstd_iHH = Math.abs(anc_ihh - der_ihh);
		}
		else {
			//No variance and thus no EHH pattern can be found
			unused_snps.add(core_snp);
			return null;
		}
		
		if(Double.isNaN(unstd_iHH) || Double.isInfinite(unstd_iHH)) {
			//Error in calculating iHH
			unused_snps.add(core_snp);
			return null;
		}
		
		//saving the successful ÆiHH SNP
		all_iHH_snp.add(core_snp);
		return unstd_iHH;
	}
}
