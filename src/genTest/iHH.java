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

	public void runStat() {
		
		System.out.println("Starting iHH Analysis");
		log.addLine("Starting iHH Analysis");
		
		//Starting iHH Analysis
		int st_index = win.getStIndex();
		for(int i = 0; i < win.getSNPs().size(); i++) {
			
			Double unstd_iHH = getUnstandardizedIHH(win.getSNPs().get(i), (st_index + i));
			
			if(unstd_iHH != null)
				all_unstd_iHH.add(unstd_iHH);

		}
		
		all_std_iHH = standardizeData(all_unstd_iHH);
		
		for(int i = 0; i < all_std_iHH.size(); i++) {
			System.out.print("iHH =\t");
			System.out.print(all_iHH_snp.get(i) + "\t");
			System.out.print(all_unstd_iHH.get(i) + "\t");
			System.out.println(all_std_iHH.get(i) + "\n");	
		}
		
		log.addLine("Out of " + win.getSNPs().size() + " SNPs, " 
				+ all_std_iHH.size() + " were successful and " + unused_snps.size() 
				+ " SNPs were unsuccessful");
		
		
	}
	
	private Double getUnstandardizedIHH(SNP win_snp, int snp_index) {
		
		double unstd_iHH = 0.0;
		
		SNP anc_snp = getAncestralSNP(win_snp, anc_types);
		
		//TODO: check if there are any reverse compliment data being thrown out here
		if(checkValidSnpComparison(win_snp, anc_snp)) {
				
			//Initial Grouping (according to ancestral or derived type)
			setHaplotypeGroups(anc_eh, der_eh, individuals, snp_index, anc_snp, win_snp);
			
			if(anc_eh.size() <= 1 || der_eh.size() <= 1) {
				//No variance and thus no EHH pattern can be found
				unused_snps.add(win_snp);
				return null;
			}
			
			//Starting EHH Analysis
			EHH anc_ehh = new EHH(win, individuals, win_snp, anc_eh, all_win);
			EHH der_ehh = new EHH(win, individuals, win_snp, der_eh, all_win);
			
			//Running Ancestral Analysis
			anc_ehh.calcSignificantEhhValues();
			double[] ehh_values_anc = anc_ehh.getEhhValues();
			int[] ehh_pos_anc = anc_ehh.getEhhPositions();
			
			//Running Derived Analysis
			der_ehh.calcSignificantEhhValues();
			double[] ehh_values_der = der_ehh.getEhhValues();
			int[] ehh_pos_der = der_ehh.getEhhPositions();
			
			double anc_ihh = integrateEhhValues(ehh_values_anc, ehh_pos_anc, win_snp, gm);
			double der_ihh = integrateEhhValues(ehh_values_der, ehh_pos_der, win_snp, gm);
			
			unstd_iHH = Math.abs(anc_ihh - der_ihh);
		}
		else {
			//No variance and thus no EHH pattern can be found
			unused_snps.add(win_snp);
			return null;
		}
		
		if(Double.isNaN(unstd_iHH) || Double.isInfinite(unstd_iHH)) {
			//Error in calculating iHH
			unused_snps.add(win_snp);
			return null;
		}
		
		all_iHH_snp.add(win_snp);
		return unstd_iHH;
	}
}
