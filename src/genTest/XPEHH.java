package genTest;

import java.util.ArrayList;
import java.util.List;

import log.Log;
import tools.ExtendedHaplotype;
import tools.GeneticMap;
import tools.Individual;
import tools.SNP;
import tools.Window;

public class XPEHH extends HaplotypeTests {
	
	private Window win;
	private Individual[] tp_individuals;//target population (tp) intersected with the cross population
	private Individual[] xp_individuals;//cross population (xp) intersected with the target population
	private GeneticMap gm;
	
	private List<Window> all_win;
	private List<Double> all_unstd_XPEHH;
	private List<SNP> valid_snps;
	
	//XPEHH statistic information
	private List<SNP> unused_snps;//make this a hash set and do a check (define this by finding the intersect of tp and xp
	private List<SNP> all_XPEHH_snps;
	private List<Double> all_XPEHH;
	
	private Log log;
	
	public XPEHH(Log log,
					Window win,
					Individual[] tp_individuals,
					Individual[] xp_individuals,
					List<Window> all_win,
					GeneticMap gm) {
		
		this.log = log;
		
		this.win = win;
		this.tp_individuals = tp_individuals;
		this.xp_individuals = xp_individuals;
		this.gm = gm;
		
		this.all_win = all_win;
		
		valid_snps = new ArrayList<SNP>();
		unused_snps = new ArrayList<SNP>();
		all_XPEHH_snps = new ArrayList<SNP>();
		all_unstd_XPEHH = new ArrayList<Double>();
		all_XPEHH = new ArrayList<Double>();
		
	}

	public void runStat() {
		
		System.out.println("Starting XPEHH Analysis");
		log.addLine("Starting XPEHH Analysis");
		
		
		//Starting XPEHH Analysis
		Individual[] all_indv = combineIndvArrays(tp_individuals, xp_individuals);
		
		List<SNP> win_snps = win.getSNPs();
		for(int i = 0; i < win_snps.size(); i++) {
			
			SNP core_snp = win_snps.get(i);
			
			EHH comb_ehh = getCombinedEHH(all_indv, core_snp);
			double last_ehh = comb_ehh.getLastEhhValue();
			
			if(last_ehh < 0.05 && last_ehh > 0.0) {
				
				SNP last_snp = comb_ehh.getLastSNP();
				
				double tp_integral = calcUnstandardEhhIntegral(core_snp, last_snp, tp_individuals);
				double xp_integral = calcUnstandardEhhIntegral(core_snp, last_snp, xp_individuals);
				
				double unstd_XPEHH = Math.log(tp_integral / xp_integral);
				all_XPEHH_snps.add(core_snp);
				all_unstd_XPEHH.add(unstd_XPEHH);
				
			}
			else {
				System.out.println("Insignificant " + last_ehh);
				unused_snps.add(core_snp);
			}	
		}
		
		all_XPEHH = standardizeData(all_unstd_XPEHH);
		
		for(int i = 0; i < all_XPEHH.size(); i++) {
			System.out.print("XPEHH =\t");
			System.out.print(all_XPEHH_snps.get(i) + "\t");
			System.out.print(all_unstd_XPEHH.get(i) + "\t");
			System.out.println(all_XPEHH.get(i) + "\n");
		}
		
		log.addLine("Out of " + win.getSNPs().size() + " SNPs, " 
				+ all_XPEHH.size() + " were successful and " + unused_snps.size() 
				+ " SNPs were unsuccessful");
	} 
	
	private double calcUnstandardEhhIntegral(SNP core_snp, SNP last_snp, Individual[] indv) {
		
		ExtendedHaplotype pop_eh = setHaplotypeGroup(indv);
		EHH pop_ehh = new EHH(win, indv, core_snp, pop_eh, all_win);
		pop_ehh.calcEhhToPosition(last_snp.getPosition());
		
		double[] ehh_vals = pop_ehh.getEhhValues();
		int[] ehh_pos = pop_ehh.getEhhPositions();
		
		return integrateEhhValues(ehh_vals, ehh_pos, core_snp, gm);
	}
	
	private EHH getCombinedEHH(Individual[] all_indv, SNP core_snp) {
		
		ExtendedHaplotype all_eh = setHaplotypeGroup(all_indv);
		
		EHH comb_ehh = new EHH(win, all_indv, core_snp, all_eh, all_win);
		comb_ehh.calcSignificantEhhValues(0.045);
		
		return comb_ehh;
	}
	
	private Individual[] combineIndvArrays(Individual[] a, Individual[] b) {
		Individual[] tot = new Individual[a.length + b.length];
		
		for(int i = 0; i < a.length; i++)
			tot[i] = a[i];
		
		for(int j = 0; j < b.length; j++) 
			tot[j + a.length] = b[j];
		
		return tot;
	}	
	
	
}
