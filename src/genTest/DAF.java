package genTest;

import java.util.ArrayList;
import java.util.List;

import log.Log;
import tools.Individual;
import tools.SNP;
import tools.Window;

public class DAF extends HaplotypeTests {
	
	private Window tp_win;
	private Individual[] tp_indv;
	
	private List<Window> xoin_wins;
	private Individual[] xp_ino_indv;//previously intersected with op
	private Individual[] op_inx_indv;//previously intersected with xp
	
	private List<SNP> anc_types;
	
	//ÆDAF statistic information
	private List<SNP> unused_snps;
	private List<SNP> all_DAF_snps;
	private List<Double> all_DAF;
	
	private Log log;
	
	public DAF(Log log, 
				Window tp_win, 
				Individual[] tp_indv,
				List<Window> xoin_wins, //I might want this to be the non-selected (intersected) window because then I can just find the index of a particular snp making life easy
				Individual[] xp_ino_indv, //this might be the 2 non-selected populations that are intersected with each other (if that's the case I'll have to start over)
				Individual[] op_inx_indv,
				List<SNP> anc_types){
		
		this.log = log;
		
		this.tp_win = tp_win;
		this.tp_indv = tp_indv;
		
		this.xoin_wins = xoin_wins;
		this.xp_ino_indv = xp_ino_indv;
		this.xp_ino_indv = op_inx_indv;
		
		this.anc_types = anc_types;
		
		unused_snps = new ArrayList<SNP>();
		all_DAF_snps = new ArrayList<SNP>();
		all_DAF = new ArrayList<Double>();
	}
	
	public void runStat() {
		
		int st_index = tp_win.getStIndex();
		List<SNP> win_snps = tp_win.getSNPs();
		for(int i = 0; i < win_snps.size(); i++) {
			
			SNP core_snp = win_snps.get(i);
			SNP anc_snp = getAncestralSNP(core_snp, anc_types);
			
			if(checkValidSnpComparison(core_snp, anc_snp)) {
				
//				Individual[] all_indv = combineIndvArrays(tp_indv, xp_indv);
				int tp_instance_der = getInstanceOfDerivedAllele(tp_indv, core_snp, anc_snp, (st_index + i));
				int xp_instance_der = getInstanceOfDerivedAllele(xp_ino_indv, core_snp, anc_snp, (st_index + i));
				
				int tot_instance_der = xp_instance_der + xp_instance_der;
				
				double xp_freq_der = (double) tp_instance_der / (double) tp_indv.length;
				double tot_freq_der = (double) tot_instance_der / ((double) tp_indv.length + (double) xp_ino_indv.length);
				
				double delta_daf = xp_freq_der - tot_freq_der;
				//TODO: save delta_daf
			}
			else
				unused_snps.add(core_snp);
			
			
			
			
		}
		
		//find instance of derived allele in the "non-selected" population (d_ns)
		//find instance of derived allele in the target population (d_s)
		//ÆDAF = d_ns - d_s
		//scores should be between -1 and 1
		
	}
	
	private int getInstanceOfDerivedAllele(Individual[] indv, SNP core_snp, SNP anc_snp, int snp_index) {
		
		int count = 0;
		
		//When core_snps's a1 = derived type
		if(core_snp.getAllele0().equals(anc_snp.getAllele0())
				&& core_snp.getAllele1().equals(anc_snp.getAllele1())) {
			for(int i = 0; i < indv.length; i++) {
				
				//adding the index of the individual with the corresponding strand (1)
				int st1_allele = indv[i].getStrand1Allele(snp_index);
				if(st1_allele == 1)
					count++;
				
				//adding the index of the individual with the corresponding strand (2)
				int st2_allele = indv[i].getStrand2Allele(snp_index);
				if(st2_allele == 1)
					count++;
			}
		}
		//When core_snps's a0 = derived type
		else if(core_snp.getAllele0().equals(anc_snp.getAllele1())
				&& core_snp.getAllele1().equals(anc_snp.getAllele0())) {
			for(int i = 0; i < indv.length; i++) {
				
				int st1_allele = indv[i].getStrand1Allele(snp_index);
				if(st1_allele == 0)
					count++;
				
				int st2_allele = indv[i].getStrand2Allele(snp_index);
				if(st2_allele == 0)
					count++;
			}
		}
	
		return count;
	}
}
