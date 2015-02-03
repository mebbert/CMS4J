package genTest;

import java.util.List;

import log.Log;
import tools.Individual;
import tools.Window;

public class Fst extends HaplotypeTests{
	
	static int NUM_OF_POPULATIONS = 3;
	
	private Window win;
	
	//All three populations intersected with each other
	private Individual[] tp_indv;
	private Individual[] xp_indv;
	private Individual[] op_indv;
	
	private Log log;
	
	
	public Fst(Log log,
				Window txin_win, 
				Individual[] tp_inx_indv,
				Individual[] xp_int_indv,
				Individual[] op_inx_indv) {
		
		this.win = txin_win;
		
		this.tp_indv = tp_inx_indv;
		this.xp_indv = xp_int_indv;
		this.op_indv = op_inx_indv;
		
		this.log = log;
	}
	
	public void runStat() {
		//calculate the freq of allele A in each of the 3 populations (CEU, YRI, JPT)
		//find weighted average of those frequencies by number of individuals (haplotypes???)
		//figure out sample variance, s^2
		//plug numbers into Fst eqn
		
		//TODO: answer question, Is there a difference in Fst if you use a0 or a1?
		//		should I just do the analysis on the derived allele type? probably not b/c then I'm limited by Ancestral data
		
	}

}
