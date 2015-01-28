package cms;

import java.util.List;

import log.Log;
import tools.GeneticMap;
import tools.Individual;
import tools.SNP;
import tools.Window;
import tools.WindowStats;
import genTest.*;

public class Stats {
	
	private iHS i;
	private iHH h;
	private XPEHH x;
	private DAF d;
	private Fst f;
	
	public Stats(Log log, 
				Window tp_win, 
				List<Window> all_win,
				Individual[] tp_indv, 
				
				Window txin_win,
				List<Window> txin_wins,
				Individual[] tpin_indv,
				Individual[] xpin_indv,
				
				List<SNP> anc_types, 
				GeneticMap gm) {
		
		i = new iHS(log, tp_win, tp_indv, anc_types, all_win, gm);
		h = new iHH(log, tp_win, tp_indv, anc_types, all_win, gm);
		x = new XPEHH(log, txin_win, txin_wins, tpin_indv, xpin_indv, gm);
//		d = new DAF(log, tp_win, tp_indv, );
		f = new Fst();
	}

	/**
	 * Runs the major 5 population genetics statistics (iHS, iHH, XPEHH, DAF,
	 * and Fst) needed for a successful CMS analysis. If one test fails the 
	 * entire window is thrown out.
	 * 
	 * @param win		The window (default size 1Mb) by which the statistics are going to be run on
	 * @return			Returns the stats for the input window wrapped in a WindowStats object
	 */
	public void getStats() {
		
//		WindowStats ws = new WindowStats(win.getStPos(), win.getEndPos());
		
		System.out.println("\n\n\t***Starting stats run***");
		
		i.runStat();
		x.runStat();
		h.runStat();
//		d.runStat();
		
		f.runStatFst();
		
//		return ws;
	}
	
}
