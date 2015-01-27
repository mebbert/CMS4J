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
				Window win, 
				Window insect_win,
				Individual[] tp_individuals, 
				Individual[] xp_individuals,
				Individual[] tp_insect_indv,
				Individual[] xp_insect_indv,
				List<SNP> anc_types, 
				List<Window> all_win,
				List<Window> insect_windows,
				GeneticMap gm) {
		
		i = new iHS(log, win, tp_individuals, anc_types, all_win, gm);
		h = new iHH(log, win, tp_individuals, anc_types, all_win, gm);
		x = new XPEHH(log, insect_win, tp_insect_indv, xp_insect_indv, insect_windows, gm);
		d = new DAF();
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
		
		i.runStat();
		x.runStat();
		h.runStat();
		
		d.runStatDAF();
		f.runStatFst();
		
//		return ws;
	}
	
}
