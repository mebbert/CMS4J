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
	
	private static int WAIT_TIME = 50;
	
	private WindowStats ws;
	
	private iHS i;
	private iHH h;
	private XPEHH x;
	private DAF d;
	private Fst f;
	
	/**
	 * By creating a Stats object you setup the environment for all population
	 * statistics to run. Threading is also setup for stats analysis.
	 * 
	 * @param log			universal log for progress and error output
	 * @param tp_win		current Window within the target population (tp)
	 * @param all_win		all Windows in the tested region, usually the chr
	 * @param tp_indv		all Individuals of the target population	
	 * @param txin_win		current Window for the intersection of target-cross populations
	 * @param txin_wins		all Windows for the intersection of target-cross populations
	 * @param tp_inx_indv	all Individuals of the target population after the intersection of target-cross populations
	 * @param xp_int_indv	all Individuals of the cross population after the intersection of target-cross populations
	 * @param xoin_wins		all Windows for the intersection of cross-outgroup populations
	 * @param xp_ino_indv	all Individuals of the cross population after the intersection of cross-outgroup populations
	 * @param op_inx_indv	all Individuals of the outgroup population after the intersection of cross-outgroup populations
	 * @param anc_types		all Ancestral types in the form of SNPs; ancestral type is a0
	 * @param gm			Genetic Map for the tested region, usually the chr
	 */
	public Stats(Log log, 
				Window tp_win, 
				List<Window> all_win,
				Individual[] tp_indv, 
				
				Window txin_win,
				List<Window> txin_wins,
				Individual[] tp_inx_indv,
				Individual[] xp_int_indv,
				
				List<Window> xoin_wins,
				Individual[] xp_ino_indv,
				Individual[] op_inx_indv,
				
				List<Window> anc_types,
				GeneticMap gm) {
		
		ws = new WindowStats(tp_win.getStPos(), tp_win.getEndPos());
		
		i = new iHS(log, tp_win, tp_indv, anc_types, all_win, gm);
		h = new iHH(log, tp_win, tp_indv, anc_types, all_win, gm);
		x = new XPEHH(log, txin_win, txin_wins, tp_inx_indv, xp_int_indv, gm);
		d = new DAF(log, tp_win, tp_indv, xoin_wins, xp_ino_indv, op_inx_indv, anc_types);
		f = new Fst(log, txin_win, tp_inx_indv, xp_int_indv, op_inx_indv);
	}

	/**
	 * Runs the major 5 population genetics statistics (iHS, iHH, XPEHH, DAF,
	 * and Fst) needed for a successful CMS analysis. If one test fails the 
	 * entire window is thrown out. It should be noted that each statistic runs
	 * on its own thread creating for a parallelized computation.
	 * 
	 * @return			Returns the stats for the input window wrapped in a WindowStats object
	 */
	public WindowStats getStats() {
		
//		System.out.println("\n\n\t***Starting stats run***");
		
		
//=====================unthreaded block==================
		
		i.runStat();
		h.runStat();
		x.runStat();
		d.runStat();
		f.runStat();
		
//========================threaded block=================
//		try {
//		
//			StatsThread i_thrd = new StatsThread(i);
//			Thread.sleep(WAIT_TIME);
//			i_thrd.start();
//			
//			StatsThread h_thrd = new StatsThread(h);
//			Thread.sleep(WAIT_TIME);
//			h_thrd.start();
//			
//			StatsThread x_thrd = new StatsThread(x);
//			Thread.sleep(WAIT_TIME);
//			x_thrd.start();
//			
//			StatsThread d_thrd = new StatsThread(d);
//			Thread.sleep(WAIT_TIME);
//			d_thrd.start();
//			
//			StatsThread f_thrd = new StatsThread(f); 
//			Thread.sleep(WAIT_TIME);
//			f_thrd.start();
//			
//			synchronize(i_thrd, h_thrd, x_thrd, d_thrd, f_thrd);
//			
//			i = (iHS) i_thrd.getTest();
//			h = (iHH) h_thrd.getTest();
//			x = (XPEHH) x_thrd.getTest();
//			d = (DAF) d_thrd.getTest();
//			f = (Fst) f_thrd.getTest();
//		
//		} catch (InterruptedException e) { //TODO: throw new threading error
//			e.printStackTrace();
//		}
		
//===================universal block=======================
		
//		i.printStats();
//		h.printStats();
//		x.printStats();
//		d.printStats();
//		f.printStats();
		
//		i.printRStats();
//		h.printRStats();
//		x.printRStats();
//		d.printRStats();
//		f.printRStats();
		
		i.logRunStats();
		h.logRunStats();
		x.logRunStats();
		d.logRunStats();
		f.logRunStats();
		
		ws.setIHS(i.getStats(), i.getSNPs());
		ws.setIHH(h.getStats(), h.getSNPs());
		ws.setXPEHH(x.getStats(), x.getSNPs());
		ws.setDAF(d.getStats(), d.getSNPs());
		ws.setFst(f.getStats(), f.getSNPs());
		
		return ws;
	}
	
	private void synchronize(StatsThread i_thrd, 
							StatsThread h_thrd, 
							StatsThread x_thrd, 
							StatsThread d_thrd, 
							StatsThread f_thrd) {
		
		for(;;) {
			if(i_thrd.isFinished()
					&& h_thrd.isFinished()
					&& x_thrd.isFinished()
					&& d_thrd.isFinished()
					&& f_thrd.isFinished())
				break;
			else
				continue;
		}
	}
}

class StatsThread extends Thread {

	private Thread thrd;
	private HaplotypeTests tst;
	
	private boolean finished;
	
	StatsThread(HaplotypeTests tst) {
		this.tst = tst;
		
		finished = false;
		
		thrd = new Thread(this);
//		thrd.start();
	}
	
	public void start() {
		thrd.start();
	}
	
	@Override
	public void run() {
		
		tst.runStat();	
		finished = true;
		
		thrd.interrupt();
	}
	
	public HaplotypeTests getTest() {
		return tst;
	}
	
	public boolean isFinished() {
		return finished;
	}
	
}


