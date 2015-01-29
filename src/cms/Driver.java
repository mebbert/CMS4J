package cms;

import io.AncestralParser;
import io.MapParser;
import io.PhasedParser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import errors.*;
import tools.GeneticMap;
import tools.Individual;
import tools.SNP;
import tools.Window;
import tools.WindowStats;
import log.Log;

public class Driver {
	
	private static int MEGABASE_CONVERSION = 1000000;
	private static String LEGEND_TYPE = "legend";
	private static String MAP_TYPE = "map";
	private static String PHASED_TYPE = "phased";
	
	private int win_size;
	private int chr_st;
	private int chr_end;
	
	private String t_pop;
	private String x_pop;
	private String o_pop;
	
	private File ph_dir;
	private File map_dir;
	private File anc_dir;
	private File out_file;
	
	//target population variables (tp)
	private List<Window> tp_wins;
	private Individual[] tp_indv;
	
	//cross population variables (xp)
	private List<Window> xp_wins;
	private Individual[] xp_indv;
	
	//out-group population variable (op)
	private List<Window> op_wins;
	private Individual[] op_indv;
	
	//intersection of target and cross populations (txin)
	private List<Window> txin_wins;//for listing the windows in the intersection
	private Individual[] tp_inx_indv;//for listing the individual alleles in the intersection
	private Individual[] xp_int_indv;
	
	//intersection of cross population and out-group population (xoin)
	private List<Window> xoin_wins;
	private Individual[] xp_ino_indv;
	private Individual[] op_inx_indv;
	
	//universal variables
	private GeneticMap gm;
	private List<SNP> anc_types;
	private List<WindowStats> win_stats;
	
	private Log log;
	
	public Driver() {
		
		win_size = 0;
		chr_st = -1;
		chr_end = -1;
		
		t_pop = "";
		x_pop = "";
		o_pop = "";
		
		ph_dir = null;
		map_dir = null;
		anc_dir = null;
		
		xp_wins = null;
		tp_wins = null;
		txin_wins = null;
		tp_indv = null;
		xp_indv = null;
		tp_inx_indv = null;
		xp_int_indv = null;
		gm = null;
		anc_types = null;
		
		log = null;
	}
	
	/**
	 * Creates a CMS Driver object that will take care of the basic functions of
	 * running the CMS algorithm. Upon creating this object, the constructor 
	 * will also check the validity of all input to ensure the rest of the 
	 * program runs as expected
	 * 
	 * @param args		Both user defined and default arguments used by CMS in String array
	 * @param log		Logger to be used by all higher level classes
	 */
	public Driver(String[] args, Log log) throws IllegalInputException {
		
		resetDataValues();
		
		this.log = log;
		
		setArgs(args);
		
	}
	
	/**
	 * Main call to starting the CMS pipeline and running the data through 
	 * both population genetic statistic and CMS analysis of those scores one
	 * chromosome at a time
	 * 
	 * @see See CMS api for complete description of CMS calculations and logic
	 */
	public void runCMS() throws Exception {
		
		for(int i = chr_st; i <= chr_end; i++) {
			parseFiles(i);
		
			intersectPopulations();
			
			getStats();
			
			runAnalysis();	//runs CMS
			
			resetDataValues();
			
			close();
		}
	}
	
	private void close() {
		System.out.println("CMS is now ending");
		log.addLine("CMS is now ending");
	}
	
	private void runAnalysis() {
		Analysis an = new Analysis();
		an.runCmsAnalysis();
	}
	
	private void getStats() {
		log.addLine("\n\t\t\t*****Starting Stats Analysis*****");
		
		for(int i = 0; i < tp_wins.size(); i++) {
			Window tp_w = tp_wins.get(i);
			Window txin_w = txin_wins.get(i);//TODO: what do I do when there aren't enough windows here...
			
			log.addLine("\nSweep between positions " + tp_w.getStPos() + " to " + tp_w.getEndPos());
			
			Stats stats = new Stats(log, 
									tp_w, 
									tp_wins, 
									tp_indv,
									
									txin_w,
									txin_wins,
									tp_inx_indv,
									xp_int_indv,
									
									anc_types, 
									gm);
			stats.getStats();
			
//			win_stats.add(stats.getStats());
			
//			=========FOR TESTING===========
			if(i == 2) //breaks at i = 3 on clean_data because of how the intersect is calculated... 
				break;
//			===============================
		}
	}
	
	private void intersectPopulations() {
		
		Individual[] tp_indv_insect = new Individual[tp_indv.length];
		Individual[] xp_indv_insect = new Individual[xp_indv.length];
		
		for(int i = 0; i < tp_indv_insect.length; i++)
			tp_indv_insect[i] = new Individual(tp_indv[i].getID(), tp_indv[i].getChr());
		for(int i = 0; i < xp_indv_insect.length; i++)
			xp_indv_insect[i] = new Individual(xp_indv[i].getID(), xp_indv[i].getChr());
		
		List<Window> wins_insect = new ArrayList<Window>();
		wins_insect.add(new Window(0, 0, 0));
		
		int last_xp_win_index = 0;
		for(int i = 0; i < tp_wins.size(); i++) {
			for(int j = last_xp_win_index; j < xp_wins.size(); j++) {
				Window tp_win = tp_wins.get(i);
				Window xp_win = xp_wins.get(j);
				
				int tp_win_st = tp_win.getStPos();
				int tp_win_end = tp_win.getEndPos();
				int xp_win_st = xp_win.getStPos();
				int xp_win_end = xp_win.getEndPos();
				
				if(tp_win_st == xp_win_st
						&& tp_win_end == xp_win_end) {			
					
					List<SNP> tp_win_snps = tp_wins.get(i).getSNPs();
					List<SNP> xp_win_snps = xp_wins.get(j).getSNPs();
					
					compareSNPs(tp_win_snps, 
								xp_win_snps, 
								wins_insect, 
								tp_win_st, 
								tp_win_end, 
								tp_win, 
								xp_win, 
								tp_indv_insect, 
								xp_indv_insect);	
				}	
			}	
		}
		
		wins_insect.remove(0);//to get rid of the initial window
		
		//set the global variables
		txin_wins = wins_insect;
		tp_inx_indv = tp_indv_insect;
		xp_int_indv = xp_indv_insect;
	}
	
	//continue to pear down this method into smaller chunks
	private void compareSNPs(List<SNP> tp_win_snps, 
								List<SNP> xp_win_snps,
								List<Window> wins_insect,
								int tp_win_st,
								int tp_win_end,
								Window tp_win,
								Window xp_win,
								Individual[] tp_indv_insect,
								Individual[] xp_indv_insect) {
		
		for(int k = 0; k < tp_win_snps.size(); k++) {
			for(int l = 0; l < xp_win_snps.size(); l++) {
				if(tp_win_snps.get(k).getPosition() == xp_win_snps.get(l).getPosition()) {
					if(!containsWindow(wins_insect, tp_win_st, tp_win_end)) {
						
						//make and put a new window window in wins_insect
						Window last_win = wins_insect.get(wins_insect.size() - 1);
						int last_win_indx = last_win.getStIndex() + last_win.getSnpListSize() - 1;
						last_win.setEndIndex(last_win_indx);
						
						Window new_win = new Window(tp_win_st, tp_win_end, (last_win_indx + 1));
						wins_insect.add(new_win);
					}
					
					Window cur_win = getCurWindow(wins_insect, tp_win_st, tp_win_end);
					int cur_win_indx = wins_insect.indexOf(cur_win);
					
					SNP tp_snp = tp_win_snps.get(k);
					SNP xp_snp = xp_win_snps.get(l);
					
					int tp_indx = tp_win.getSnpIndex(tp_snp);
					int xp_indx = xp_win.getSnpIndex(xp_snp);
					
					addAllelesToIndividuals(tp_indx, 
											xp_indx, 
											tp_snp, 
											xp_snp, 
											tp_indv_insect, 
											xp_indv_insect);
					
					cur_win.addSNP(tp_snp);
					cur_win.setEndIndex(cur_win.getStIndex() + cur_win.getSnpListSize() - 1);
					
					wins_insect.set(cur_win_indx, cur_win);
					
				}
			}	
		}
	}
	
	private void addAllelesToIndividuals(int tp_indx, 
											int xp_indx, 
											SNP tp_snp, 
											SNP xp_snp,
											Individual[] tp_indv_insect,
											Individual[] xp_indv_insect) {
		
		//Adding alleles to target population's individuals
		for(int m = 0; m < tp_indv_insect.length; m++) {
			Integer str_1 = tp_indv[m].getStrand1Allele(tp_indx);
			Integer str_2 = tp_indv[m].getStrand2Allele(tp_indx);
			
			tp_indv_insect[m].addAlleleToStrand1(str_1.toString());
			tp_indv_insect[m].addAlleleToStrand2(str_2.toString());
		}
		
		//Adding alleles to cross population's individuals
		for(int i = 0; i < xp_indv_insect.length; i++) {
			
			Integer str_1 = xp_indv[i].getStrand1Allele(xp_indx);
			Integer str_2 = xp_indv[i].getStrand2Allele(xp_indx);
			
			//switch allele types because they are reported on opposite a0 or a1 column
			if(tp_snp.getAllele0().equals(xp_snp.getAllele1())) {
				
				if(str_1 == 0)
					str_1 = 1;
				else
					str_1 = 0;
				
				if(str_2 == 0)
					str_2 = 1;
				else
					str_2 = 0;
			}
			
			xp_indv_insect[i].addAlleleToStrand1(str_1.toString());
			xp_indv_insect[i].addAlleleToStrand2(str_2.toString());
		}
	}
	
	private Window getCurWindow(List<Window> wins, int st, int end) {
		
		for(Window w : wins) {
			if(w.getStPos() == st && w.getEndPos() == end) 
				return w;
		}
		
		return null;
	}
	
	private boolean containsWindow(List<Window> wins, int st, int end) {
		
		
		for(Window w : wins) {
			if(w.getStPos() == st && w.getEndPos() == end)
				return true;
		}
		
		return false;
	}
	
	private void parseFiles(int chr) throws Exception {
		log.addLine("\nLoading referenced data into memory for chromosome " + chr);
		
		String lg_tp_path = getPhasedPath(ph_dir, LEGEND_TYPE, chr, t_pop);
		String ph_tp_path = getPhasedPath(ph_dir, PHASED_TYPE, chr, t_pop);//for target population
		
		String lg_xp_path = getPhasedPath(ph_dir, LEGEND_TYPE, chr, x_pop);
		String ph_xp_path = getPhasedPath(ph_dir, PHASED_TYPE, chr, x_pop);//for cross population
		
		String lg_op_path = getPhasedPath(ph_dir, LEGEND_TYPE, chr, o_pop);
		String ph_op_path = getPhasedPath(ph_dir, PHASED_TYPE, chr, o_pop);
		
		String map_path = getPath(map_dir, MAP_TYPE, chr);
		String anc_path = getPath(anc_dir, LEGEND_TYPE, chr);
		
		PhasedParser tp_pp = new PhasedParser(lg_tp_path, ph_tp_path, chr, log);
		PhasedParser xp_pp = new PhasedParser(lg_xp_path, ph_xp_path, chr, log);
		PhasedParser op_pp = new PhasedParser(lg_op_path, ph_op_path, chr, log);
		MapParser mp = new MapParser(map_path, log);
		AncestralParser ap = new AncestralParser(anc_path, log);
		
		tp_wins = tp_pp.parseLegend(win_size);
		tp_indv = tp_pp.parsePhased(chr);
		
		xp_wins = xp_pp.parseLegend(win_size);
		xp_indv = xp_pp.parsePhased(chr);
		
		op_wins = op_pp.parseLegend(win_size);
		op_indv = op_pp.parsePhased(chr);
		
		gm = mp.parseGeneMap();
		anc_types = ap.parseAncestralTypes();
	}
	
	private String getPath(File dir, String type, int chr) 
			throws UnknownFileException {
		
		String chr_check = "chr" + chr;
		
		String[] all_files = dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains(type) 
					&& file_name.contains(chr_check)
					&& file_name.charAt(0) != '.')
				return dir.getAbsolutePath() + File.separator + file_name;
		}
		
		throw new UnknownFileException(log, dir);
	}
	
	private String getPhasedPath(File dir, String type, int chr, String pop) 
			throws UnknownFileException {
		
		String chr_check = "chr" + chr;
		
		String[] all_files = dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			
			if(file_name.contains(type)
					&& file_name.contains(chr_check)
					&& file_name.contains(pop)
					&& file_name.charAt(0) != '.') {
				return dir.getAbsolutePath() + File.separator + file_name;
			}
		}
		
		throw new UnknownFileException(log, dir);
	}
	
	private void resetDataValues() {
		
		tp_wins = new ArrayList<Window>();
		tp_indv = new Individual[1];
		xp_indv = new Individual[1];
		gm = new GeneticMap();
		anc_types = new ArrayList<SNP>();
		win_stats = new ArrayList<WindowStats>();
	}
	
	private void setArgs(String[] args) throws IllegalInputException {
		
		log.add("\nParameter Check");
		
		ph_dir = new File(args[0]);
		if(!ph_dir.isDirectory()) {
			String msg = "Error: Phased directory path does not exist";
			throw new IllegalInputException(log, msg);
		}
		log.add(".");
		
		map_dir = new File(args[1]);
		if(!map_dir.isDirectory()) {
			String msg = "Error: Map directory path does not exist";
			throw new IllegalInputException(log, msg);
		}
		log.add(".");
		
		anc_dir = new File(args[2]);
		if(!anc_dir.isDirectory()) {
			String msg = "Error: Ancestor directory path does not exist";
			throw new IllegalInputException(log, msg);
		}
		log.add(".");
		
		try {
			out_file = new File(args[3]);
			out_file.createNewFile();
		} catch (IOException e) {
			String msg = "Error: In creating out file; check path";
			throw new IllegalInputException(log, msg);
		}
		log.add(".");
		
		t_pop = args[4];
		if(!t_pop.equals("CEU") && !t_pop.equals("YRI") 
				&& !t_pop.equals("JPT") && !t_pop.equals("CHB")) {
			String msg = "Error: Target population declaration not recognized";
			throw new IllegalInputException(log, msg);
		}
		log.add(".");
		
		x_pop = args[5];
		if(!x_pop.equals("CEU") && !x_pop.equals("YRI") 
				&& !x_pop.equals("JPT") && !x_pop.equals("CHB")) {
			String msg = "Error: Cross population declaration not recognized";
			throw new IllegalInputException(log, msg);
		}
		if(x_pop.equals(t_pop)) {
			String msg = "Error: Cross population declaration cannont be " 
					+ "the same as target population declaration";
			throw new IllegalInputException(log, msg);
		}
		log.add(".");
		
		o_pop = args[6];
		if(!o_pop.equals("CEU") && !o_pop.equals("YRI") 
				&& !o_pop.equals("JPT") && !o_pop.equals("CHB")) {
			String msg = "Error: Out-group population declaration not recognized";
			throw new IllegalInputException(log, msg);
		}
		if(o_pop.equals(t_pop) || o_pop.equals(x_pop)) {
			String msg = "Error: Out-group population declaration cannont be " 
					+ "the same as target population declaration";
			throw new IllegalInputException(log, msg);
		}
		log.add(".");
		
		chr_st = getChrSt(args[7]);
		log.add(".");
		
		chr_end = getChrEnd(args[7]);
		log.add(".");
		
		win_size = getWindowSize(args[8]);
		
		log.addLine(" complete");
	}
	
	private int getWindowSize(String in) throws IllegalInputException {
		
		int in_size = -1;
		
		try {
			double win_size_in = Double.parseDouble(in) * MEGABASE_CONVERSION;
			
			in_size = (int) win_size_in;
			
		} catch (NumberFormatException e) {
			String msg = "Error: Window size invalid format";
			throw new IllegalInputException(log, msg);
		}
		
		if(in_size <= 0 || in_size > (100 * MEGABASE_CONVERSION)) {
			String msg = "Error: Window size declaration invalid";
			throw new IllegalInputException(log, msg);
		}
		
		return in_size;
	}
	
	private int getChrSt(String chr_range) throws IllegalInputException {
		
		String[] st_end = chr_range.split("-");
		
		if(st_end.length != 2) {
			String msg = "Error: Chromosome declaration " + chr_range + " is an invalid format";
			throw new IllegalInputException(log, msg);
		}
		
		int st = -1;
		try {
			st = Integer.parseInt(st_end[0]);
		} catch (NumberFormatException e) {
			String msg = "Error: Chromosome number format is incorrect";
			throw new IllegalInputException(log, msg);
		}
		
		if(st < 1 || st > 22) {
			String msg = "Error: Start chromosome declaration " + st + " is out of bounds";
			throw new IllegalInputException(log, msg);
		}
		
		return st;
	}
	
	private int getChrEnd(String chr_range) throws IllegalInputException {
		
		String[] st_end = chr_range.split("-");
		
		int end = -1;
		try {
			end = Integer.parseInt(st_end[1]);
		} catch (NumberFormatException e) {
			String msg = "Error: Chromosome number format is incorrect";
			throw new IllegalInputException(log, msg);
		}
		
		if(end < 1 || end > 22) {
			String msg = "Error: End chromosome declaration out of bounds";
			throw new IllegalInputException(log, msg);
		}
		
		if(end < chr_st) {
			String msg = "Error: End chromosome and start chromosome invalid order";
			throw new IllegalInputException(log, msg);
		}
		
		return end;
		
	}

}
