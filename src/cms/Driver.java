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
	private static String EMF_TYPE = "emf";
	private static String MAP_TYPE = "map";
	private static String PHASED_TYPE = "phased";
	
	private int win_size;
	private int chr_st;
	private int chr_end;
	
	//population declarations
	private String t_pop;
	private String x_pop;
	private String o_pop;
	
	//directories and files for accessing and writing data
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
	private List<Window> anc_types;
	private List<WindowStats> win_stats;
	
	//log for progress and error checking
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
		op_wins = null;
		
		txin_wins = null;
		xoin_wins = null;
		
		tp_indv = null;
		xp_indv = null;
		op_indv = null;
		
		tp_inx_indv = null;
		xp_int_indv = null;
		xp_ino_indv = null;
		op_inx_indv = null;
		
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
	
	/*
	 * Runs the CMS analysis of the 5 population statistics
	 * Not quite sure what exactly this does yet. Still working on that part
	 */
	private void runAnalysis() {
		
		System.out.println("position\tiHS\tXPEHH\tiHH\tDAF\tFst");
		for(int i = 0; i < win_stats.size(); i++)
			System.out.print(win_stats.get(i));
		
		Analysis an = new Analysis();
		an.runCmsAnalysis(win_stats);
	}
	
	/*
	 * Loops through all the windows in a given chromosome calculating all the
	 * populations statistics necessary for CMS calculation. Requires there to 
	 * be both a target population window and a intersecting window to overlap
	 * otherwise it skips that particular window and goes to the next.
	 */
	private void getStats() {
		log.addLine("\n\t\t\t*****Starting Stats Analysis*****");
		
		System.out.println("Population Stats Progress:");
		int progress = 0;
		for(int i = 0; i < tp_wins.size() && i < txin_wins.size(); i++) {
			Window tp_w = tp_wins.get(i);
			Window txin_w = txin_wins.get(i);
			
			log.addLine("\nSweep between positions " + tp_w.getStPos() + " to " + tp_w.getEndPos());
			
			Stats stats = new Stats(log, 
									tp_w, 
									tp_wins, 
									tp_indv,
									
									txin_w,
									txin_wins,
									tp_inx_indv,
									xp_int_indv,
									
									xoin_wins,
									xp_ino_indv,
									op_inx_indv,
									
									anc_types, 
									gm);
			
			win_stats.add(stats.getStats());
			
			double percent = ((double)(win_stats.size()) / tp_wins.size()) * 1000;
//			int progress_dif = progress - (int) percent;
//			for(int j = 0; j < progress_dif; j++)
//				System.out.print("=");
			progress = (int) percent;
			System.out.print(progress + "% ");
		}
		System.out.println("100%");
	}
	
	/*
	 * In order to eliminate the variation that comes from population to population
	 * the intersection of each population needs to be made. This series of 
	 * functions gets rid of the SNPs not found in the other populations. It is
	 * important to note that all populations are intersected with the Cross
	 * Population (xp) and therefore all relate back to that .legend file setup.
	 * This is important for statistic calculations, particularly DAF and Fst.
	 */
	private void intersectPopulations() {
		
		intersectCrossWithTargetPopulations();
		
		intersectCrossWithOutgroupPopulations();
	}
	
	private void intersectCrossWithOutgroupPopulations() {
		
		Individual[] xp_indv_insect = new Individual[xp_indv.length];
		Individual[] op_indv_insect = new Individual[op_indv.length];
		
		for(int i = 0; i < xp_indv_insect.length; i++)
			xp_indv_insect[i] = new Individual(xp_indv[i].getID(), xp_indv[i].getChr());
		for(int i = 0; i < op_indv_insect.length; i++)
			op_indv_insect[i] = new Individual(op_indv[i].getID(), op_indv[i].getChr());
		
		List<Window> wins_insect = new ArrayList<Window>();
		
		compareWindows(wins_insect, xp_wins, xp_indv, xp_indv_insect, op_wins, op_indv, op_indv_insect);

		//set the global variables
		xoin_wins = wins_insect;
		xp_ino_indv = xp_indv_insect;
		op_inx_indv = op_indv_insect;
	}
	
	private void intersectCrossWithTargetPopulations() {
		
		Individual[] xp_indv_insect = new Individual[xp_indv.length];
		Individual[] tp_indv_insect = new Individual[tp_indv.length];
		
		for(int i = 0; i < xp_indv_insect.length; i++)
			xp_indv_insect[i] = new Individual(xp_indv[i].getID(), xp_indv[i].getChr());
		for(int i = 0; i < tp_indv_insect.length; i++)
			tp_indv_insect[i] = new Individual(tp_indv[i].getID(), tp_indv[i].getChr());
		
		List<Window> wins_insect = new ArrayList<Window>();
		
		compareWindows(wins_insect, xp_wins, xp_indv, xp_indv_insect, tp_wins, tp_indv, tp_indv_insect);

		//set the global variables
		txin_wins = wins_insect;
		xp_int_indv = xp_indv_insect;
		tp_inx_indv = tp_indv_insect;
	}
	
	private void compareWindows(List<Window> wins_insect,
								List<Window> p1_wins,
								Individual[] p1_indv,
								Individual[] p1_indv_insect,
								List<Window> p2_wins,
								Individual[] p2_indv,
								Individual[] p2_indv_insect) {
		
		wins_insect.add(new Window(0, 0, 0));
		
		for(int i = 0; i < p1_wins.size(); i++) {
			for(int j = 0; j < p2_wins.size(); j++) {
				Window p1_win = p1_wins.get(i);
				Window p2_win = p2_wins.get(j);
				
				int p1_win_st = p1_win.getStPos();
				int p1_win_end = p1_win.getEndPos();
				int p2_win_st = p2_win.getStPos();
				int p2_win_end = p2_win.getEndPos();
				
				if(p1_win_st == p2_win_st
						&& p1_win_end == p2_win_end) {			
					
					List<SNP> p1_win_snps = p1_wins.get(i).getSNPs();
					List<SNP> p2_win_snps = p2_wins.get(j).getSNPs();
					
					compareSNPs(p1_win_snps, 
								p2_win_snps, 
								wins_insect, 
								p1_win_st, 
								p2_win_end, 
								p1_win, 
								p2_win, 
								p1_indv,
								p2_indv,
								p1_indv_insect, 
								p2_indv_insect);	
				}
			}	
		}
		
		wins_insect.remove(0);//to get rid of the initial window
	}
	
	//continue to pear down this method into smaller chunks
	private void compareSNPs(List<SNP> p1_win_snps, 
								List<SNP> p2_win_snps,
								List<Window> wins_insect,
								int p1_win_st,
								int p1_win_end,
								Window p1_win,
								Window p2_win,
								Individual[] p1_indv,
								Individual[] p2_indv,
								Individual[] p1_indv_insect,
								Individual[] p2_indv_insect) {
		
		for(int k = 0; k < p1_win_snps.size(); k++) {
			for(int l = 0; l < p2_win_snps.size(); l++) {
				if(p1_win_snps.get(k).sameAs(p2_win_snps.get(l))) { 
					if(!containsWindow(wins_insect, p1_win_st, p1_win_end)) {
						
						//make and put a new window window in wins_insect
						Window last_win = wins_insect.get(wins_insect.size() - 1);
						int last_win_indx = last_win.getStIndex() + last_win.getSnpListSize() - 1;
						last_win.setEndIndex(last_win_indx);
						
						Window new_win = new Window(p1_win_st, p1_win_end, (last_win_indx + 1));
						wins_insect.add(new_win);
					}
					
					Window cur_win = getCurWindow(wins_insect, p1_win_st, p1_win_end);
					int cur_win_indx = wins_insect.indexOf(cur_win);
					
					SNP p1_snp = p1_win_snps.get(k);
					SNP p2_snp = p2_win_snps.get(l);
					
					int p1_indx = p1_win.getSnpIndex(p1_snp);
					int p2_indx = p2_win.getSnpIndex(p2_snp);
					
					addAllelesToIndividuals(p1_indx, 
											p2_indx, 
											p1_snp, 
											p2_snp, 
											p1_indv,
											p2_indv,
											p1_indv_insect, 
											p2_indv_insect);
					
					cur_win.addSNP(p1_snp);
					cur_win.setEndIndex(cur_win.getStIndex() + cur_win.getSnpListSize() - 1);
					
					wins_insect.set(cur_win_indx, cur_win);
					
				}
			}	
		}
	}
	
	private void addAllelesToIndividuals(int p1_indx, 
											int p2_indx, 
											SNP p1_snp, 
											SNP p2_snp,
											Individual[] p1_indv,
											Individual[] p2_indv,
											Individual[] p1_indv_insect,
											Individual[] p2_indv_insect) {
		
		//Adding alleles to p1 population's individuals
		for(int m = 0; m < p1_indv_insect.length; m++) {
			Integer str_1 = p1_indv[m].getStrand1Allele(p1_indx);
			Integer str_2 = p1_indv[m].getStrand2Allele(p1_indx);
			
			p1_indv_insect[m].addAlleleToStrand1(str_1.toString());
			p1_indv_insect[m].addAlleleToStrand2(str_2.toString());
		}
		
		//Adding alleles to p2 population's individuals
		for(int i = 0; i < p2_indv_insect.length; i++) {
			
			Integer str_1 = p2_indv[i].getStrand1Allele(p2_indx);
			Integer str_2 = p2_indv[i].getStrand2Allele(p2_indx);
			
			//switch allele types because they are reported on opposite a0 or a1 column
			if(p1_snp.getAllele0().equals(p2_snp.getAllele1())) {
				
				if(str_1 == 0)
					str_1 = 1;
				else
					str_1 = 0;
				
				if(str_2 == 0)
					str_2 = 1;
				else
					str_2 = 0;
			}
			
			p2_indv_insect[i].addAlleleToStrand1(str_1.toString());
			p2_indv_insect[i].addAlleleToStrand2(str_2.toString());
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
	
	/*
	 * Used to parse the .legend .phased and genetic map files for Phased, Map,
	 * and Ancestral data types. In order for it to work properly the following
	 * must be in you file system
	 * 		-Phased directory path passed as argument
	 * 			-Phased data for 3 different populations (both .legend and .phased)
	 * 			-Population type (CEU, YRI, JPT, or CHB) must be in the file name
	 * 			-File type ("phased" or "legend") must be in the file name
	 * 			-Chromosome ("chr" + [1-22]) must be in the file name
	 * 		-Ancestral directory path passed as argument
	 * 			-File type of "legend" must be in the file name
	 * 			-Chromosome ("chr" + [1-22]) must be in the file name 
	 * 		-Map directory path passed as argument
	 * 			-Chromosome ("chr" + [1-22]) must be in the file name
	 */
	private void parseFiles(int chr) throws Exception {
		log.addLine("\nLoading referenced data into memory for chromosome " + chr);
		
		String lg_tp_path = getPhasedPath(ph_dir, LEGEND_TYPE, chr, t_pop);
		String ph_tp_path = getPhasedPath(ph_dir, PHASED_TYPE, chr, t_pop);//for target population
		
		String lg_xp_path = getPhasedPath(ph_dir, LEGEND_TYPE, chr, x_pop);
		String ph_xp_path = getPhasedPath(ph_dir, PHASED_TYPE, chr, x_pop);//for cross population
		
		String lg_op_path = getPhasedPath(ph_dir, LEGEND_TYPE, chr, o_pop);
		String ph_op_path = getPhasedPath(ph_dir, PHASED_TYPE, chr, o_pop);
		
		String map_path = getMapPath(map_dir, chr);
		String anc_path = getAncestralPath(anc_dir, chr);
		
		PhasedParser tp_pp = new PhasedParser(lg_tp_path, ph_tp_path, chr, log);
		PhasedParser xp_pp = new PhasedParser(lg_xp_path, ph_xp_path, chr, log);
		PhasedParser op_pp = new PhasedParser(lg_op_path, ph_op_path, chr, log);
		MapParser mp = new MapParser(map_path, log);
		AncestralParser ap = new AncestralParser(anc_path, chr, log);
		
		tp_wins = tp_pp.parseLegend(win_size);
		tp_indv = tp_pp.parsePhased(chr);
		
		xp_wins = xp_pp.parseLegend(win_size);
		xp_indv = xp_pp.parsePhased(chr);
		
		op_wins = op_pp.parseLegend(win_size);
		op_indv = op_pp.parsePhased(chr);
		
		gm = mp.parseGeneMap();
		anc_types = ap.parseAncestralTypes(out_file);
	}
	
	private String getAncestralPath(File dir, int chr) 
			throws UnknownFileException {
		
		String chr_check = "chr" + chr;
		
		String[] all_files = dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains(LEGEND_TYPE) 
					&& file_name.contains(chr_check)
					&& file_name.charAt(0) != '.')
				return dir.getAbsolutePath() + File.separator + file_name;
			if(file_name.contains(EMF_TYPE) 
					&& file_name.contains(chr_check)
					&& file_name.charAt(0) != '.')
				return dir.getAbsolutePath() + File.separator + file_name;
		}
		
		throw new UnknownFileException(log, dir);
	}
	
	private String getMapPath(File dir, int chr) 
			throws UnknownFileException {
		
		String chr_check = "chr" + chr;
		
		String[] all_files = dir.list();
		for(int i = 0; i < all_files.length; i++) {
			
			String file_name = all_files[i];
			if(file_name.contains(MAP_TYPE) 
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
		anc_types = new ArrayList<Window>();
		win_stats = new ArrayList<WindowStats>();
	}
	
	/*
	 * Runs a check on all the parameters to make sure they are exactly what is
	 * needed to run CMS4J. Check Log output for specific error report. Once 
	 * one error is found it breaks out of function and closes CMS.
	 */
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
