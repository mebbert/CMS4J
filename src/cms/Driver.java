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
	
	private String pop;
	private String x_pop;
	private File ph_dir;
	private File map_dir;
	private File anc_dir;
	private File out_file;
	
	private List<Window> windows;
	private List<Window> xp_windows;
	private Individual[] individuals;
	private Individual[] xp_individuals;
	private GeneticMap gm;
	private List<SNP> anc_types;
	private List<WindowStats> win_stats;
	
	private Log log;
	
	public Driver() {
		
		win_size = 0;
		chr_st = -1;
		chr_end = -1;
		
		pop = "";
		x_pop = "";
		ph_dir = null;
		map_dir = null;
		anc_dir = null;
		
		xp_windows = null;
		windows = null;
		individuals = null;
		xp_individuals = null;
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
		
		for(int i = 0; i < windows.size(); i++) {
			Window w = windows.get(i);
			log.addLine("\nSweep between positions " + w.getStPos() + " to " + w.getEndPos());
			
			Stats stats = new Stats(log, w, individuals, xp_individuals, anc_types, windows, xp_windows, gm);
			stats.getStats();
			
//			win_stats.add(stats.getStats());
			
//			=========FOR TESTING===========
			if(i == 3)
				break;
//			===============================
		}
	}
	
	private void parseFiles(int chr) throws Exception {
		log.addLine("\nLoading referenced data into memory for chromosome " + chr);
		
		String lg_path = getPhasedPath(ph_dir, LEGEND_TYPE, chr, pop);
		String ph_tp_path = getPhasedPath(ph_dir, PHASED_TYPE, chr, pop);//for target population
		String lg_xp_path = getPhasedPath(ph_dir, LEGEND_TYPE, chr, x_pop);
		String ph_xp_path = getPhasedPath(ph_dir, PHASED_TYPE, chr, x_pop);//for cross population
		String map_path = getPath(map_dir, MAP_TYPE, chr);
		String anc_path = getPath(anc_dir, LEGEND_TYPE, chr);
		
		PhasedParser pp = new PhasedParser(lg_path, ph_tp_path, chr, log);
		PhasedParser xp_pp = new PhasedParser(lg_xp_path, ph_xp_path, chr, log);
		MapParser mp = new MapParser(map_path, log);
		AncestralParser ap = new AncestralParser(anc_path, log);
		
		windows = pp.parseLegend(win_size);
		individuals = pp.parsePhased(chr);
		xp_windows = xp_pp.parseLegend(win_size);
		xp_individuals = xp_pp.parsePhased(chr);
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
		
		windows = new ArrayList<Window>();
		individuals = new Individual[1];
		xp_individuals = new Individual[1];
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
		
		pop = args[4];
		if(!pop.equals("CEU") && !pop.equals("YRI") 
				&& !pop.equals("JPT") && !pop.equals("CHB")) {
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
		if(x_pop.equals(pop)) {
			String msg = "Error: Cross population declaration cannont be " 
					+ "the same as target population declaration";
			throw new IllegalInputException(log, msg);
		}
			
		log.add(".");
		
		chr_st = getChrSt(args[6]);
		log.add(".");
		
		chr_end = getChrEnd(args[6]);
		log.add(".");
		
		win_size = getWindowSize(args[7]);
		
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
