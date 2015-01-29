package io;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import errors.FileParsingException;
import log.Log;
import tools.Individual;
import tools.Window;

public class PhasedParser {
	
	private static int TEST_LINES = 10; //for testing the first 10 lines of a file to ensure it is the right format
	
	private boolean skip_first_line;
	private int chr;
	
	private String lg_path;
	private String ph_path;
	
	private Scanner lg_scan;
	private Scanner ph_scan;
	
	private Log log;
	
	public PhasedParser() {
		
		lg_scan = null;
		ph_scan = null;
		
		log = null;
	}
	
	public PhasedParser(String lg_path, String ph_path, int chr, Log log) throws FileParsingException {
		
		this.lg_path = lg_path;
		this.ph_path = ph_path;
		this.chr = chr;
		this.log = log;
		
		skip_first_line = false;
		
		try {
			lg_scan = new Scanner(new File(lg_path));
			ph_scan = new Scanner(new File(ph_path));
		} catch (FileNotFoundException e) {
			String msg = "File not found in the /Phased directory";
			throw new FileParsingException(log, msg);
		}
	}
	
	/**
	 * For use on phased legend files only where allele types can be both SNPs
	 * and indels
	 * 
	 * @param win_size		This size of the analysis window on chromosome
	 * @return				Returns all the windows as defined by _legend file
	 */
	public List<Window> parseLegend(int win_size) throws FileParsingException {
		log.addLine("Importing legend data from " + lg_path);
		
		checkLegendFile();
		
		List<Window> all_win = new ArrayList<Window>();
		
		String line = "";
		if(skip_first_line)
			line = lg_scan.nextLine(); //skips the first line header
		
		int st_pos = 0;
		int end_pos = win_size;
		int rs_pos = 0;
		int index = 0;
		Window cur_win = new Window();
		
		do{
			if(lg_scan.hasNextLine())
				line = lg_scan.nextLine();
			else {
				cur_win.setEndIndex(index - 1);
				all_win.add(cur_win);
				break;
			}
			
			String[] line_arr = line.split("\\s+");
			rs_pos = Integer.parseInt(line_arr[1]);
			
			if(rs_pos > (st_pos + win_size)) {
				if(!cur_win.equals(new Window()))
					all_win.add(cur_win);
				
				while (rs_pos > (st_pos + win_size)) {
					st_pos += win_size;
					end_pos += win_size;
				}
				cur_win.setEndIndex(index - 1);
				cur_win = new Window(st_pos, end_pos, index);
			}
			
			cur_win.addSNP(rs_pos, line_arr[2], line_arr[3], line_arr[0]);
			index++;
				
		} while(st_pos <= rs_pos);
		
//		for(int i = 0; i < all_win.size(); i++) {
//			System.out.println(all_win.get(i).toString());
//			log.addLine(all_win.get(i).toString());
//		}
		
		return all_win;
	}
	
	/**
	 * For used on phased files where allele types are binary: a0 or a1
	 * 
	 * @param chr_in		Chromosme where phased data was taken from
	 * @return				Returns a set of Individual containing all individuals with phased data for each strand and unique id
	 */
	public Individual[] parsePhased(int chr_in) throws FileParsingException {
		log.addLine("Importing phased data from " + ph_path);
		
		checkPhasedFile();
		
		ArrayList<Individual> all_indv = new ArrayList<Individual>();
		
		byte chr = (byte) chr_in;
		int id_index = 0;
		while(ph_scan.hasNextLine()) {
			
			Individual indv = new Individual(id_index, chr);
			String strand1 = ph_scan.nextLine();
			String strand2 = ph_scan.nextLine();
			
			if(strand1.length() == strand2.length()) {
				Scanner str1 = new Scanner(strand1);
				Scanner str2 = new Scanner(strand2);
				
				boolean check = false;
				while(str1.hasNext()) {
					check = indv.addAlleleToStrand1(str1.next());
					if(!check) {
						String msg = "Error: Phased file " + ph_path
								+ " has an problem with an allele on line "
								+ (((1 + id_index)*2) - 1);
						throw new FileParsingException(log, msg);
					}
				}
				while(str2.hasNext()) {
					check = indv.addAlleleToStrand2(str2.next());
					if(!check) {
						String msg = "Error: Phased file " + ph_path
								+ " has an problem with an allele on line "
								+ (1 + id_index)*2;
						throw new FileParsingException(log, msg);
					}
				}
			} 
			else {
				String msg = "Error: Phased file " + ph_path 
						+ " has unequal Individual strand lengths at lines " 
						+ (((1 + id_index)*2) - 1) + " and " + (1 + id_index)*2
						+"\n\t*Make sure there isn't excess white space at the end of the line";
				throw new FileParsingException(log, msg);
			}
			
			id_index++;
			all_indv.add(indv);
		}
		
		Individual[] all_indv_arr = new Individual[all_indv.size()];
		for(int i = 0; i < all_indv.size(); i++) {
			all_indv_arr[i] = all_indv.get(i);
		}
		
		//FOR PRINTING ONLY
//		for(int i = 0; i < all_indv_arr.length; i++) {
//			System.out.println(i + "\n" + all_indv_arr[i]);
//			log.addLine(i.toString());
//		}
		
		return all_indv_arr;
	}
	
	private void checkLegendFile() throws FileParsingException {
		
		Scanner temp_scan = null;
		try {
			temp_scan = new Scanner(new File(lg_path));
			
			//to check the first line; to skip or not to skip, that is the question
			String first_line = temp_scan.nextLine();
			String[] first_line_arr = first_line.split("\\s+");
			if(first_line_arr[1].contains("pos"))
				skip_first_line = true;
			if(first_line_arr[0].contains("rs") && first_line_arr[0].length() < 3)
				skip_first_line = true;
			
			for(int i = 0; i < TEST_LINES; i++) {
				
				String line = temp_scan.nextLine();
				String[] line_arr = line.split("\\s+");
				
				if(line_arr.length < 4 || line_arr.length > 7) {
					String msg = "Error: Legend file " + lg_path + " has invalid number of columns";
					throw new FileParsingException(log, msg);
				}
				
				Integer.parseInt(line_arr[1]);//check for position in position column
			}
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Legend File not found in the /Phased directory";
			throw new FileParsingException(log, msg);
		} catch (NumberFormatException e) {
			String msg = "Error: Legend File " + lg_path + " has incorrect colum formatting";
			throw new FileParsingException(log, msg);
		}
	}
	
	private void checkPhasedFile() throws FileParsingException {
		
		Scanner temp_scan = null;
		try {
			temp_scan = new Scanner(new File(ph_path));
			
			int line_length = 0;
			for(int i = 0; i < TEST_LINES; i++) {
				String line = temp_scan.nextLine();
				String[] line_arr = line.split("\\s+");
				
				if(i != 0 && line_length != line_arr.length) {
					String msg = "Error: Phased file " + ph_path + " line length not equal";
					throw new FileParsingException(log, msg);
				}
				
				for(int j = 0; j < line_arr.length; j++) {
					if(!line_arr[j].equals("0") && !line_arr[j].equals("1")) {
						String msg = "Error: Phased file " + ph_path + "is incorrect format, see line " 
								+ (i + 1) + " position " + (j + 1);
						throw new FileParsingException(log, msg);
					}
				}
				
				line_length = line_arr.length;
			}
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Phased File not found in the /Phased directory";
			throw new FileParsingException(log, msg);
		}
		
	}
}


//------------------------------------------------------------------------------
//---------------------------Package-only class---------------------------------	
//------------------------------------------------------------------------------
//class rsNum {
//	private int num;
//	private String modifyer;
//	
//	public rsNum(Log log, String rs) {
//		
//		StringBuilder sb_mod = new StringBuilder();
//		StringBuilder sb_num = new StringBuilder();
//		
//		boolean is_rs_modifyer = true;
//		for(int i = 0; i < rs.length(); i++) {
//			
//			char index = rs.charAt(i);
//			if(!Character.isDigit(index)) {
//				is_rs_modifyer = false;
//				sb_mod.append(index);
//			}
//			else if(is_rs_modifyer) {
//				sb_mod.append(index);
//			}
//			else {
//				sb_num.append(index);
//			}
//		}
//		
//		num = Integer.parseInt(sb_num.toString());
//		modifyer = sb_mod.toString();
//		
//	}
//	
//	public rsNum(Log log, int chr, int line_num, int pos, String rs) {
//		
//		StringBuilder sb_mod = new StringBuilder();
//		StringBuilder sb_num = new StringBuilder();
//		
//		boolean is_rs_modifyer = true;
//		for(int i = 0; i < rs.length(); i++) {
//			
//			char index = rs.charAt(i);
//			if(!Character.isDigit(index)) {
//				is_rs_modifyer = false;
//				sb_mod.append(index);
//			}
//			else if(is_rs_modifyer) {
//				sb_mod.append(index);
//			}
//			else {
//				sb_num.append(index);
//			}
//		}
//		
//		try {
//			
//			num = Integer.parseInt(sb_num.toString());
//			modifyer = sb_mod.toString();
//			
//		} catch (NumberFormatException e) {
//			num = pos;
//			modifyer = chr + ":";
//			log.add("\tWARNING: changing SNP id around line " + line_num + " to " + modifyer + num);
//			log.addLine(" for computational purposes only");
//			log.addLine("\t\t*original id " + rs + " is not valid");
//		}
//	}
//	
//	public int getNum() {
//		return num;
//	}
//	
//	public String getModifyer() {
//		return modifyer;
//	}
//}
//------------------------------------------------------------------------------	
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

