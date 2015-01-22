package io;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

import errors.FileParsingException;
import tools.SNP;
import log.Log;

public class AncestralParser {
	
	private static int TEST_LINES = 10; //for testing the first 10 lines of a file to ensure it is the right format
	
	private boolean skip_first_line;
	
	private String anc_path;
	private Scanner anc_scan;
	private Log log;

	public AncestralParser(String anc_path, Log log) throws FileParsingException {
		
		this.anc_path = anc_path;
		this.log = log;
		
		skip_first_line = false;
		
		try {
			anc_scan = new Scanner(new File(anc_path));
		} catch (FileNotFoundException e) {
			String msg = "Error: Legend File not found in the /Ancestral directory";
			throw new FileParsingException(log, msg);
		}
	}
	
	/**
	 * For use on ancestral legend files; only the first four columns are parsed
	 * Ancestral type is defined as a0 allele
	 * 
	 * @return		Returns a list of SNP objects to be referenced later on when calculating population genetic statistics
	 */
	public List<SNP> parseAncestralTypes() throws FileParsingException {
		log.addLine("Importing Ancestral data from " + anc_path);
		
		List<SNP> anc_types = new ArrayList<SNP>();
		
		checkLegendFile();
		
		if(skip_first_line)
			anc_scan.nextLine(); //skips the first line header
		
		while(anc_scan.hasNextLine()) {
			
			String[] ln_arr = anc_scan.nextLine().split("\\s+");
			
			rsNum rs = new rsNum(log, ln_arr[0]);
			int pos = Integer.parseInt(ln_arr[1]);
			anc_types.add(new SNP(pos, rs.getNum(), ln_arr[2], ln_arr[3], rs.getModifyer()));
		}
		
//		for(SNP s : anc_types) {
//			System.out.println(s.toString());
//			log.addLine(s.toString());
//		}
		
		return anc_types;
	}
	
	private void checkLegendFile() throws FileParsingException {
		
		Scanner temp_scan = null;
		try {
			temp_scan = new Scanner(new File(anc_path));
			
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
					String msg = "Error: Legend file " + anc_path + " has invalid number of columns";
					throw new FileParsingException(log, msg);
				}
				
				Integer.parseInt(line_arr[1]);//check for position in position column
			}
			
		} catch (FileNotFoundException e) {
			String msg = "Error: Legend File not found in the /Ancestral directory";
			throw new FileParsingException(log, msg);
		} catch (NumberFormatException e) {
			String msg = "Error: Legend File " + anc_path + " has incorrect colum formatting";
			throw new FileParsingException(log, msg);
		}
	}
}
