package io;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;
import java.util.Scanner;

import errors.FileParsingException;
import tools.SimDist;
import log.Log;

public class SimulationParser {
	
	private final int NUM_TESTS = 5;
	private final String neut_path = "sim_data" + File.separator + "neutral_simulation.tsv";
	private final String sel_path = "sim_data" + File.separator + "selection_simulation.tsv";
	
	private Log log;
	
	public SimulationParser(Log log) {
		
		this.log = log;
	}
	
	public SimDist[] getNeutralSimulations() throws FileParsingException {
		
		return parseSimulatedData(neut_path);
		
	}
	
	public SimDist[] getSelectedSimulations() throws FileParsingException {
		
		return parseSimulatedData(sel_path);
	}
	
	/*
	 * Simulations and scores need to be stored into arrays where:
	 * 		[0] = iHS data
	 * 		[1] = iHH data
	 * 		[2] = Fst data
	 * 		[3] = DAF data
	 * 		[4] = XPEHH data
	 */
	private SimDist[] parseSimulatedData(String file_path) throws FileParsingException {
		
		SimDist[] dists = new SimDist[NUM_TESTS];
		for(int i = 0; i < NUM_TESTS; i++) 
			dists[i] = new SimDist();
		
		try {
			Scanner scan = new Scanner(new File(file_path));
			scan.nextLine();//skip header line
			
			while(scan.hasNext()) {
				
				String[] vals = scan.nextLine().split("\\s+");
				
				if(vals.length != NUM_TESTS) {
					String msg = "Error: Irregularities in number of simulated data columns";
					throw new FileParsingException(log, msg);
				}
				
				for(int i = 0; i < NUM_TESTS; i++) {
					double val = Double.parseDouble(vals[i]);
					dists[i].addSimValue(val);
				}
			}
		
		} catch (FileNotFoundException e) {
			String msg = "Error: Could not find simulated data file";
			throw new FileParsingException(log, msg);
		} catch (NumberFormatException e) {
			String msg = "Error: Invalid simulated data values";
			throw new FileParsingException(log, msg);
		}
		
		return dists;
	}
	
	private void printSimDistArr(SimDist[] dists) {
		
		for(int i = 0; i < NUM_TESTS; i++) {
			List<Double> vals = dists[i].getSimVals();
			for(int j = 0; j < vals.size(); j++)
				System.out.print(vals.get(j) + "  ");
			
			System.out.println();
		}
	}

}
