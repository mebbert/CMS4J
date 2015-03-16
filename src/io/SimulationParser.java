package io;

import java.util.List;

import tools.SimDist;
import log.Log;

public class SimulationParser {
	
	private final String neut_path = "";
	private final String sel_path = "";
	
	private Log log;
	
	public SimulationParser(Log log) {
		
		this.log = log;
	}
	
	public SimDist[] getNeutralSimulations(int num_tests) {
		
		//goes through the neutral file and gets all the bin values
		//puts the values into the SimDist obj
		//does that for all 5 tests
		
		return null;
	}
	
	public SimDist[] getSelectedSimulations(int num_tests) {
		
		return null;
	}

}
