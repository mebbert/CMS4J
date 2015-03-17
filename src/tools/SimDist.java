package tools;

import java.util.ArrayList;
import java.util.List;

public class SimDist {
	
	private List<Double> sim_vals;
	
	public SimDist() {
		
		sim_vals = new ArrayList<Double>();
	}

	public void addSimValue(double val) {
		
		sim_vals.add(val);
	}
	
	public Double getProb(Double score) {
		
		//TODO: This.
		
		
		return null;
	}
	
	public List<Double> getSimVals() {
		return sim_vals;
	}
	
}
