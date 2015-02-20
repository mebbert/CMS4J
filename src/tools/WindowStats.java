package tools;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class WindowStats {
	
	private int st_pos;
	private int end_pos;
	
	private List<SNP> ihs_snps;
	private List<SNP> xpehh_snps;
	private List<SNP> ihh_snps;
	private List<SNP> daf_snps;
	private List<SNP> fst_snps;
	
	private List<Double> ihs_stats;
	private List<Double> xpehh_stats;
	private List<Double> ihh_stats;
	private List<Double> daf_stats;
	private List<Double> fst_stats;
	
	
	public WindowStats(int st_pos, int end_pos) {
		
		this.st_pos = st_pos;
		this.end_pos = end_pos;
		
		ihs_snps = null;
		xpehh_snps = null;
		ihh_snps = null;
		daf_snps = null;
		fst_snps = null;
		
		ihs_stats = null;
		xpehh_stats = null;
		ihh_stats = null;
		daf_stats = null;
		fst_stats = null;
	}


	public int getTotNumSNPs() {
		
		Set<Integer> all_pos = new HashSet<Integer>();
		
		for(SNP s : ihs_snps)
			all_pos.add(s.getPosition());
		for(SNP s : xpehh_snps)
			all_pos.add(s.getPosition());
		for(SNP s : ihh_snps)
			all_pos.add(s.getPosition());
		for(SNP s : daf_snps)
			all_pos.add(s.getPosition());
		for(SNP s : fst_snps)
			all_pos.add(s.getPosition());
		
		return all_pos.size();
	}
	
	//TODO: test what happens when you try to access a value that is right at the boarder of the end_pos
	public int getNextPosition(int prev_pos) {
		
		int nxt_pos = end_pos;
		
		nxt_pos = comparePositions(nxt_pos, prev_pos, ihs_snps);
		nxt_pos = comparePositions(nxt_pos, prev_pos, xpehh_snps);
		nxt_pos = comparePositions(nxt_pos, prev_pos, ihh_snps);
		nxt_pos = comparePositions(nxt_pos, prev_pos, daf_snps);
		nxt_pos = comparePositions(nxt_pos, prev_pos, fst_snps);
		
		if(nxt_pos == prev_pos) {//check this...
			//at end
			return -1;
		}
		
		return nxt_pos;
	}
	
	public int getStPos() {
		return st_pos;
	}
	
	public int getEndPos() {
		return end_pos;
	}

	public List<SNP> getIHSsnps() {
		return ihs_snps;
	}

	public List<Double> getIHSstats() {
		return ihs_stats;
	}

	public void setIHS(List<Double> ihs_stats, List<SNP> ihs_snps) {
		this.ihs_stats = ihs_stats;
		this.ihs_snps = ihs_snps;
	}

	public List<SNP> getXPEHHsnps() {
		return xpehh_snps;
	}

	public List<Double> getXPEHHstats() {
		return xpehh_stats;
	}

	public void setXPEHH(List<Double> xpehh_stats, List<SNP> xpehh_snps) {
		this.xpehh_stats = xpehh_stats;
		this.xpehh_snps = xpehh_snps;
	}

	public List<SNP> getIHHsnps() {
		return ihh_snps;
	}
	
	public List<Double> getIHHstats() {
		return ihh_stats;
	}

	public void setIHH(List<Double> ihh_stats, List<SNP> ihh_snps) {
		this.ihh_stats = ihh_stats;
		this.ihh_snps = ihh_snps;
	}

	public List<SNP> getDAFsnps() {
		return daf_snps;
	}

	public List<Double> getDAFstats() {
		return daf_stats;
	}

	public void setDAF(List<Double> daf_stats, List<SNP> daf_snps) {
		this.daf_stats = daf_stats;
		this.daf_snps = daf_snps;
	}

	public List<SNP> getFSTsnps() {
		return fst_snps;
	}
	
	public List<Double> getFSTstats() {
		return fst_stats;
	}

	public void setFst(List<Double> fst_stats, List<SNP> fst_snps) {
		this.fst_stats = fst_stats;
		this.fst_snps = fst_snps;
	}
	
	private int comparePositions(int nxt_pos, int prev_pos, List<SNP> snps) {
		
		for(int i = 0; i < snps.size(); i++) {
			SNP s = snps.get(i);
			if(s.getPosition() > prev_pos && s.getPosition() <= nxt_pos) {
				return s.getPosition();
			}
		}
		
		return nxt_pos;
	}

}
