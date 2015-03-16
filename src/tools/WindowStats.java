package tools;

import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
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
	
	public List<SNP> getAllSNPs() {
		
		List<SNP> all_snps = new LinkedList<SNP>();
		
		all_snps = buildAllSNPs(all_snps, ihs_snps);
		all_snps = buildAllSNPs(all_snps, xpehh_snps);
		all_snps = buildAllSNPs(all_snps, ihh_snps);
		all_snps = buildAllSNPs(all_snps, daf_snps);
		all_snps = buildAllSNPs(all_snps, fst_snps);
		
		Collections.sort(all_snps);
			
		return all_snps;
	}
	
	public int getTotNumSNPs() {
		
		List<SNP> all_snps = getAllSNPs();
		
		return all_snps.size();
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
	
	public Double getScore(List<SNP> snps, List<Double> stats, SNP snp) {
		
		for(int i = 0; i < snps.size(); i++) {
			if(snps.get(i).sameAs(snp))
				return stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getIhsScore(SNP snp) {
		
		for(int i = 0; i < ihs_snps.size(); i++) {
			if(ihs_snps.get(i).sameAs(snp))
				return ihs_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getIhhScore(SNP snp) {
		
		for(int i = 0; i < ihh_snps.size(); i++) {
			if(ihh_snps.get(i).sameAs(snp))
				return ihh_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getXpehhScore(SNP snp) {
		
		for(int i = 0; i < xpehh_snps.size(); i++) {
			if(xpehh_snps.get(i).sameAs(snp))
				return xpehh_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getDafScore(SNP snp) {
		
		for(int i = 0; i < daf_snps.size(); i++) {
			if(daf_snps.get(i).sameAs(snp))
				return daf_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	public Double getFstScore(SNP snp) {
		
		for(int i = 0; i < fst_snps.size(); i++) {
			if(fst_snps.get(i).sameAs(snp))
				return fst_stats.get(i);
		}
		
		return Double.NaN;
	}
	
	private List<SNP> buildAllSNPs(List<SNP> all_snps, List<SNP> snps) {
		
		for(int i = 0; i < snps.size(); i++) {
			if(!containsSNP(all_snps, snps.get(i)))
					all_snps.add(snps.get(i));
		}
		
		return all_snps;
	}
	
	private boolean containsSNP(List<SNP> all_snps, SNP snp) {
		
		for(SNP s : all_snps) {
			if(s.sameAs(snp))
				return true;
		}
		
		return false;
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
	
	/*
	 * For testing the output of the different WindowStats objects
	 * Not called at all hence the SuppressWarning
	 */
	@SuppressWarnings("unused")
	private String printLists(List<SNP> snps, List<Double> stats) {
		
		StringBuilder sb = new StringBuilder();
		
		sb.append("SNPS:\t" + snps.size() + "\n");
		for(int i = 0; i < snps.size(); i++) {
			sb.append(snps.get(i) + "\n");
		}
		
		sb.append("Scores:\t" + stats.size() + "\n");
		for(int i = 0; i < stats.size(); i++) {
			sb.append(stats.get(i) + "\n");
		}
		
		return sb.toString();
	}
	
	@Override
	public String toString() {
		
		StringBuilder sb = new StringBuilder();
		
		List<SNP> all_snps = getAllSNPs();
		
		for(int i = 0; i < all_snps.size(); i++) {
			
			SNP cur_snp = all_snps.get(i);

			Double iHS_score = getScore(ihs_snps, ihs_stats, cur_snp);
			Double XPEHH_score = getScore(xpehh_snps, xpehh_stats, cur_snp);
			Double iHH_score = getScore(ihh_snps, ihh_stats, cur_snp);
			Double DAF_score = getScore(daf_snps, daf_stats, cur_snp);
			Double Fst_score = getScore(fst_snps, fst_stats, cur_snp);
			
			sb.append(cur_snp.getSnpID() + "\t");
			sb.append(cur_snp.getPosition() + "\t");
			sb.append(iHS_score + "\t");
			sb.append(XPEHH_score + "\t");
			sb.append(iHH_score + "\t");
			sb.append(DAF_score + "\t");
			sb.append(Fst_score + "\n");
		}
		
		return sb.toString();
	}
}
