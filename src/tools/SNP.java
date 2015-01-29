package tools;

public class SNP {
	
	private int pos;
	
	private String snp_id;
	private String a0;
	private String a1;
	
	
	public SNP() {
		pos = -1;
		a0 = "";
		a1 = "";
		snp_id = "";
	}
	
	public SNP(int pos, String a0, String a1, String snp_id) {

		this.pos = pos;
		this.a0 = a0;
		this.a1 = a1;
		this.snp_id = snp_id;
	}
	
	public String getAllele0() {
		return a0;
	}
	
	public String getAllele1() {
		return a1;
	}
	
	public int getPosition() {
		return pos;
	}

	@Override
	public String toString() {
		return "SNP [pos=" + pos + ", snp_id=" + snp_id + ", a0=" + a0 
				+ ", a1=" + a1 + "]";
	}
}
