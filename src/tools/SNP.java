package tools;

public class SNP {
	
	private int pos;
	private int rs_num;
	
	private String a0;
	private String a1;
	private String rs_modifyer;
	
	public SNP() {
		pos = -1;
		rs_num = -1;
		a0 = "";
		a1 = "";
		rs_modifyer = "";
	}
	
	public SNP(int pos, int rs_num, String a0, String a1, String rs_modifyer) {

		this.pos = pos;
		this.rs_num = rs_num;
		this.a0 = a0;
		this.a1 = a1;
		this.rs_modifyer = rs_modifyer;
	}
	
	public boolean equalTo(SNP s) {
		
		if(pos != s.getPosition())
			return false;
		if(rs_num != s.getRsNum())
			return false;
		if(rs_modifyer.equals(s.getRsModifyer()))
			return true;
		
		return false;
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
	
	public int getRsNum() {
		return rs_num;
	}
	
	public String getRsModifyer() {
		return rs_modifyer;
	}

	@Override
	public String toString() {
		return "SNP [pos=" + pos + ", a0=" + a0 + ", a1=" + a1 + ", rs_number=" 
				+ rs_modifyer + rs_num + "]";
	}
}
