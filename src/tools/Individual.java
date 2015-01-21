package tools;

import java.util.Arrays;

public class Individual {
	
	private final int ARRAY_SIZE = 256; // 32 bytes at a time
	
	private int id;
	private int str1_index;
	private int str2_index;
	private byte chr;
	private boolean[] strand1;
	private boolean[] strand2;
	
	public Individual() {
		
		id = 0;
		chr = 0;
		str1_index = 0;
		str2_index = 0;
		strand1 = new boolean[ARRAY_SIZE];
		strand2 = new boolean[ARRAY_SIZE];
	}
	
	public Individual(int id, byte chr) {
		
		this.id = id;
		this.chr = chr;
		str1_index = 0;
		str2_index = 0;
		strand1 = new boolean[ARRAY_SIZE];
		strand2 = new boolean[ARRAY_SIZE];
	}
	
	public boolean addAlleleToStrand1(String a) {
		
		boolean valid = true;
		if(!isValidAllele(a))
			valid = false;
		
		char allele = a.charAt(0);
		
		if((str1_index + 1) == strand1.length)
			strand1 = increaseSize(strand1);
		
		if(!valid)
			strand1[str1_index] = false;
		if(allele == '0')
			strand1[str1_index] = false;
		else if(allele == '1')
			strand1[str1_index] = true;
		
		str1_index++;
		return valid;
	}
	
	public boolean addAlleleToStrand2(String a) {
		
		boolean valid = true;
		if(!isValidAllele(a))
			valid = false;
		
		char allele = a.charAt(0);
		
		if((str2_index + 1) == strand2.length)
			strand2 = increaseSize(strand2);
		
		if(!valid)
			strand2[str2_index] = false;
		else if(allele == '0')
			strand2[str2_index] = false;
		else if(allele == '1')
			strand2[str2_index] = true;
		
		str2_index++;
		return valid;
	}
	
	public int getStrand1Allele(int index) {
		
		boolean a = strand1[index];
		if(a == true)
			return 1;
		else
			return 0;
	}
	
//	public boolean getStrand1Allele(int index) {
//		return strand1[index];
//	}
	
	public int getStrand2Allele(int index) {
		
		boolean a = strand2[index];
		if(a == true)
			return 1;
		else
			return 0;
	}
	
	private boolean[] increaseSize(boolean[] str) {
		boolean[] new_str = new boolean[str.length * 2];
		
		for(int i = 0; i < str.length; i++) 
			new_str[i] = str[i];
		
		return new_str;
	}
	
	private boolean isValidAllele(String a) {
		
		if(a.length() != 1)
			return false;
		
		char ch = a.charAt(0);
		if(!Character.isDigit(ch))
			return false;
		
		if(ch != '1' && ch != '0')
			return false;
		
		return true;
	}
	
	public String toString(int st, int end) throws IndexOutOfBoundsException{
		StringBuilder sb = new StringBuilder();
		sb.append("Individual [id=" + id + ", chr=" + chr + "]\n");
		
		sb.append("\tStrand1 [");
		if(st > strand1.length)
			throw new IndexOutOfBoundsException();
		if(st < 0)
			st = 0;
		if(end > strand1.length)
			end = strand1.length;
		for(int i = st; i < end; i++) 
			sb.append(strand1[i]);
		sb.append("]\n");
		
		sb.append("\tStrand2 [");
		if(st > strand2.length)
			throw new IndexOutOfBoundsException();
		if(end > strand2.length)
			end = strand2.length;
		for(int i = st; i < end; i++) 
			sb.append(strand2[i]);
		sb.append("]\n");
		
		return sb.toString();
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("Individual [id=" + id + ", chr=" + chr + "]\n");
		
		sb.append("\tStrand1 [");
		for(int i = 0; i < str1_index; i++) 
			sb.append(getStrand1Allele(i));
		sb.append("]\n");
		
		sb.append("\tStrand2 [");
		for(int i = 0; i < str2_index; i++) 
			sb.append(getStrand2Allele(i));
		sb.append("]\n");
		
		return sb.toString();
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + chr;
		result = prime * result + id;
		result = prime * result + Arrays.hashCode(strand1);
		result = prime * result + Arrays.hashCode(strand2);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Individual other = (Individual) obj;
		if (id != other.id)
			return false;
		if (chr != other.chr)
			return false;
		if (!Arrays.equals(strand1, other.strand1))
			return false;
		if (!Arrays.equals(strand2, other.strand2))
			return false;
		return true;
	}
}
