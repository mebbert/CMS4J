package cms;

import errors.IllegalInputException;
import log.Log;

public class CMS {

	/** Composite of Multiple Scores (CMS) Java implementation GWAS study version 1.0
	 * @author Hayden Smith
	 * 
	 * For memory management make sure to allocate at least 2x the size of your Phased directory
	 * For optimal speed allocate 4x the size of your Phased directory
	 * 
	 * @param args0		Phased Files Directory; unrelated individual dataset (non-trios)
	 * @param args1		Genetic Map Directory
	 * @param args2		Ancestral States Directory
	 * @param args3		Output Directory Name
	 * @param args4		Target Population (CEU, YRI, JPT, CHB)
	 * @param args5		Cross Population (CEU, YRI, JPT, CHB)
	 * @param args6		Outgroup Population (CEU, YRI, JPT, CHB)
	 * @param args7		Chromosome Declaration (i.e. 1-22)
	 * @param args8		Window Size in Mb
	 */
	public static void main(String[] args) {
		
		Log log = new Log();
		log.addLine("\t\t\t*****Starting CMS*****\n");
		
		try {
			args = setupArgs(args, log);
		
			Driver dv = new Driver(args, log);
			dv.runCMS();
		} catch (Exception e) {
			System.out.println("CMS Died Prematurely." 
					+ " Check log output for troubleshooting.");
			
			log.addLine("\n\nCMS Died Prematurely. Error in computation.");
			
			e.printStackTrace();
		}
		
		cleanup(log);
	}
	
	private static void cleanup(Log log) {
		
		log.close();
	}
	
	private static String[] setupArgs(String[] args, Log log) 
				throws IllegalInputException{
		
		if(args.length != 9) {
			String msg = "Error: Parameter length incorrect";
			throw new IllegalInputException(log, msg);
		}
		
		log.addLine("Working Parameters");
		log.addLine("Phased Dir:\t" + args[0]);
		log.addLine("Map Dir:\t" + args[1]);
		log.addLine("Ancestor Dir:\t" + args[2]);
		log.addLine("Out Dir:\t" + args[3]);
		log.addLine("Target Pop:\t" + args[4]);
		log.addLine("Cross Pop:\t" + args[5]);
		log.addLine("Outgroup Pop:\t" + args[6]);
		log.addLine("Chr Range:\t" + args[7]);
		log.addLine("Window Size:\t" + args[8] + "Mb");
		
//		String[] def_args = createDefaultArgs(log);
//		String[] args_in = args;
//		args = new String[def_args.length + args_in.length];
//		
//		for(int i = 0; i < args_in.length; i++) {
//			args[i] = args_in[i];
//		}
//		for(int i = 0; i < def_args.length; i++) {
//			args[i + args_in.length] = def_args[i];
//		}
		
		return args;
	}
	
//	private static String[] createDefaultArgs(Log log) {
//		String[] def_args = new String[4];
//		
//		def_args[0] = "CEU";
//		def_args[1] = "21-21";
//		def_args[2] = "hg19";
//		def_args[3] = ".5";
//		
//		log.addLine("\nDefault Parameters");
//		log.addLine("Population:\tCEU");
//		log.addLine("Chr Range:\t21-21");
//		log.addLine("Genome Build:\thg19");
//		log.addLine("Window Size:\t0.5Mb");
//		
//		return def_args;
//	}
	
//	This is a test

}
