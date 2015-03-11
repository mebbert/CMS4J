package errors;

import tools.SNP;
import log.Log;

@SuppressWarnings("serial")
public class EhhComputationException extends Exception {
	
	public EhhComputationException() {}
	
	public EhhComputationException(Log log, SNP core_snp, SNP end_snp) {
		log.addLine("\n");
		log.addLine("Error in computing EHH values at CORE_" + core_snp);
		log.addLine("\t*Found when looking at " + end_snp);
		log.addLine("\t*This cannot be resolved computationally; consider revising data");
	}

}
