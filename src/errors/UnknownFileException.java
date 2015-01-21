package errors;

import java.io.File;

import log.Log;

@SuppressWarnings("serial")
public class UnknownFileException extends Exception {

	public UnknownFileException(Log log, File dir) {
		
		log.addLine("There is an error with reading files from " 
				+ dir.getAbsolutePath());
		log.addLine("\t*check that you have the correct flags in your file names");
		log.addLine("\t*and go to api for parameter descriptions");
	}
}
