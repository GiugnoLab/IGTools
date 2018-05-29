package igtools.common.exceptions;

public class DataVersioningException extends Exception{
	
	public DataVersioningException(long expected, long found){
		super("Wrong data versioning, expected["+expected+"] found["+found+"]");
	}

}
