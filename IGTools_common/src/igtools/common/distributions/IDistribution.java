package igtools.common.distributions;


/**
 * General Distribution interface
 * 
 * @author vbonnici
 *
 */
public interface IDistribution {
	
	/*
	/**
	 * @param x
	 * @return the value for x from the estimated distribution 
	 *
	public double geValue(double x);
	*/
	
	/**
	 * @param x
	 * @return the value for x from the estimated distribution 
	 */
	public Double getValue(Double x);
}
