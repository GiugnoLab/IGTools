package igtools.common.distributions;


import java.util.Map;
import java.util.TreeMap;

import flanagan.analysis.Regression;


/**
 * Theoretic Distribution estimation from real distributions
 * 
 * @author vbonnici
 *
 */
public interface EstimatedDistribution extends IDistribution{

	/**
	 * Estimate, parameters and/or values of the distribution by real data points.
	 * 
	 * @param multDistr, a real distribution in the for  multDistr[i][0] = key, multDistr[i][1] = value
	 */
	public void estimateDistrParameter(double[][] multDistr);
	
	
	
	
	
	
	/**
	 * Geometric distribution
	 * @author vbonnici
	 *
	 */
	public class GeometricBased implements EstimatedDistribution{
		public double ep = 0.0;
		
		@Override
		public void estimateDistrParameter(double[][] multDistr){
			ep = 0.0;
			double n = 0.0;
			
	        for(int i=0; i<multDistr.length; i++){
	            ep += multDistr[i][0] * multDistr[i][1];
	            n += multDistr[i][1];
	            
	        }
                
                
//           for(int i=0; i<multDistr.length; i++){
//	            ep += multDistr[i][1];
//	            n += 1;
//	            
//	        }
                
	        ep = n/ep;
	        
		}
		
		@Override
		public Double getValue(Double x){
			return Math.pow(1.0 - ep, x) * ep;
		}
	}
	
	
	/*public class ExponentialBased implements EstimatedDistribution{
		public double estimateDistrParameter(double[][] multDistr){
                        double sum = 0.0;
                        for(int i=0; i<multDistr.length; i++){
                            sum += multDistr[i][1];
                        }
                    
			double avg = 0.0;
			double n = 0.0;
			for(int i=0; i<multDistr.length; i++){
				n += multDistr[i][1];
				avg += multDistr[i][0] * multDistr[i][1];
                                
                                //avg += multDistr[i][0] * (multDistr[i][1] / sum); 
                                //avg += multDistr[i][0] ;
			}
			return 1.0 / (avg / n);
                        
                        //return 1.0 / (avg / ((double)multDistr.length));
		}
		public double getEstValue(double x, double ep){
			return ep * Math.exp(-(ep * x)); 
		}
	}*/
        
        
        
        /*
        public class ExponentialBased implements EstimatedDistribution{
		public double estimateDistrParameter(double[][] multDistr){
                        double sum_d = 0.0;
                        for(int i=0; i<multDistr.length; i++){
                            sum_d += multDistr[i][1];
                        }
                        
                        double x_avg =  multDistr[multDistr.length -1 ][0] - multDistr[0][0];
                        
                    
			double avg = 0.0;
			for(int i=0; i<multDistr.length; i++){
                                avg += multDistr[i][0];
			}
                        
                        
                        double h_sum = 0.0;
                        for(int i=0; i<multDistr.length; i++){
                            h_sum += (multDistr[i][0]-x_avg) * (multDistr[i][1] - avg);
                        }
                        
                        double l_sum = 0.0;
                        for(int i=0; i<multDistr.length; i++){
                            l_sum += (multDistr[i][0]-x_avg) * (multDistr[i][1] - avg) * (multDistr[i][0]-x_avg) * (multDistr[i][1] - avg);
                        }
                        l_sum = Math.sqrt(l_sum);
                         return h_sum / l_sum;
			
		}
		public double getEstValue(double x, double ep){
			return ep * Math.exp(-(ep * x)); 
		}
	}
        */
        
        /*public class ExponentialBased implements EstimatedDistribution{
		public double estimateDistrParameter(double[][] multDistr){
                        double sum = 0.0;
                        for(int i=0; i<multDistr.length; i++){
                            sum += multDistr[i][1];
                        }
                    
			double avg = 0.0;
			for(int i=0; i<multDistr.length; i++){
				avg += multDistr[i][0] * multDistr[i][1];
			}
                        avg = avg / sum;
                        
                        double sd = 0.0;
                        for(int i=0; i<multDistr.length; i++){
                            sd += multDistr[i][1] * ((multDistr[i][0] - avg ) * (multDistr[i][0] - avg ) );
			}
                        sd = Math.sqrt(sd / sum);
                        
                        
			return sd / avg;
                        
                        //return 1.0 / (avg / ((double)multDistr.length));
		}
		public double getEstValue(double x, double ep){
			return ep * Math.exp(-(ep * x)); 
		}
	}*/
        
	/**
	 * Exponential distribution estimated by maximum likelihood method
	 * 
	 * @author vbonnici
	 *
	 */
        public class ExponentialBased implements EstimatedDistribution{
        	public double ep = 0.0;
        	public ExponentialBased(){}
        	@Override
				public void estimateDistrParameter(double[][] multDistr){
		                        //double sum = (double)(multDistr[multDistr.length-1][0]);
		                    double n = (double)(multDistr.length);
		                    
		                    
		                    
		                    
		                    double m_avg = 0.0;
		                        for(int i=1; i<multDistr.length; i++){
		                            m_avg += Math.abs(multDistr[i][1] - multDistr[i-1][1]);
		                        }
		                        m_avg = m_avg / n;
		                        System.out.println(m_avg);
		                    
		                    /*double sum = 0.0;
		                    for(int i=0; i<multDistr.length; i++){
		                        if(multDistr[i][1] > sum)
		                                sum = multDistr[i][1];
		                    }*/
		                    
		                    
		                        double avg = 0.0;
		                        for(int i=0; i<multDistr.length; i++){
		                            avg += multDistr[i][1];
		                        }
		                        avg = avg / n;
		
		                        
		                        this.ep =  1 / avg;
		                        //return avg;
				}
				
				@Override
				public Double getValue(Double x){
					return ep * Math.exp(-(ep * x)); 
                } 
	}
        
    /**
     * Exponential distribution estimated by a good method, see Flanagan (UCL) library.
     * 
     * @author vbonnici
     *
     */
    public class FExponentialBased implements EstimatedDistribution{

    	public enum TYPE{
    		TWO_PARAMETERS,
    		ONE_PARAMETER,
    		STANDARD,
    		SIMPLE,
    		ONE_MINUS,
    		MULTIPLE
    	}
    	
    	public TYPE regType = TYPE.TWO_PARAMETERS;
    	public int nofExp = 1;
    	
    	public Map<Double,Double> yCalc = new TreeMap<Double,Double>();
    	
    	public FExponentialBased(){
    	}
    	
    	
		@Override
		public void estimateDistrParameter(double[][] multDistr) {
			double[] x = new double[multDistr.length];
			double[] y = new double[multDistr.length];
			
			double sum = 0.0;
			for(int i=0; i<multDistr.length; i++)
				sum += multDistr[i][1];
			
			for(int i=0; i<multDistr.length; i++){
				x[i] = multDistr[i][0];
				y[i] = multDistr[i][1] / sum;
			}
			
			Regression reg = new Regression(x,y);
			
			switch(regType){
			case MULTIPLE:
				reg.exponentialMultiple(nofExp);
				break;
			case ONE_MINUS:
				reg.oneMinusExponential();
				break;
			case ONE_PARAMETER:
				reg.exponentialOnePar();
				break;
			case SIMPLE:
				reg.exponentialSimple();
				break;
			case STANDARD:
				reg.exponentialStandard();
				break;
			case TWO_PARAMETERS:
				reg.exponential();
				break;
			default:
				reg.exponential();
				break;
			
			}
			
			double[] ycalc = reg.getYcalc();
			
			yCalc.clear();
			
			if(x.length>1){
				yCalc.put(x[0], ycalc[0]);
			}
			for(int i=1; i<x.length; i++){
				yCalc.put(x[i], ycalc[i]);
			}
		}

		@Override
		public Double getValue(Double x) {
			if(yCalc.get(x) == null)return 0.0;
			return yCalc.get(x);
		}
    	
    }
}
