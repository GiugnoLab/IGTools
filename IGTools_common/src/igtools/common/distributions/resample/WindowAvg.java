package igtools.common.distributions.resample;

import java.util.Map;

public class WindowAvg extends DistributionResampler{

	public double p = 0.01; 
	
	public void resample(double[][] distr){
		int wlength = (int)Math.ceil(((double)distr.length)*p/2.0);
		//System.out.println("wlength "+wlength);
		if(wlength > 1){
			double v = 0.0;
			double count = 0.0;
			int j;
			
			double[] t = new double[distr.length];
			for(int i=0; i<distr.length; i++){
				t[i] = distr[i][1];
			}
			
			
			for(int i=0; i<distr.length; i++){
				v = 0.0;
				count = 0.0;
				j = i-wlength;
				if(j<0) j=0;
				for(;j<=i+wlength && j<distr.length; j++){
					v += t[j];
					count++;
				}
				distr[i][1] = v/count;
//				System.out.println(i+" "+distr[i][1]);
			}
		}
	}
	
	
	
	
	/**
	 * remove noise by averaging windows
	 * 
	 * @param distr
	 * @param p
	 */
	public void resample(Map<Double,Double> distr){
		int wlength = (int)Math.ceil(((double)distr.size())*p/2.0);
//		System.out.println("wlength "+wlength);
		if(wlength > 1){
			double v = 0.0;
			double count = 0.0;
			int j;
			
			double[] t = new double[distr.size()];
			
			int i = 0;
			for(Map.Entry<Double,Double> en : distr.entrySet()){
				t[i] = en.getValue();
				i++;
			}
			
			i = 0;
			for(Map.Entry<Double,Double> en : distr.entrySet()){
				v = 0.0;
				count = 0.0;
				j = i-wlength;
				if(j<0) j=0;
				for(;j<=i+wlength && j<t.length; j++){
					v += t[j];
					count++;
				}
				en.setValue(v/count);
				i++;
//				System.out.println(i+" "+distr[i][1]);
			}
		}
	}

}
