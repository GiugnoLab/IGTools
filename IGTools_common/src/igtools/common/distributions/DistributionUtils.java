package igtools.common.distributions;

import java.util.Map;
import java.util.TreeMap;

public class DistributionUtils {

	
	/**
	 * convert the distribution in a frequency distribution
	 * 
	 * @param distr, must be in the form distr[i][0]=key, distr[i][1]=value
	 */
	public static void normalize(double[][] distr){
		double sum = 0.0;
		for(int i=0;i<distr.length; i++){
			sum += distr[i][1];
		}
		for(int i=0;i<distr.length; i++){
			distr[i][1] /= sum;
		}
	}
	
	
	/**
	 * convert the distribution in a frequency distribution
	 * 
	 * @param distr
	 */
	public static void normalize(Map<Double,Double> distr){
		Double sum = 0.0;
		for(Map.Entry<Double, Double> en : distr.entrySet()){
			sum += en.getValue();
		}
		for(Map.Entry<Double, Double> en : distr.entrySet()){
			en.setValue(en.getValue() / sum);
		}
	}
	
	
	
	
	public static double[][] toArray(Map<Double,Double> distr){
		double[][] ret = new double[distr.size()][2];
		int i=0;
		for(Map.Entry<Double, Double> en : distr.entrySet()){
			ret[i][0] = en.getKey();
			ret[i][1] = en.getValue();
			i++;
		}
		return ret;
	}
	
	public static TreeMap<Double,Double> toTreeMap(double[][] distr){
		TreeMap<Double,Double> dmap = new TreeMap<Double,Double>();
		for(int i=0; i<distr.length; i++){
			dmap.put(distr[i][0], distr[i][1]);
		}
		return dmap;
	}
	
	
	
	
	/**
	 * remove noise by averaging windows
	 * 
	 * @param distr
	 * @param p
	 */
	public static void resample(double[][] distr, double p){
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
	public static void resample(Map<Double,Double> distr, double p){
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
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	 * remove noise by appling a linear kernel 
	 * 
	 * @param distr
	 * @param p
	 */
	public static void linearResample(double[][] distr, double p){
		int wlength = (int)Math.ceil(((double)distr.length)*p/2.0);
		//System.out.println("wlength "+wlength);
		if(wlength > 1){
			double v = 0.0;
			double count = 0.0;
			int j;
			
			
			double[] kernel = new double[wlength];
			for(int i=0; i<wlength; i++){
				kernel[i] = 1.0 / ((double)i+1);
			}
			for(int i=0; i<kernel.length; i++)
				System.out.print(kernel[i]+" ");
			System.out.println();
				
			
			double[] t = new double[distr.length];
			for(int i=0; i<distr.length; i++){
				t[i] = distr[i][1];
			}
			
			
			for(int i=0; i<distr.length; i++){
				v = 0.0;
				count = 0.0;
				
				v += t[i];
				count++;
				
				for(j=0; j<kernel.length; j++){
					if(i-j-1>=0){
						v += t[i-j-1] * kernel[j];
						count++;
					}
				}
				
				for(j=0; j<kernel.length; j++){
					if(i+j+1<distr.length){
						v += t[i+j+1] * kernel[j];
						count++;
					}
				}
				
				distr[i][1] = v/count;
//				System.out.println(i+" "+distr[i][1]+" "+v+" "+count);
			}
			
			
//			for(int i=0; i<distr.length; i++){
//				v = 0.0;
//				count = 0.0;
//				j = i-wlength;
//				if(j<0) j=0;
//				for(;j<=i+wlength && j<distr.length; j++){
//					v += t[j];
//					count++;
//				}
//				distr[i][1] = v/count;
////				System.out.println(i+" "+distr[i][1]);
//			}
		}
	}
	
	
	
	
	/**
	 * remove noise by appling a linear kernel
	 * 
	 * @param distr
	 * @param p
	 */
	public static void linearResample(Map<Double,Double> distr, double p){
		int wlength = (int)Math.ceil(((double)distr.size())*p/2.0);
		//System.out.println("wlength "+wlength);
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
			
			
			double[] kernel = new double[wlength];
			for(i=0; i<wlength; i++){
				kernel[i] = 1.0 / ((double)i+1);
			}
			
			
			i = 0;
			for(Map.Entry<Double,Double> en : distr.entrySet()){
				v = 0.0;
				count = 0.0;
				
				v += t[i];
				count++;
				
				for(j=0; j<kernel.length; j++){
					if(i-j-1>=0){
						v += t[i-j-1] * kernel[j];
						count++;
					}
				}
				
				for(j=0; j<kernel.length; j++){
					if(i+j+1<t.length){
						v += t[i+j+1] * kernel[j];
						count++;
					}
				}
				
				en.setValue(v/count);
				i++;
			}
			
			
//			i = 0;
//			for(Map.Entry<Double,Double> en : distr.entrySet()){
//				v = 0.0;
//				count = 0.0;
//				j = i-wlength;
//				if(j<0) j=0;
//				for(;j<=i+wlength && j<t.length; j++){
//					v += t[j];
//					count++;
//				}
//				en.setValue(v/count);
//				i++;
////				System.out.println(i+" "+distr[i][1]);
//			}
		}
	}
	
	
	
	
	
	
	public static void lowPass(double[][] distr, double alpha){
		double[] out = new double[distr.length];
		//out[0] = distr[0][1];
		out[0] = 0;
		for(int i=1; i<distr.length; i++){
			out[i] = out[i-1] + alpha *  (distr[i][1] - out[i-1]);
		}
		for(int i=distr.length-2; i>=0; i--){
			out[i] = out[i+1] + alpha *  (distr[i][1] - out[i+1]);
		}
		
		for(int i=0; i<distr.length; i++){
			distr[i][1] = out[i];
		}
	}
	
	
	public static void lowPass(double[][] distr, double edistr, double alpha){
		double[] out = new double[distr.length];
		//out[0] = distr[0][1];
		out[0] = edistr;
		for(int i=1; i<distr.length; i++){
			out[i] = out[i-1] + alpha *  (distr[i][1] - out[i-1]);
		}
		for(int i=distr.length-2; i>=0; i--){
			out[i] = out[i+1] + alpha *  (distr[i][1] - out[i+1]);
		}
		
		for(int i=0; i<distr.length; i++){
			distr[i][1] = out[i];
		}
	}
}
