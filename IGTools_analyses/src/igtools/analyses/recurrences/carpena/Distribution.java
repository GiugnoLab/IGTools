//package igtools.analyses.toRemove;
package igtools.analyses.recurrences.carpena;

import java.util.Map;
import java.util.TreeMap;

public class Distribution {

	
	public TreeMap<Double, Double> distribution = new TreeMap<Double,Double>();
	
	private double nofDistrInside = 0.0;
	
	public Distribution(){
	}
	
//	public void addDistribution(double[] values){
//		TreeMap<Double, Double> distr = new TreeMap<Double, Double>();
//		for(int i=0; i<values.length; i++){
//			Double count = distr.get(values[i]);
//			if(count == null){
//				distr.put(values[i], 1.0);
//			}
//			else{
//				distr.put(values[i], count +1.0);
//			}
//		}
//	}
//	public void addDistribution(TreeMap<Double, Double> distr){
//		
//	}
	
	public double getNofDistrInside(){
		return this.nofDistrInside;
	}
	
	
	public void lazyAddDistribution(double[] values){
		TreeMap<Double, Double> distr = new TreeMap<Double, Double>();
		for(int i=0; i<values.length; i++){
			Double count = distr.get(values[i]);
			if(count == null){
				distr.put(values[i], 1.0);
			}
			else{
				distr.put(values[i], count +1.0);
			}
		}
		lazyAddDistribution(distr);
	}
	public void lazyAddDistribution(TreeMap<Double, Double> distr){
		for(Map.Entry<Double, Double> entry : distr.entrySet()){
			Double count = distribution.get(entry.getKey());
			if(count == null){
				distribution.put(entry.getKey(), entry.getValue());
			}
			else{
				distribution.put(entry.getKey(), count + entry.getValue());
			}
		}
		nofDistrInside++;
	}
	public void flushLazyOperations(){
		for(Map.Entry<Double, Double> entry : distribution.entrySet()){
			entry.setValue(entry.getValue() / nofDistrInside);
		}
	}
	
	public void normalizeByMax(){
		Double max = Double.NEGATIVE_INFINITY;
		for(Map.Entry<Double, Double> entry : distribution.entrySet()){
			if(entry.getValue() > max)
				max = entry.getValue();
		}
		for(Map.Entry<Double, Double> entry : distribution.entrySet()){
			entry.setValue(entry.getValue() / max);
		}
	}
	
	public void normalizeByMean(){
		Double mean = 0.0;
		for(Map.Entry<Double, Double> entry : distribution.entrySet()){
			mean += entry.getValue();
		}
		if(mean != 0.0){
			mean /= (double)distribution.size();
			for(Map.Entry<Double, Double> entry : distribution.entrySet()){
				entry.setValue(entry.getValue() / mean);
			}
		}
	}
        
        public void convertInFrequency(){
                Double sum = 0.0;
                for(Map.Entry<Double, Double> entry : distribution.entrySet()){
                     sum += entry.getValue();
                }
                if(sum != 0.0){
                    for(Map.Entry<Double, Double> entry : distribution.entrySet()){
                            entry.setValue(entry.getValue() / sum);
                    }
                }
            
        }
        
        public void interpolate(){
        	TreeMap<Double, Double> toadd = new TreeMap<Double,Double>();
        	
        	boolean doit = false;
        	Double prev_v = 0.0;
        	Double prev_k = 0.0;
        	Double v,k;
        	for(Map.Entry<Double, Double> entry : distribution.entrySet()){
        		k = entry.getKey();
        		v = entry.getValue();
        		
        		if(doit){
        			if(k > prev_k+1){
        				for(double i=prev_k+1; i<k; i++){
        					toadd.put(i, (((v-prev_v)*(i-prev_k))/(k-prev_k)) + prev_v );
        				}
        			}
        		}
        		else{
        			doit = true;
        		}
        		prev_v = v;
        		prev_k = k;
           }
        	
        	distribution.putAll(toadd);
        }
	
	
	public static void normalizeByMax(double[] values){
		double max = values[0];
		for(int i=1; i<values.length; i++){
			if(values[i] > max)	max = values[i];
		}
		for(int i=0; i<values.length; i++){
			values[i] /= max;
		}
	}
	
	
}
