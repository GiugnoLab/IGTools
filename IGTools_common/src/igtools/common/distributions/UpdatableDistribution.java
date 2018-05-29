package igtools.common.distributions;

import java.util.Map;
import java.util.TreeMap;

public class UpdatableDistribution {

	
	private TreeMap<Double,Double> distr = new TreeMap<Double,Double>();
	private double nofDistrInside = 0.0;
	
	
	public UpdatableDistribution(){
	}
	
	public TreeMap<Double,Double> getDistribution(){
		return distr;
	}
	
	public void increaseNofDistribtuonsInside(){
		nofDistrInside++;
	}
	
	
	private void update(Double key, Double value, double nof){
		Double v = distr.get(key);
		if(v == null){
			distr.put(key, value / nof);
		}
		else{
			distr.put(key, ((v/nof-1.0) + value) / nof);
		}
	}
	
	public void update(Map<Double,Double> d){
		nofDistrInside++;
		for(Map.Entry<Double, Double> entry : d.entrySet()){
			update(entry.getKey(), entry.getValue(), nofDistrInside);
		}
	}
	
	
	
	
	
	private void updateLazy(Double key, Double value){
		Double v = distr.get(key);
		if(v == null){
			distr.put(key, value);
		}
		else{
			distr.put(key, v + value);
		}
	}
	public void updateLazy(Map<Double,Double> d){
		nofDistrInside++;
		for(Map.Entry<Double, Double> entry : d.entrySet()){
			updateLazy(entry.getKey(), entry.getValue());
		}
	}
	public void updateLazy(double[][] d){
		nofDistrInside++;
		for(int i=0; i< d.length; i++){
			updateLazy(d[i][0], d[i][1]);
		}
	}
	public void flushLazy(){
		for(Map.Entry<Double, Double> entry : distr.entrySet()){
			entry.setValue(entry.getValue() / nofDistrInside);
		}
	}
	
	
	public void normalize(){
		Double max = Double.NEGATIVE_INFINITY;
		for(Double value : distr.values()){
			if (value > max) max = value;
		}
		for(Map.Entry<Double, Double> entry : distr.entrySet()){
			entry.setValue(entry.getValue() / max);
		}
	}
	
}
