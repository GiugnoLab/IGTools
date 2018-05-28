package igtools.analyses;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;


/**
 * Multiplicity/CoMultiplicity and Frequency/CoFrequency Distributions
 * 
 * 
 * @author vbonnici
 *
 */
public class Frequency {

	
	private INELSA nelsa;
	private Map<Integer, Double> nof = new HashMap<Integer,Double>();
	private Map<Integer, Double> totalMult = new HashMap<Integer,Double>();
	
	public Frequency(INELSA nelsa){
		this.nelsa = nelsa;
	}
	
	private double nof(int k){
		Double n = nof.get(k);
		if(n == null){
			n = (double)nelsa.nof(k);
			nof.put(k,n);
		}
		return n;
	}
	
	private double totalMult(int k){
		Double n = totalMult.get(k);
		if(n == null){
			n = 0.0;
			IELSAIterator it = nelsa.begin(k);
			while(it.next()){
				n += it.multiplicity(); 
			}
			totalMult.put(k, n);
		}
		return n;
	}
	
	public double freq(IELSAIterator it){
		double nof = nof(it.k());
		return ((double)it.multiplicity()) / nof;
	}
	
	
	public TreeMap<Integer,Integer> comult_map(int k){
		TreeMap<Integer,Integer> comult = new TreeMap<Integer,Integer>();
		comult_map(k,comult);
		return comult;
	}
	public void comult_map(int k, TreeMap<Integer,Integer> comult){
		Integer cm;
		int mult;
		IELSAIterator it = nelsa.begin(k);
		while(it.next()){
			mult = it.multiplicity();
			cm = comult.get(mult);
			if(cm == null){
				comult.put(mult, 1);
			}
			else{
				comult.put(mult, cm+1);
			}
		}
	}
	
	/**
	 * 
	 * ret[i][0] = multiplicity
	 * ret[i][1] = comultiplicity
	 * 
	 * @param k
	 * @return
	 */
	public int[][] comult_array(int k){
		TreeMap<Integer,Integer> comultsMap = comult_map(k);
		int[][] comult = new int[comultsMap.size()][2];
		int i=0;
		for(Map.Entry<Integer,Integer> entry : comultsMap.entrySet()){
			comult[i][0] = entry.getKey();
			comult[i][1] = entry.getValue();
			i++;
		}
		return comult;
	}
	
	
	
	
	public TreeMap<Double,Integer> cofreq_map(int k){
		TreeMap<Double,Integer> comult = new TreeMap<Double,Integer>();
		cofreq_map(k,comult);
		return comult;
	}
	public void cofreq_map(int k, TreeMap<Double,Integer> comult){
		Integer cm;
		double nof_k = totalMult(k);//nof(k);
		double freq;
		IELSAIterator it = nelsa.begin(k);
		while(it.next()){
			freq = ((double)it.multiplicity()) / nof_k;
			cm = comult.get(freq);
			if(cm == null){
				comult.put(freq, 1);
			}
			else{
				comult.put(freq, cm+1);
			}
		}
	}
	
	/**
	 * ret[i][0] = frequency
	 * ret[i][1] = cofrequency
	 * 
	 * @param k
	 * @return
	 */
	public double[][] cofreq_array(int k){
		TreeMap<Double,Integer> cofreqMap = cofreq_map(k);
		double[][] cofreq = new double[cofreqMap.size()][2];
		int i=0;
		for(Map.Entry<Double,Integer> entry : cofreqMap.entrySet()){
			cofreq[i][0] = entry.getKey();
			cofreq[i][1] = entry.getValue();
			i++;
		}
		return cofreq;
	}
	
}
