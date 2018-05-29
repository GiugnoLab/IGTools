package igtools.analyses.recurrences.carpena;
//package igtools.analyses.toRemove;


import java.util.Map;
import java.util.TreeMap;

import igtools.analyses.random.Random;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3Sequence;
import igtools.common.util.Pair;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.INELSA;
import igtools.dictionaries.elsa.NELSA;

public class ExtractByDistrComp {

	
	private B3Sequence iSeq;
	private INELSA iNelsa;
	
	private B3Sequence randSeq;
	private INELSA randNelsa;
	
	
	public  ExtractByDistrComp(B3Sequence iSeq, INELSA iNelsa){
		this.iSeq = iSeq;
		this.iNelsa = iNelsa;
	}
	
	public void generateRandomSequence(){
		if(randSeq == null){
			System.out.println("Generating random sequence...");
			randSeq = Random.generateUniform(iSeq.length());
			System.out.println("Reputting Ns on it...");
			for(int i=0; i<iSeq.length(); i++){
				if(iSeq.getB3(i) == B3Nucleotide.N_CODE)
					randSeq.setB3(i, B3Nucleotide.N_CODE);
			}
		}
		if(randNelsa == null){
			System.out.println("Computing its NELSA...");
			randNelsa = new NELSA(randSeq);
			System.out.println("done.");
		}
	}
	
	
	
	public interface WordListener{
		public void selected(B3Nucleotide[] word);
	}
	
	
	
	public void releaseResources(){
		this.byLengthDistr.clear();
	}
	
	
	private TreeMap<Integer,Distribution> byLengthDistr = new TreeMap<Integer, Distribution>();
	
	
	public Distribution getCachedByLengthDistr(INELSA nelsa, Integer k, boolean normalizeByMax){
		Distribution distr = byLengthDistr.get(k);
		if(distr == null){
			distr = getByLengthDistr(nelsa, k, normalizeByMax);
			byLengthDistr.put(k, distr);
		}
		return distr;
	}
	public Distribution getByLengthDistr(INELSA nelsa, Integer k, boolean normalizeByMax){
		//Distribution distr = byLengthDistr.get(k);
		//if(distr == null){
			Distribution distr = new Distribution();
			IELSAIterator it = nelsa.begin(k);
			while(it.next()){
				double[] dists = DistalRecurrence.double_proper_codistances(it);
				distr.lazyAddDistribution(dists);
			}
			distr.flushLazyOperations();
			if(normalizeByMax)
				distr.normalizeByMax();
			//byLengthDistr.put(k, distr);
		//}
		return distr;
	}
	
	
	
	
	
	public Distribution getDistribution(IELSAIterator it, boolean normalizeByMax){
		double[] dists = DistalRecurrence.double_proper_codistances(it);
		Distribution itDistr = new Distribution();
		itDistr.lazyAddDistribution(dists);
		if(normalizeByMax)
			itDistr.normalizeByMax();
		return itDistr;
	}
	
	public void run(WordListener listener, TreeMap<Integer,Integer> k_minMult, boolean normalizeByMax){
		generateRandomSequence();
		
		byLengthDistr.clear();
		
		
		for(Map.Entry<Integer,Integer> entry : k_minMult.entrySet()){
			int k = entry.getKey();
			int minMult = entry.getValue();
			
			
			int m;
			IELSAIterator it = iNelsa.begin(k);
		
			while(it.next()){
				m = it.multiplicity();
				if(m >= minMult){
					//get it distance distribution
					Distribution itDistr = getDistribution(it, normalizeByMax);
					Distribution lengthDistr = getCachedByLengthDistr(randNelsa, k,  normalizeByMax);
					//compare
					Double compValue = simpleCompare(itDistr, lengthDistr);
					
				}
			}
		}
		
	}
	
	
	public double simpleCompare(Distribution a, Distribution b){
		double ret = 0.0;
		
		Double bValue;
		for(Map.Entry<Double,Double> entry : a.distribution.entrySet()){
			bValue = b.distribution.get(entry.getKey());
			if(bValue == null){
				//check it
				ret += entry.getValue() * Math.log(entry.getValue());
			}
			else{
				ret += entry.getValue() * Math.log(entry.getValue() / bValue);
			}
		}
		
		return ret;
	}
	
	
	
	public static double KullbackLeibler(Distribution a, Distribution b){
		double ret = 0.0;
		Double bValue;
		for(Map.Entry<Double,Double> entry : a.distribution.entrySet()){
                    
                   if(entry.getValue() != 0.0){
						bValue = b.distribution.get(entry.getKey());
			                        
			                        //System.out.println(entry.getKey() +" "+ entry.getValue() +" "+ bValue);
			                        
						if(bValue == null  ||  bValue == 0.0){
			                            //System.out.println("_"+ (entry.getValue() * Math.log(entry.getValue())));
							//check it, theoretically
							ret += entry.getValue() * Math.log(entry.getValue());
			                            //ret += 0.0;
						}
						else{
			                            if(Double.isNaN(Math.log(entry.getValue()))  ||   Double.isInfinite(Math.log(entry.getValue() / bValue))){
			                                //System.out.println("-"+ (entry.getValue() * Math.log(entry.getValue())));
			                                ret += entry.getValue() * Math.log(entry.getValue());
			                            }
			                            else{
			                                //System.out.println("+"+ (entry.getValue() * Math.log(entry.getValue() / bValue)) );
			                            	ret += entry.getValue() * Math.log(entry.getValue() / bValue);
			                            }
						}
                        //System.out.println("-");
                   }
        }
		//return ret * -1.0;
		return ret;
	}
	
	
	
	
	public static double interplolatedKullbackLeibler(Distribution a, Distribution b){
		double ret = 0.0;
		
		
		
		Double bValue;
		for(Map.Entry<Double,Double> entry : a.distribution.entrySet()){
           if(entry.getValue() != 0.0){
				bValue = b.distribution.get(entry.getKey());   
				if(bValue == null  ||  bValue == 0.0 || Double.isNaN(Math.log(entry.getValue()))  ||   Double.isInfinite(Math.log(entry.getValue() / bValue))){
					ret += entry.getValue() * Math.log(entry.getValue());
				}
				else{
                    ret += entry.getValue() * Math.log(entry.getValue() / bValue);
				}
           }
        }
		return ret;
	}
	
	
	private static double KL(Double a_v, Double b_v){
		if(b_v == null  ||  b_v == 0.0 || Double.isNaN(Math.log( a_v))  ||   Double.isInfinite(Math.log( a_v / b_v))){
			return  a_v * Math.log( a_v);
		}
		else{
            return  a_v * Math.log( a_v / b_v);
		}
	}
	public static double KullbackLeiblerToEstimatedGeometric(Distribution a_distr, int[] codists, ExtractByDistrComp obj){
		double ret = 0.0;
		Distribution g_distr =  obj.estimateGeometric(codists);
		
		
		Double prev_a_v = 0.0, prev_g_v = 0.0;
    	Double prev_d = 0.0;
    	Double a_v, g_v ,d;
    	
		Double nof = 0.0;
		
		
    	boolean doit = false;
    	for(Map.Entry<Double, Double> entry : a_distr.distribution.entrySet()){
    		d = entry.getKey();
    		a_v = entry.getValue();
    		g_v = g_distr.distribution.get(d);
    		
    		
    		if(doit){
    			if(d > prev_d+1){
    				for(double i=prev_d+1; i<d; i++){
    					nof++;
    					KL(
    							(((a_v-prev_a_v)*(i-prev_d))/(d-prev_d)) + prev_a_v,
    							(((g_v-prev_g_v)*(i-prev_d))/(d-prev_d)) + prev_g_v
    					);
    				}
    			}
    			else{
    				ret += KL(a_v, g_v);
    			}
    		}
    		else{
    			doit = true;
    			
    			ret += KL(a_v, g_v);
    		}
    		nof++;
    		
    		prev_a_v = a_v;
    		prev_g_v = g_v;
    		prev_d = d;
       }
    	
//		Double bValue;
//		for(Map.Entry<Double,Double> entry : a.distribution.entrySet()){
//           if(entry.getValue() != 0.0){
//				bValue = b.distribution.get(entry.getKey());   
//				if(bValue == null  ||  bValue == 0.0 || Double.isNaN(Math.log(entry.getValue()))  ||   Double.isInfinite(Math.log(entry.getValue() / bValue))){
//					ret += entry.getValue() * Math.log(entry.getValue());
//				}
//				else{
//                    ret += entry.getValue() * Math.log(entry.getValue() / bValue);
//				}
//           }
//        }
		
		
		//return ret / ((double)a_distr.distribution.size());
    	return ret / nof;
	}
	
	
	
	public Distribution estimateGeometric(int[] distances){

        int[][] codists = DistalRecurrence.co_recurrences_array(distances, false);
        
        double ep = 0.0;
		double n = 0.0;
		
        for(int i=0; i<codists[0].length; i++){
            ep += (double)codists[0][i] * (double)codists[1][i];
            n += (double)codists[1][i];
        }
        ep = n/ep;
        
        TreeMap<Double, Double> distr = new TreeMap<Double, Double>();
        for(int i=0; i<codists[0].length; i++){
        	distr.put((double)codists[0][i], Math.pow(1.0-ep, (double)codists[0][i]) * ep);
        }
        
        Distribution d = new Distribution();
        d.lazyAddDistribution(distr);
        
        return d;
	}
}
