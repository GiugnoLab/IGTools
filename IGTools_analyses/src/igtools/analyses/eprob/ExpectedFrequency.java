package igtools.analyses.eprob;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3Sequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

import java.util.Map;
import java.util.TreeMap;

public abstract class ExpectedFrequency {

	
	protected NELSA nelsa;
	protected B3Sequence b3seq;
	protected Map<Integer, Double> totalCounts = new TreeMap<Integer, Double>();
	
	
	public ExpectedFrequency(NELSA nelsa){
		this.nelsa = nelsa;
		this.b3seq = nelsa.b3seq();
	}
	
	
	public double totalCount(int k){
		Double ret = totalCounts.get(k);
		if(ret == null){
			int nof = 0;
			
			if(k == 1){
				int[] sa = nelsa.sa();
				for(int i=0; i<nelsa.length(); i++){
					if(b3seq.getB3(sa[i]) != B3Nucleotide.N_CODE)	nof++;
					//else nof++;
				}
			}
			else{
				IELSAIterator it = nelsa.begin(k);
				while(it.next()){
					nof += it.multiplicity();
				}
			}
			
			ret = (double)nof;
			totalCounts.put(k, ret);
		}
		return ret;
	}
	
//	public double realFrequencyOf(IELSAIterator it){
//		return ((double)it.multiplicity()) / totalCount(it.k());
//	}
	public double toFrequency(int m, int k){
		return ((double)m) / totalCount(k);
	}
	public int toMultiplicity(double f, int k){
		//return (int)Math.floor(f * totalCount(k));
		return (int)Math.ceil(f * totalCount(k));
	}
	
	public abstract double expFreq(IELSAIterator wit, int order);
}
