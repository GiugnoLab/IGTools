package igtools.analyses.eprob;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3ArraySequence;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;


/**
 * Expected frequency by prefix-suffix like in a Markov Chain.
 * 
 * The prefix is ever k-1.
 * The suffix depends on the chosen order.
 * The order is the overlap length between the prefix and the suffix.
 * 
 * @author vbonnici
 *
 */
public class ExpectedFrequencyByMM extends ExpectedFrequency{

	public ExpectedFrequencyByMM(NELSA nelsa) {
		super(nelsa);
	}

	@Override
	public double expFreq(IELSAIterator wit, int order) {
		
		int k = wit.k();
		
		if( order<0 ) return Double.NaN;
		
        double prob = 0.0;
        
        if(k == 1){
			return 0.25;
			//return Math.abs( (mult / totalCount(k)) - 0.25d );
		}
        else{
        	int r_size = order + 1;
            int l_size = k - 1;
            //int m_size = r_size - 1;
            
            //double sum_k = totalCount(k);
            double sum_r = totalCount(r_size); 
            double sum_l = totalCount(l_size);
            //double sum_m = totalCount(m_size);
            
            IELSAIterator it;
            //IELSAIterator realmer;
            
//            double prob;
            double iprob;
            
	        if(order == 0){
	        	B3Nucleotide[] wit_kmer = wit.kmer();
	        	//realmer = nelsa.find(new B3ArraySequence(wit.kmer(), 0, k));
	           	 
	           	prob = 1.0;
	           	 
	           	it = nelsa.find(new B3ArraySequence(wit_kmer, 0, k-1));
	            if(it!=null)    iprob = ((double)it.multiplicity()) / sum_l;
	            else            iprob = 0;
	            prob *= iprob;
	            
	            it = nelsa.find(new B3ArraySequence(wit_kmer, k-1, k));
	            if(it!=null)    iprob = ((double)it.multiplicity()) / sum_r;
	            else            iprob = 0;
	            prob *= iprob;
	        }
	        else{
	        	
	        	
	        	B3Nucleotide[] wit_kmer = wit.kmer();
	           	//realmer = nelsa.find(new B3ArraySequence(wit.kmer(), 0, k));
	           	 
	           	prob = 1.0;
	            
	            it = nelsa.find(new B3ArraySequence(wit_kmer, k-order-1, k-1));
	            if(it!=null)    iprob = ((double)it.multiplicity()) / sum_r;
	            else            iprob = 0;
	            if(iprob == 0)
	            	prob = 0;
	            else
	            	prob /= iprob;
	            
	            if(prob != 0){
	            	it = nelsa.find(new B3ArraySequence(wit_kmer, 0, k-1));
	                if(it!=null)    iprob = ((double)it.multiplicity()) / sum_l;
	                else            iprob = 0;
	                prob *= iprob;
	                
	                it = nelsa.find(new B3ArraySequence(wit_kmer, k-order-1, k));
	                if(it!=null)    iprob = ((double)it.multiplicity()) / sum_r;
	                else            iprob = 0;
	                prob *= iprob;
	            }
	        }
        }
		
		return prob;
	}

}
