package igtools.analyses.eprob;


import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

/**
 * Expected frequency by sub-factors
 * Order = sub-factor overlapping length
 * 
 * @author vbonnici
 *
 */
public class ExpectedFrequencyBySubFactors extends ExpectedFrequency{

	
	public ExpectedFrequencyBySubFactors(NELSA nelsa) {
		super(nelsa);
	}

	public double expFreq(IELSAIterator wit, int order){
		int k = wit.k();
		
		if( order<0 || order>k-1) return Double.NaN;
		//if( order>k-1 || (k==2 && order<1)) return -1.0;
		
		//double mult = (double)wit.multiplicity();
		if(k == 1){
			return 0.25d;
			//return Math.abs( (mult / totalCount(k)) - 0.25d );
		}
		else{
			if(order == 0){
				double sum_hk = totalCount(1);
				double prob = 1.0;
				double iprob = 0.0;
				
				String skmer = B3Nucleotide.toString(wit.kmer());
				
				IELSAIterator it;
				for(int i=0; i<skmer.length(); i++){
					it = nelsa.find(new B3LLSequence(skmer.substring(i, i+1)));
                    if(it!=null)    iprob = ((double)it.multiplicity()) / sum_hk;
                    else            iprob = 0;
                    prob *= iprob;
				}
				
				return prob;
			}
			else if(order == k-2){
				double prob = 0.0;
				
				//double sum_k = totalCount(k);
				double sum_pk = totalCount(k-1);
				double sum_ppk = totalCount(k-2);
				
                double prob_l = 0.0;
                double prob_c = 0.0;
                double prob_r = 0.0;
                
                String skmer = B3Nucleotide.toString(wit.kmer());

                String tkmer_l = skmer.substring(0, skmer.length() -1);
                String tkmer_c = skmer.substring(1, skmer.length() -1);
                String tkmer_r = skmer.substring(1, skmer.length());
                
                IELSAIterator it = nelsa.find(new B3LLSequence(tkmer_l));
	                if(it!=null)    prob_l = ((double)it.multiplicity()) / sum_pk;
	                else            prob_l = 0;
	            it = nelsa.find(new B3LLSequence(tkmer_c));
	                if(it!=null)    prob_c = ((double)it.multiplicity()) / sum_ppk;
	                else            prob_c = 0;
	            it = nelsa.find(new B3LLSequence(tkmer_r));
	                if(it!=null)    prob_r = ((double)it.multiplicity()) / sum_pk;
	                else            prob_r = 0;
	                
	            if(prob_c!=0)   prob = (prob_l * prob_r) / prob_c;
	            else prob = (prob_l * prob_r);
	            
	            return prob;
			}
			else{
				double sum_hk = totalCount(order+1);
				double sum_lk = totalCount(order);
				
				double prob = 1.0;
				double iprob = 0.0;
				
				String skmer = B3Nucleotide.toString(wit.kmer());
				
				IELSAIterator it;
				
				for(int i=0; i<skmer.length()-order; i++){
                    it = nelsa.find(new B3LLSequence(skmer.substring(i, i+order+1)));
                    if(it!=null)    iprob = ((double)it.multiplicity()) / sum_hk;
                    else            iprob = 0;
                    prob *= iprob;
                }
                
                
                for(int i=1; i<skmer.length()-order; i++){
                    it = nelsa.find(new B3LLSequence(skmer.substring(i, i+order)));
                    if(it!=null)    iprob = ((double)it.multiplicity()) / sum_lk;
                    else            iprob = 0;
                    prob /= iprob;
                }
                
                return prob;
			}
		}
	}
	
	
}
