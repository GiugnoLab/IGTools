package igtools.dictionaries.intersection;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;


/**
 * Performs only intersection.
 * Comparisons are done by B3Nucleotide.toLexicoOrder 
 *  
 * @author vbonnici
 *
 */
public class NumericIntersection {
	
	public interface IntersectionListerner{
		/**
		 * The reported multiplicity is the sum.
		 * 
		 * @param kmer
		 * @param multiplicity
		 */
		public void intersection(B3Nucleotide[] kmer, int multiplicity);
	}
	
	B3LLSequence[] seqs;
	IELSAIterator[] its;
	long[] orders;
	int k;
	
	IntersectionListerner listener = null;
	
	
	public long intersectionSize = 0;
	public long unionSize = 0;
	
	
	public NumericIntersection(B3LLSequence[] seqs, IELSAIterator[] its, int k, IntersectionListerner listener){
		this.seqs = seqs;
		this.its = its;
		this.orders = new long[its.length];
		this.k = k;
		this.listener = listener;
	}
	
	public void intersection(){
		
		int maxkmer_i = -1;
		B3Nucleotide[][] kmers = new B3Nucleotide[its.length][];
		
		boolean isint = true;
		boolean isok = true;
		int mult;
		int i;
		
		for(i=0; i<its.length; i++){
			if(its[i].next()){
				kmers[i] = new B3Nucleotide[k];
				its[i].kmer(kmers[i]);
				orders[i] = B3Nucleotide.toLexicoOrder(kmers[i]);
				
				if(maxkmer_i == -1){
					maxkmer_i = i;
				}
				else{
					isint &= (orders[i] == orders[maxkmer_i]);
					if(orders[i] > orders[maxkmer_i]){
						maxkmer_i = i;
					}
				}
			}
			else{
				return;
			}
		}
		
		if(isint){
			mult = 0;
			for(i=0; i<its.length; i++){
				mult += its[i].multiplicity();
			}
			listener.intersection(kmers[0], mult);
		}
		
		
		
		while(true){
			if(its[maxkmer_i].next()){
				its[maxkmer_i].kmer(kmers[maxkmer_i]);
				orders[maxkmer_i] = B3Nucleotide.toLexicoOrder(kmers[maxkmer_i]);
				
				
				for(i=0; i<its.length; i++){
					
					if(i != maxkmer_i){
						isok = true;
						while((isok = its[i].next())){
							its[i].kmer(kmers[i]);
							orders[i] = B3Nucleotide.toLexicoOrder(kmers[i]);
							if(orders[i] == orders[maxkmer_i]){
								break;
							}
							else if(orders[i] > orders[maxkmer_i]){
								maxkmer_i = i;
								i = -1;
								break;
							}
						}
						if(!isok)
							return;
					}
				}
				
			}
			else{
				return;
			}
			
			mult = 0;
			for(i=0; i<its.length; i++){
				if(orders[i] != orders[maxkmer_i]){
					System.out.println(i+" "+maxkmer_i+" "+orders[i]+" "+orders[maxkmer_i]+" "+its[i].multiplicity()+" "+its[maxkmer_i].multiplicity()+" "+B3Nucleotide.toString(kmers[i])+" "+B3Nucleotide.toString(kmers[maxkmer_i]));
					return;
				}
				mult += its[i].multiplicity();
			}
			listener.intersection(kmers[0], mult);
		}

	}
	
	
}
