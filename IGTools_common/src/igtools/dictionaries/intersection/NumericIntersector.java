package igtools.dictionaries.intersection;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;


/**
 * Performs only intersection.
 * Comparisons are done by B3Nucleotide.toLexicoOrder
 * 
 * It does not report multiplicity, like NumericIntersection, but I'm not sure which one is correct.
 *  
 * @author vbonnici
 *
 */
public class NumericIntersector {
	public interface IntersectorListerner{
		public void intersection(B3Nucleotide[] kmer);
	}
	
	
	B3LLSequence[] seqs;
	long[] orders;
	IELSAIterator[] its;
	int k;
	IntersectorListerner listener;
	
	
	public int intersectionSize = 0;
	
	
	public NumericIntersector(B3LLSequence[] seqs, IELSAIterator[] its, int k, IntersectorListerner listener){
		this.seqs = seqs;
		this.its = its;
		this.orders = new long[its.length];
		this.k = k;
		this.listener = listener;
	}
	
	public void intersect(){
		boolean cicle = true;
		for(int i=0; i<its.length; i++){
			if(!its[i].next()){
				cicle = false;
				break;
			}
		}
		
		if(cicle){
			B3Nucleotide[][] kmers = new B3Nucleotide[its.length][];
			for(int i=0; i<its.length; i++){
				kmers[i] = new B3Nucleotide[k];
				its[i].kmer(kmers[i]);
				orders[i] = B3Nucleotide.toLexicoOrder(kmers[i]);
			}
			
			int kmer_i = 0;
			
			long compare;
			
			
			while(cicle){
				
				for(int i=0; i<its.length; i++){
						//compare = B3Nucleotide.compare(kmers[kmer_i], kmers[i]);
						compare = orders[kmer_i] - orders[i]; 
						
						if(compare == 0){
						}
						else{
							if(compare > 0){
								while(compare >0){
									if(!its[i].next()){
										cicle = false;
										break;
									}
									else{
										its[i].kmer(kmers[i]);
										orders[i] = B3Nucleotide.toLexicoOrder(kmers[i]);
										//compare = B3Nucleotide.compare(kmers[kmer_i], kmers[i]);
										compare = orders[kmer_i] - orders[i];
									}
								}
								
								if(cicle){
									if(compare == 0){
									}
									else if(compare<0){
										kmer_i = i;
										i = -1;
									}
									else{
										System.out.println("OPS");
									}
								}
							}
							else{
								kmer_i = i;
								i = -1;
							}
						}
				}	
				if(cicle){
					intersectionSize++;
					if(listener != null)
						listener.intersection(kmers[kmer_i]);
					
					if(!its[kmer_i].next()){
						cicle = false;
						break;
					}
					else{
						its[kmer_i].kmer(kmers[kmer_i]);
						orders[kmer_i] = B3Nucleotide.toLexicoOrder(kmers[kmer_i]);
					}
				}
			}
		}
	}
	
}
