package igtools.dictionaries.intersection;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.OnDiskNELSAIteratorV2;

public class NumericUnion {
	public interface UnionListerner{
		public void intersection(B3Nucleotide[] kmer);
		public void union(B3Nucleotide[] kmer, int multiplicity);
	}
	
	protected B3LLSequence[] seqs;
	protected IELSAIterator[] its;
	protected int k;
	
	protected UnionListerner listener = null;
	
	
	public long intersectionSize = 0;
	public long unionSize = 0;
	
	
	public NumericUnion(B3LLSequence[] seqs, IELSAIterator[] its, int k, UnionListerner listener){
		this.seqs = seqs;
		this.its = its;
		this.k = k;
		this.listener = listener;
	}
	
	public void union(){
		long[] orders = new long[its.length];
		
		int minkmer_i = -1;
		//int n_minkmer_i = -1;
//		int o_minkmer_i = -1;
		//int compare;
		int ii;
		//boolean isint = true;
		
		B3Nucleotide[][] kmers = new B3Nucleotide[its.length][];
		
		for(int i=0; i<its.length; i++){
			if(its[i].next()){
				kmers[i] = new B3Nucleotide[k];
				its[i].kmer(kmers[i]);
				orders[i] = B3Nucleotide.toLexicoOrder(kmers[i]);
				
				if(minkmer_i == -1){
					minkmer_i = i;
				}
				else{
					//if(B3Nucleotide.compare(kmers[i],kmers[minkmer_i])<0){
					if(orders[i] < orders[minkmer_i]){
						minkmer_i = i;
					}
				}
			}
			else{
				kmers[i] = null;
				orders[i] = Long.MAX_VALUE;
			
				if(its[i] instanceof OnDiskNELSAIteratorV2)
					try{((OnDiskNELSAIteratorV2)its[i]).close();}catch(Exception e){};
				its[i] = null;
			}
		}
		
		//System.out.println("first min: "+minkmer_i+" "+B3Nucleotide.toString(kmers[minkmer_i]));
		
		
		boolean cicle = (minkmer_i != -1);
		long minorder;
		int mult;
		
		while(cicle){
			
			minorder = orders[minkmer_i];
			
			mult = 0;
			for(ii=0; ii<its.length; ii++){
				if(its[ii] != null && (orders[ii] == minorder)){
					mult += its[ii].multiplicity();
				}
			}
				
			listener.union(kmers[minkmer_i], mult);
			
			for(ii=0; ii<its.length; ii++){
				if(its[ii]!=null && (orders[ii] == minorder)){
					
					if(its[ii].next()){
						its[ii].kmer(kmers[ii]);
						orders[ii] = B3Nucleotide.toLexicoOrder(kmers[ii]);
					}
					else{
						//its[ii] = null;
						orders[ii] = Long.MAX_VALUE;
						if(its[ii] instanceof OnDiskNELSAIteratorV2)
							try{((OnDiskNELSAIteratorV2)its[ii]).close();}catch(Exception e){};
						its[ii] = null;
					}
				}
			}
			
			minkmer_i = -1;
			minorder = Long.MAX_VALUE;
			for(ii=0; ii<its.length; ii++){
				if(its[ii]!=null && (orders[ii]< minorder)){
					minkmer_i = ii;
					minorder = orders[ii];
				}
			}
			
			cicle = (minkmer_i != -1);
		}
	}
	
	
}
