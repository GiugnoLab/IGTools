package igtools.dictionaries.intersection;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.OnDiskNELSAIteratorV2;

public class NumericUnionIntersection {
	public interface UnionListerner{
		public void intersection(B3Nucleotide[] kmer, int multiplicity, int nofSeqs);
		public void union(B3Nucleotide[] kmer, int multiplicity, int nofSeqs);
	}
	
	B3LLSequence[] seqs;
	IELSAIterator[] its;
	long[] orders;
	int k;
	
	UnionListerner listener = null;
	
	
	public long intersectionSize = 0;
	public long unionSize = 0;
	
	
	public NumericUnionIntersection(B3LLSequence[] seqs, IELSAIterator[] its, int k, UnionListerner listener){
		this.seqs = seqs;
		this.its = its;
		this.orders = new long[its.length];
		this.k = k;
		this.listener = listener;
	}
	
	public void union(){
		
		int minkmer_i = -1;
		//int n_minkmer_i = -1;
//		int o_minkmer_i = -1;
		//int compare;
		int ii;
		boolean isint = true;
		
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
		int nofSeqs;
		
		//B3Nucleotide[] ok= new B3Nucleotide[k];
		//long oo = 0L;
		
		while(cicle){
			
			minorder = orders[minkmer_i];
			
			isint = true;
			mult = 0;
			nofSeqs = 0;
			for(ii=0; ii<its.length; ii++){
				if(its[ii] != null && (orders[ii] == minorder)){
					mult += its[ii].multiplicity();
					nofSeqs++;
				}
				else{
					isint = false;
				}
			}
				
			listener.union(kmers[minkmer_i], mult, nofSeqs);
			if(isint){
				listener.intersection(kmers[minkmer_i], mult, nofSeqs);
			}
			
			for(ii=0; ii<its.length; ii++){
				if(its[ii]!=null && (orders[ii] == minorder)){
					
					if(its[ii].next()){
//						copy(kmers[ii], ok);
//						oo = orders[ii];
						
						its[ii].kmer(kmers[ii]);
						orders[ii] = B3Nucleotide.toLexicoOrder(kmers[ii]);
						
//						if(equal(ok,kmers[ii]) || orders[ii] <= oo){
//							System.out.println("# "+ii+" "+B3Nucleotide.toString(ok)+" "+oo+" "+B3Nucleotide.toString(kmers[ii])+" "+orders[ii]);
//							return;
//						}
					}
					else{
//						System.out.println("# "+ii+" is at the end");
//						System.out.println("# "+B3Nucleotide.toString(kmers[ii])+" "+orders[ii]);
						
						orders[ii] = Long.MAX_VALUE;
						if(its[ii] instanceof OnDiskNELSAIteratorV2)
							try{((OnDiskNELSAIteratorV2)its[ii]).close();}catch(Exception e){};
						its[ii] = null;
						
//						return;
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
