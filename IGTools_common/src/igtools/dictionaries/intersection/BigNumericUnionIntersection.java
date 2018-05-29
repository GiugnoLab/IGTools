package igtools.dictionaries.intersection;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.OnDiskNELSAIteratorV2;

public class BigNumericUnionIntersection {
	public interface UnionListerner{
		public void intersection(B3Nucleotide[] kmer, int multiplicity, int nofSeqs);
		public void union(B3Nucleotide[] kmer, int multiplicity, int nofSeqs);
	}
	
	B3LLSequence[] seqs;
	IELSAIterator[] its;
	int k;
	
	UnionListerner listener = null;
	
	
	public long intersectionSize = 0;
	public long unionSize = 0;
	
	
	public BigNumericUnionIntersection(B3LLSequence[] seqs, IELSAIterator[] its, int k, UnionListerner listener){
		this.seqs = seqs;
		this.its = its;
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
				
				if(minkmer_i == -1){
					minkmer_i = i;
				}
				else{
					if(B3Nucleotide.compare(kmers[i],kmers[minkmer_i])<0){
					//if(orders[i] < orders[minkmer_i]){
						minkmer_i = i;
					}
				}
			}
			else{
				kmers[i] = null;
			
				if(its[i] instanceof OnDiskNELSAIteratorV2)
					try{((OnDiskNELSAIteratorV2)its[i]).close();}catch(Exception e){};
				its[i] = null;
			}
		}
		
		boolean cicle = (minkmer_i != -1);
		int mult;
		int nofSeqs;
		
//		B3Nucleotide[] ok= new B3Nucleotide[k];
//		long oo = 0L;
		
		while(cicle){
			
			//minorder = orders[minkmer_i];
			
			isint = true;
			mult = 0;
			nofSeqs = 0;
			for(ii=0; ii<its.length; ii++){
				if(its[ii] != null && (B3Nucleotide.areEqual(kmers[ii], kmers[minkmer_i]))){
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
				if(ii != minkmer_i){
					if(its[ii]!=null && (B3Nucleotide.areEqual(kmers[ii], kmers[minkmer_i]))){
						if(its[ii].next()){
							its[ii].kmer(kmers[ii]);
						}
						else{
							if(its[ii] instanceof OnDiskNELSAIteratorV2)
								try{((OnDiskNELSAIteratorV2)its[ii]).close();}catch(Exception e){};
							its[ii] = null;
						}
					}
				}
			}
			if(its[minkmer_i].next()){
				its[minkmer_i].kmer(kmers[minkmer_i]);
			}
			else{
				if(its[minkmer_i] instanceof OnDiskNELSAIteratorV2)
					try{((OnDiskNELSAIteratorV2)its[minkmer_i]).close();}catch(Exception e){};
				its[minkmer_i] = null;
			}
			
			
			
			
			minkmer_i = -1;
			for(ii=0; ii<its.length; ii++){
				if(its[ii]!=null &&   ((minkmer_i == -1) ||  (B3Nucleotide.compare(kmers[ii], kmers[minkmer_i])<0) )){
					minkmer_i = ii;
				}
			}
			
			cicle = (minkmer_i != -1);
		}
	}
	
}
