package igtools.dictionaries.intersection;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.OnDiskNELSAIteratorV2;

public class BigNumericUnion  extends NumericUnion{
	
//	B3LLSequence[] seqs;
//	IELSAIterator[] its;
//	//long[] orders;
//	int k;
//	
//	UnionListerner listener = null;
//	
//	
//	public long intersectionSize = 0;
//	public long unionSize = 0;
	
	//protected B3Nucleotide[] minkmer;
	
	
	public BigNumericUnion(B3LLSequence[] seqs, IELSAIterator[] its, int k, UnionListerner listener){
		super(seqs,its,k,listener);
		
		//minkmer = new B3Nucleotide[k];
		
		//this.seqs = seqs;
		//this.its = its;
		//this.orders = new long[its.length];
		//this.k = k;
		//this.listener = listener;
	}
	
	public void union(){
		B3Nucleotide[][] orders = new B3Nucleotide[its.length][];
		for(int i=0; i<orders.length; i++){
			orders[i] = new B3Nucleotide[k];
		}
		
		int minkmer_i = -1;
		int ii;
		
		B3Nucleotide[][] kmers = new B3Nucleotide[its.length][];
		
		for(int i=0; i<its.length; i++){
			if(its[i].next()){
				kmers[i] = new B3Nucleotide[k];
				its[i].kmer(kmers[i]);
				//orders[i] = B3Nucleotide.toLexicoOrder(kmers[i]);
				B3Nucleotide.copy(kmers[i], orders[i]);
				
				if(minkmer_i == -1){
					minkmer_i = i;
				}
				else{
					//if(orders[i] < orders[minkmer_i]){
					if(B3Nucleotide.compare(orders[i],orders[minkmer_i]) <0){
						minkmer_i = i;
					}
				}
			}
			else{
				kmers[i] = null;
				//orders[i] = Long.MAX_VALUE;
				orders[i] = null;
			
				if(its[i] instanceof OnDiskNELSAIteratorV2)
					try{((OnDiskNELSAIteratorV2)its[i]).close();}catch(Exception e){};
				its[i] = null;
			}
		}
		
		System.out.println("first min: "+minkmer_i+" "+B3Nucleotide.toString(kmers[minkmer_i]));
		
		
		boolean cicle = (minkmer_i != -1);
		//long minorder;
		B3Nucleotide[] minorder = new B3Nucleotide[k];
		int mult;
		
		while(cicle){
			
			//minorder = orders[minkmer_i];
			B3Nucleotide.copy(orders[minkmer_i], minorder);
			
			mult = 0;
			for(ii=0; ii<its.length; ii++){
				//if(its[ii] != null && (orders[ii] == minorder)){
				if(its[ii] != null && (B3Nucleotide.compare(orders[ii],minorder) == 0)){
					mult += its[ii].multiplicity();
				}
			}
				
			listener.union(kmers[minkmer_i], mult);
			
			for(ii=0; ii<its.length; ii++){
				//if(its[ii]!=null && (orders[ii] == minorder)){
				if(its[ii] != null && (B3Nucleotide.compare(orders[ii],minorder) == 0)){
					
					if(its[ii].next()){
						its[ii].kmer(kmers[ii]);
						//orders[ii] = B3Nucleotide.toLexicoOrder(kmers[ii]);
						B3Nucleotide.copy(kmers[ii],orders[ii]);
					}
					else{
						//orders[ii] = Long.MAX_VALUE;
						orders[ii] = null;
						if(its[ii] instanceof OnDiskNELSAIteratorV2)
							try{((OnDiskNELSAIteratorV2)its[ii]).close();}catch(Exception e){};
						its[ii] = null;
					}
				}
			}
			
			minkmer_i = -1;
			//minorder = Long.MAX_VALUE;
			for(ii=0; ii<its.length; ii++){
				//if(its[ii]!=null && (orders[ii]< minorder)){
				if(its[ii]!=null && ((minkmer_i == -1) || (B3Nucleotide.compare(orders[ii],minorder)<0))){
					minkmer_i = ii;
					//minorder = orders[ii];
					B3Nucleotide.copy(orders[ii], minorder);
				}
			}
			
			cicle = (minkmer_i != -1);
		}
	}
	
//	public void union(){
//		
//		int minkmer_i = -1;
//		//int compare;
//		int ii;
//	//	boolean isint = true;
//		
//		B3Nucleotide[][] kmers = new B3Nucleotide[its.length][];
//		
//		for(int i=0; i<its.length; i++){
//			if(its[i].next()){
//				kmers[i] = new B3Nucleotide[k];
//				its[i].kmer(kmers[i]);
//				//orders[i] = B3Nucleotide.toLexicoOrder(kmers[i]);
//				
//				if(minkmer_i == -1){
//					minkmer_i = i;
//				}
//				else{
//					if(B3Nucleotide.compare(kmers[i],kmers[minkmer_i])<0){
//					//if(orders[i] < orders[minkmer_i]){
//						minkmer_i = i;
//					}
//				}
//			}
//			else{
//				kmers[i] = null;
//				//orders[i] = Long.MAX_VALUE;
//			
//				if(its[i] instanceof OnDiskNELSAIteratorV2)
//					try{((OnDiskNELSAIteratorV2)its[i]).close();}catch(Exception e){};
//				its[i] = null;
//			}
//		}
//		
//		System.out.println("first min: "+minkmer_i+" "+B3Nucleotide.toString(kmers[minkmer_i]));
//		
//		
//		boolean cicle = (minkmer_i != -1);
//		//long minorder;
//		int mult;
//		
//		while(cicle){
//			
//			//minorder = orders[minkmer_i];
//			
//			for(ii=0; ii<k; ii++){
//				minkmer[ii] = kmers[minkmer_i][ii];
//			}
//			
//			mult = 0;
////			for(ii=0; ii<its.length; ii++){
////				if(its[ii] != null && (B3Nucleotide.areEqual(kmers[ii], kmers[minkmer_i]))){
////				//if(its[ii] != null && (orders[ii] == minorder)){
////					mult += its[ii].multiplicity();
////				}
////			}
////				
////			listener.union(kmers[minkmer_i], mult);
//			
//			for(ii=0; ii<its.length; ii++){
//				//if(its[ii]!=null && (orders[ii] == minorder)){
//				//if(its[ii]!=null && (B3Nucleotide.areEqual(kmers[ii], kmers[minkmer_i]))){
//				if(its[ii]!=null && (B3Nucleotide.areEqual(kmers[ii], minkmer))){
//					mult += its[ii].multiplicity();
//					
//					if(its[ii].next()){
//						its[ii].kmer(kmers[ii]);
//						//orders[ii] = B3Nucleotide.toLexicoOrder(kmers[ii]);
//					}
//					else{
//						//its[ii] = null;
//						//orders[ii] = Long.MAX_VALUE;
//						if(its[ii] instanceof OnDiskNELSAIteratorV2)
//							try{((OnDiskNELSAIteratorV2)its[ii]).close();}catch(Exception e){};
//						its[ii] = null;
//						kmers[ii] = null;
//					}
//				}
//			}
//			
//			listener.union(minkmer, mult);
//			
//			minkmer_i = -1;
//			//minorder = Long.MAX_VALUE;
//			for(ii=0; ii<its.length; ii++){
//				//if(its[ii]!=null && (orders[ii]< minorder)){
//				//if(its[ii]!=null &&   ((minkmer_i == -1) ||  (B3Nucleotide.compare(kmers[ii], kmers[minkmer_i])<0) )){
//				if(its[ii]!=null &&   ((minkmer_i == -1) ||  (B3Nucleotide.compare(kmers[ii], minkmer)<0) )){
//					minkmer_i = ii;
//					//minorder = orders[ii];
//				}
//			}
//			
//			cicle = (minkmer_i != -1);
//		}
//	}
	
	
}