package igtools.dictionaries.intersection;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSA;
import igtools.dictionaries.elsa.IELSAIterator;

public class NumericUnionIterator implements IELSAIterator {
	
	private B3LLSequence[] seqs;
	private IELSAIterator[] its;
	private int k;
	
	private long[] orders;
	private B3Nucleotide[][] kmers;
	private boolean first;
	private int minkmer_i;
	private long minorder;
	
	private long multiplicity;
	private long[] multiplicities;
	
	public NumericUnionIterator(B3LLSequence[] seqs, IELSAIterator[] its, int k){
		this.seqs = seqs;
		this.its = its;
		this.k = k;

		this.orders = new long[its.length];
		this.kmers = new B3Nucleotide[its.length][];
		this.multiplicity = 0;
		this.multiplicities = new long[its.length];
		this.minorder = Long.MAX_VALUE;
		this.first = true;
	}
	
	
	public void print(){
		System.out.print("orders [ ");
		for(int i=0; i<orders.length; i++){
			if(its[i] != null) System.out.print(orders[i]+" ");
			else System.out.print("null ");
		}
		System.out.println("]");
		System.out.print("kmers [ ");
		for(int i=0; i<kmers.length; i++){
			if(its[i] != null) System.out.print(B3Nucleotide.toString(kmers[i])+" ");
			else System.out.print("null ");
		}
		System.out.println("]");
		System.out.print("multiplicities [ ");
		for(int i=0; i<multiplicities.length; i++){
			if(its[i] != null) System.out.print(multiplicities[i]+" ");
			else System.out.print("null ");
		}
		System.out.println("]");
		System.out.println("minkmer_i["+minkmer_i+"] minorder["+minorder+"]");
	}
	
	
	public long oder(){
		return this.minorder;
	}
	

	@Override
	public IELSA elsa() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}

	@Override
	public int istart() {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

	@Override
	public int iend() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}

	@Override
	public int k() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public B3Nucleotide[] kmer() {
		return this.kmers[minkmer_i];
	}

	@Override
	public void kmer(B3Nucleotide[] ns) {
		for(int i=0; i<ns.length && i<this.kmers[minkmer_i].length; i++){
			ns[i] = this.kmers[minkmer_i][i];
		}
	}

	@Override
	public int multiplicity() {
		return (int)this.multiplicity;
	}

	public long longMultiplicity() {
		return this.multiplicity;
	}
	public long[] multiplicities(){
		return this.multiplicities;
	}
	
	@Override
	public int[] positions() {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

	@Override
	public int[] sortedPositions() {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

	@Override
	public boolean isMinimalHapax() {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

	@Override
	public boolean isGlobalMaximalRepeat() {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

	@Override
	public boolean next() {
		
		if(this.first){
			this.first = false;
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
					its[i] = null;
				}
			}
			
			if(minkmer_i != -1){
				minorder = orders[minkmer_i];
				multiplicity = 0;
				for(int ii=0; ii<its.length; ii++){
					if(its[ii] != null && (orders[ii] == minorder)){
						multiplicity += its[ii].multiplicity();
						multiplicities[ii] = its[ii].multiplicity();
					}
					else{
						multiplicities[ii] = 0;
					}
				}
				return true;
			}
			else{
				return false;
			}
		}
		else{
			if(minkmer_i != -1){
//				System.out.println("# "+minorder+" "+orders[minkmer_i]);
				//print();
				
				minorder =  orders[minkmer_i];
				
				for(int ii=0; ii<its.length; ii++){
					if(its[ii]!=null && (orders[ii] == minorder)){
//						System.out.println("+ "+ii);
						if(its[ii].next()){
							its[ii].kmer(kmers[ii]);
							orders[ii] = B3Nucleotide.toLexicoOrder(kmers[ii]);
						}
						else{
							its[ii] = null;
							orders[ii] = Long.MAX_VALUE;
							its[ii] = null;
						}
					}
					else{
//						System.out.println("- "+ii);
					}
				}
				//print();
//				System.out.println("...");
				
				
				minkmer_i = -1;
				minorder = Long.MAX_VALUE;
				for(int ii=0; ii<its.length; ii++){
					if(its[ii]!=null && (orders[ii]< minorder)){
						minkmer_i = ii;
						minorder = orders[ii];
					}
				}
				
				if(minkmer_i != -1){
					multiplicity = 0;
					for(int ii=0; ii<its.length; ii++){
						if(its[ii] != null && (orders[ii] == minorder)){
							multiplicity += its[ii].multiplicity();
							multiplicities[ii] = its[ii].multiplicity();
						}
						else{
							multiplicities[ii] = 0;
						}
					}
					return true;
				}
				//return false;
			}
			//return false;
		}
		return false;
	}

	@Override
	public boolean prev() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean good() {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

	@Override
	public boolean hasNext() {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

	@Override
	public boolean hasPrev() {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

	@Override
	public IELSAIterator clone() {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

	@Override
	public int compare(IELSAIterator it) {
		// TODO Auto-generated method stub
				throw new UnsupportedOperationException();
	}

}
