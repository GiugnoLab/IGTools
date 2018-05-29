package igtools.common.sequence;

import igtools.common.exceptions.StreamException;
import igtools.common.kmer.b2.array.B2LLKmer;
import igtools.common.kmer.b2.array.B2LLKmerFiller;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.streams.INucleotideStream;

//public class B3LLSequence implements ISequence<B2LLKmer>, CharSequence{
public class B2LLSequence implements Comparable<B2LLSequence>{
	
	public final static int BITS_PER_BLOCK = 32;
	public final static int BITS_PER_MER = 2;
	public final static int MERS_PER_BLOCK = BITS_PER_BLOCK / BITS_PER_MER;
	
	public final static int FIRST_SHIFT = BITS_PER_BLOCK - BITS_PER_MER;//32 - 2 = 30 
	public final static int MASK = 0x3;//111 
	
	public final static int ALL_ONES = 0xFFFFFFFF;
	
	
	public int[] blocks;
	public int length;
	
	/**
	 * 
	 * @param length number of nucleotides
	 */
	public B2LLSequence(int length){
		blocks = new int[(length / MERS_PER_BLOCK) + ( (length % MERS_PER_BLOCK) > 0 ? 1 : 0)];//ceil(length / MER_PER_BLOCK)
		this.length = length;
	}
	
	public B2LLSequence(String s){
		clear(s.length());
		for(int i=0; i<s.length(); i++){
			setB3(i, B3Nucleotide.by(s.charAt(i)).code());
		}
	}
	
	public B2LLSequence(B3Nucleotide[] s){
		clear(s.length);
		for(int i=0; i<s.length; i++){
			setB3(i, s[i].code());
		}
	}
	
	public void warp(int length, int[] blocks){
		this.blocks = blocks;
		this.length = length;
		if(blocks.length < (length / MERS_PER_BLOCK) + ( (length % MERS_PER_BLOCK) > 0 ? 1 : 0)){
			throw new Error("expected blocks length is less that the given length");
		}
	}
	
	public void clear(int length){
		blocks = new int[(length / MERS_PER_BLOCK) + ( (length % MERS_PER_BLOCK) > 0 ? 1 : 0)];//ceil(length / MER_PER_BLOCK)
		this.length = length;
	}
	
	public int length(){
		return length;
	}
	
	public B2LLSequence clone(){
		B2LLSequence seq = new B2LLSequence(this.length);
		for(int i=0; i<seq.blocks.length; i++){
			seq.blocks[i] = this.blocks[i];
		}
		return seq;
	}
	
	
	public void fillBy(INucleotideStream<B3Nucleotide> gstream) throws StreamException{
		int shift = FIRST_SHIFT;
		int iblock = 0;
//		B3Nucleotide n;
//		while((n = gstream.next()) != null){
//			blocks[iblock] |= n.code() << shift;
		int n;
		while((n = gstream.nextCode()) != INucleotideStream.NULL_CODE){
			blocks[iblock] |= n << shift;
			
			shift -= BITS_PER_MER;
			if(shift < 0){
				shift = FIRST_SHIFT;
				iblock++;
			}
		}
	}
	
	
	public B2LLKmer get(int pos, int length){
		if(pos + length > this.length) return null;
		
		B2LLKmerFiller filler = new B2LLKmerFiller(length);
		
		int iblock = (pos / MERS_PER_BLOCK);
		int mask = (MASK << FIRST_SHIFT ) >>> ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		int shift = FIRST_SHIFT - ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		int n;
		
		for(int i=0;i<length; i++){
			n = (blocks[iblock] & mask) >>> shift;
			
//			System.out.println(n+"	"+Nucleotide.by(n));
			
			if(n > 3){
				return null;
			}
			
			filler.push(n);
			
			mask = mask >>> BITS_PER_MER;
			shift -= BITS_PER_MER;
			
			if(mask < 7){
				mask = (MASK << FIRST_SHIFT);
				shift = FIRST_SHIFT;
				iblock++;
			}
		}
		
		return new B2LLKmer(length, filler.code());
	}
	
	public int getB2(int pos){
		if(pos >= this.length) return -1;
		
		int iblock = (pos / MERS_PER_BLOCK);
		int mask = (MASK << FIRST_SHIFT ) >>> ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		int shift = FIRST_SHIFT - ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		
		int n = (blocks[iblock] & mask) >>> shift;
		if(n>3)
			return -1;
		return n;
	}
	
	public int getB3(int pos){
		if(pos >= this.length) return -1;
		
		int iblock = (pos / MERS_PER_BLOCK);
		int mask = (MASK << FIRST_SHIFT ) >>> ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		int shift = FIRST_SHIFT - ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		
		int n = (blocks[iblock] & mask) >>> shift;
		return n;
	}
	
	
	
	public int getB3(int pos,B3Nucleotide[] ns){
//		int l = ns.length - ( this.length - (pos + ns.length));
		int l = 0;
		
		if(pos >=this.length) return 0;
		if(pos + ns.length >= this.length) 
			l = this.length - pos;
		else 
			l = ns.length;
		
		
		int iblock = (pos / MERS_PER_BLOCK);
		int mask = (MASK << FIRST_SHIFT ) >>> ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		int shift = FIRST_SHIFT - ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		int n;
		
		for(int i=0;i<l; i++){
			n = (blocks[iblock] & mask) >>> shift;
			
			ns[i] = B3Nucleotide.by(n);
			
			mask = mask >>> BITS_PER_MER;
			shift -= BITS_PER_MER;
			
			if(mask < 7){
				mask = (MASK << FIRST_SHIFT);
				shift = FIRST_SHIFT;
				iblock++;
			}
		}
		
		return l;
	}
	
	
	public void setB3(int pos, int n){
		//TODO
		int iblock = (pos / MERS_PER_BLOCK);
		
//		System.out.println("@"+iblock+"#"+blocks.length);
		
		if(pos>= this.length){
//			int iblock = (pos / MERS_PER_BLOCK);
//			System.out.println("[---------------]");
			this.length = pos+1;
			
			if(iblock>=blocks.length){
//				System.out.println("#---------------#");
				int[] tblocks = new int[iblock+1];
				for(int i=0; i<blocks.length; i++)
					tblocks[i] = blocks[i];
				blocks = tblocks;
			}
		}
		
//		int iblock = (pos / MERS_PER_BLOCK);
		int mask = (MASK << FIRST_SHIFT ) >>> ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		int shift = FIRST_SHIFT - ((pos % MERS_PER_BLOCK) * BITS_PER_MER);
		
		mask = ~mask;
		blocks[iblock] = (blocks[iblock] & mask);
		blocks[iblock] = (blocks[iblock] | (n << shift));
	}
	
	
	
	
	
	public int countBads(){
		int count = 0;
		
		int iblock = (0 / MERS_PER_BLOCK);
		int mask = (MASK << FIRST_SHIFT ) >>> ((0 % MERS_PER_BLOCK) * BITS_PER_MER);
		int shift = FIRST_SHIFT - ((0 % MERS_PER_BLOCK) * BITS_PER_MER);
		int n;
		
		for(int i=0;i<this.length; i++){
			n = (blocks[iblock] & mask) >>> shift;

			//ns[i] = C2Nucleotide.by(n);
			if(B3Nucleotide.isBAD(n)){
			//if(B3Nucleotide.by(n) == C2Nucleotide.BAD){
				count++;
			}
			
			
			mask = mask >>> BITS_PER_MER;
			shift -= BITS_PER_MER;
			
			if(mask < 7){
				mask = (MASK << FIRST_SHIFT);
				shift = FIRST_SHIFT;
				iblock++;
			}
		}
		
		return count;
	}
	
	
	
	/* CharSequence */
	
	public char charAt(int index) {
		//return B3Nucleotide.by(getB3(index)).getChar();
		return B3Nucleotide.charFor(getB2(index));
	}

	public byte[] toB2ByteArray(){
		byte[] a = new byte[this.length];
		for(int i=0; i<this.length; i++){
			a[i] = (byte)(this.getB2(i) & 0xff);
		}
		return a;
	}
	
	
	public int indexOf(int code, int from, int to){
//		System.out.println("indexOf("+code+", "+from+", "+to+")");
		
		int i = from;
		
		int iblock = (i / MERS_PER_BLOCK);
		int mask = (MASK << FIRST_SHIFT ) >>> ((i % MERS_PER_BLOCK) * BITS_PER_MER);
		int shift = FIRST_SHIFT - ((i % MERS_PER_BLOCK) * BITS_PER_MER);
		int n;
		
//		System.out.print("\t");
		for( ; i<to ; i++){
			n = (blocks[iblock] & mask) >>> shift;
//			System.out.print(B3Nucleotide.charFor(n));
			if(n == code)
				break;
			
			mask = mask >>> BITS_PER_MER;
			shift -= BITS_PER_MER;
			
			if(mask < 7){
				mask = (MASK << FIRST_SHIFT);
				shift = FIRST_SHIFT;
				iblock++;
			}
		}
		
//		System.out.println("\t<"+i);
		
		return i;
	}
	
	
	//@Override
	public int lcp(int pos_a, int pos_b) {
		int lcp = 0;
		
		int i = pos_a;
		int iblock_i = (i / MERS_PER_BLOCK);
		int mask_i = (MASK << FIRST_SHIFT ) >>> ((i % MERS_PER_BLOCK) * BITS_PER_MER);
		int shift_i = FIRST_SHIFT - ((i % MERS_PER_BLOCK) * BITS_PER_MER);
		int n_i;
		
		int j = pos_b;
		int iblock_j = (j / MERS_PER_BLOCK);
		int mask_j = (MASK << FIRST_SHIFT ) >>> ((j % MERS_PER_BLOCK) * BITS_PER_MER);
		int shift_j = FIRST_SHIFT - ((j % MERS_PER_BLOCK) * BITS_PER_MER);
		int n_j;
		
		
		while(i<this.length && j<this.length){
			n_i = (blocks[iblock_i] & mask_i) >>> shift_i;
			n_j = (blocks[iblock_j] & mask_j) >>> shift_j;
			
			if(n_i != n_j)
				break;
		
			lcp++;
		
			
			mask_i = mask_i >>> BITS_PER_MER;
			shift_i -= BITS_PER_MER;	
			if(mask_i < 7){
				mask_i = (MASK << FIRST_SHIFT);
				shift_i = FIRST_SHIFT;
				iblock_i++;
			}
			mask_j = mask_j >>> BITS_PER_MER;
			shift_j -= BITS_PER_MER;	
			if(mask_j < 7){
				mask_j = (MASK << FIRST_SHIFT);
				shift_j = FIRST_SHIFT;
				iblock_j++;
			}
			
			i++;
			j++;
		}
		
		return lcp;
	}
	
	
	public B2LLSequence reverseComplementCopy(){
		B2LLSequence rc = new B2LLSequence(this.length());
		int to = (length / 2);
    	if(length%2 != 0)
    		to++;
    	
    	int tmp;
    	for(int i=0; i<to; i++){    		
    		tmp = B3Nucleotide.rcs[getB3(i)].code();
    		rc.setB3(i, B3Nucleotide.rcs[getB3(length - i - 1)].code()  );
    		//ns[i] = rcs[ ns[ns.length - i - 1].code ];
    		rc.setB3(length -i - 1, tmp);
    		//ns[ns.length - i - 1] = tmp;
    	}
    	return rc;
	}
	
	/* ToString */
	
	@Override
	public String toString(){
		String s = "";
		for(int i=0; i<this.length(); i++)
			s += this.charAt(i);
		return s;
	}
        
        
        public String getString(){
		String s = "";
		for(int i=0; i<this.length(); i++)
			s += this.charAt(i);
		return s;
	}

		@Override
		public int compareTo(B2LLSequence o) {
			int l = length() <= o.length() ? length() : o.length();
			for(int i=0; i<l; i++){
				if(getB2(i) != o.getB2(i)){
					return getB2(i) - o.getB2(i);
				}
			}
			return 0;
		}
}
