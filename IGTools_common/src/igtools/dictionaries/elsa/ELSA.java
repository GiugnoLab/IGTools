package igtools.dictionaries.elsa;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.util.Arrays;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3Sequence;
import igtools.common.versioning.DataVersioning;

public class ELSA implements IELSA{
	public static final long serialVersionUID = DataVersioning.ELSA |	0x0000000000000001L;
	
	protected B3Sequence s = null;
	protected int r_limit;
	protected int s_length;
	
	protected int[] sa;
	protected int[] lcp;
	
	public ELSA(B3Sequence s){
		init(s);
	}
	
	protected void init(B3Sequence s){
		this.s = s;
		r_limit = s.length()-1;
		s_length = s.length();
		
		byte maxValue = B3Nucleotide.A_CODE;
		if(B3Nucleotide.C_CODE > maxValue) maxValue = B3Nucleotide.C_CODE;
		if(B3Nucleotide.G_CODE > maxValue) maxValue = B3Nucleotide.G_CODE;
		if(B3Nucleotide.T_CODE > maxValue) maxValue = B3Nucleotide.T_CODE;
		if(B3Nucleotide.N_CODE > maxValue) maxValue = B3Nucleotide.N_CODE;
		
		//Timer timer = new Timer();
		//timer.reset();
		//System.out.print("sa\t");
		
		sa = ByteSAIS.buildSuffixArray(s, 0, s.length(), maxValue);
		
		//System.out.println(timer.getElapsedSecs()+" sec.");
		//System.out.print("lcp\t");
		//timer.reset();
		
		lcp = computeLCP(s, 0, s.length(), sa);
	}
	
	public ELSA(){
		s = null;
		r_limit = 0;
		s_length = 0;
		sa = null;
		lcp = null;
	}
	public void setSequence(B3Sequence s){
		this.s = s;
	}

	@Override
	public int length() {
		return s_length;
	}
	@Override
	public B3Sequence b3seq() {
		return s;
	}
	@Override
	public int[] sa() {
		return sa;
	}
	@Override
	public int[] lcp() {
		return lcp;
	}
	
	
	/**
   * Calculate longest prefix (LCP) array for an existing suffix array and input. Index
   * <code>i</code> of the returned array indicates the length of the common prefix
   * between suffix <code>i</code> and <code>i-1<code>. The 0-th
   * index has a constant value of <code>-1</code>.
   * <p>
   * The algorithm used to compute the LCP comes from
   * <tt>T. Kasai, G. Lee, H. Arimura, S. Arikawa, and K. Park. Linear-time longest-common-prefix
   * computation in suffix arrays and its applications. In Proc. 12th Symposium on Combinatorial
   * Pattern Matching (CPM ’01), pages 181–192. Springer-Verlag LNCS n. 2089, 2001.</tt>
   */
	protected static int[] computeLCP(B3Sequence input, final int start, final int length, int[] sa)
	  {
	      final int [] rank = new int [length];
	      for (int i = 0; i < length; i++)
	          rank[sa[i]] = i;
	      int h = 0;
	      final int [] lcp = new int [length];
	      for (int i = 0; i < length; i++)
	      {
	          int k = rank[i];
	          if (k == 0)
	          {
	              lcp[k] = -1;
	          }
	          else
	          {
	              final int j = sa[k - 1];
	              while (i + h < length && j + h < length
	                  && input.getB3(start + i + h) == input.getB3(start + j + h))
	              {
	                  h++;
	              }
	              lcp[k] = h;
	          }
	          if (h > 0) h--;
	      }

	      return lcp;
	  }
	
	
	private class _Iterator implements IELSAIterator{
		private int k;
		protected int istart;//inclusive
		protected int iend;//exclusive
		protected int limit;
                private ELSA elsa;
		_Iterator(int _k, int _istart, int _iend, ELSA _elsa){
			k = _k;
			limit = s_length - k;
			istart = _istart;
			iend = _iend;
                        elsa = _elsa;
		}
		@Override
		public int istart() {
			return istart;
		}
		@Override
		public int iend() {
			return iend;
		}
		@Override
		public B3Nucleotide[] kmer() {
			B3Nucleotide[] seq = new B3Nucleotide[k];
			s.getB3(sa[istart], seq);
			return seq;
		}
                @Override
		public void kmer(B3Nucleotide[] ns) {
			s.getB3(sa[istart], ns);
		}
		@Override
		public int multiplicity() {
			return iend - istart;
		}
		@Override
		public int[] positions() {
			return Arrays.copyOfRange(sa, istart, iend);
		}
		@Override
		public int[] sortedPositions() {
			int[] pos =  Arrays.copyOfRange(sa, istart, iend);
			Arrays.sort(pos);
			return pos;
		}
		
		@Override
		public boolean next() {
			int l = iend;
			while(l<sa.length && (sa[l] > limit)){
				l++;
			}
			
			if(l >= sa.length){
				istart = -1;
				iend = 0;
				return false;
			}
			
			istart = l;
			l++;
			while(l<sa.length && lcp[l] >= k && sa[l] < limit){
				l++;
			}
			iend = l;
			return true;
		}
		@Override
		public boolean prev() {
			int l = istart;
			if(l == sa.length)
				l = sa.length-1;
			while(l>=0 && sa[l]>limit)
				l--;
			
			if(l<0){
				istart = -1;
				iend = 0;
				return false;
			}
			
			iend = l+1;
			while(l>=0 && lcp[l]>=k && sa[l]<limit){
				l--;
			}
			istart = l;
			return true;
		}
		
		@Override
		public boolean isMinimalHapax() {
			return ((iend - istart) == 1) &&
					(
						(k<2) ||
						(
							(istart==0 || lcp[istart]==k-1) 		&& 
							(iend==r_limit || lcp[iend+1]==k-1)
						)
					);
//			 if((iend - istart) == 1){
//				if(k<2){//k = 1
//					return true;
//				}
//				else{
//					return  (istart==0 || lcp[istart-1]==k-1) && (iend==r_limit || lcp[iend+1]==k-1); 
//				}
//				
//			 }
//			 //else is not an hapax
//			 return false;
			
		}
		@Override
		public boolean isGlobalMaximalRepeat() {
			return ((iend - istart) > 1) &&
					(
						(istart==0 || lcp[istart]==k) 		&&
						(iend==r_limit || lcp[iend+1]==k)
					);
		}
		@Override
		public boolean hasNext() {
			// TODO Auto-generated method stub
			return false;
		}
		@Override
		public boolean hasPrev() {
			// TODO Auto-generated method stub
			return false;
		}

        @Override
        public IELSAIterator clone() {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public int compare(IELSAIterator it) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public int k() {
            return k;
        }

        @Override
        public ELSA elsa() {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public boolean good() {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
	}

	@Override
	public IELSAIterator begin(int k) {
		_Iterator it = new _Iterator(k, -1, 0, this);
		return it;
	}

//	@Override
//	public IELSAIterator end(int k) {
//		//TODO verify limits
//		return new _Iterator(this, k, s_length, s_length+1);
//	}

	@Override
	public IELSAIterator rbegin(int k) {
		return new _Iterator(k, s_length, s_length+1,this);
	}

//	@Override
//	public IELSAIterator rend(int k) {
//		//TODO verify limits
//		return new _Iterator(this, k, -1, 0);
//	}

	@Override
	public void load(String file) throws Exception {
		DataInputStream reader = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
		long version = reader.readLong();
		if(version !=  ELSA.serialVersionUID){
			reader.close();
			throw new Exception("Wrong ELSA file version. \nExpected: "+ Long.toBinaryString(ELSA.serialVersionUID) +"\nfound: \n"+Long.toBinaryString(version));
		}
		s_length = reader.readInt();
		r_limit = s_length-1;
		sa = new int[s_length];
		lcp = new int[s_length];
		for(int i=0; i<s_length; i++){
			sa[i] = reader.readInt();
			lcp[i] = reader.readInt();
		}
		reader.close();
	}

	@Override
	public void save(String file) throws Exception {
		DataOutputStream writer = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file,false)));
		writer.writeLong(ELSA.serialVersionUID);
		writer.writeInt(s_length);
		for(int i=0; i<s_length; i++){
			writer.writeInt(sa[i]);
			writer.writeInt(lcp[i]);
		}
		writer.flush();
		writer.close();
	}

	@Override
	public void print(int k) {
		System.out.println("i\tSA\tLCP\tsuff");
		for(int i=0; i< sa.length; i++){
			System.out.print(i +"\t"+ sa[i] +"\t"+ lcp[i] +"\t");
			  
			int nn = s.length() - sa[i] < k ? s.length() - sa[i] : k;
			B3Nucleotide[] suff = new B3Nucleotide[nn];
			s.getB3(sa[i], suff);
			for(int j=0; j<suff.length; j++)
				System.out.print(suff[j]);  
			System.out.println();
		}
	}

	@Override
	public int nof(int k) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}
	@Override
	public int nof_mults(int k) {
		throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
	}

	
}
