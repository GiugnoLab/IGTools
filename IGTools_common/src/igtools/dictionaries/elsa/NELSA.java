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

public class NELSA extends ELSA implements INELSA, IQuerableELSA{
	public static final long p_serialVersionUID = DataVersioning.NELSA |	0x0000000000000001L;
	public static final long serialVersionUID = DataVersioning.NELSA |	0x0000000000000002L;
	
	protected int[] ns;
	
	public NELSA(){
		super();
		ns = null;
	}
	
	public NELSA(B3Sequence b3seq){
		init(b3seq);
	}
	
	protected void init(B3Sequence b3seq){
		super.init(b3seq);
		ns = computeNS(s,sa);
	}
	
	
//	/**
//	 * ns[i] > 0 ever
//	 * if seq[sa[i]]==N  =>  ns[i] = 1
//	 * if seq[sa[i]...sa[i]+1]==AN  =>  ns[i] = 2
//	 */
	//new version
	/**
	 * ns[i] >= 0
	 * if seq[sa[i]]==N  =>  ns[i] = 0
	 * if seq[sa[i]...sa[i]+1]==AN  =>  ns[i] = 1
	 */
	public int[] ns(){
		return ns;
	}
	
	
	private static int[] computeNS(B3Sequence input, int[] sa){
		  int[] nn = new int[input.length()];
		  
		  int pn = input.length();
		  for(int i=input.length() - 1; i>=0; i--){
			  if(input.getB3(i) == B3Nucleotide.N_CODE){
				  nn[i] = 0;
				  pn = i;
			  }
			  else{
				  nn[i] = pn - i;
			  }
		  }
		  
		  int[] fn = new int[input.length()];
		  for(int i=0; i<fn.length; i++){
			  //fn[i] = nn[sa[i]] + 1;
			  fn[i] = nn[sa[i]];//new version
		  }
		  return fn;
	  }
	
	@Override
	public void print(int k){
		  System.out.println("i\tSA\tLCP\tFN\tsuff");
		  for(int i=0; i< sa.length; i++){
			  System.out.print(i +"\t"+ sa[i] +"\t"+ lcp[i] +"\t"+ ns[i] +"\t");
			  
			  int nn = s.length() - sa[i] < k ? s.length() - sa[i] : k;
			  B3Nucleotide[] suff = new B3Nucleotide[nn];
			  s.getB3(sa[i], suff);
			  for(int j=0; j<suff.length; j++)
				  System.out.print(suff[j]);
			  
			  System.out.println();
		  }
        }
	
	private class _Iterator implements IELSAIterator{
        private NELSA nelsa;
		private int k;
		protected int istart;//inclusive
		protected int iend;//exclusive
		protected int limit;
		_Iterator(NELSA _nelsa, int _k, int _istart, int _iend){
            nelsa = _nelsa;
			k = _k;
			limit = s_length - k;
			istart = _istart;
			iend = _iend;
		}
		
		void reset(int _k, int _istart, int _iend){
			k = _k;
			limit = s_length - k;
			istart = _istart;
			iend = _iend;
		}
		
        @Override
		public int k() {
			return k;
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
//			System.out.println("------+");
//			System.out.println(istart+" "+iend);
			
			int l = iend;
			
//			System.out.println(l+"\t"+sa.length+"\t"+ns[l]+"\t"+k+"\t"+sa[l]+"\t"+limit);
			//while(l<sa.length && (ns[l]<=k || sa[l]>limit)){
			while(l<sa.length && (ns[l]<k || sa[l]>limit)){//new version
//				System.out.println(".."+l+" "+ns[l]);
				l++;
			}
			
//			System.out.println("l "+l);
			
			if(l >= sa.length){
				istart = -1;
				iend = 0;
				return false;
			}
			
			istart = l;
			l++;
			
//			System.out.println("l "+l);
			//while(l<sa.length && lcp[l]>=k && sa[l]<limit && ns[l]>k){
			while(l<sa.length && lcp[l]>=k && sa[l]<limit && ns[l]>=k){//new version
//				System.out.println(",,"+l+" "+ns[l]);
				l++;
			}
//			System.out.println("l "+l);
			
			iend = l;
			
//			System.out.println("------");
			return true;
		}
		@Override
		public boolean prev() {
			int l = istart;
                        l--;
			if(l == sa.length)
				l = sa.length-1;
			while(l>=0 && (ns[l]<k || sa[l]>limit))//new version
			//while(l>=0 && (ns[l]<=k || sa[l]>limit))
				l--;
			
			if(l<0){
				istart = -1;
				iend = 0;
				return false;
			}
//			
//			iend = l+1;
			//while(l>=0 && lcp[l]>=k && sa[l]<limit && ns[l]>k){
			while(l>=0 && lcp[l]>=k && sa[l]<limit && ns[l]>=k){//new version
				l--;
			}
//			istart = l;
                        iend = istart;
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
			//return iend < r_limit;
                    return (istart+1 < r_limit && iend < r_limit);
		}
		@Override
		public boolean hasPrev() {
			return istart > 0;
		}

        @Override
        public IELSAIterator clone() {
            return new _Iterator(nelsa, k,istart,iend);
        }

        @Override
        public int compare(IELSAIterator it) {
            return istart - it.istart();
        }

        @Override
        public ELSA elsa() {
            return nelsa;
        }

        @Override
        public boolean good() {
            //TODO
//                    System.out.println("["+istart+","+iend+","+r_limit+"]");
            return (istart<=iend && istart >=0 && iend <= r_limit+1);
        }
	}

	
	/*
	 * IELSAIterator it = nelsa.begin(k);
		while(it.next()){
		}
	 */
	@Override
	public IELSAIterator begin(int k) {
		_Iterator it = new _Iterator(this, k, -1, 0);
		return it;
	}
	@Override
	public IELSAIterator rbegin(int k) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
		//return new _Iterator(k, s_length, s_length+1);
	}
        
	@Override
	public int nof(int k) {
		int nof = 0;
		IELSAIterator it = begin(k);
		while(it.next()){
			nof++;
		}
		return nof;
	}
	
	@Override
	public int nof_mults(int k) {
		int nof = 0;
		IELSAIterator it = begin(k);
		while(it.next()){
			nof += it.multiplicity();
		}
		return nof;
	}
        
        
        
        public int[] sortedNsByPositions(IELSAIterator it){
            int[][]  data = new int[2][it.multiplicity()];
            data[0] = Arrays.copyOfRange(sa, it.istart(), it.iend());
            data[1] = Arrays.copyOfRange(ns, it.istart(), it.iend());
            
            igtools.common.util.Arrays.qsortOnFirstRow(data);
            data[0] = null;
            
            return data[1];
        }
        
        
        
        
    private int cc(int i){
		if(i>r_limit)
			return -1;
		return s.getB3(i);
	}
	private int cc(B3Sequence seq, int i){
		//if(i>seq.length())
		if(i>=seq.length())
			return -1;
		return seq.getB3(i);
	}
	
	private int cc(B3Nucleotide[] seq, int i){
		//if(i>seq.length())
		if(i>=seq.length)
			return -1;
		return seq[i].code();
	}
	
	private int rc(B3Nucleotide[] seq, int i){
		//if(i>seq.length())
		if(i>=seq.length)
			return -1;
		return B3Nucleotide.compl(seq[seq.length - (i+1)].code());
	}
        
        
        @Override
        public boolean contains(B3Sequence kmer) {
            return (find(kmer)!=null);
        }

        @Override
        public IELSAIterator find(B3Sequence q) {
            if(q.length()<1)
                return null;
            
            int l = 0;
            int r = r_limit;
            int ll,lr, rl, rr;
            int m;
            int cq;

            for(int i=0; i<q.length(); i++){
                    cq = cc(q,i);

                    if(cq < cc(sa[l] +i))
                            return null;
                    if(cq > cc(sa[r] +i))
                            return null;


                    ll = l;
                    lr = r;
                    while(true){
                            if(ll==lr)
                                    break;
                            m = ll + ((lr-ll)/2);
                            if(cq <= cc(sa[m] +i)){
                                    lr = m;
                            }
                            else if(cq > cc(sa[m] +i)){
                                    ll = m+1;
                            }
                    }
                    l = ll;
                    if(cq != cc(sa[l] +i))
                            return null;


                    rl = l;
                    rr = r;
                    while(true){
                            if(rl==rr)
                                    break;
                            m = rl + ((rr-rl)/2);
                            if((rr-rl)%2 != 0)
                                    m++;
                            if(cq < cc(sa[m] +i)){
                                    rr = m-1;
                            }
                            else if(cq >= cc(sa[m] +i)){
                                    rl = m;
                            }
                    }
                    r = rr;
                    if(cq != cc(sa[r] +i))
                            return null;

            }
            return new _Iterator(this, q.length(), l, r+1);
        }
        
        
        
        
        public IELSAIterator find(B3Nucleotide[] q) {
            if(q.length<1)
                return null;
            
            int l = 0;
            int r = r_limit;
            int ll,lr, rl, rr;
            int m;
            int cq;

            for(int i=0; i<q.length; i++){
                    cq = cc(q,i);

                    if(cq < cc(sa[l] +i))
                            return null;
                    if(cq > cc(sa[r] +i))
                            return null;


                    ll = l;
                    lr = r;
                    while(true){
                            if(ll==lr)
                                    break;
                            m = ll + ((lr-ll)/2);
                            if(cq <= cc(sa[m] +i)){
                                    lr = m;
                            }
                            else if(cq > cc(sa[m] +i)){
                                    ll = m+1;
                            }
                    }
                    l = ll;
                    if(cq != cc(sa[l] +i))
                            return null;


                    rl = l;
                    rr = r;
                    while(true){
                            if(rl==rr)
                                    break;
                            m = rl + ((rr-rl)/2);
                            if((rr-rl)%2 != 0)
                                    m++;
                            if(cq < cc(sa[m] +i)){
                                    rr = m-1;
                            }
                            else if(cq >= cc(sa[m] +i)){
                                    rl = m;
                            }
                    }
                    r = rr;
                    if(cq != cc(sa[r] +i))
                            return null;

            }
            return new _Iterator(this, q.length, l, r+1);
        }
        
        public IELSAIterator find_rc(B3Nucleotide[] q) {
            if(q.length<1)
                return null;
            
            int l = 0;
            int r = r_limit;
            int ll,lr, rl, rr;
            int m;
            int cq;

            for(int i=0; i<q.length; i++){
                    cq = rc(q,i);

                    if(cq < cc(sa[l] +i))
                            return null;
                    if(cq > cc(sa[r] +i))
                            return null;


                    ll = l;
                    lr = r;
                    while(true){
                            if(ll==lr)
                                    break;
                            m = ll + ((lr-ll)/2);
                            if(cq <= cc(sa[m] +i)){
                                    lr = m;
                            }
                            else if(cq > cc(sa[m] +i)){
                                    ll = m+1;
                            }
                    }
                    l = ll;
                    if(cq != cc(sa[l] +i))
                            return null;


                    rl = l;
                    rr = r;
                    while(true){
                            if(rl==rr)
                                    break;
                            m = rl + ((rr-rl)/2);
                            if((rr-rl)%2 != 0)
                                    m++;
                            if(cq < cc(sa[m] +i)){
                                    rr = m-1;
                            }
                            else if(cq >= cc(sa[m] +i)){
                                    rl = m;
                            }
                    }
                    r = rr;
                    if(cq != cc(sa[r] +i))
                            return null;

            }
            return new _Iterator(this, q.length, l, r+1);
        }
        
        
        public int first(B3Nucleotide[] q) {
            if(q.length<1)
                return -1;
            
            int l = 0;
            int r = r_limit;
            int ll,lr, rl, rr;
            int m;
            int cq;

            for(int i=0; i<q.length; i++){
                    cq = cc(q,i);

                    if(cq < cc(sa[l] +i))
                            return -1;
                    if(cq > cc(sa[r] +i))
                            return -1;


                    ll = l;
                    lr = r;
                    while(true){
                            if(ll==lr)
                                    break;
                            m = ll + ((lr-ll)/2);
                            if(cq <= cc(sa[m] +i)){
                                    lr = m;
                            }
                            else if(cq > cc(sa[m] +i)){
                                    ll = m+1;
                            }
                    }
                    l = ll;
                    if(cq != cc(sa[l] +i))
                            return -1;


                    rl = l;
                    rr = r;
                    while(true){
                            if(rl==rr)
                                    break;
                            m = rl + ((rr-rl)/2);
                            if((rr-rl)%2 != 0)
                                    m++;
                            if(cq < cc(sa[m] +i)){
                                    rr = m-1;
                            }
                            else if(cq >= cc(sa[m] +i)){
                                    rl = m;
                            }
                    }
                    r = rr;
                    if(cq != cc(sa[r] +i))
                            return -1;

            }
            //return new _Iterator(this, q.length, l, r+1);
            return l;
        }
        
        
        public int first_rc(B3Nucleotide[] q) {
            if(q.length<1)
                return -1;
            
            int l = 0;
            int r = r_limit;
            int ll,lr, rl, rr;
            int m;
            int cq;

            for(int i=0; i<q.length; i++){
                    cq = rc(q,i);

                    if(cq < cc(sa[l] +i))
                            return -1;
                    if(cq > cc(sa[r] +i))
                            return -1;


                    ll = l;
                    lr = r;
                    while(true){
                            if(ll==lr)
                                    break;
                            m = ll + ((lr-ll)/2);
                            if(cq <= cc(sa[m] +i)){
                                    lr = m;
                            }
                            else if(cq > cc(sa[m] +i)){
                                    ll = m+1;
                            }
                    }
                    l = ll;
                    if(cq != cc(sa[l] +i))
                            return -1;


                    rl = l;
                    rr = r;
                    while(true){
                            if(rl==rr)
                                    break;
                            m = rl + ((rr-rl)/2);
                            if((rr-rl)%2 != 0)
                                    m++;
                            if(cq < cc(sa[m] +i)){
                                    rr = m-1;
                            }
                            else if(cq >= cc(sa[m] +i)){
                                    rl = m;
                            }
                    }
                    r = rr;
                    if(cq != cc(sa[r] +i))
                            return -1;

            }
            //return new _Iterator(this, q.length, l, r+1);
            return l;
        }
        
        
        
        
        

        @Override
        public IELSAIterator findOrNear(B3Sequence kmer, boolean preferForward) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public int firstSAIndex(int pl, int pr) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public int lastSAIndex(int pl, int pr) {
        	int l = 0;
    		int r = r_limit;
    		int ll,lr, rl, rr;
    		int m;
    		int cq;
    		for(int i=0; i<pr-pl+1; i++){  			
    			cq = cc(pl + i);
    			if(cq < cc(sa[l] +i)){
    				return -1;
    			}
    			if(cq > cc(sa[r] +i)){
    				return -1;
    			}
    			
    			ll = l;
    			lr = r;
    			while(true){
    				if(ll==lr)
    					break;
    				m = ll + ((lr-ll)/2);
    				if(cq <= cc(sa[m] +i)){
    					lr = m;
    				}
    				else if(cq > cc(sa[m] +i)){
    					ll = m+1;
    				}
    			}
    			l = ll;
    			if(cq != cc(sa[l] +i)){
    				return -1;
    			}

    			rl = l;
    			rr = r;
    			while(true){
    				if(rl==rr)
    					break;
    				m = rl + ((rr-rl)/2);
    				if((rr-rl)%2 != 0)
    					m++;
    				if(cq < cc(sa[m] +i)){
    					rr = m-1;
    				}
    				else if(cq >= cc(sa[m] +i)){
    					rl = m;
    				}
    			}
    			r = rr;
    			if(cq != cc(sa[r] +i)){
    				return -1;
    			}
    			
    		}
    		return r;
        }

        @Override
        public int toUnique(int pl) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }

        @Override
        public IELSAIterator toUniqueIterator(int pl) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        
        @Override
        public int toUnique(int pl, int pr) {
        	int l = 0;
    		int r = r_limit;
    		int ll,lr, rl, rr;
    		int m;
    		int cq;
    		
    		int i=0;
    		while(true){
    		//for(int i=0; i<q.length(); i++){
    			cq = cc(pl + i);

    			if(cq < cc(sa[l] +i))
    				return -1;
    			if(cq > cc(sa[r] +i))
    				return -1;
    			

    			ll = l;
    			lr = r;
    			while(true){
    				if(ll==lr)
    					break;
    				m = ll + ((lr-ll)/2);
    				if(cq <= cc(sa[m] +i)){
    					lr = m;
    				}
    				else if(cq > cc(sa[m] +i)){
    					ll = m+1;
    				}
    			}
    			l = ll;
    			if(cq != cc(sa[l] +i))
    				return -1;
    			
    			
    			rl = l;
    			rr = r;
    			while(true){
    				if(rl==rr)
    					break;
    				m = rl + ((rr-rl)/2);
    				if((rr-rl)%2 != 0)
    					m++;
    				if(cq < cc(sa[m] +i)){
    					rr = m-1;
    				}
    				else if(cq >= cc(sa[m] +i)){
    					rl = m;
    				}
    			}
    			r = rr;
    			if(cq != cc(sa[r] +i))
    				return -1;
    			
    			if((pl+i >= pr) && r == l)
    				return i;
    			i++;
    		}
    		//return r-l+1;
    		
        }
	
	
	@Override
	public void load(String file) throws Exception {
		DataInputStream reader = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
		long version = reader.readLong();
		
		
		
		if(version == NELSA.p_serialVersionUID){
			s_length = reader.readInt();
            System.out.println("s length: "+s_length);
			r_limit = s_length-1;
			sa = new int[s_length];
			lcp = new int[s_length];
			ns = new int[s_length];
			for(int i=0; i<s_length; i++){
				sa[i] = reader.readInt();
				lcp[i] = reader.readInt();
				ns[i] = reader.readInt() -1;
			}
		}
		else if(version == NELSA.serialVersionUID){
			s_length = reader.readInt();
            System.out.println("s length: "+s_length);
			r_limit = s_length-1;
			sa = new int[s_length];
			lcp = new int[s_length];
			ns = new int[s_length];
			for(int i=0; i<s_length; i++){
				sa[i] = reader.readInt();
				lcp[i] = reader.readInt();
				ns[i] = reader.readInt();
			}
		}
		else{
			reader.close();
			throw new Exception("Wrong NELSA file version. \nExpected: "+ Long.toBinaryString(NELSA.serialVersionUID) +"\nfound: \n"+Long.toBinaryString(version));
		}
		
		reader.close();
	}

	@Override
	public void save(String file) throws Exception {
		DataOutputStream writer = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file,false)));
		writer.writeLong(NELSA.serialVersionUID);
		writer.writeInt(s_length);
		for(int i=0; i<s_length; i++){
			writer.writeInt(sa[i]);
			writer.writeInt(lcp[i]);
			writer.writeInt(ns[i]);
		}
		writer.flush();
		writer.close();
	}
	
	
	
	
	

	
	public interface LeafListener{
		public void isLeaf(int pos, int k);
	}
	public void getLeafs(LeafListener listener){
	    throw new UnsupportedOperationException("Not supported yet.");
	}

}
