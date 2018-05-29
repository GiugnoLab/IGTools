package igtools.dictionaries.elsa;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3Sequence;

import java.io.RandomAccessFile;

public class DiskNELSATreeNavigator {

	private B3Sequence b3seq;
	private RandomAccessFile raf = null;
	
	private static long p_start = 12l;
	private int sa_start = 0;//inclusive
	private int sa_end = 0;//inclusive
	private int k;
	
	public DiskNELSATreeNavigator(B3Sequence b3seq, RandomAccessFile raf, int start, int end, int k){
		this.b3seq = b3seq;
		this.raf = raf;
		this.sa_start = start;
		this.sa_end = end;
		this.k = k;
	}
	
	public DiskNELSATreeNavigator clone(){
		return new DiskNELSATreeNavigator(b3seq, raf, sa_start, sa_end, k);
	}
	
	
	public void close() throws Exception{
		if(this.raf != null){
			raf.close();
		}
	}
	
	
	public long multiplicity(){
		return sa_end - sa_start + 1;
	}
	public int k(){
		return k;
	}
	public int sa_start(){
		return sa_start;
	}
	public int sa_end(){
		return sa_end;
	}
	
	private int getInnerSA(int p) throws Exception{
		//raf.seek(p_start + (this.sa_start * 12l) + (p * 12l));
		raf.seek(p_start + (p * 12l));
		return raf.readInt();
	}
	private int getInnerChildCode(int p) throws Exception{
		return b3seq.getB3(getInnerSA(p) + k);
	}
	
	
	
	public DiskNELSATreeNavigator getChild(int code) throws Exception{
//		public IELSAIterator find(B3Sequence q) {
//		int l = 0;
//		int r = r_limit;
//		int ll,lr, rl, rr;
//		int m;
//		int cq;
//		cq = cc(q,i);
//	    if(cq < cc(sa[l] +i))return null;
//	    if(cq > cc(sa[r] +i))return null;
		
		int l = sa_start;
		int r = sa_end;
		//long m = (long)Math.floor((r-l)/2);
		
		int c_l = getInnerChildCode(l); if(c_l > code) return null;
		int c_r = getInnerChildCode(r); if(c_r < code) return null;
		
		//int l_m = getInnerChildCode(m);
		int c_m;
		
		
		int m;
		int ll = l;
		int rr = r;
		while(true){
			if(ll == rr) break;
			m = ll + ((rr-ll)/2);
			c_m = getInnerChildCode(m);
			if(code <= c_m) rr = m;
			//else if(code > c_m) ll = m+1;
			else ll = m+1;
		}
		l = ll;
		c_l = getInnerChildCode(l);
		if(code != c_l) return null;
		
//		ll = l;
//	    lr = r;
//	    while(true){
//	            if(ll==lr)
//	                    break;
//	            m = ll + ((lr-ll)/2);
//	            if(cq <= cc(sa[m] +i)){
//	                    lr = m;
//	            }
//	            else if(cq > cc(sa[m] +i)){
//	                    ll = m+1;
//	            }
//	    }
//	    l = ll;
//	    if(cq != cc(sa[l] +i)) return null;
		
		
		
//		long rl = l;
//		long rr = r;
		ll = l;
		rr = r;
		while(true){
			if(ll == rr) break;
			m = ll + ((rr-ll)/2);
			if((rr-ll)%2 != 0) m++;
			c_m = getInnerChildCode(m);
			if(code < c_m) rr = m-1;
			else ll = m;
		}
		r = rr;
		c_r = getInnerChildCode(r);
		if(code != c_r) return null;
		
//		rl = l;
//	    rr = r;
//	    while(true){
//	            if(rl==rr)
//	                    break;
//	            m = rl + ((rr-rl)/2);
//	            if((rr-rl)%2 != 0)
//	                    m++;
//	            if(cq < cc(sa[m] +i)){
//	                    rr = m-1;
//	            }
//	            else if(cq >= cc(sa[m] +i)){
//	                    rl = m;
//	            }
//	    }
//	    r = rr;
//	    if(cq != cc(sa[r] +i)) return null;
		
//		return new _Iterator(this, q.length(), l, r+1);
		
//		System.out.println("["+B3Nucleotide.charFor(getInnerChildCode(l))+","+B3Nucleotide.charFor(getInnerChildCode(r))+"]");
		return new DiskNELSATreeNavigator(b3seq, raf, l, r, k+1);
	}
	
	public static DiskNELSATreeNavigator begin(B3Sequence b3seq, String file) throws Exception{
		RandomAccessFile nelsaRAF = new RandomAccessFile(file, "r");
//		System.out.println("nof bytes: "+nelsaRAF.length());
		long version = nelsaRAF.readLong();
		if(version !=  NELSA.serialVersionUID){
			nelsaRAF.close();
			throw new Exception("Wrong NELSA file version. \nExpected: "+ Long.toBinaryString(NELSA.serialVersionUID) +"\nfound: \n"+Long.toBinaryString(version));
		}
		int nelsaLength = nelsaRAF.readInt();
		return new DiskNELSATreeNavigator(b3seq, nelsaRAF, 0, nelsaLength - 1, 0);
	}
	
}
