package igtools.dictionaries.elsa;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.RandomAccessFile;
import java.util.Arrays;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.util.Timer;

public class OnDiskNELSAIteratorV2 implements IELSAIterator {

	public int k;
	
//	private boolean onDiskSeq = true;
	private B3LLSequence b3seq = null;
//	private DataInputStream seqStream = null;
	
	private DataInputStream nelsaStream = null;
	private RandomAccessFile nelsaRAF = null;
	
	private int seqLength = 0;
	private int nelsaLength = 0;
	
	public B3Nucleotide[] kmer;
	
	protected int sa_istart;
	protected int istart;//inclusive
	protected int iend;//exclusive
	protected int limit;
	protected int r_limit;
	
	private long seek_pos = 0;
	
	private int c_sa=-1, c_lcp=-1, c_ns=-1;

	public OnDiskNELSAIteratorV2(B3LLSequence seq, String nelsaFile, int k) throws Exception{
		this.k = k;
		b3seq = seq;
		seqLength = b3seq.length();
		initNELSA(nelsaFile);
		
		if(seqLength != nelsaLength)
			throw new Exception("Discarding seq and nelsa length: "+seqLength+" "+nelsaLength);
		
		r_limit = seqLength - 1;
		limit = seqLength - k;
		istart = -1;
		sa_istart = -1;
		iend = 0;
	}
	
	public void close() throws Exception{
		if(nelsaStream != null)
			nelsaStream.close();
		if(nelsaRAF != null)
			nelsaRAF.close();
	}
	
	public B3LLSequence b3seq(){
		return b3seq;
	}
	
	

	private void initNELSA(String nelsaFile) throws Exception{
		nelsaRAF = new RandomAccessFile(nelsaFile, "r");
		long version = nelsaRAF.readLong();
		if(version !=  NELSA.serialVersionUID){
			nelsaRAF.close();
			throw new Exception("Wrong NELSA file version. \nExpected: "+ Long.toBinaryString(NELSA.serialVersionUID) +"\nfound: \n"+Long.toBinaryString(version));
		}
		nelsaLength = nelsaRAF.readInt();
		seek_pos = nelsaRAF.getFilePointer();
		nelsaStream = new DataInputStream(new BufferedInputStream(new FileInputStream(nelsaRAF.getFD())));
		if(nelsaLength>0)
			nextRow();
	}

	@Override
	public IELSA elsa() {
		throw new UnsupportedOperationException("Not supported yet.");
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
	public int k() {
		return k;
	}

	@Override
	public B3Nucleotide[] kmer() {
		B3Nucleotide[] seq = new B3Nucleotide[k];
		b3seq.getB3(sa_istart, seq);
		return seq;
	}

	@Override
	public void kmer(B3Nucleotide[] ns) {
		b3seq.getB3(sa_istart, ns);
	}

	@Override
	public int multiplicity() {
		return iend - istart;
	}

	@Override
	public int[] positions() {
		if(multiplicity() > 1){
			try{
				nelsaRAF.seek(seek_pos + (3*4*istart));
				nelsaStream = new DataInputStream(new BufferedInputStream(new FileInputStream(nelsaRAF.getFD())));
				
				int[] poss = new int[multiplicity()];
				
				for(int i=0; i<poss.length; i++){
					nextRow();
					poss[i] = c_sa;
				}
				
				if(iend < r_limit)
					nextRow();
				
				return poss;
			}catch(Exception e){
				e.printStackTrace();
				System.out.println(e);
				return null;
			}
		}
		else{
			int[] poss = new int[1];
			poss[0] = sa_istart;
			return poss;
		}
	}

	@Override
	public int[] sortedPositions() {
		int[] poss = positions();
		Arrays.sort(poss);
		return poss;
	}

	@Override
	public boolean isMinimalHapax() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public boolean isGlobalMaximalRepeat() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Not supported yet.");
	}
	
	
	public void nextRow() throws Exception{		
		c_sa = nelsaStream.readInt();
		c_lcp = nelsaStream.readInt();
		c_ns = nelsaStream.readInt();
	}
	public void prevRow() throws Exception{
		if(istart > 0){
			nelsaRAF.seek(seek_pos + (3*4*(istart-1)));
			nelsaStream = new DataInputStream(new BufferedInputStream(new FileInputStream(nelsaRAF.getFD())));
			c_sa = nelsaStream.readInt();
			c_lcp = nelsaStream.readInt();
			c_ns = nelsaStream.readInt();
		}
	}
	
	@Override
	public boolean next(){
		
		if(iend < nelsaLength){
		try{
			int l = iend;
			while(l<nelsaLength && (c_ns<k || c_sa>limit)){//new version
				l++;
				if(l<nelsaLength) nextRow();
			}
			if(l >= nelsaLength){
				istart = -1;
				iend = 0;
				return false;
			}
			sa_istart = c_sa;
			istart = l;
			l++;
			if(l<nelsaLength) nextRow();
			while(l<nelsaLength && c_lcp>=k && c_sa<limit && c_ns>=k){//new version
				l++;
				if(l<nelsaLength) nextRow();
			}
			iend = l;
			
			//if(c_lcp >= k && c_ns < k)
//			if(c_ns == 0)
//				return false;
			return true;
			
//			int l = iend;
//			if(l<nelsaLength){
//				//TODO
//				//while(l<nelsaLength && (c_ns<=k || c_sa>limit)){
//				while(l<nelsaLength && (c_ns<k || c_sa>limit)){
//					l++;
//					if(l<nelsaLength)
//						nextRow();
//				}
//			}
//			
//			if(l >= nelsaLength){
//				istart = -1;
//				iend = 0;
//				return false;
//			}
//			
//			istart = l;
//			sa_istart = c_sa;
//			l++;
//			if(l<nelsaLength){
//				nextRow();
//				//TODO
//				//while(l<nelsaLength && (c_lcp>=k && c_sa<limit && c_ns>k)){
//				while(l<nelsaLength && (c_lcp>=k && c_sa<limit && c_ns>=k)){
//					l++;
//					if(l<nelsaLength)
//						nextRow();
//				}
//			}
//			iend = l;
//			return true;
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
			return false;
		}
		}
		return false;
	}

	@Override
	public boolean prev() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public boolean good() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public boolean hasNext() {
		 return (istart+1 < r_limit && iend < r_limit);
		// TODO Auto-generated method stub
		//throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public boolean hasPrev() {
		// TODO Auto-generated method stub
		throw new UnsupportedOperationException("Not supported yet.");
	}

	@Override
	public int compare(IELSAIterator it) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	@Override
	public OnDiskNELSAIteratorV2 clone(){
		throw new UnsupportedOperationException("Not supported yet.");
	}
	
	
	
	public static boolean equals(B3Nucleotide[] ns1, B3Nucleotide[] ns2){
		for(int i=0; i<ns1.length; i++){
			if(ns1[i].code() != ns2[i].code())
				return false;
		}
		return true;
	}
	public static boolean equals(int[] pos1, int[] pos2){
		for(int i=0; i<pos1.length; i++){
			if(pos1[i]!= pos2[i])
				return false;
		}
		return true;
	}
	public static void println(int[] pos){
		for(int i=0; i<pos.length; i++)
			System.out.print(pos[i]+" ");
		System.out.println();
	}
	
	
	public static void main(String[] args){
		String iseq = null;
		String inelsa = null;
		int k = 0;
		try{
			iseq =args[0];
			inelsa = args[1];
			k = Integer.parseInt(args[2]);
		}catch(Exception e){
		}
		
		try{
			Timer timer = new Timer();
			
			System.out.println(iseq);
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(iseq);
			System.out.println("done "+timer.getElapsedSecs() +"sec.\n");
			
			timer.reset();
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(inelsa);
			System.out.println("done "+timer.getElapsedSecs()+" sec.");
			nelsa.setSequence(b3seq);
			//nelsa.print(k);
			
			IELSAIterator nit = nelsa.begin(k);
			OnDiskNELSAIteratorV2 dit = new OnDiskNELSAIteratorV2(b3seq, inelsa, k);
			
			//int count = 0;
			
			B3Nucleotide[] nns = new B3Nucleotide[k];
			B3Nucleotide[] dns = new B3Nucleotide[k];
			int[] npos;
			int[] dpos;
			
			while(nit.next()){
				nit.kmer(nns);
				npos=nit.positions();
				
				//System.out.println("nit ["+nit.istart()+","+nit.iend()+"] "+ B3Nucleotide.toString(nit.kmer()));
				//println(npos);
				
				if(!dit.next()){
					System.out.println("RAM NELSA is longer than Disk NELSA");
					System.out.println("last kmer: "+ B3Nucleotide.toString(nit.kmer()));
					System.out.println("RAM istart: "+nit.istart());
					System.out.println("RAM iend: "+nit.iend());
					System.out.println("Disk istart: "+dit.istart());
					System.out.println("Disk iend: "+dit.iend());
					break;
				}
				dit.kmer(dns);
				dpos=dit.positions();
				if(!equals(nns, dns)  || (nit.istart()!=dit.istart()) || (nit.iend()!=dit.iend())){
					System.out.println("RAM NELSA is different from Disk NELSA");
					System.out.println("last kmer: "+ B3Nucleotide.toString(nit.kmer()));
					System.out.println("RAM istart: "+nit.istart());
					System.out.println("RAM iend: "+nit.iend());
					System.out.println("last kmer: "+ B3Nucleotide.toString(dit.kmer()));
					System.out.println("Disk istart: "+dit.istart());
					System.out.println("Disk iend: "+dit.iend());
					break;
				}
				else if(!equals(npos, dpos)){
					System.out.println("RAM NELSA pos are different from Disk NELSA");
					System.out.println("last kmer: "+ B3Nucleotide.toString(nit.kmer()));
					System.out.println("RAM istart: "+nit.istart());
					System.out.println("RAM iend: "+nit.iend());
					println(npos);
					System.out.println("last kmer: "+ B3Nucleotide.toString(dit.kmer()));
					System.out.println("Disk istart: "+dit.istart());
					System.out.println("Disk iend: "+dit.iend());
					println(dpos);
					break;
				}
				else{
					//System.out.println("ok ["+nit.istart()+","+nit.iend()+"] ["+dit.istart()+","+dit.iend()+"]");
				}
				//System.out.println("-");
			}
			
			if(dit.next()){
				System.out.println("RAM NELSA is shorter than Disk NELSA");
				System.out.println("last kmer: "+ B3Nucleotide.toString(dit.kmer()));
				System.out.println("RAM istart: "+nit.istart());
				System.out.println("RAM iend: "+nit.iend());
				System.out.println("Disk istart: "+dit.istart());
				System.out.println("Disk iend: "+dit.iend());
			}
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}
