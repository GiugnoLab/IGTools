package igtools.dictionaries.elsa;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;

import igtools.common.sequence.B3Sequence;
import igtools.common.versioning.DataVersioning;

//TODO
public class OnDiskELSA implements IELSA{
	public static final long serialVersionUID = DataVersioning.OnDiskELSA |	0x0000000000000001L;
	
	//private final int COLUMNS_PER_ROW = 2;
	//private int block_size = 128 * 1024 * 1024;//128 MB in Bytes
	
	private String file;
	private FileInputStream fis;
	private BufferedInputStream bis;
	private DataInputStream dis;
	//DataInputStream reader = new DataInputStream(new BufferedInputStream(new FileInputStream(file)));
	
	private B3Sequence b3seq;
	private int s_length;
	
	private int[] sa_block;
	private int[] lcp_block;
	
	public OnDiskELSA(B3Sequence b3seq, String ELSAFile) throws Exception{
		this.b3seq = b3seq;
		
		this.file = ELSAFile;
		fis = new FileInputStream(this.file);
		bis = new BufferedInputStream(fis);
		dis = new DataInputStream(bis);
		
		long v = dis.readLong();
		if(v !=  ELSA.serialVersionUID){
			dis.close();
			throw new Exception("Wrong ELSA file version, readed by OnDiskELSA. \nExpected: "+ Long.toBinaryString(ELSA.serialVersionUID) +"\nfound: \n"+Long.toBinaryString(v));
		}
		
		s_length = dis.readInt();
		if(s_length != b3seq.length()){
			dis.close();
			throw new Exception("Found wrong lengths between ELSA and B3Sequence.\nELSA: "+s_length+"\nB3Sequence: "+b3seq.length());
		}
	}

	@Override
	public int length() {
		return s_length;
	}

	@Override
	public B3Sequence b3seq() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int[] sa() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int[] lcp() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public IELSAIterator begin(int k) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public IELSAIterator rbegin(int k) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void load(String file) throws Exception {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void save(String file) throws Exception {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void print(int k) {
		// TODO Auto-generated method stub
		
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
