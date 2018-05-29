package igtools.common.streams;

import java.io.BufferedInputStream;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;

import igtools.common.exceptions.*;
import igtools.common.kmer.b3.unit.B3UnitLLKmer;
import igtools.common.sequence.B3LLSequence;
import igtools.common.nucleotide.B3Nucleotide;

public class B3LLStream implements INucleotideStream<B3Nucleotide>{
	
	//public static int NULL_CODE = B3LLKmer.ALL_ONES;
	
	private DataInputStream reader = null;
	private int iblock = -1;
	private int ielement = -1;
	private int ipos = 0;
	private int slength = 0;
	
	private int icode = NULL_CODE;
	
	
	private static int[] masks;
	private static int[] shifts;
	static{
		masks = new int[B3UnitLLKmer.MERS_PER_BLOCK];
		shifts = new int[B3UnitLLKmer.MERS_PER_BLOCK];
		int mask = B3UnitLLKmer.MASK << B3UnitLLKmer.FIRST_SHIFT;
		int shift = B3UnitLLKmer.FIRST_SHIFT;
		for(int i=0;i<masks.length;i++){
			masks[i] = mask;
			shifts[i] = shift;
			mask = mask >>> B3UnitLLKmer.BITS_PER_MER;
			shift -= B3UnitLLKmer.BITS_PER_MER;
		}
	}
	
	
	public B3LLStream(String ifile) throws StreamException, DataVersioningException, IOException{
		this.reader = new DataInputStream(new BufferedInputStream(new FileInputStream(ifile)));
		
		long version = reader.readLong();
		if(version !=  B3LLSequence.serialVersionUID){
			reader.close();
			throw new DataVersioningException(B3LLSequence.serialVersionUID, version);
		}
		
		slength = reader.readInt();
		iblock = reader.readInt();
		ielement = -1;
		ipos = -1;
	}
	
	private void nextICode() throws StreamException{
		if(ipos < slength-1){
			ipos++;
			try{
				ielement++;
				if(ielement >= masks.length){
					ielement = 0;
					iblock = reader.readInt();
				}
				
				icode = (iblock & masks[ielement])>>>shifts[ielement];
				
			}catch(Exception e){
				StreamException se = new StreamException("Exception on reading from stream ["+e.getMessage()+"]");
				se.addSuppressed(e);
				throw se;
			}
		}
		else
			icode = NULL_CODE;
	}
	
	@Override
	public int nextCode() throws StreamException {
		nextICode();
		return icode;
	}

	@Override
	public B3Nucleotide next() throws StreamException {
		nextICode();
		if(icode != NULL_CODE)
			return B3Nucleotide.by(icode);
		return null;
	}

	@Override
	public int position(){
		return ipos;
	}

	@Override
	public void seek(int position) throws StreamException {
		//throw new StreamException("Unimplemented function");
		if(position >= slength)
			throw new StreamException("Position is greater than sequence length");
		while(ipos < position && ipos < slength){
			//TODO
			nextICode();
		}
	}

	@Override
	public int length(){
		return slength;
	}

	@Override
	public void close() throws StreamException {
		try{
			reader.close();
		}catch(IOException e){
			StreamException se = new StreamException("Exception on closing stream ["+e.getMessage()+"]");
			se.addSuppressed(e);
			throw se;
		}
	}

}
