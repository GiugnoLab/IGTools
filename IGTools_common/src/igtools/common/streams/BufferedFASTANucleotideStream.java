package igtools.common.streams;

import igtools.common.exceptions.StreamException;
import igtools.common.nucleotide.B3Nucleotide;

import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;

public class BufferedFASTANucleotideStream implements INucleotideStream<B3Nucleotide>{
	
	private final static int BUFFER_SIZE = 1024;
	
	private InputStreamReader ireader;
	
	private char[] buffer = new char[BUFFER_SIZE];
	private char read;
	private int bufferReaded = -1;
	private int bufferPointer = BUFFER_SIZE;
	
	private int position = -1;
	
	public BufferedFASTANucleotideStream(String file) throws Exception{
		ireader = new InputStreamReader(new FileInputStream(file));
	}
	public BufferedFASTANucleotideStream(InputStreamReader istream) throws Exception{
		ireader = istream;
	}

	@Override
	public B3Nucleotide next() throws StreamException {
		B3Nucleotide n;
		while(nextChar()){
			
			//skip comments' line
			if(read=='>' || read==';'){
				while(nextChar() && read!='\n' && read!='\r');
				if(bufferReaded == -1)
					return null;
			}
			
			n = B3Nucleotide.by(read);
			
			if(n != null){
//				System.out.println("+"+read+">"+n+"-");
				position++;
				return n;
			}
			//else skip current char (nor ATCG or other FAST code) and carry on
		}
		return null;
	}

	@Override
	public int nextCode() throws StreamException {
		int n;
		while(nextChar()){
			//skip comments' line
			if(read=='>' || read==';'){
				while(nextChar() && read!='\n' && read!='\r');
				if(bufferReaded == -1)
					return NULL_CODE;
					//return B3Nucleotide.NULL_CODE;
			}
			
//			System.out.println("r["+read+"]");
			
			n = B3Nucleotide.codeFor(read);
			
			if(n != B3Nucleotide.NULL_CODE){
//				System.out.println("+"+read+">"+B3Nucleotide.by(n)+"-");
				position++;
				return n;
			}
			//else skip current char (nor ATCG or other FAST code) and carry on
		}
		return NULL_CODE;
		//return  B3Nucleotide.NULL_CODE;
	}

	@Override
	public int position() {
		return position;
	}

	@Override
	public void seek(int position) throws StreamException {
		while(this.position != position  &&  nextCode()!=NULL_CODE){
		}
	}

	@Override
	public int length() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public void close() throws StreamException {
		try{
			ireader.close();
		}catch(IOException e){
			StreamException se = new StreamException("Exception on closing stream ["+e.getMessage()+"]");
			se.addSuppressed(e);
			throw se;
		}
	}
	
	
	private boolean nextChar(){
		if(bufferPointer == BUFFER_SIZE){
//			System.out.println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@òò");
			try{
				bufferPointer = 0;
				bufferReaded = ireader.read(buffer);
				if(bufferReaded == -1)
					return false;
				read = buffer[0];
			}catch(Exception e){
				e.printStackTrace();
				System.out.println(e);
			}
		}
		else if(bufferPointer == bufferReaded){
			bufferReaded = -1;
			return false;
		}
		
		read = buffer[bufferPointer];
		bufferPointer++;
		return true;
	}
}
