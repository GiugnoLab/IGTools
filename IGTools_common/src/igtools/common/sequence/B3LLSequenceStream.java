package igtools.common.sequence;

import igtools.common.exceptions.StreamException;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.streams.INucleotideStream;

public class B3LLSequenceStream implements INucleotideStream<B3Nucleotide>{

	private B3LLSequence b3seq = null;
	private int pos;
	private int length;
	
	public B3LLSequenceStream(B3LLSequence b3seq){
		this.b3seq = b3seq;
		this.length = b3seq.length();
		this.pos = 0;
	}
	
	
	
	@Override
	public B3Nucleotide next() throws StreamException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int nextCode() throws StreamException {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int position() {
		return pos;
	}

	@Override
	public void seek(int position) throws StreamException {
		if(position > pos){
			while(pos < position && pos < length)
				pos++;
		}	
	}

	@Override
	public int length() {
		return length;
	}

	@Override
	public void close() throws StreamException {
		// TODO Auto-generated method stub
		
	}

	
}
