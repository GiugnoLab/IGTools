package igtools.common.util;

import igtools.common.nucleotide.B2Nucleotide;

import java.util.Random;

public class NucleotideGenerator {

	public static enum Distribution{
		UNIFORM, 
		DUO
	}
	
	private static int[] duomap = {
		B2Nucleotide.A.code(),
		B2Nucleotide.A.code(),
		B2Nucleotide.T.code(),
		B2Nucleotide.T.code(),
		B2Nucleotide.C.code(),
		B2Nucleotide.G.code()
		};//A,A,T,T,C,G
	
	private static int[] unifmap = {
		B2Nucleotide.A.code(),//<-A<-0 
		B2Nucleotide.C.code(),//<-C<-1 
		B2Nucleotide.G.code(),//<-G<-2 
		B2Nucleotide.T.code() //<-T<-3
	};
	
	
	private Distribution dist;
	private Random r;
	
	public NucleotideGenerator(Distribution dist){
		this.dist = dist;
		r = new Random();
	}
	
	
	public int nextInt(int upper){
		return r.nextInt(upper);
	}
	public int nextInt(){
		return nextInt(this.dist);
	}
	
	public int nextInt(Distribution dist){
		switch(dist){
		case UNIFORM:
			return nextIntUniform();
		case DUO:
			return nextIntDuo();
		}
		return -1;
	}
	
	public int nextIntUniform(){
		return r.nextInt(4);
	}
	
	public int nextIntDuo(){
		return duomap[r.nextInt(6)];
	}
	
	public B2Nucleotide nextN(){
		return B2Nucleotide.by(unifmap[nextInt()]);
	}
	
}
