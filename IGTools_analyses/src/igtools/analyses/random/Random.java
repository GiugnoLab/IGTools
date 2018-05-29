package igtools.analyses.random;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;



/**
 * Random sequences generators.
 * 
 * @author vbonnici
 *
 */
public class Random {
	
	
	public static B3Sequence generatePseudoUniform(int length){
		B3Sequence b3seq = new B3LLSequence(length);
		int c;
		for(int i=0; i<length; i++){
			c = (int)(Math.random()*4);
			b3seq.setB3(i, c);
		}
		return b3seq;
	}
	
	
	public static B3Sequence generateUniform(int length){
		B3Sequence b3seq = new B3LLSequence(length);
		
		//that's to make sure that the a-posterior probability of each nucleobase is exactly length/4
		B3Sequence urn = new B3LLSequence(length);
		
		int i;
		for(i=0; i<length; i++)
			urn.setB3(i, i%4);
		
		/*int i = 0; int l = length/4; 
		for(; i<l; i++)			urn.setB3(i, B3Nucleotide.A_CODE);		
		l = length/2;
		for(; i<l; i++)			urn.setB3(i, B3Nucleotide.C_CODE);		
		l = length*3/4;
		for(; i<l; i++)			urn.setB3(i, B3Nucleotide.G_CODE);
		for(; i<length; i++)	urn.setB3(i, B3Nucleotide.T_CODE);
		*/
		
		int c;
		for(i=0; i<length; i++){
			c = (int)(Math.random()*(length-i));
			b3seq.setB3(i, urn.getB3(c));
			urn.setB3(c, urn.getB3(length -i -1));
		}
		return b3seq;
	}
	
	public static B3Sequence generateUniformShuffled(int length){
		B3Sequence b3seq = new B3LLSequence(length);
		
		//that's to make sure that the a-posterior probability of each nucleobase is exactly length/4
		B3Sequence urn = new B3LLSequence(length);
		int i = 0; int l = length/4; 
		for(; i<l; i++)			urn.setB3(i, B3Nucleotide.A_CODE);		l += length/4;
		for(; i<l; i++)			urn.setB3(i, B3Nucleotide.C_CODE);		l += length/4;
		for(; i<l; i++)			urn.setB3(i, B3Nucleotide.G_CODE);
		for(; i<length; i++)	urn.setB3(i, B3Nucleotide.T_CODE);
		
		int c;
		int p;
		
		for(int j=0; j<10; j++){
			for(i=0; i<length; i++){
				c = urn.getB3(i);
				p = (int)(Math.random()*(length));
				urn.setB3(i, urn.getB3(p));
				urn.setB3(p, c);
			}
		}
		
		for(i=0; i<length; i++){
			c = (int)(Math.random()*(length-i));
			b3seq.setB3(i, urn.getB3(c));
			urn.setB3(c, urn.getB3(length -i -1));
		}
		return b3seq;
	}
	
	
	public static void replaceNsUniformly(B3Sequence b3seq){
		int nCount = b3seq.countBads();
		if(nCount > 0){
			B3Sequence urn = new B3LLSequence(nCount);
			int i = 0; int l = nCount/4; 
			for(; i<l; i++)			urn.setB3(i, B3Nucleotide.A_CODE);		l += nCount/4;
			for(; i<l; i++)			urn.setB3(i, B3Nucleotide.C_CODE);		l += nCount/4;
			for(; i<l; i++)			urn.setB3(i, B3Nucleotide.G_CODE);
			for(; i<nCount; i++)	urn.setB3(i, B3Nucleotide.T_CODE);
			
			int count = 0;
			int nSub = 0;
			int c;
			for(i=0; i<b3seq.length(); i++){
				if(b3seq.getB3(i) == B3Nucleotide.N_CODE){
					count++;
					c = (int)(Math.random()*(nCount - nSub));
					b3seq.setB3(i, urn.getB3(c));
					urn.setB3(c, urn.getB3(nCount -nSub -1));
				}
			}
			System.out.println("count" + count);
		}
		
	}
}
