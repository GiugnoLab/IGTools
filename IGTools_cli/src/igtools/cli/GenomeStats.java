package igtools.cli;

import igtools.analyses.Inform;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.NELSA;


/**
 * output:   sequence length nof_N effective_lengt length_ratio minimal_hapax_length maximal_repeat_length MF(G)  entropy(MF(G))
 * @author vbonnici
 *
 */
public class GenomeStats {

	
	public static void usage(){
		System.out.println("Usage: cmd iseq.3bit iseq.nelsa");
	}
	
	public static void main(String[] args){
		
		String a_iseq = "";
		String a_inelsa = "";
		
		try{
			a_iseq = args[0];
			a_inelsa = args[1];
		}catch(Exception  e){
			usage();
			System.exit(0);
		}
		
		
		try{

			System.out.println("seq");
			B3LLSequence b3seq = B3LLSequence.load(a_iseq);
			System.out.println("nelsa");
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			nelsa.setSequence(b3seq);
			
			System.out.println("lengths");
			int slength = b3seq.length();
			int nbads = b3seq.countBads();
			
			System.out.println("mhl");
			int mhl = Inform.mhl(nelsa);
			System.out.println("mrl");
			int mrl = Inform.Mrl(nelsa);
			System.out.println("mg");
			int mf = Inform.kCompleteness(nelsa) + 1;
			System.out.println("e1");
			double e1 = Inform.entropy(nelsa, mf);
			System.out.println("e2");
			double e2 = Inform.entropy(nelsa, mf-1);
			
			System.out.print("# ALL " + a_iseq +"\t"+ slength +"\t"+nbads +"\t"+ (slength-nbads) +"\t"+ (((double) (slength-nbads) )/((double) slength )) +"\t");
			System.out.print(mhl+"\t"+mrl+"\t"+mf+"\t"+e1+"\t"+e2);
			System.out.println();
			
			System.out.println("|G| "+slength);
			System.out.println("|N| "+nbads);
			System.out.println("|G-N| "+(slength - nbads));
			System.out.println("|G-N| / |G| "+((((double) (slength-nbads) )/((double) slength ))));
			System.out.println("MHL "+mhl);
			System.out.println("MRL "+mrl);
			System.out.println("MFL "+mf);
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
		
	}
}
