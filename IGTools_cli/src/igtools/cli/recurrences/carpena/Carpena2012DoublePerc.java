package igtools.cli.recurrences.carpena;

import java.util.Arrays;

import igtools.analyses.recurrences.carpena.ExtractByDistrComp;
import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;


/**
 * Run the Caprena1012 algorithm to assign a clustering coefficient to words in D_k
 * 
 * @author vbonnici
 *
 */
public class Carpena2012DoublePerc {

	
	public static void usage(){
		System.out.println("Usage: cmd sequence.3bit nelsa.nelsa high_p_thr low_p_thr k");
	}
	
	
	
	
	public static double getPercentile(double percentile, double[] ovalues){
		//double[] ovalues = Arrays.copyOf(values, values.length);
		//Arrays.sort(ovalues);
		
		int n = ovalues.length;
		for(int i=ovalues.length-1; i>=0; i--){
			if(!Double.isNaN(ovalues[i]))
				break;
			n--;
		}
		
		int p = (int)((percentile * ((double)n)) + (0.5));

		return ovalues[n -p -1];
	}
	
	public static void main(String[] args){
		String a_iseq = "";
		String a_inelsa = "";
		String a_type = "";
		double a_h_thr = 0.0;
		double a_l_thr = 0.0;
		int k = 0;
		
		
		try{
			a_iseq = args[0];
			a_inelsa = args[1];
			
			a_h_thr = Double.parseDouble(args[2]);
			if(a_h_thr < 0 || a_h_thr>1)
				throw new Exception("percentile must be in ]0,1[");
			
			a_l_thr = Double.parseDouble(args[3]);
			if(a_l_thr < 0 || a_l_thr>1)
				throw new Exception("percentile must be in ]0,1[");
			
			
			k = Integer.parseInt(args[4]);
			if(k<1)
				throw new Exception("k must be > 0");
			
			
		}catch(Exception e){
			e.printStackTrace();
			usage();
			System.exit(1);
		}
		
		try{
			
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(a_iseq);
			System.out.println("done");
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(a_inelsa);
			System.out.println("done");
			nelsa.setSequence(b3seq);
			
			
			ExtractByDistrComp ebdc = new ExtractByDistrComp(b3seq, nelsa);
			
			IELSAIterator it = nelsa.begin(k);
			int nofk = 0;
			while(it.next())	if(it.multiplicity() > 1)	nofk++;
			
			igtools.analyses.recurrences.carpena.Carpena2009 c2009 = new igtools.analyses.recurrences.carpena.Carpena2009(nelsa);
			double[] values = c2009.get_sigma_nors(k);
			
			
			int it_i = 0;
			
			
			
			double[] ovalues = Arrays.copyOf(values, values.length);
			Arrays.sort(ovalues);
			
			if(ovalues.length > 1){
				
				if(ovalues.length == 2){
					it_i = 0;
					it = nelsa.begin(k);
					while(it.next()){
						if(it.multiplicity() > 1){
								System.out.println(B3Nucleotide.toString(it.kmer()) +"\t"+ values[it_i]);
								System.out.println("S: "+ B3Nucleotide.toString(it.kmer()));
							it_i++;
						}
					}
				}
				else{
					double hp = getPercentile(a_h_thr, ovalues);
					double lp = getPercentile(a_l_thr, ovalues);
					
					System.out.println("hp "+hp);
					System.out.println("lp "+lp);
					
					it_i = 0;
					it = nelsa.begin(k);
					while(it.next()){
						if(it.multiplicity() > 1){
							//System.out.println(B3Nucleotide.toString(it.kmer()) +"\t"+ values[it_i]);
							if(values[it_i] >= lp  &&  values[it_i] <= hp){
								System.out.println(B3Nucleotide.toString(it.kmer()) +"\t"+ values[it_i]);
								System.out.println("S: "+ B3Nucleotide.toString(it.kmer()));
							}
							it_i++;
						}
					}
				}
			}
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}
