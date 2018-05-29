package igtools.cli.coverage;

import igtools.common.sequence.B3LLSequence;
import igtools.common.util.AvgSD;
import igtools.common.util.Timer;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

import java.util.Arrays;
import java.util.Set;
import java.util.TreeSet;


/**
 * get the coverage of H_k(G) and R_k(G) for several values of k.
 * 
 * @author vbonnici
 *
 */
public class HRCoverage {

	public static void usage(){
		System.out.println("Usage:  cmd  isequence.3bit isequence.nelsa <<k>>");
	}
	
	public static void main(String[] args){
			
	//		DecimalFormat df = new DecimalFormat("#.##");
			
			String a_iseq = null;
			String a_inelsa = null;
			
			Set<Integer> kk = new TreeSet<Integer>();
			
			try{
				a_iseq = args[0];
				a_inelsa = args[1];
				
				if(args.length == 2){
					usage();
					System.exit(1);
				}
				
				for(int i=2; i<args.length; i++){
					kk.add(Integer.parseInt(args[i]));
				}
				
			}catch(Exception e){
				usage();
				System.exit(1);
			}
			
			
			try{
				
				Timer timer = new Timer();
				
				System.out.println(a_iseq);
				System.out.println("Loading sequence...");
				B3LLSequence b3seq = B3LLSequence.load(a_iseq);
				System.out.println("done "+timer.getElapsedSecs() +"sec.\n");
				
				int s_length = b3seq.length();
				System.out.println("Sequence length: " + s_length);
				
				timer.reset();
				System.out.println("Loading NELSA...");
				NELSA nelsa = new NELSA();
				nelsa.load(a_inelsa);
				System.out.println("done "+timer.getElapsedSecs()+" sec.");
				nelsa.setSequence(b3seq);
				
				
//				double elength = b3seq.length() - b3seq.countBads();
				int nof_r, nof_h;
				double[] r_cov = new double[b3seq.length()];
				double[] h_cov = new double[b3seq.length()];
				double[] cov;
				int[] sa = nelsa.sa();
				
				System.out.println("k\tnof_r\tnof_h\tr_scov\tH_scov\tR_pavg\tR_psd\tH_pavg\tH_psd");
				
				
				for(Integer k  : kk){
					Arrays.fill(r_cov, 0.0);
					Arrays.fill(h_cov, 0.0);
					nof_r = nof_h = 0;
					
					
					IELSAIterator it = nelsa.begin(k);
					while(it.next()){
						if(it.multiplicity() == 1){
							cov = h_cov;
							nof_h++;
						}
						else{
							cov = r_cov;
							nof_r++;
						}
						
						for(int i=it.istart(); i<it.iend(); i++){
							for(int j=sa[i]; j<sa[i]+k; j++){
								cov[j]++;
							}
						}
					}
					
					
					System.out.print(k+"\t"+nof_r+"\t"+nof_h+"\t");
					System.out.print(sequenceCoverage(r_cov)+"\t");
					System.out.print(sequenceCoverage(h_cov)+"\t");
					double avg = AvgSD.avg(r_cov, 0.0);
					System.out.print(avg+"\t");
					System.out.print(AvgSD.sd(r_cov, avg, 0.0)+"\t");
					avg = AvgSD.avg(h_cov, 0.0);
					System.out.print(avg+"\t");
					System.out.print(AvgSD.sd(h_cov, avg, 0.0)+"\n");
				}
				
				
				
			}catch(Exception e){
				e.printStackTrace();
				System.out.println(e);
			}
	}
	
	private static int sequenceCoverage(double[] cov){
		int count = 0;
		for(int i=0; i<cov.length; i++){
			if(cov[i] != 0){
				count++;
			}
		}
		return count;
	}
	
//	private double positionalCoverageAVG(double[] cov){
//		double count = 0.0;
//		double avg = 0.0;
//		for(int i=0; i<cov.length; i++){
//			if(cov[i] != 0.0){
//				count++;
//				avg += cov[i];
//			}
//		}
//		return avg;
//	}
//	private double positionalCoverageSD(double[] cov, double avg){
//		double count = 0.0;
//		double sd = 0.0;
//		for(int i=0; i<cov.length; i++){
//			sd += 
//		}
//	}
}
