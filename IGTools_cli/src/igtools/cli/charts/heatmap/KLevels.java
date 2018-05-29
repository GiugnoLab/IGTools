package igtools.cli.charts.heatmap;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import javax.imageio.ImageIO;

import igtools.cli.util.TwoColumnsConfigFile;
import igtools.common.charts.KHeatMap;
import igtools.common.kmer.b2.unit.B2UnitRLKmer;
import igtools.common.sequence.B3LLSequence;
import igtools.common.util.Timer;


/**
 * Print an heatmap of kmer multiplicities (over the input sequences) for k=[1,6].
 * The multiplicity of a kmer is the sum of the multiplicity over all the sequences.
 * One row, one k, with kmers lexicographically ordered on the x axis.
 *  
 * @author vbonnici
 *
 */
public class KLevels {

	
	
	public static void usage(){
		System.out.println("Usage: cmd  <config_file> <<iseq.3bit>> ofile.png");
		
	}
	
	public static void main(String[] args){
		String a_cfile = null;
		String[] a_iseqs = null;
		String a_ofile = null;
		
		try{
			a_cfile = args[0];
			a_iseqs = new String[args.length-2];
			
			int i=1;
			for(; i<args.length-1; i++)
				a_iseqs[i-1] = args[i];
			
			a_ofile = args[i];
			
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		
		try{
		
			int MAX_K = 6;
			
			double[][] values = new double[MAX_K][];
//			System.out.println("values.length: "+values.length);
			for(int i=0; i<MAX_K; i++){
				values[i] = new double[(int)Math.pow(4, i+1)];
//				System.out.println("values["+i+"].length: "+values[i].length);
				for(int j=0; j<values[i].length; j++)
					values[i][j] = 0;
			}
			
			int omask = B2UnitRLKmer.MAX_CODE[1];
			int amask = B2UnitRLKmer.MAX_CODE[MAX_K];
			int code = 0 & amask;
			int nn = MAX_K+1;
			
//			System.out.println("omask: "+ Integer.toBinaryString(omask));
//			System.out.println("amask: "+ Integer.toBinaryString(amask));
//			System.out.println("code: "+ Integer.toBinaryString(code));
			
			
			Timer timer = new Timer();

			for(int s=0; s<a_iseqs.length; s++){
				System.out.println(a_iseqs[s]);
				System.out.println("Loading sequence...");
				B3LLSequence b3seq = B3LLSequence.load(a_iseqs[s]);
				System.out.println("done "+timer.getElapsedSecs() +"sec.\n");
				
				int c;
				for(int i=0; i<b3seq.length(); i++){
					c = b3seq.getB3(i);
					if(c>3){//N or NULL
						nn = 0;
					}
					else{
						if(nn<MAX_K+1)
							nn++;
						code = code << B2UnitRLKmer.BITS_PER_MER;
						code |= (c & omask);
						code &= amask;
						
						for(int j=0; j<MAX_K; j++){
							if(i>=j && nn>j){
								values[j][code & B2UnitRLKmer.MAX_CODE[j+1]]++;
							}
						}
					}
				}
			}
			
			for(int i=0; i<4; i++)
				System.out.println(values[0][i]);
			
			String[] labels = new String[MAX_K];
			for(int i=0; i<labels.length; i++)
				labels[i] = "k= "+(i+1);
			KHeatMap heatmap = new KHeatMap();
			Map<String,String> hmconfig = new HashMap<String,String>();
			TwoColumnsConfigFile.parse(a_cfile, hmconfig);
			heatmap.configByMap(hmconfig);
			heatmap.values = values;
			heatmap.labels = labels;
			heatmap.draw();
			
			if(heatmap.normalizeByRow){
				for(int i=0; i<values.length; i++){
					System.out.println("["+(i+1)+"]max_value: "+heatmap.max[i]);
				}
			}
			else{
				for(int i=0; i<values.length; i++){
					System.out.println("["+(i+1)+"]max_value: "+heatmap.max[0]);
				}
			}
			
			ImageIO.write(heatmap.img, "png", new File(a_ofile));
			
		
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
			
		}
	}
}
