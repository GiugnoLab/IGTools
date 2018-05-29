package igtools.cli.charts.heatmap;

import igtools.cli.util.TwoColumnsConfigFile;
import igtools.common.charts.KHeatMap;
import igtools.common.kmer.b2.unit.B2UnitRLKmer;
import igtools.common.sequence.B3LLSequence;
import igtools.common.util.Timer;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import javax.imageio.ImageIO;


/**
 * Build and heatmap of k-mer multiplicities for each of the input sequence.
 * One row, one sequence, with kmers lexicographically ordered on the x axis.
 * 
 * @author vbonnici
 *
 */
public class MultsHeatMaps {
	public static void usage(){
		System.out.println("Usage: cmd <config_file> <<iseq.3bit>> <k> ofile.png");
		
	}
	
	public static void main(String[] args){
		String a_cfile = null;
		String[] a_iseqs = null;
		String a_ofile = null;
		int k = 0;
		
		try{
			a_cfile = args[0];
			System.out.println("a_cfile: "+a_cfile);
			a_iseqs = new String[args.length-3];
			
			int i=1;
			for(; i<args.length-2; i++){
				a_iseqs[i-1] = args[i];
				System.out.println("a_iseqs[i]: "+a_iseqs[i-1]);
			}
			
			k = Integer.parseInt(args[i]);
			System.out.println("k: "+k);
			if(k<1 || k>6)
				throw new Exception();
			a_ofile = args[i+1];
			System.out.println("a_ofile: "+a_ofile);
			
			
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		
		try{
			Map<String,String> hmconfig = new HashMap<String,String>();
			TwoColumnsConfigFile.parse(a_cfile, hmconfig);
			KHeatMap heatmap = new KHeatMap();
			heatmap.configByMap(hmconfig);
			
			double[][] values = new double[a_iseqs.length][(int)Math.pow(4, k)];
			for(int i=0; i<values.length; i++){
				for(int j=0; j<values[i].length; j++)
					values[i][j] = 0;
			}
			
			int omask = B2UnitRLKmer.MAX_CODE[1];
			int amask = B2UnitRLKmer.MAX_CODE[k];
			int code = 0 & amask;
			int nn = k;
			int k_1 = k-1;
			
			
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
						if(nn<k)
							nn++;
						code = code << B2UnitRLKmer.BITS_PER_MER;
						code |= (c & omask);
						code &= amask;
						
						if(i>=k_1 && nn>k_1)
							values[s][code & B2UnitRLKmer.MAX_CODE[k]]++;
					}
				}
			}
			
			String[] labels = new String[a_iseqs.length];
			String name="";
			int ext_i;
			File f;
			for(int i=0; i<labels.length; i++){
				f = new File(a_iseqs[i]);
				name = f.getName();
				ext_i = name.lastIndexOf('.');
				if(ext_i >0)
					name = name.substring(0, ext_i);
				labels[i] = name;
			}
			
			
			heatmap.values = values;
			heatmap.labels = labels;
			heatmap.draw();
			
			ImageIO.write(heatmap.img, "png", new File(a_ofile));
			
		
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
			
		}
	}
}
