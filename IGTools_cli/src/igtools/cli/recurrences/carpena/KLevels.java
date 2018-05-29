package igtools.cli.recurrences.carpena;

import igtools.analyses.recurrences.carpena.Carpena2009;
import igtools.analyses.recurrences.carpena.Carpena2009Factory;
import igtools.analyses.recurrences.distances.ProperMinimalRecurrenceDistancesExtractor;
import igtools.cli.util.TwoColumnsConfigFile;
import igtools.common.charts.KHeatMap;
import igtools.common.distributions.DistributionUtils;
import igtools.common.distributions.distance.DistributionDistance;
import igtools.common.distributions.distance.KLDistance;
import igtools.common.distributions.distance.KSDistance;
import igtools.common.distributions.distance.PointDiffSUMDistance;
import igtools.common.sequence.B3LLSequence;
import igtools.dictionaries.elsa.CompleteUnitIterator;
import igtools.dictionaries.elsa.IELSAIterator;
import igtools.dictionaries.elsa.NELSA;

import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import javax.imageio.ImageIO;

public class KLevels {
	public static void usage(){
		System.out.println("Usage: cmd b3seq nelsa <sigma|sigma_nor|C> <confgi_file> ofile.png");
		System.out.println("OUTPUT\t <S: kmer value>");
	}
	
	
	enum DistType {SIGMA, SIGMA_NOR, C};
	
	public static void main(String[] args){
		String a_b3seq = "";
		String a_nelsa = "";
		DistType a_dtype = DistType.SIGMA; 

		String a_cfile = "";
		String a_ofile = "";
		
		try{
			a_b3seq = args[0];
			a_nelsa = args[1];
			
			
			if(args[2].compareTo("sigma") == 0){
				a_dtype = DistType.SIGMA;
			}
			else if(args[2].compareTo("sigma_nor") == 0){
				a_dtype = DistType.SIGMA_NOR;
			} 
			else if(args[2].compareTo("C") == 0){
				a_dtype = DistType.C;
			}
			else
				throw new Exception();
			
			a_cfile = args[3];
			a_ofile = args[4];
			
		}catch(Exception e){
			usage();
			System.exit(1);
		}
		
		try{
			
			System.out.println("Loading sequence...");
			B3LLSequence b3seq = B3LLSequence.load(a_b3seq);
			System.out.println("done");
			System.out.println("Loading NELSA...");
			NELSA nelsa = new NELSA();
			nelsa.load(a_nelsa);
			System.out.println("done");
			nelsa.setSequence(b3seq);
			
			
			/*====================================================================================================================*/
			//set hetamp value arrays
			
			int MAX_K = 6;
			
			double[][] values = new double[MAX_K][];
			for(int i=0; i<MAX_K; i++){
				values[i] = new double[(int)Math.pow(4, i+1)];
				for(int j=0; j<values[i].length; j++)
					values[i][j] = 0.0;
			}
			
			
			/*====================================================================================================================*/
			//putting values
			
			Carpena2009 c2009 = Carpena2009Factory.factory(nelsa, true, false);
			
			for(int k=1; k<=MAX_K; k++){
				System.out.println("k "+k);
				
				final CompleteUnitIterator cit = new CompleteUnitIterator(k);
                int cit_code;
                IELSAIterator it;
                double value;
                do{
                	value = 0.0;
                    cit_code = cit.code();
                    it = nelsa.find(new B3LLSequence(cit.kmer()));
                    if(it != null && it.multiplicity()>1){
                    	                   	
    					
    					switch(a_dtype){
						case C:
							value = c2009.C(it);
							break;
						case SIGMA:
							value = c2009.sigma(it);
							break;
						case SIGMA_NOR:
							value = c2009.sigma_nor(it);
							break;
						default:
							break;
    					
    					}
    						
                    }
                    
                    values[k -1][cit_code]  =value;
                    
                }while(cit.next());
				
				
			}

			
			
			
			/*====================================================================================================================*/
			//create and save hetmap
			
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
			
			
			/*====================================================================================================================*/
			
			
		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e);
		}
	}
}
