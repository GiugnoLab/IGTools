package igtools.cli.dictionaries;

import java.text.DecimalFormat;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.streams.BufferedFASTANucleotideStream;
import igtools.common.streams.INucleotideStream;

public class Frequencies {

	
	public static void main(String[] args){
		String iseq=null;
		
		try{
			iseq = args[0];
		}catch(Exception e){
			System.exit(1);
		}
		
		
		try{
			DecimalFormat df = new DecimalFormat("#.##########");
			INucleotideStream<B3Nucleotide> fastaStream = new BufferedFASTANucleotideStream(iseq);
			
			int[] k1mers = new int[4];
			int[][] k2mers = new int[4][4];
			int[][][] k3mers = new int[4][4][4];
			
			int code = 0;
			int p_code = -1;
			int pp_code = -1;
			
			while((code = fastaStream.nextCode()) != INucleotideStream.NULL_CODE){
				if(code < 4){
					k1mers[code]++;
					
					if(p_code > -1 && p_code < 4){
						k2mers[p_code][code]++;
						
						if(pp_code > -1 && pp_code < 4){
							k3mers[pp_code][p_code][code]++;
						}
					}
				}
				
				pp_code = p_code;
				p_code = code;
			}
			
//		for(int i=0; i<k3mers.length; i++){
//			for(int j=0; j<k3mers[i].length; j++){
//				for(int k=0; k<k3mers[i][j].length; k++){
//					System.out.print(k3mers[i][j][k] + "\t");
//				}
//			}
//		}
//		System.out.println();
			
			double[] k1ratio = new double[4];
			double[][] k2ratio = new double[4][4];
			double[][][] k3ex = new double[4][4][4];
			double[][][] k3ratio= new double[4][4][4];
			
			double sum = 0;
			
			
			for(int i=0; i<k1mers.length; i++)
				sum += k1mers[i];
			for(int i=0; i<k1mers.length; i++){
				k1ratio[i] = k1mers[i] / sum;
			}
			
			sum = 0;
			for(int i=0; i<k2mers.length; i++){
				for(int j=0; j<k2mers[i].length; j++){
					sum += k2mers[i][j];
				}
			}
			for(int i=0; i<k2mers.length; i++){
				for(int j=0; j<k2mers[i].length; j++){
					k2ratio[i][j] = k2mers[i][j] / sum;
				}
			}
			
			
			sum = 0;
			for(int i=0; i<k3mers.length; i++){
				for(int j=0; j<k3mers[i].length; j++){
					for(int k=0; k<k3mers[i][j].length; k++){
						sum += k3mers[i][j][k];
					}
				}
			}
			for(int i=0; i<k3mers.length; i++){
				for(int j=0; j<k3mers[i].length; j++){
					for(int k=0; k<k3mers[i][j].length; k++){
						k3ratio[i][j][k] = k3mers[i][j][k] / sum;
					}
				}
			}
			for(int i=0; i<k3mers.length; i++){
				for(int j=0; j<k3mers[i].length; j++){
					for(int k=0; k<k3mers[i][j].length; k++){
						//k3ex[i][j][k] = k2ratio[i][j] * k2ratio[j][k];
						k3ex[i][j][k] = (k2ratio[i][j] * k2ratio[j][k]) / k1ratio[j];
					}
				}
			}
			
			
			
			for(int i=0; i<k3mers.length; i++){
				for(int j=0; j<k3mers[i].length; j++){
					for(int k=0; k<k3mers[i][j].length; k++){
						System.out.print(k3mers[i][j][k] + "\t");
					}
				}
			}
			System.out.println();
			
			
			for(int i=0; i<k3mers.length; i++){
				for(int j=0; j<k3mers[i].length; j++){
					for(int k=0; k<k3mers[i][j].length; k++){
						System.out.print(k3ratio[i][j][k] + "\t");
					}
				}
			}
			System.out.println();
			
			sum = 0;
			for(int i=0; i<k3mers.length; i++){
				for(int j=0; j<k3mers[i].length; j++){
					for(int k=0; k<k3mers[i][j].length; k++){
						System.out.print(df.format(k3ex[i][j][k]) + "\t");
						sum += k3ex[i][j][k];
					}
				}
			}
			System.out.print(sum);
			System.out.println();
			
			
			
			for(int i=0; i<k3mers.length; i++){
				for(int j=0; j<k3mers[i].length; j++){
					for(int k=0; k<k3mers[i][j].length; k++){
						System.out.print((k3ratio[i][j][k] - k3ex[i][j][k]) + "\t");
					}
				}
			}
			System.out.println();
			
			
			
//			int goods = 0;
//			int bads = 0;
//			
//			code = fastaStream.nextCode();
//			if(code != INucleotideStream.NULL_CODE){
//				if(code < 4){
//					k1mers[code]++;
//					goods++;
//				}
//				else{
//					bads++;
//				}
//				p_code = code;
//			}
//			
//			
//			
//			while((code = fastaStream.nextCode()) != INucleotideStream.NULL_CODE){
//				if(code < 4){
//					goods++;
//					k1mers[code]++;
//					if(p_code < 4)
//						k2mers[p_code][code]++;
//				}
//				else{
//					bads++;
//				}
//				p_code = code;
//				
//			}
//			
//			System.out.println(goods+"\t"+bads);
//			
//			for(int i=0; i<k1mers.length; i++){
//				System.out.print(k1mers[i] + "\t");
//			}
//			System.out.println();
//			
//			for(int i=0; i<k2mers.length; i++){
//				for(int j=0; j<k2mers[i].length; j++){
//					System.out.print(k2mers[i][j] + "\t");
//				}
//			}
//			System.out.println();
//			
//			
//			double sum = 0;
//			
//			double[] k1mersp = new double[4];
//			
//			for(int i=0; i<k1mers.length; i++)
//				sum += k1mers[i];
//			for(int i=0; i<k1mers.length; i++){
//				k1mersp[i] = k1mers[i] / sum;
//				System.out.print(( k1mersp[i]  )+"\t");
//			}
//			System.out.println();
//			
//			for(int i=0; i<k2mers.length; i++){
//				for(int j=0; j<k2mers[i].length; j++){
//					System.out.print((k1mersp[i] * k1mersp[j]) + "\t");
//				}
//			}
//			System.out.println();
//			
//			sum = 0;
//			for(int i=0; i<k2mers.length; i++){
//				for(int j=0; j<k2mers[i].length; j++){
//					sum += k2mers[i][j];
//				}
//			}
//			
//			for(int i=0; i<k2mers.length; i++){
//				for(int j=0; j<k2mers[i].length; j++){
//					System.out.print((k2mers[i][j] / sum)+"\t");
//				}
//			}
//			System.out.println();
			
			
		}catch(Exception e){
			System.out.println(e);
			e.printStackTrace();
			System.exit(2);
		}
		System.exit(0);
		
		
	}
}
