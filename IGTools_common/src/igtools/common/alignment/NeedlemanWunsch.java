/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package igtools.common.alignment;

import igtools.common.nucleotide.B3Nucleotide;
import igtools.common.sequence.B3LLSequence;
import igtools.common.sequence.B3Sequence;
import igtools.common.sequence.ReverseB3Sequence;

/**
 *
 * @author vbonnici
 */
public class NeedlemanWunsch {
    
    private final static int gapD = -1;
    
    private final static int[][] simMatrix = 
    {  
        {1,-1,-1,-1,-1},
        {-1,1,-1,-1,-1},
        {-1,-1,1,-1,-1},
        {-1,-1,-1,1,-1},
        {-1,-1,-1,-1,1},
    };
    
    
    private B3Sequence s1;
    private B3Sequence s2;
    
    private int[][] mF;
    
    private B3Sequence as1;
    private B3Sequence as2;
    
    public NeedlemanWunsch(B3Sequence s1, B3Sequence s2){
        this.s1 = s1;
        this.s2 = s2;
        
        fillMf();
        getAligedSequences();
    }
    
    
    public int[][] getScoreMatrix(){
        return mF;
    }
    public B3Sequence getAlignedS1(){
        return as1;
    }
    public B3Sequence getAlignedS2(){
        return as2;
    }
    
    
    private void fillMf(){
        mF = new int[s1.length() +1][s2.length() +1];
        
        for(int i=0; i<mF.length; i++)
            mF[i][0] = gapD * i;
        for(int j=0; j<mF[0].length; j++)
            mF[0][j] = gapD * j;
        
        for(int i=0; i<mF.length; i++){
            for(int j=0; j<mF[i].length; j++){
                if (i == 0) {
                        mF[i][j] = -j;
                } else if (j == 0) {
                        mF[i][j] = -i;
                } else {
                        mF[i][j] = 0;
                }
            }
        }
        
        
        int match, delete, insert;
        
        for(int i=1; i<mF.length; i++){
            for(int j=1; j<mF[i].length; j++){
                
                match = mF[i-1][j-1] + simMatrix[s1.getB3(i-1)] [s2.getB3(j-1)];
                delete = mF[i-1][j] + gapD;
                insert = mF[i][j-1] + gapD;
                
                mF[i][j] = Math.max(match, Math.max(delete, insert));
                
            }
        }
    }
    
    private void getAligedSequences(){
        as1 = new B3LLSequence(s1.length());//presumed length
        as2 = new B3LLSequence(s2.length());//presumed length
        
        
        int i = s1.length();
        int j = s2.length();
        
        int a_i = Math.max(s1.length() - 1, s2.length() -1);
        
        while(i>0 && j>0){
            
            //System.out.println("("+i+","+j+")["+mF[i][j]+"]["+(mF[i-1][j-1]+simMatrix[s1.getB3(i)][s2.getB3(j)])+"]");
            
            
            if(mF[i][j]==mF[i-1][j-1]+simMatrix[s1.getB3(i-1)][s2.getB3(j-1)]){
                as1.setB3(a_i, s1.getB3(i-1));
                as2.setB3(a_i, s2.getB3(j-1));
                i--;
                j--;
            }
            else if(mF[i][j]==mF[i-1][j]+gapD){
                as1.setB3(a_i, s1.getB3(i-1));
                as2.setB3(a_i, B3Nucleotide.N_CODE);//TODO
                i--;
            }
            else{
                as1.setB3(a_i, B3Nucleotide.N_CODE);//TODO
                as2.setB3(a_i, s2.getB3(j-1));
                j--;
            }
            a_i--;
        }
        
        if(i>0){
            while(i>0){
                as1.setB3(a_i, s1.getB3(i-1));
                as2.setB3(a_i, B3Nucleotide.N_CODE);//TODO
                i--;
                a_i--;
            }
                
        }
        if(j>0){
            while(j>0){
                as1.setB3(a_i, B3Nucleotide.N_CODE);//TODO
                as2.setB3(a_i, s2.getB3(j-1));
                j--;
                a_i--;
            }
        }
        
    }
    
    
    
    
    
    public static void main(String[] args){
        B3LLSequence s1 = new B3LLSequence("TAATC");
        B3LLSequence s2 = new B3LLSequence("TACGAGTATGGAGC");
        //NeedlemanWunsch alignment = new NeedlemanWunsch(s1, s2);
        NeedlemanWunsch alignment = new NeedlemanWunsch(new ReverseB3Sequence(s1), new ReverseB3Sequence(s2));
        
        System.out.println(new ReverseB3Sequence(alignment.getAlignedS1()));
        System.out.println(new ReverseB3Sequence(alignment.getAlignedS2()));
        
        int[][] mat = alignment.getScoreMatrix();
        for(int i=0; i<mat.length; i++){
            for(int j=0; j<mat[i].length; j++){
                System.out.print(mat[i][j]+"\t");
            }
            System.out.println();
        }
    }
    
    
}
