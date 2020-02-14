/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package needleman_wunsch;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

/**
 *
 * @author Hilal KAYNARCA
 */
public class NeedlemanWunsch {
    // scoring variables
    public double match, mismatch, gap;

    /**
     * This creates a new NeedlemanWunsch object which will perform string
     * comparisons using the provided scoring system. The three values represent
     * the cost of performing each operation.
     * 
     * @param match
     *            is the (typically positive) change in score when the ith
     *            character of the two strings match
     * @param mismatch
     *            is the (typically negative) change in score when the ith
     *            character of the two strings do not match but are forced to
     *            align
     * @param gap
     *            is the (typically negative) change in score when the ith
     *            character of the two strings do not match, and a "gap" is
     *            added
     */
    public NeedlemanWunsch(double match, double mismatch, double gap) {
        this.match = match;
        this.mismatch = mismatch;
        this.gap = gap;
    }

    /**
     * This computes the optimal global alignment using the Needleman-Wunsch
     * algorithm (dynamic programming). It uses the scoring system determined
     * when this object was first initialized. See in-line comments for more
     * details about the operation of this algorithm.
     * 
     * This implementation only finds the edit-distance. To backtrack and
     * recover the optimal alignment, look at the bottom-right corner of the
     * `dp` table and work back up.
     * 
     * @param a
     *            is the first string
     * @param b
     *            is the second string
     * @return the edit distance
     */
    public double align(String a, String b) {
        // create a dynamic programming table such that:
        // - dp[0][0] = align(a.substring(0,0), b.substring(0,0))
        // - dp[1][0] = align(a.substring(0,1), b.substring(0,0))
        // - dp[1][1] = align(a.substring(0,1), b.substring(0,1))
        // ...
        // - dp[a.length()+1][b.length()+1] = align(a, b)
        double[][] dp = new double [a.length() + 1][b.length() + 1];

        // initialize dp table with empty string alignment values
        // - align("","") = 0;
        dp[0][0] = 0;
        // - align("","...") = "...".length() * gap
        for (int i = 1; i < dp.length; i++)
            dp[i][0] = dp[i - 1][0] + gap;
        // - align("...","") = "...".length() * gap
        for (int j = 1; j < dp[0].length; j++)
            dp[0][j] = dp[0][j - 1] + gap;

        // compute all the alignments
        for (int i = 1; i < dp.length; i++) {
            for (int j = 1; j < dp[i].length; j++) {
                // initialize best_score
                double best_score = Double.MIN_VALUE;
                // perfect match?
                if (a.charAt(i - 1) == b.charAt(j - 1))
                    if (best_score < dp[i - 1][j - 1] + match)
                        best_score = dp[i - 1][j - 1] + match;
                // force match?
                if (a.charAt(i - 1) != b.charAt(j - 1))
                    if (best_score < dp[i - 1][j - 1] + mismatch)
                        best_score = dp[i - 1][j - 1] + mismatch;
                // insert gap in `a`?
                if (best_score < dp[i - 1][j] + gap)
                    best_score = dp[i - 1][j] + gap;
                // insert gap in `b`?
                if (best_score < dp[i][j - 1] + gap)
                    best_score = dp[i][j - 1] + gap;
                // assign best_score
                dp[i][j] = best_score;
                //System.out.print("Skor:"+dp[i][j] +"i:"+i+"j:"+j);
            }
        }


        // print out dp table
        for (int i = 0; i < dp.length; i++)
            System.out.println(java.util.Arrays.toString(dp[i]));

        // return the bottom-right corner of the dp table
        return dp[dp.length - 1][dp[0].length - 1];
    }

    /**
     * Test the Needleman-Wunsch algorithm using two DNA fragments.
     */
    public static void main(String[] args) throws FileNotFoundException {
        
        ArrayList<Double> skor= new ArrayList<Double>();
        
        String seq_1,seq_2;
            boolean first = true;
String[] sekans = new String[100000];
        try (Scanner sc = new Scanner(new File("1K_Sequence.fasta"))) {
            
            int i=0;
            while (sc.hasNextLine()) {
                String line =sc.nextLine().trim(); 
                               
                if (line.length()>0&&line.charAt(0) == '>') {                
                       if (first)                      
                        first = false; 
                      // System.out.println(line);
                } 
                else if(line.length()>0) {                 
                   sekans[i]=line;
                   System.out.println(sekans[i]);
                   i++;
                   
               
                                                                 
                }   
                
            }
        }
        for(int k=0;k<1000;k++){
            seq_1=sekans[k];
            for(int l=0;l<1000;l++){
            seq_2=sekans[l];
            NeedlemanWunsch nw = new NeedlemanWunsch(3.621354295, -2.451795405,-1.832482334);
                 double score = nw.align(seq_1,seq_2);
                 //System.out.println("Expected: 0, Actual: " + score);
                 skor.add(score);
                // System.out.println("k="+k+"l="+l+"skor= " + score);
                 
                 
            }}
        
        for(int i=0;i<1000;i++){
        for(int j=0;j<1000;j++){
        System.out.println("Skor:"+skor.get(j) +"i:"+i+"j:"+j);
        }
        }
        
       for(int a=0;a<1000000;a++){
          // System.out.println(skor.get(a));
         Collections.sort(skor,Collections.reverseOrder());
       }
       
      //System.out.println(skor);
       
     for(int i=0; i<20;i++ )
        {
            System.out.println(skor.get(i));
           
        }
 
       
       // double score = nw.align(,);
        
    }
}