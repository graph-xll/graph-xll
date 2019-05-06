/*
 * RedistDeg.java
 * Used sequentially after SortGraphAscBg/GGtAscBg. After Sort+Bg 
 *   the effective degrees are accumulated in the middle 
 *   (a hill shape distribution). This program redistributes 
 *   the degrees by controlled permutation of the labels.
 *   (It cannot be random shuffling because it has to be 
 *   invertible).  
 * G and Gt do not need to be done simultaneously as long as the parameters used are the same.
 *   However, for Gt, might need to rename the file after.
 *           - can use a t-flag for this.
 * Usage: 
 *       java RedistDeg basename  tflag  step
 *           where 
 *              basename is the WebGraph basename - after Sort+Bg/GGtAscBg.
 *              tflag = "t" for transpose (e.g., "g" for non-transpose).
 *              step is the stepsize (positive integer)
 * Output files: 
 *        basename-RD.graph
 *      or
 *        basename-RD-t.graph // for the transpose graph.
 * -
 * Requires: it.unimi.dsi.webgraph libraries.
 * -
 * Version 1.00 - first version 
 *       - 16 Nov 2018 - Yudi Santoso
 */ 
import java.io.IOException;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.Arrays;

import it.unimi.dsi.webgraph.BVGraph;
import it.unimi.dsi.webgraph.ImmutableGraph;
import it.unimi.dsi.webgraph.IncrementalImmutableSequentialGraph;

public class RedistDeg {

   String basename;
   ImmutableGraph G;
   String tflag;
   int step;
   int n;
//   int[] deg;
	
   public RedistDeg(String basename, String tflag, int step) throws Exception {
      this.basename = basename;
      this.tflag = tflag;
      this.step = step;

      String inname;
      if (tflag.equals("t")){
          inname = basename+"-t";
          System.out.println("Transpose");
      } else {
          inname = basename;
      }

      G = ImmutableGraph.loadMapped(inname);
      n = G.numNodes();
//      deg = new int[n];
   }

   public void RAndSave() throws Exception {

      String outname;
      if (tflag.equals("t")){
          outname = basename+"-RD-t";
      } else {
          outname = basename+"-RD";
      }
      		
      final IncrementalImmutableSequentialGraph g = new IncrementalImmutableSequentialGraph();
      ExecutorService executor = Executors.newSingleThreadExecutor();
      final Future<Void> future = executor.submit( new Callable<Void>() {
         public Void call() throws IOException {
                        BVGraph.store( g, outname );
                        return null;
         }
      } );

// newId translation:
      int[] Bx = new int[n];   // new to old
      int[] Cx = new int[n];   // old to new
      int offset = 0;   
//      step = 1000;    // step size - taken as an input parameter
      int j = 0;          // step multiplier
      for (int v=0; v<n; v++){
         int x = offset + j*step;
         if (x < n){
             Bx[v] = x;
             j++;
         } else {
             offset++;
             j = 0;
             Bx[v] = offset;
             j++;
         }
         Cx[Bx[v]] = v;
      }

      for(int v=0; v<n; v++) {
         if (v%1_000_000 == 0) System.out.println(v);
         int v_deg = G.outdegree(Bx[v]);
         int[] v_succ = G.successorArray(Bx[v]);  // translate this into new order
			
         int[] v_succ2 = new int[v_deg];    // use upper bound in size

         for(int i=0; i<v_deg; i++) {
	      v_succ2[i] = Cx[v_succ[i]];   // the translation
         }
// No need to trim

// Sort each neighbors list:
         Arrays.sort(v_succ2);     // webgraph compression requires this

         g.add(v_succ2, 0, v_deg);  
      }

      g.add( IncrementalImmutableSequentialGraph.END_OF_GRAPH );
      future.get();
      executor.shutdown();
   }
	
   public static void main(String[] args) throws Exception {
      long startTime = System.currentTimeMillis();
		
      String basename = args[0];
      String tflag = args[1];
      int step = Integer.parseInt(args[2]); 
		
      RedistDeg t = new RedistDeg(basename, tflag, step);

      t.RAndSave();
		
      System.out.println("Total time elapsed = " + (System.currentTimeMillis() - startTime) / 1000.0 + " seconds");
   }
}


