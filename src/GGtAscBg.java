/*
 * GGtAscBg.java
 *      - derived from SortGraphAscBg.java
 * Takes a directed WebGraph and its transpose. 
 * Compare their max-out-degree.
 * Call the one with larger degree X. 
 * Sorts the nodes of both according to the out-degree of X 
 *   in ascending order. If two nodes have same degree, 
 *   sorts on the id. Relabels the nodes, of G and Gt 
 *   simultaneously, then, filters the result and saves 
 *   only the neighbors with higher (new) node id.
 * Usage: java GGtAscBg basename
 *         where basename is the WebGraph basename
 * Note: both G and Gt must be present.
 * Output files: 
 *        basename-GGt.graph
 *        basename-GGt-t.graph
 * -
 * Requires: net.mintern.primitive and it.unimi.dsi.webgraph
 *           libraries.
 * -
 * Version 1.00 - first version 
 *       - 6 Aug, 2018 - Yudi Santoso
 */ 
import java.io.IOException;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.Arrays;
import net.mintern.primitive.Primitive;

import it.unimi.dsi.webgraph.BVGraph;
import it.unimi.dsi.webgraph.ImmutableGraph;
import it.unimi.dsi.webgraph.IncrementalImmutableSequentialGraph;

public class GGtAscBg {

    String basename;
    ImmutableGraph G;
    ImmutableGraph Gt;  // transpose of G
    int n;
    int nG;
    int nGt;
    long mG;
    int maxdegG;
    long mGt;
    int maxdegGt;
    int[] deg;
	
    public GGtAscBg(String basename) throws Exception {
        this.basename = basename;
		
        G = ImmutableGraph.loadMapped(basename);
        Gt = ImmutableGraph.loadMapped(basename+"-t");
        nG = G.numNodes();
        nGt = Gt.numNodes();
        if (nG != nGt) System.out.println("WARNING!! Number of nodes do not match:" + nG + " != " + nGt); 
        n = nG;
        deg = new int[n];
        maxdegG = 0; mG = 0; maxdegGt = 0; mGt = 0;
        for(int v=0; v<n; v++) {
            int v_degG = G.outdegree(v);
            mG += v_degG;
            if(v_degG > maxdegG) maxdegG = v_degG;
            int v_degGt = Gt.outdegree(v);
            mGt += v_degGt;
            if(v_degGt > maxdegGt) maxdegGt = v_degGt;
        }
        System.out.println("n=" + n + ", mG=" + mG + ", maxdegG=" + maxdegG  + ", mGt=" + mGt + ", maxdegGt=" + maxdegGt);

    }

    public void SortAndSave() throws Exception {
	
        int[] idx = new int[n];   // new node labels
// Set the degrees to be sorted:
        for (int v=0; v<n; v++){
            if (maxdegG > maxdegGt){
                deg[v] = G.outdegree(v);
            } else {
                deg[v] = Gt.outdegree(v);
            }
            idx[v] = v;
        }

// Sort ascending:
        Primitive.sort(idx, (o1,o2) -> Integer.compare(deg[o1], deg[o2]));

        int[] vtx = new int[n];
        for(int i = 0; i < n; i++) vtx[idx[i]] = i;   // the new labels

// Create new WebGraph for G:
        final IncrementalImmutableSequentialGraph g = new IncrementalImmutableSequentialGraph();
        ExecutorService executor = Executors.newSingleThreadExecutor();
        final Future<Void> future = executor.submit( new Callable<Void>() {
            public Void call() throws IOException {
                        BVGraph.store( g, basename+"-GGt" );
                        return null;
            }
        } );
        
// Go through the adjacency list of G:
        for(int v=0; v<n; v++) {
            if (v%1_000_000 == 0) System.out.println(v);
            int v_deg = G.outdegree(idx[v]);
            int[] v_succ = G.successorArray(idx[v]);  // translate this into new order
            int[] v_succ2 = new int[v_deg];    // use upper bound in size

            int j=0;
            for(int i=0; i<v_deg; i++) {
// Filter out the smaller neighbours 
                if (v<vtx[v_succ[i]]) {
			        v_succ2[j]=vtx[v_succ[i]];
                    j++;
                }
            }
            int v_deg2B = j;
            int[] v_succ2B = new int[v_deg2B];    
            for(int i=0; i<v_deg2B; i++) {
                v_succ2B[i]=v_succ2[i];
            }
            Arrays.sort(v_succ2B);   // webgraph compression requires this

            g.add(v_succ2B, 0, v_deg2B);  
        }

        g.add( IncrementalImmutableSequentialGraph.END_OF_GRAPH );
        future.get();
        executor.shutdown();   // Done with G

// Create new WebGraph for Gt:
        final IncrementalImmutableSequentialGraph gt = new IncrementalImmutableSequentialGraph();
        ExecutorService executor2 = Executors.newSingleThreadExecutor();
        final Future<Void> future2 = executor2.submit( new Callable<Void>() {
            public Void call() throws IOException {
                        BVGraph.store( gt, basename+"-GGt-t" );
                        return null;
            }
        } );
        
// Go through the adjacency list of Gt:
        for(int v=0; v<n; v++) {
            if (v%1_000_000 == 0) System.out.println(v);
            int v_deg = Gt.outdegree(idx[v]);
            int[] v_succ = Gt.successorArray(idx[v]);  // translate this into new order
            int[] v_succ2 = new int[v_deg];    // use upper bound in size

            int j=0;
            for(int i=0; i<v_deg; i++) {
// Filter out the smaller neighbours 
                if (v<vtx[v_succ[i]]) {
			        v_succ2[j]=vtx[v_succ[i]];
                    j++;
                }
            }
            int v_deg2B = j;
            int[] v_succ2B = new int[v_deg2B];    
            for(int i=0; i<v_deg2B; i++) {
                v_succ2B[i]=v_succ2[i];
            }
            Arrays.sort(v_succ2B);   // webgraph compression requires this

            gt.add(v_succ2B, 0, v_deg2B);  
        }

        gt.add( IncrementalImmutableSequentialGraph.END_OF_GRAPH );
        future2.get();
        executor2.shutdown();   // Done with Gt

    }
	
    public static void main(String[] args) throws Exception {
        long startTime = System.currentTimeMillis();
		
        String basename = args[0]; 
		
        GGtAscBg t = new GGtAscBg(basename);

        t.SortAndSave();
		
        System.out.println("Total time elapsed = " + (System.currentTimeMillis() - startTime) / 1000.0 + " seconds");
    }
}


