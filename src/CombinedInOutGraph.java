/*
 * CombinedInOutGraph.java
 * This program builds graph data structure which combine ingoing and outgoing edges 
 * in the adjacency lists, as used by Parimalarangan. The connection type is encoded 
 * using 2 bits out of the 32 bits as either 01 (outgoing), 10 (ingoing), or 11 (both).
 *
 * Dependency: 
 *       - WebGraph library.
 *       - java.io
 * Input: A directed graph and its transpose in webgraph format.
 * Usage: java CombinedInOutGraph basename
 *      Note: both the graph and the transpose graph must be present.
 * Output: Combined adjacency list in plain text format.
 *            basename-CDt.txt
 *
 * Algorithm: Use two pointers, one each for G and Gt.
 *
 * Version 1.00 
 *      - 4 Nove 2018 - Yudi Santoso 
 */
import it.unimi.dsi.webgraph.ImmutableGraph;
//import java.util.stream.IntStream;
import java.io.File;
import java.io.IOException;
import java.io.FileWriter;
import java.io.Writer;
//import java.util.Scanner;


public class CombinedInOutGraph {

    String basename;
    ImmutableGraph G;   // a directed graph
    ImmutableGraph Gt;  // transpose of G
    int n;
    int nG;
    int nGt;
    long mG;
    int maxdegG;
    long mGt;
    int maxdegGt;
	
    public CombinedInOutGraph(String basename) throws Exception {
        this.basename = basename;
        G = ImmutableGraph.loadMapped(basename);
        Gt = ImmutableGraph.loadMapped(basename+"-t");
        nG = G.numNodes();
        nGt = Gt.numNodes();
        if (nG != nGt) System.out.println("WARNING!! Number of nodes do not match:" + nG + " != " + nGt); 
        n = nG;
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


    public void createCompactDt() throws Exception {

        String outName;
        outName = basename+"-CDt.txt";
        File outTxtFile = new File(outName);
        FileWriter writer = new FileWriter(outTxtFile);
 
        writer.write(n + "\n");
        int edge_cnt = 0;       

        for(int u=0; u<n; u++) {
            if(u%1_000_000 == 0)
                System.out.println(u);
        	
            int[] u_neighbors = G.successorArray(u);
            int[] ut_neighbors = Gt.successorArray(u);
            int u_deg = G.outdegree(u);
            int ut_deg = Gt.outdegree(u);

            writer.write(u + "\t");

            for(int i=0, j=0; i<u_deg || j<ut_deg; ) {
                int v=0;
                int e1 = 0;   // edge uv
                if (u_deg != 0 && ut_deg != 0){
                    if (i == u_deg){       // i done (cannot be both done)
                        v = ut_neighbors[j];
                        e1 = 2;
                        v = (v << 2) | e1;
                        writer.write(v + "\t");
                        j++;
                    }
                    else if (j == ut_deg){  // j done (cannot be both done)
                        v = u_neighbors[i];
                        e1 = 1;
                        v = (v << 2) | e1;
                        writer.write(v + "\t");
                        i++;
                    }
                    else if (u_neighbors[i] < ut_neighbors[j]){
                        v = u_neighbors[i];
                        e1 = 1;
                        v = (v << 2) | e1;
                        writer.write(v + "\t");
                        i++;
                    }
                    else if (u_neighbors[i] > ut_neighbors[j]){
                        v = ut_neighbors[j];
                        e1 = 2;
                        v = (v << 2) | e1;
                        writer.write(v + "\t");
                        j++;
                    }
                    else if (u_neighbors[i] == ut_neighbors[j]){
                        v = u_neighbors[i];
                        e1 = 3;
                        v = (v << 2) | e1;
                        writer.write(v + "\t");
                        i++; j++;
                    }
                } else if (u_deg != 0 ){   // only u_neighbors
                    v = u_neighbors[i];
                    e1 = 1;
                    v = (v << 2) | e1;
                    writer.write(v + "\t");
                    i++;
                } else if (ut_deg != 0 ){  // only ut_neighbors
                    v = ut_neighbors[j];
                    e1 = 2;
                    v = (v << 2) | e1;
                    writer.write(v + "\t");
                    j++;
                }

                edge_cnt++;

            }            

            writer.write("\n");

        }

        System.out.println("Number of edges: " + edge_cnt);	
        writer.write(edge_cnt + "\n");   // to check that all edges is counted.

        writer.close();
    }


    public static void main(String[] args) throws Exception {

        long startTime = System.currentTimeMillis();
        String basename = args[0]; 
		
        CombinedInOutGraph t = new CombinedInOutGraph(basename);
        t.createCompactDt();
		
        System.out.println("Total time elapsed = " + (System.currentTimeMillis() - startTime) / 1000.0 + " seconds");
    }

}
