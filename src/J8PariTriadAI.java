/*
 * J8PariTriadAI.java
 * Counting/Enumerating Connected Triads in a Directed Graph
 *      using Intersection algorithm
 *      with Java8 parallel stream
 *      and input Compact Adjacency List in plain text.
 * Dependency:
 *       - java.util.stream.IntStream
 * Input: A directed graph with compact adjacency list in plain text format.
 * Usage: java J8PariTriadAI basename
 * Output: Count of connected triads for each type.
 *       Can be modified to enumerate the triads.
 * Algorithm: Use Edge Iteration
 * Version 1.0  -- only total count
 *      - Nov 7, 2018 - Yudi Santoso
 * Version 1.1  -- counts of each type
 *      - Nov 8, 2018 - Yudi Santoso
 * Version 1.2  -- use array for each neighbours
 *      - Nov 8, 2018 - Yudi Santoso
 * Version 1.3  -- use ArrayList<ArrayList>
 *      - Nov 9, 2018 - Yudi Santoso
 * Version 1.4  -- use trove for ArrayList
 *      - Nov 10, 2018 - Yudi Santoso
 * Version 1.5  -- clean up and bypass TriType to 7types
 *      - Nov 10, 2018 - Yudi Santoso
 */

import java.util.stream.IntStream;
import java.util.ArrayList;
import java.util.List;
import java.io.*;
import gnu.trove.list.array.*;

public class J8PariTriadAI {

   String filename;
   int n;
   long m;
   int maxdeg;
   int[] degs;
   int[] intersection;
   ArrayList<TIntArrayList> adjList;
   int[] Tri7Type;

   public J8PariTriadAI(String basename) throws Exception {

      Tri7Type = new int[64];
      Tri7Type[0] = 0;
      Tri7Type[1] = 0;
      Tri7Type[2] = 0;
      Tri7Type[3] = 0;
      Tri7Type[4] = 0;
      Tri7Type[5] = 0;
      Tri7Type[6] = 0;
      Tri7Type[7] = 0;
      Tri7Type[8] = 0;
      Tri7Type[9] = 0;
      Tri7Type[10] = 0;
      Tri7Type[11] = 0;
      Tri7Type[12] = 0;
      Tri7Type[13] = 0;
      Tri7Type[14] = 0;
      Tri7Type[15] = 0;
      Tri7Type[16] = 0;
      Tri7Type[17] = 0;
      Tri7Type[18] = 0;
      Tri7Type[19] = 0;
      Tri7Type[20] = 0;
      Tri7Type[21] = 2;
      Tri7Type[22] = 2;
      Tri7Type[23] = 4;
      Tri7Type[24] = 0;
      Tri7Type[25] = 1;
      Tri7Type[26] = 2;
      Tri7Type[27] = 3;
      Tri7Type[28] = 0;
      Tri7Type[29] = 3;
      Tri7Type[30] = 5;
      Tri7Type[31] = 6;
      Tri7Type[32] = 0;
      Tri7Type[33] = 0;
      Tri7Type[34] = 0;
      Tri7Type[35] = 0;
      Tri7Type[36] = 0;
      Tri7Type[37] = 2;
      Tri7Type[38] = 1;
      Tri7Type[39] = 3;
      Tri7Type[40] = 0;
      Tri7Type[41] = 2;
      Tri7Type[42] = 2;
      Tri7Type[43] = 5;
      Tri7Type[44] = 0;
      Tri7Type[45] = 4;
      Tri7Type[46] = 3;
      Tri7Type[47] = 6;
      Tri7Type[48] = 0;
      Tri7Type[49] = 0;
      Tri7Type[50] = 0;
      Tri7Type[51] = 0;
      Tri7Type[52] = 0;
      Tri7Type[53] = 5;
      Tri7Type[54] = 3;
      Tri7Type[55] = 6;
      Tri7Type[56] = 0;
      Tri7Type[57] = 3;
      Tri7Type[58] = 4;
      Tri7Type[59] = 6;
      Tri7Type[60] = 0;
      Tri7Type[61] = 6;
      Tri7Type[62] = 6;
      Tri7Type[63] = 7;

      this.filename = basename+"-CDt.txt";

      adjList = new ArrayList<TIntArrayList>();
      BufferedReader br = new BufferedReader(new FileReader(filename));

// Read in the adjacency list. Should put these in a try:
      String line;
      line = br.readLine();   // first line contains n only
      n = Integer.parseInt(line);
      for (int i=0; i<n; i++){
          adjList.add(new TIntArrayList());
      }
      degs = new int[n];
      maxdeg = 0;
      int mx = 0;    // check count of edges
      for (int i=0; i<n; i++){
          if(i%1_000_000 == 0) System.out.println(i);   // to check progress of loading
          line = br.readLine();
          String[] items = line.split("\t");
          int v = Integer.parseInt(items[0]);    // the node
          if (v != i) System.out.println("WARNING: missing node!");
          int deg = items.length - 1;   // the degree
          degs[i] = deg;
          for (int j=1; j<=deg; j++){
              adjList.get(v).add(Integer.parseInt(items[j]));    // populate the adjacency list
          }
          mx += deg;
          if(deg > maxdeg)
            maxdeg = deg;
      }

      line = br.readLine();   // last line contains m only
      m = Integer.parseInt(line);
      if (mx != m) System.out.println("WARNING: number of edges is not right!");

      br.close();        // Always close files.

      intersection = new int[maxdeg];
      System.out.println("n=" + n + ", m=" + m + ", maxdeg=" + maxdeg);
   }


   public void computeTriangles() throws Exception {

      long[] t_count = IntStream.range(0,n).parallel().mapToObj(u -> {
          if(u%1_000_000 == 0) System.out.println(u);

          int u_deg = degs[u];
          int[] u_N = adjList.get(u).toArray();     //  u_neighbors


          long[] tp_count = new long[8];
          for(int i=0; i<u_deg; i++) {
             int v = u_N[i];     //  u_neighbors[i];
             int uv = 3 & v;     // get the 2-bits code
             v = v >>> 2;         // shift 2-bits to the right to get the node label
             if (u < v){       // to avoid double counting and self-loop
                int v_deg = degs[v];
                int[] v_N = adjList.get(v).toArray();     //  v_neighbors
                int[] ucount = intersection(u,v,u_deg,v_deg, u_N, v_N, uv);
                   for (int t=1; t<8; t++) {
                       tp_count[t] += ucount[t];
                   }
             }
          }
          return tp_count;
      }).reduce(new long[]{0,0,0,0,0,0,0,0}, (a,b)->sumArrays(a,b,8));

      long t_total_cnt = 0;
      for (int t=1; t<8; t++){
          System.out.println("Number of triangles type " + t + ": " + t_count[t]);
          t_total_cnt += t_count[t];
      }
      System.out.println("Total number of triad triangles: " + t_total_cnt);
      System.out.println("Number of trust triangles: " + (t_count[2] + t_count[3] + 2*t_count[4]+ 2*t_count[5]+ 3*t_count[6]+ 6*t_count[7]));
      System.out.println("Number of cycle triangles: " + (t_count[1] + t_count[3] + t_count[6]+ 2*t_count[7]));

   }


   int[] intersection(int u, int v,  int u_deg, int v_deg, int[] u_N, int[] v_N, int uv) {

      int[] ucount = new int[8];  // triangle counts for each types
      for (int t=0; t<8; t++) ucount[t] = 0;
      int tricode = 0;
      int stype = 0;    // seventype

      for(int i=0,j=0; i<u_deg && j<v_deg; ) {

         int uni = u_N[i];
         int uw = 3 & uni;
         uni = uni >>> 2;

         int vnj = v_N[j];
         int vw = 3 & vnj;
         vnj = vnj >>> 2;

         if(uni == vnj) {
         //  the if condition below is to avoid double counting and self-loop:
             if(uni > v){
// triad types:
                tricode = TriCode(uv, vw, uw);
                stype = Tri7Type[tricode];   // 7 types
                ucount[stype]++;
// For enumeration:
//                System.out.println("u,v,w, tricode, tritype, stype: " + u + "\t" + v + "\t" + uni + "\t" + tricode + "\t" + tritype + "\t" + stype);
             }
             i++; j++;
             continue;
         }

         if(uni < vnj) {
            i++;
            continue;
         }

         if(uni > vnj) {
            j++;
            continue;
         }
      }

//      return index;
      return ucount;
   }




    int TriCode(int uv, int vw, int uw){
        // Note: we have uvw here, compared to Pajek and Pari vuw. (i.e., u<->v).
        int Luv = uv & 1;
        int Lvu = (uv & 2) >>> 1;
        int Lvw = vw & 1;
        int Lwv = (vw & 2) >>> 1;
        int Luw = uw & 1;
        int Lwu = (uw & 2) >>> 1;
        int tricode = Luv + 2*(Lvu + 2*(Luw + 2*(Lwu + 2*(Lvw + 2*Lwv))));
        return tricode;
    }

    static long[] sumArrays(long[] A, long[] B, int N){
        long[] C = new long[N];
        for (int i=0;i<N;i++){
            C[i] = A[i] + B[i];
        }
        return C;
    }


   public static void main(String[] args) throws Exception {

      long startTime = System.currentTimeMillis();

      String basename = args[0];

      J8PariTriadAI t = new J8PariTriadAI(basename);

      System.out.println("Loading time = " + (System.currentTimeMillis() - startTime) / 1000.0 + " seconds");

      t.computeTriangles();

      System.out.println("Total time elapsed = " + (System.currentTimeMillis() - startTime) / 1000.0 + " seconds");
   }

}
