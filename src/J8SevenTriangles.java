/*
 * J8SevenTriangles.java
 * Counting/Enumerating Triangles in a Directed Graph,
 *      classified into 7 types of triangles (i.e., 7 counters).
 *      Suppose triangle (u,v,w) : edges e1=uv, e2=vw, e3=wu (notice this!)       
 *      - Type 1:  ->, ->, ->    (1 cycle)
 *      - Type 2:  ->, ->, <-    (1 trust)
 *      - Type 3:  ->, ->, <->   (1 trust + 1 cycle)
 *      - Type 4:  ->, <-, <->   (2 trust merge)
 *      - Type 5:  <-, ->, <->   (2 trust parting)
 *      - Type 6:  <->, <->, ->  (4trust + 1 cycle)
 *      - Type 7:  <->, <->, <-> (6 trusts + 2 cycles)
 *      Edges types:
 *      .  .     0
 *      .->.     1
 *      .<-.     2
 *      .<->.    3
 *      There are 3*3*3 = 27 combinations of 3 non-zero edges.
 *      Triangle types in edge types combinations:
 *      - Type 1:  111, 212, 122, 222, 121, 211 
 *      - Type 2:  112, 221
 *      - Type 3:  113, 312, 132, 223, 321, 231
 *      - Type 4:  123, 311, 232 
 *      - Type 5:  213, 322, 131
 *      - Type 6:  133, 233, 313, 323, 331, 332
 *      - Type 7:  333
 * Use parallel stream (Java 8).
 * Dependency: 
 *       - WebGraph library.
 *       - java.util.stream.IntStream
 * Input: A directed graph and its transpose in webgraph format.
 * Usage: java J8SevenTriangles basename
 *      Note: both the graph and the transpose graph must be present.
 * Output: Seven counts of triangles, one for each type.
 *       Can be modified to enumerate the triangles.
 * Algorithm: Use Edge Iteration, on both G and Gt simultaneously.
 *            Use Four Pointers iteration to find common neighbors.
 *            Use flags for edge type/direction.
 *                 None  = 0
 *                 ->    = 1
 *                 <-    = 2
 *                 <->   = 3
 * Notes:
 *      - use long for mG and mGt 
 * Version 1.00 - derived from SevenTriangles 1.30
 *      - use Java 8 stream Parallel 
 *      - Aug 06, 2018 - Yudi Santoso 
 * Version 1.10 - add trust and cycle counts
 * Version 1.20 - clean up code
 *      - Nov 08, 2018 - Yudi Santoso 
 * Version 1.30 - use lookup table for triad type to avoid multiple 
 *                computation.
 *      - Nov 10, 2018 - Yudi Santoso 
 */
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.stream.IntStream;

public class J8SevenTriangles {

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
    int[][][] de7Type;
	
    public J8SevenTriangles(String basename) throws Exception {

        de7Type = new int[4][4][4];    // to avoid computing detype again and again
        for (int i=0; i<4; i++){
            for (int j=0; j<4; j++){
                for (int k=0; k<4; k++){
                    de7Type[i][j][k] = detype(i, j, k);
                }
            }
        }

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


    public void computeTriangles() throws Exception {

        long[] t_count = IntStream.range(0,n).parallel().mapToObj(u -> {
            if(u%1_000_000 == 0) System.out.println(u);
            ImmutableGraph H = G.copy();
            ImmutableGraph Ht = Gt.copy();
            int[] u_neighbors = H.successorArray(u);
            int[] ut_neighbors = Ht.successorArray(u);
            int u_deg = H.outdegree(u);
            int ut_deg = Ht.outdegree(u);
            long[] tp_count = new long[8];
            for(int i=0, j=0; i<u_deg || j<ut_deg; ) {
                int v=0;
                int e1 = 0;   // edge uv
                if (u_deg != 0 && ut_deg != 0){
                    if (i == u_deg){       // i done (cannot be both done)
                        v = ut_neighbors[j];
                        e1 = 2;
                        j++;
                    }
                    else if (j == ut_deg){  // j done (cannot be both done)
                        v = u_neighbors[i];
                        e1 = 1;
                        i++;
                    }
                    else if (u_neighbors[i] < ut_neighbors[j]){
                        v = u_neighbors[i];
                        e1 = 1;
                        i++;
                    }
                    else if (u_neighbors[i] > ut_neighbors[j]){
                        v = ut_neighbors[j];
                        e1 = 2;
                        j++;
                    }
                    else if (u_neighbors[i] == ut_neighbors[j]){
                        v = u_neighbors[i];
                        e1 = 3;
                        i++; j++;
                    }
                } else if (u_deg != 0 ){   // only u_neighbors
                    v = u_neighbors[i];
                    e1 = 1;
                    i++;
                } else if (ut_deg != 0 ){  // only ut_neighbors
                    v = ut_neighbors[j];
                    e1 = 2;
                    j++;
                }

                if (u < v){  // to avoid double counting
                             // also, when v=0 it means that e1=0, so do not iterate
                    int[] v_neighbors = H.successorArray(v);
                    int[] vt_neighbors = Ht.successorArray(v);
                    int v_deg = H.outdegree(v);
                    int vt_deg = Ht.outdegree(v);

                    int[] ucount = new int[8];   // temp counters

                    ucount = census4(u, v, e1, u_neighbors, ut_neighbors, v_neighbors, vt_neighbors, u_deg, ut_deg, v_deg, vt_deg);

                    for (int t=1; t<8; t++) {
                        tp_count[t] += ucount[t];
                    }
                }
                continue;
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
	
	
    int[] census4(int u, int v, int e1, int[] u_N, int[] ut_N, int[] v_N, int[] vt_N, int u_deg, int ut_deg, int v_deg, int vt_deg) {
		
        int[] ucount = new int[8];  // triangle counts for each types
        for (int t=0; t<8; t++) ucount[t] = 0;

        for(int i=0,j=0,k=0,l=0; (i<u_deg || j<ut_deg) && (k<v_deg || l < vt_deg); ) {
            int A = (i == u_deg ? n : u_N[i]);
            int B = (j == ut_deg ? n : ut_N[j]);
            int C = (k == v_deg ? n : v_N[k]);
            int D = (l == vt_deg ? n : vt_N[l]);
            if(A < v) {
                i++;
                continue;
            }
            if(B < v) {
                j++;
                continue;
            }
            if(C < v) {
                k++;
                continue;
            }
            if(D < v) {
                l++;
                continue;
            }
            int w = 0;
            int e2 = 0;   // edge vw
            int e3 = 0;   // edge uw
            int type = 0;
// Check the edges to the smallest neighbor.
// Note that each case is mutually disjoint to others.
            if(A < B && A < C && A < D) {
                i++;
            }
            else if(B < A && B < C && B < D) {
                j++;
            }
            else if(C < A && C < B && C < D) {
                k++;
            }
            else if(D < A && D < B && D < C) {
                l++;
            }
            else if(A == B && A < C && A < D) {
                i++; j++;
            }
            else if(C == D && C < A && C < B) {
                k++; l++;
            }
            else if(A == C && A < B && A < D) {
                e2 = 1;  e3 = 2;
                w = A;
                i++; k++;
            }
            else if(A == D && A < B && A < C) {
                e2 = 2;  e3 = 2;
                w = A;
                i++; l++;
            }
            else if(B == C && B < A && B < D) {
                e2 = 1;  e3 = 1;
                w = B;
                j++; k++;
            }
            else if(B == D && B < A && B < C) {
                e2 = 2;  e3 = 1;
                w = B;
                j++; l++;
            }
            else if(A == B && A == C && A < D) {
                e2 = 1;  e3 = 3;
                w = A;
                i++; j++; k++;
            }
            else if(A == B && A == D && A < C) {
                e2 = 2;  e3 = 3;
                w = A;
                i++; j++; l++;
            }
            else if(A == C && A == D && A < B) {
                e2 = 3;  e3 = 2;
                w = A;
                i++; k++; l++;
            }
            else if(B == C && B == D && B < A) {
                e2 = 3;  e3 = 1;
                w = B;
                j++; k++; l++;
            }
            else if(A == B && A == C && A == D ) {
                e2 = 3;  e3 = 3;
                w = A;
                i++; j++; k++; l++;
            }

            if (v < w) {    // to avoid double counting
// Now determine the type:
             //   type = detype(e1,e2,e3);
                type = de7Type[e1][e2][e3];    // use lookup array
// For enumeration:
//			    System.out.println("u,v,w, [e1e2e3], type: " + u + "\t" + v + "\t" + w + "\t[" + e1 + e2 + e3 + "] " + type);
// Add to the count:
                ucount[type]++; 
            }
		    continue;
        }
        return ucount;
    }
	

    int detype(int e1, int e2, int e3){
        int type = 0;
        int T = ((e1 << 2) | e2) << 2 | e3;   
        int t1 = T >> 3;  // the left 3 bits
        int t2 = T & 7;   // the right 3 bits      
        int Sigma = Hamming(T, 6);
        if (Sigma == 6){
            type = 7;
            return type;
        } else if (Sigma == 5){
            type = 6;
            return type;
        } else if (Sigma == 4){
            int Lambda = Hamming(e1 & e2 & e3, 2);
            if (Lambda == 1){
                type = 3;
                return type;
            } else if (Lambda == 0){
                int Lambda3 = Hamming(t1 & t2, 3);
                if (Lambda3 == 2){
                    type = 4;
                    return type;
                } else if (Lambda3 == 1){
                    type = 5;
                    return type;
                }
            } 
        } else if (Sigma == 3){
            int Lambda = Hamming(e1 & e2 & e3, 2);
            if (Lambda == 1){
                type = 1;
                return type;
            } else if (Lambda == 0){
                type = 2;
                return type;
            } 
        } 
        return type;
	}


    int Hamming(int N, int M){
        int ham = 0;
        for (int i=0; i<M; i++){
            if ((N >>> i) % 2 == 1) ++ham;
        }
        return ham;
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
		
        J8SevenTriangles t = new J8SevenTriangles(basename);

        t.computeTriangles();
		
        System.out.println("Total time elapsed = " + (System.currentTimeMillis() - startTime) / 1000.0 + " seconds");
    }

}
