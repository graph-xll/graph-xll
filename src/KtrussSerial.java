/*
*  K-truss decomposition
*  Serial algorithm with optimized implementation
*/
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.Arrays;
import java.util.BitSet;
import java.io.File;
import java.io.PrintStream;
import java.io.FileOutputStream;
import java.util.Scanner;

public class KtrussSerial {
  final ImmutableGraph G; // WebGraph object
  final int n; // number of nodes
  final int m; // number of edges
  final int[] degree; // degree for each node
  final int[] edgeTail;
  final int[] offset;  // edge offset for edgeHead array
  AtomicIntegerArray support;
  BitSet removed; // to record if an edge is removed or not
  int[] edgePos;
  int[] sortedEdge;

  public KtrussSerial(String basename) throws Exception {
    G = ImmutableGraph.loadMapped(basename);
    n = G.numNodes();
    degree = new int[n];
    offset = new int[n+1];
    // Since we have preprocessed the graph, there are no self-loops anymore.
    // m is the number of edges in the undirected graph.
    m = (int) (G.numArcs()/2);
    edgeTail = new int[m];
    // Populate edgeTail such that u < v.
    System.out.println("Start populating edgeTail.");
    int startPos = 0;
    for (int u = 0; u < n; u++) {
      offset[u] = startPos;
      // Get the neighbors of u
      degree[u] = G.outdegree(u);
      int[] u_successor = G.successorArray(u);
      Arrays.sort(u_successor);
      int count = 0;
      for (int j = 0; j < u_successor.length; j++) {
        int v = u_successor[j];
        if (u <  v) {
          edgeTail[startPos+count] = v;
          count++;
        }
      }
      startPos += count;
    }
    offset[n] = offset[n-1];
    System.out.println("Finished.");
    System.out.println("Total number of edges = " + m);
    removed = new BitSet(m);
    edgePos = new int[m];
    sortedEdge = new int[m];
    support = new AtomicIntegerArray(m);
  }

  public void supportCompute() {
    // Step 1: compute support for each edge.
    System.out.println();
    System.out.println("Starting step 1: support computing for each edge.");
    long startTime1 = System.currentTimeMillis();
    // Since m is huge, we divide the edge list into segments.
    // We parallelize each segment at a time.
    int numOfSeg = 10000;
    int segLength = m/numOfSeg;
    int residual = m%numOfSeg;
    int[] index = new int[segLength];
    for (int segCount = 0; segCount < numOfSeg; segCount++) {
      int start = segCount * segLength;
      // Each index is the index in edgeHead (or edgeTail) array.
      Arrays.parallelSetAll(index, (i) -> (start+i));
      Arrays.stream(index)
            .parallel()
            .boxed()
            .mapToInt(e -> {
              int u = findHead(e);
              int v = edgeTail[e];
              // We need to identify triangles containing edge uv.
              // We define a triangle to be uvw such that u < v < w.
              // We just need to find all w.
              // Do the intersection of w1 and w2
              int w1 = e + 1;
              int w2 = offset[v];
              while (w1 < offset[u+1] && w2 < offset[v+1]) {
                if (edgeTail[w1] < edgeTail[w2]) {
                  w1++;
                } else if (edgeTail[w1] > edgeTail[w2]) {
                  w2++;
                } else {
                  // edgeTail[w1] = edgeTail[w2]
                  // Now we find a triangle.
                  // We add 1 to edge uv, vw, and uw's support
                  support.getAndIncrement(e);
                  support.getAndIncrement(w1);
                  support.getAndIncrement(w2);
                  w1++;
                  w2++;
                }
              }
              return 0;
      }).sum();
      // System.out.println("Calculated support for " + (segCount+1)*100.0/numOfSeg + "% edges.");
    }
    // Calculate support for the residual edges
    for (int k = 0; k < residual; k++) {
      int start = m - residual;
      int e = start + k;
      int u = findHead(e);
      int v = edgeTail[e];
      // We need to identify triangles containing edge uv.
      // We define a triangle to be uvw such that u < v < w.
      // We just need to find all w.
      // Do the intersection of w1 and w2
      int w1 = e + 1;
      int w2 = offset[v];
      while (w1 < offset[u+1] && w2 < offset[v+1]) {
        if (edgeTail[w1] < edgeTail[w2]) {
          w1++;
        } else if (edgeTail[w1] > edgeTail[w2]) {
          w2++;
        } else {
          // edgeTail[w1] = edgeTail[w2]
          // Now we find a triangle.
          // We add 1 to edge uv, vw, and uw's support
          support.getAndIncrement(e);
          support.getAndIncrement(w1);
          support.getAndIncrement(w2);
          w1++;
          w2++;
        }
      }
    }
    System.out.println("maximum support: " + findMaxSupport(support));
    long endTime1 = System.currentTimeMillis();
    System.out.println("Time elapsed (sec) for support computing = " + (endTime1 - startTime1)/1000.0);
  }

  private int findMaxSupport(AtomicIntegerArray support) {
    int max = 0;
    for (int i = 0; i < support.length(); i++) {
      int temp = support.get(i);
      if (temp > max) {
        max = temp;
      }
    }
    return max;
  }

  private int findHead(int target) {
    // find index such that offset[index] <= target
    // and target < offset[index+1]
    int start = 0;
    int end = offset.length - 1;
    while (start + 1 < end) {
      int mid = start + (end - start)/2;
      if (offset[mid] <= target) {
        start = mid;
      } else {
        end = mid;
      }
    }
    if (offset[end] <= target) {
      return end;
    }
    if (offset[start] <= target) {
      return start;
    }
    return 0;
  }

  public void ktrussCompute() {
    // Step 2: sort the edges in ascending order of their support
    // Store the position information in the edgePos array.
    System.out.println();
    System.out.println("Starting step 2: sorting the edges in ascending order of their support.");
    long startTime2 = System.currentTimeMillis();
    int maxSupport = findMaxSupport(support);
    int[] bin = new int[maxSupport+2];
    // Find the number of edges in each bin
    for (int i = 0; i < m; i++) {
      bin[support.get(i)]++;
    }
    // Modify the bin array such that
    // bin[i] is the total number of edges with suport < i.
    // So the number of edges with support = i will be bin[i+1] - bin[i].
    // bin[0] = 0
    int sum = 0;
    for (int i = 0; i <= maxSupport + 1; i++) {
      int temp = bin[i];
      bin[i] = sum;
      sum += temp;
    }
    // Do bin-sort of the edges
    // sortedEdge stores the edges in the ascending order of support.
    // edgePos stores the position of an edge in the sortedEdge array.
    for (int i = 0; i < m; i++) {
      edgePos[i] = bin[support.get(i)];
      sortedEdge[edgePos[i]] = i;
      bin[support.get(i)]++;
    }
    // Restore the bin array
    for (int i = maxSupport; i >= 1; i--) {
      bin[i] = bin[i-1];
    }
    bin[0] = 0;
    long endTime2 = System.currentTimeMillis();
    System.out.println("Time elapsed (sec) for edge sorting = " + (endTime2 - startTime2)/1000.0);

    // Step 3: k-truss computation
    System.out.println();
    System.out.println("Starting step 3: k-truss computation.");
    long startTime3 = System.currentTimeMillis();
    int k = 2; // (k >= 2)
    int edge = 0; // index in the sortedEdge array
    int edgeRemoved = 0; // number of edges removed
    while (edgeRemoved < m) {
      System.out.print("k = " + k + " ");
      // Iterate over the edges in sortedEdge with suport <= k-2
      while (edge < m && support.get(sortedEdge[edge]) <= k-2) {
        int id_uv = sortedEdge[edge];
        int u = findHead(id_uv); // edge's head
        int v = edgeTail[id_uv]; // edge's tail
        // Make sure that degree(u) < degree(v)
        if (degree[u] > degree[v]) {
          // Swap u and v
          int temp = u;
          u = v;
          v = temp;
        }
        // For each neighbor of u
        for (int w : G.successorArray(u)) {
          // Check if <u, w> is an edge or not
          // Check if <u, w> has been removed or not
          // If it is an edge, we should get the egde id.
          // Let's first find the index range in edgeHead array using offset.
          // Then use binary search in edgeTail array within the range.
          int id_uw = Arrays.binarySearch(edgeTail, offset[Math.min(u,w)], offset[Math.min(u,w) + 1], Math.max(u,w));
          // If <u,w> is an edge and has not been removed
          if (id_uw >= 0 && !removed.get(id_uw)) {
            // Check if <v, w> is an edge or not
            // Check if <v, w> has been removed or not
            if (w != v) {
              int id_vw = Arrays.binarySearch(edgeTail, offset[Math.min(v,w)], offset[Math.min(v,w) + 1], Math.max(v,w));
              // If <v, w> is an edge and not removed
              if (id_vw >= 0 && !removed.get(id_vw)) {
                if (support.get(id_vw) > k - 2) {
                  // Now we need to decrease the support of edge vw by 1
                  // And also we need to modify vw's position in the sortedEdge array
                  // as well as the bin array.
                  int support_vw = support.get(id_vw); // support for edge vw
                  int edgePos_vw = edgePos[id_vw]; // vw's position in the sortedEdge array
                  int startPos_vw  = bin[support_vw];// start position for support_vw
                  int id_first_vw = sortedEdge[startPos_vw];// edge id for the start position for the support_vw bin
                  // Swap id_vw and id_first_vw
                  // Modify the sortedEdge and the edgePos array
                  if (id_vw != id_first_vw) {
                    edgePos[id_first_vw] = edgePos_vw;
                    edgePos[id_vw] = startPos_vw;
                    sortedEdge[startPos_vw] = id_vw;
                    sortedEdge[edgePos_vw] = id_first_vw;
                  }
                  // Right shift the start position of the support_vw bin by 1
                  bin[support_vw]++;
                  // Decrease the support of edge vw by 1
                  support.getAndDecrement(id_vw);
                }

                // Do the same process on edge uw
                if (support.get(id_uw) > k - 2) {
                  // Now we need to decrease the support of edge uw by 1
                  // And also we need to modify uw's position in the sortedEdge array
                  // as well as the bin array.
                  int support_uw = support.get(id_uw); // support for edge uw
                  int edgePos_uw = edgePos[id_uw]; // uw's position in the sortedEdge array
                  int startPos_uw  = bin[support_uw];// start position for support_uw
                  int id_first_uw = sortedEdge[startPos_uw];// edge id for the start position for the support_uw bin
                  // Swap id_uw and id_first_uw
                  // Modify the sortedEdge and the edgePos array
                  if (id_uw != id_first_uw) {
                    edgePos[id_first_uw] = edgePos_uw;
                    edgePos[id_uw] = startPos_uw;
                    sortedEdge[startPos_uw] = id_uw;
                    sortedEdge[edgePos_uw] = id_first_uw;
                  }
                  // Right shift the start position of the support_uw bin by 1
                  bin[support_uw]++;
                  // Decrease the support of edge uw by 1
                  support.getAndDecrement(id_uw);
                }
              }
            }
          }
        }
        // Remove "edge" from G
        removed.set(id_uv);
        edgeRemoved++;
        edge++;
      }
      System.out.println("Edge removed = " + edgeRemoved);
      k++;
    }
    long endTime3 = System.currentTimeMillis();
    System.out.println("Time elapsed (sec) for k-truss computation = " + (endTime3 - startTime3)/1000.0);

    // Find the maximum truss
    // kmax is the maximum truss.
    // Or we can find the maximum support in the current support array.
    // Then max truss = max support + 2
    int maxTruss = findMaxSupport(support) + 2;
    System.out.println("Max truss = " + maxTruss);
  }

  public static void main(String[] args) throws Exception{
    System.out.println("Please enter graph's base name");
    Scanner sc = new Scanner(System.in);
    String basename = sc.next();

    KtrussSerial kt = new KtrussSerial(basename+"-noself-sym");
    kt.supportCompute();
    kt.ktrussCompute();

    // System.out.println("Writing results to file. Be patient...");
    // PrintStream ps = new PrintStream(new File(basename+"_ktruss_sequential.txt"));
    // for (int i = 0; i < kt.m; i++) {
    //   ps.println(kt.support.get(i)+2);
    // }
  }
}
