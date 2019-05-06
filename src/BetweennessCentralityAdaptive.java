import java.util.List;
import java.io.File;
import java.io.PrintStream;
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.Arrays;
import java.util.ArrayList;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.random.RandomDataGenerator;
import java.util.Scanner;

// adaptive sampling
public class BetweennessCentralityAdaptive {
  ImmutableGraph G; // WebGraph object
  ImmutableGraph G_T;
  int n; // number of nodes
  double[] score; // Betweenness centrality score array
  Integer[] index;
  int k; // actual number of pivots
  double C;
  int[] pivots; // pivots index, randomly chosen
  int CORE; // number of cores of your CPU has
  int cutoff; // maximum number of pivots
  int topK;  // number of top-score vertices to look at


  public BetweennessCentralityAdaptive(String basename, int core, double c) throws Exception {
    G = ImmutableGraph.loadMapped(basename);
    G_T = ImmutableGraph.loadMapped(basename+"-t");
    n = G.numNodes();
    score= new double[n];
    CORE = core;
    Arrays.parallelSetAll(score, (i) -> 0);
    index = new Integer[CORE];
    Arrays.parallelSetAll(index, (i) -> i);
    RandomDataGenerator rand = new RandomDataGenerator();
    pivots = rand.nextPermutation(n,n);
    this.C = c;
    k = 0; // this number should be CORE * some integer
    cutoff = n / 50;
    topK = 1000; // we focus on the top 1000 high-score vertices
  }

 // compute betweenness centrality score
  public double[] scoreCompute() {
    while (k < cutoff) {
      double[] result = Arrays.stream(index)
                    .parallel()
                    .map( u -> {
                    double[] scoreNow = new double[n];
                    int[] p = pivots.clone();
                    int source = p[u+k];
                    System.out.println("Start visiting node " + source);
                    ImmutableGraph H = G.copy();
                    ImmutableGraph H_T = G_T.copy();
                    // queue stores the visited nodes in BFS order
                    List<Integer> queue = new ArrayList<Integer>();
                    // breakPoint stores start positions (index in queue) of each BFS layer
                    // For example:
                    // breakPoint[0]=0  -->     *         depth = 0
                    // breakPoint[1]=1  -->    * *        depth = 1
                    // breakPoint[2]=3  -->   * * *       depth = 2
                    // breakPoint[3]=6  -->  * * * *      depth = 3
                    // breakPoint[4]=10 --> * * * * *     depth = 4
                    //
                    // breakPoint[0] = 0 (index (in queue) of the root node),
                    // breakPoint[1] = 1 (index (in queue) of the first node of layer 1 (depth = 1))
                    // We can think this number as the number of nodes above this layer.
                    // breakPoint[2] = 3 (index (in queue) of the first node of layer 2 (depth = 2))
                    // or the number of nodes above this layer is 3.
                    // breakPoint[k] = (index (in queue) of the first node of layer k (depth = k))
                    // The last element of breakPoint will be the size of queue which is n.
                    List<Integer> breakPoint = new ArrayList<Integer>();
                    // sigma[v] stores the number of shortest paths from source to v
                    long[] sigma = new long[n];
                    sigma[source] = 1;
                    // distance[v] stores the distance from source to v
                    long[] distance = new long[n];
                    Arrays.parallelSetAll(distance, (i) -> -1);
                    distance[source] = 0;
                    // use BFS to find shortest paths
                    queue.clear();
                    queue.add(source);
                    breakPoint.clear();
                    breakPoint.add(0);

                    int depth; // depth of the BFS spanning tree
            				for(depth = 0; queue.size() != breakPoint.get(breakPoint.size() - 1); depth++) {
            					breakPoint.add(queue.size());
            					int start = breakPoint.get(depth);
            					int end = breakPoint.get(depth + 1);

            					for(int i = start; i < end; i++) {
            						int v = queue.get(i);
                        long sigma_v = sigma[v];
            						for(int w : H.successorArray(v)) {
            							if (distance[w] == -1) {
            								distance[w] = depth + 1;
            								queue.add(w);
            							}
                          if (distance[w] == depth + 1) {
            								sigma[w] += sigma_v;
            							}
            						}
            					}
            				}

                   System.out.println("Ending BFS in node " + source);
                    // accumulation: back-propagation of dependencies
                    double[] delta = new double[n];
                    while(--depth > 0) {
            					int start = breakPoint.get(depth);
            					int end = breakPoint.get(depth + 1);

            					for(int i = start; i < end; i++) {
            						int w = queue.get(i);
            						for(int v : H_T.successorArray(w)) {
            							if (distance[v] == depth - 1) {
                            delta[v] += (1 + delta[w]) * sigma[v] / sigma[w];
                          }
                        }
                        if (w != source) {
                          scoreNow[w] += delta[w];
                        }

            					}
            				}
                    System.out.println("node " + source + " visited.");
                    //System.gc();
                    return scoreNow;})
                    .reduce(new double[n],(x,y) -> MathArrays.ebeAdd(x,y));

      k += CORE;
      score = MathArrays.ebeAdd(result, score);
      int count = Arrays.stream(score).parallel().mapToInt(u -> (u > n*C)?1:0).sum();
      if (count > topK) {
        break;
      }
    }

    // rescale score
    double scale = 1.0*n/(n-1)/(n-2)/k;
    score = MathArrays.scale(scale, score);

    return score;
  }


  public static void main(String[] args) throws Exception {
    System.out.println("Please enter graph's basename.");
    Scanner sc = new Scanner(System.in);
    String basename = sc.next();
    System.out.println("Please enter the number of cores your CPU has.");
    int core = sc.nextInt();
    System.out.println("Please enter the constant C (2,3,...).");
    double c = sc.nextDouble();

		long startTime = System.currentTimeMillis();

    PrintStream ps2 = new PrintStream(new File(basename+"_betweennessCentrality_WG_"+c+"_timing.txt"));

		System.out.println("Starting " + basename);
		BetweennessCentralityAdaptive bc = new BetweennessCentralityAdaptive(basename,core,c);

    System.out.println("Time elapsed (sec) for webgraph loading = " + (System.currentTimeMillis() - startTime)/1000.0);
    ps2.println("Time elapsed (sec) for webgraph loading = " + (System.currentTimeMillis() - startTime)/1000.0);
		//storing the score value for each node in a file.
		PrintStream ps = new PrintStream(new File(basename+"_betweennessCentrality_WG_"+c+".txt"));

    long t2 = System.currentTimeMillis();
	  double[] res = bc.scoreCompute();
    System.out.println("Time elapsed (sec) for score computing = " + (System.currentTimeMillis() - t2)/1000.0);
    ps2.println("Time elapsed (sec) for score computing = " + (System.currentTimeMillis() - t2)/1000.0);
    System.out.println("The number of samples k = " + bc.k);
    ps2.println("The number of samples k = " + bc.k);
    long t3 = System.currentTimeMillis();

		for(int i=0; i<res.length; i++) {
			ps.println(i + "\t" + res[i]);
		}

    System.out.println("Time elapsed (sec) for file writing = " + (System.currentTimeMillis() - t3)/1000.0);
    ps2.println("Time elapsed (sec) for file writing = " + (System.currentTimeMillis() - t3)/1000.0);

    System.out.println("Total time elapsed (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);
    ps2.println("Total time elapsed (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);

    long t4 = System.currentTimeMillis();
    System.out.println("Sorting... be patient");
    Runtime.getRuntime().exec("sort -g -r -k2,2 -k1,1n "+basename+"_betweennessCentrality_WG_"+c+".txt -o "+basename+"_betweennessCentrality_WG_"+c+"_sorted.txt").waitFor();
    System.out.println("Time elapsed (sec) for sorting  = " + (System.currentTimeMillis() - t4) / 1000.0);


	}
}
