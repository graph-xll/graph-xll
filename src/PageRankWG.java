import java.io.File;
import java.io.PrintStream;
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.stream.DoubleStream;
import java.util.Arrays;
// import java.util.BitSet;
import java.util.Scanner;
import org.apache.commons.math3.util.MathArrays;

public class PageRankWG {
  ImmutableGraph G; // WebGraph object
  ImmutableGraph G_T; // the transpose graph of G
  int n; // number of nodes
  int E; // number of edges
  double[] score; // page-rank score array
  int iteration;
  int MAX_ITER; // max iteration time ~100
  double D; // ~ 0.85
  Integer[] index;
  double eps; // tolerance
  double residual;// residual euclidean distance between PR(t+1) and PR(t)
  // BitSet stop; // if the score change for every node is smaller than the tolerance, stop the program.

  public PageRankWG(String basename, int niter, double damping, double eps) throws Exception {
    G = ImmutableGraph.loadMapped(basename);
    G_T = ImmutableGraph.loadMapped(basename+"-t");
    n = G.numNodes();
    score = new double[n];
    index = new Integer[n];
    Arrays.parallelSetAll(index, (i) -> i);
    Arrays.parallelSetAll(score, (value) -> 1.0/n);
    MAX_ITER = niter;
    D = damping;
    this.eps = eps;
    // stop = new BitSet(n);
    // stop.set(0, n);
    iteration = 0;
  }

  public double[] scoreCompute() {
    while (true) {
      iteration++;
      System.out.println("iteration" + iteration);
      double[]  scoreNow = Arrays.stream(index)
                                 .parallel()
                                 .mapToDouble( u -> {
                                      ImmutableGraph H = G.copy(); // soft copy
                                      ImmutableGraph H_T = G_T.copy();

                                      int[] P_u = H_T.successorArray(u); // get the incoming (linked to u) nodes of u
                                      int[] d_u = new int[P_u.length]; // D_u array stores the outdegree of u's incoming nodes
                                      for (int i = 0; i < d_u.length; i++) {
                                        d_u[i] = H.outdegree(P_u[i]);
                                      }

                                      double sum = 0;
                                      for (int i = 0; i < d_u.length; i++) {
                                        sum += score[P_u[i]] / d_u[i];
                                      }
                                      double u_score = (1 - D) / n + D * sum;
                                      return u_score;})
                                      .toArray();


      residual = MathArrays.distance(scoreNow, score);
      score = scoreNow.clone();

      if (iteration >= MAX_ITER) {
        System.out.println("You have reached the iteration limit.");
        break;
      }

      if (residual < eps) {
        System.out.println("You have reached the tolerance.");
        break;
      }
    }
    return score;
  }

  public static void main(String[] args) throws Exception {
    System.out.println("Please enter graph's basename.");
    Scanner sc = new Scanner(System.in);
    String basename = sc.next();
    System.out.println("Please enter the maximum number of iteration (~100).");
    int niter = sc.nextInt();
    System.out.println("Please enter the damping factor (~0.85).");
    double damping = sc.nextDouble();
    System.out.println("Please enter the tolerance (1E-10, 1E-12, etc.).");
    double eps = Double.parseDouble(sc.next());

    long startTime = System.currentTimeMillis();

    PrintStream ps2 = new PrintStream(new File(basename+"_pageRank_WG_"+niter+"_"+eps+"_timing.txt"));

		System.out.println("Starting " + basename);

    PageRankWG pr = new PageRankWG(basename, niter, damping, eps);

    System.out.println("Time elapsed (sec) for webgraph loading = " + (System.currentTimeMillis() - startTime)/1000.0);
    ps2.println("Time elapsed (sec) for webgraph loading = " + (System.currentTimeMillis() - startTime)/1000.0);

    //storing the score value for each node in a file.
		PrintStream ps = new PrintStream(new File(basename+"_pageRank_WG_"+niter+"_"+eps+".txt"));

    long t2 = System.currentTimeMillis();
	  double[] res = pr.scoreCompute();
    System.out.println("Time elapsed (sec) for score computing = " + (System.currentTimeMillis() - t2)/1000.0);
    ps2.println("Time elapsed (sec) for score computing = " + (System.currentTimeMillis() - t2)/1000.0);
    long t3 = System.currentTimeMillis();

		for(int i=0; i<res.length; i++) {
			ps.println(i + "\t" + res[i]);
		}

    System.out.println("Time elapsed (sec) for file writing = " + (System.currentTimeMillis() - t3)/1000.0);
    ps2.println("Time elapsed (sec) for file writing = " + (System.currentTimeMillis() - t3)/1000.0);

    System.out.println("Total time elapsed (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);
    ps2.println("Total time elapsed (sec) = " + (System.currentTimeMillis() - startTime)/1000.0);

    System.out.println("Total number of iterations = " + pr.iteration);
    ps2.println("Total number of iterations = " + pr.iteration);

    System.out.println("The residual euclidean distance = " + pr.residual);
    ps2.println("The redisual euclidean distance = " + pr.residual);

    long t4 = System.currentTimeMillis();
    System.out.println("Sorting... be patient");
    Runtime.getRuntime().exec("sort -g -r -k2,2 "+basename+"_pageRank_WG_"+niter+"_"+eps+".txt -o "+basename+"_pageRank_WG_"+niter+"_"+eps+".txt").waitFor();
    System.out.println("Time elapsed (sec) for sorting  = " + (System.currentTimeMillis() - t4) / 1000.0);
	}
}
