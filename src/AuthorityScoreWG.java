import java.io.File;
import java.io.PrintStream;
import it.unimi.dsi.webgraph.ImmutableGraph;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.Arrays;
import org.apache.commons.math3.stat.StatUtils;
import java.util.Scanner;
import org.apache.commons.math3.util.MathArrays;

/*
* t(A)*A
*/

public class AuthorityScoreWG {
  ImmutableGraph G; // WebGraph object
  ImmutableGraph G_T; // transpose
  int n; // number of nodes
  double[] score; // hub score array
  double lambda; // eigen value
  int iteration;
  int MAX_ITER; // max iteration time ~100
  Integer[] index;
  double eps; // tolerance
  double residual;// residual euclidean distance between score(t+1) and score(t)

  public AuthorityScoreWG(String basename,int niter,double eps) throws Exception {
    G = ImmutableGraph.loadMapped(basename);
    G_T = ImmutableGraph.loadMapped(basename+"-t");
    n = G.numNodes();
    score = new double[n];
    index = new Integer[n];
    Arrays.parallelSetAll(index, (i) -> i);
    Arrays.parallelSetAll(score, (value) -> 1.0);
    MAX_ITER = niter;
    this.eps = eps;
    lambda = 1.0;
    iteration = 0;
  }

  // the scoreCompute method
  public double[] scoreCompute() {
    while (true) {
      iteration++;
      System.out.println("iteration" + iteration);
      // A*x
      double[] temp1 = Arrays.stream(index)
                    .parallel()
                    .mapToDouble( u -> {
                    ImmutableGraph H = G.copy(); // soft copy
                    int[] u_seccessor = H.successorArray(u); // get the successor nodes of u
                    double sum = 0;
                    for (int i = 0; i < u_seccessor.length; i++) {
                      sum += score[u_seccessor[i]];
                    }
                    return sum;})
                    .toArray();
      // t(A)*A*x
      double[] temp2 = Arrays.stream(index)
                    .parallel()
                    .mapToDouble( u -> {
                    ImmutableGraph H_T = G_T.copy(); // soft copy
                    int[] u_seccessor = H_T.successorArray(u); // get the successor nodes of u
                    double sum = 0;
                    for (int i = 0; i < u_seccessor.length; i++) {
                      sum += temp1[u_seccessor[i]];
                    }
                    return sum;})
                    .toArray();

      // find the max value in score which will be the eigen value
      double eigenValue = StatUtils.max(temp2);
      // normalize the vector
      double norm = MathArrays.safeNorm(temp2);
      temp2 = MathArrays.scale(1.0/norm, temp2);
      // calculate the residual euclidean distance
      residual = MathArrays.distance(temp2, score);

      lambda = eigenValue;
      score = temp2.clone();

      // conditions to stop
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
    System.out.println("Please enter the tolerance (1E-10, 1E-12, etc.).");
    double eps = Double.parseDouble(sc.next());

    long startTime = System.currentTimeMillis();

    PrintStream ps2 = new PrintStream(new File(basename+"_authorityScore_WG_"+niter+ "_" + eps+"_timing.txt"));

		System.out.println("Starting " + basename);
		AuthorityScoreWG authority = new AuthorityScoreWG(basename,niter,eps);

    System.out.println("Time elapsed (sec) for webgraph loading = " + (System.currentTimeMillis() - startTime)/1000.0);
    ps2.println("Time elapsed (sec) for webgraph loading = " + (System.currentTimeMillis() - startTime)/1000.0);
		//storing the score value for each node in a file.
		PrintStream ps = new PrintStream(new File(basename+"_authorityScore_WG_"+niter+ "_"+ eps+".txt"));

    long t2 = System.currentTimeMillis();
	  double[] res = authority.scoreCompute();
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

    System.out.println("Total number of iterations = " + authority.iteration);
    ps2.println("Total number of iterations = " + authority.iteration);

    System.out.println("The residual euclidean distance = " + authority.residual);
    ps2.println("The redisual euclidean distance = " + authority.residual);

    System.out.println("The eigen value = " + authority.lambda);
    ps2.println("The eigen value = " + authority.lambda);

    long t4 = System.currentTimeMillis();
    System.out.println("Sorting... be patient");
    Runtime.getRuntime().exec("sort -g -r -k2,2 -k1,1n "+basename+"_authorityScore_WG_"+niter+"_"+eps+".txt -o "+basename+"_authorityScore_WG_"+niter+"_"+eps+".txt").waitFor();
    System.out.println("Time elapsed (sec) for sorting  = " + (System.currentTimeMillis() - t4) / 1000.0);


	}
}
