package se.kth.jabeja;

import org.apache.log4j.Logger;
import se.kth.jabeja.config.Config;
import se.kth.jabeja.config.NodeSelectionPolicy;
import se.kth.jabeja.io.FileIO;
import se.kth.jabeja.rand.RandNoGenerator;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class Jabeja {
  final static Logger logger = Logger.getLogger(Jabeja.class);
  private final Config config;
  private final HashMap<Integer/* id */, Node/* neighbors */> entireGraph;
  private final List<Integer> nodeIds;
  private int numberOfSwaps;
  private int round;
  private float T;
  private boolean resultFileCreated = false;
  private boolean differentAnnealing = true;
  private float alpha;
  private float t_min = 0.00001f;
  private int annealingCount = 0;
  private int resetMax = 30;
  private int reset = resetMax;
  private int numOfResets = 0;
  private int minEdgeCut = Integer.MAX_VALUE;
  private int linearAnnealingReset = 400;
  private enum AnnealingType {
    EXPONENTIAL, LOGISTIC, HEAVY_TAILED
  }
  private AnnealingType anneal = AnnealingType.HEAVY_TAILED;
  private Random random = new Random();

  // -------------------------------------------------------------------
  public Jabeja(HashMap<Integer, Node> graph, Config config) {
    this.entireGraph = graph;
    this.nodeIds = new ArrayList(entireGraph.keySet());
    this.round = 0;
    this.numberOfSwaps = 0;
    this.config = config;
    this.T = config.getTemperature();
    if (differentAnnealing) {
      this.T = 1.0f;
      this.alpha = 0.93f;
    }
  }

  // -------------------------------------------------------------------
  public void startJabeja() throws IOException {
    for (round = 0; round < config.getRounds(); round++) {
      for (int id : entireGraph.keySet()) {
        sampleAndSwap(id);
      }

      // one cycle for all nodes have completed.
      // reduce the temperature
      saCoolDown();
      report();
    }
    System.out.println("Annealing type: " + anneal);
    System.out.println("Total resets: " + numOfResets);
    System.out.println("Minimum Edge Cut: " + minEdgeCut);
  }

  /**
   * Simulated analealing cooling function
   */
  private void saCoolDown() {
    if (differentAnnealing) {
      T = T * alpha;
      if (T < t_min) {
        T = t_min;
        reset -= 1;
        if (reset == 0) {
          numOfResets++;
          T = 1 - (numOfResets * 0.1f);
          reset = resetMax;
        }
      }
    } else {
      if (round == linearAnnealingReset) {
        T = config.getTemperature();
      }
      if (T > 1)
        T -= config.getDelta();
      if (T < 1)
        T = 1;
    }
  }

  /**
   * Sample and swap algorith at node p
   * 
   * @param nodeId
   */
  private void sampleAndSwap(int nodeId) {
    Node partner = null;
    Node nodep = entireGraph.get(nodeId);

    if (config.getNodeSelectionPolicy() == NodeSelectionPolicy.HYBRID
        || config.getNodeSelectionPolicy() == NodeSelectionPolicy.LOCAL) {
      Integer[] neighborNodes = getNeighbors(nodep);
      Node bestPartner = findPartner(nodeId, neighborNodes);
      partner = bestPartner;
    }

    if (config.getNodeSelectionPolicy() == NodeSelectionPolicy.HYBRID
        || config.getNodeSelectionPolicy() == NodeSelectionPolicy.RANDOM) {
      if (partner != null) {

      } else {
        Integer[] neighborNodes = getSample(nodeId);
        Node bestPartner = findPartner(nodeId, neighborNodes);
        partner = bestPartner;
      }
    }
    if (partner != null) {
      numberOfSwaps++;
      int tempColor = nodep.getColor();
      nodep.setColor(partner.getColor());
      partner.setColor(tempColor);
    }
  }

  public Node findPartner(int nodeId, Integer[] nodes) {

    Node nodep = entireGraph.get(nodeId);

    Node bestPartner = null;
    double highestBenefit = 0;

    Node[] neighborNodes = new Node[nodes.length];
    int idx = 0;

    for (Integer qId : nodes) {
      Node nodeq = entireGraph.get(qId);

      neighborNodes[idx++] = nodeq;

      if (nodep.getColor() == nodeq.getColor()) {
        continue;
      }

      int dp = getDegree(nodep, nodep.getColor());
      int dq = getDegree(nodeq, nodeq.getColor());
      int dpNew = getDegree(nodep, nodeq.getColor());
      int dqNew = getDegree(nodeq, nodep.getColor());

      double benefitCurrent = Math.pow(dp, config.getAlpha()) + Math.pow(dq, config.getAlpha());
      double benefitNew = Math.pow(dpNew, config.getAlpha()) + Math.pow(dqNew, config.getAlpha());

      if (differentAnnealing) {
        if (benefitNew > benefitCurrent) {
          if (benefitNew > highestBenefit) {
            highestBenefit = benefitNew;
            bestPartner = nodeq;
          }
        }
      } else {
        if (benefitNew * T > benefitCurrent) {
          if (benefitNew > highestBenefit) {
            highestBenefit = benefitNew;
            bestPartner = nodeq;
          }
        }
      }
    }

    if (bestPartner == null && differentAnnealing && T > t_min) {
      annealingCount++;
      int iterations = 100;
      for (int i = 0; i < iterations; i++) {
        int randomIndex = random.nextInt();
        Node randNeighbor = neighborNodes[Math.abs(randomIndex) % neighborNodes.length];
        float acceptance = acceptanceProbability(nodep, randNeighbor);
        if (acceptance > random.nextFloat()) {
          bestPartner = randNeighbor;
          break;
        }
      }
    }
    return bestPartner;
  }

  private float acceptanceProbability(Node nodep, Node nodeq) {
    int dp = getDegree(nodep, nodep.getColor());
    int dq = getDegree(nodeq, nodeq.getColor());
    int dpNew = getDegree(nodep, nodeq.getColor());
    int dqNew = getDegree(nodeq, nodep.getColor());

    double benefitCurrent = Math.pow(dp, config.getAlpha()) + Math.pow(dq, config.getAlpha());
    double benefitNew = Math.pow(dpNew, config.getAlpha()) + Math.pow(dqNew, config.getAlpha());

    if (benefitNew + 5 < benefitCurrent) {
      return 0.0f;
    }

    if (anneal == AnnealingType.EXPONENTIAL) {
      return (float) Math.pow(Math.E, (benefitNew - benefitCurrent) / T);
    } else if (anneal == AnnealingType.LOGISTIC) {
      return (float) (1.0 / (1.0 + Math.exp(5*(benefitCurrent - benefitNew) / T)));
    } else {
      return (float) Math.pow(1 + ((benefitCurrent - benefitNew) / (3*T)), -3);
    }
  }

  /**
   * The the degreee on the node based on color
   * 
   * @param node
   * @param colorId
   * @return how many neighbors of the node have color == colorId
   */
  private int getDegree(Node node, int colorId) {
    int degree = 0;
    for (int neighborId : node.getNeighbours()) {
      Node neighbor = entireGraph.get(neighborId);
      if (neighbor.getColor() == colorId) {
        degree++;
      }
    }
    return degree;
  }

  /**
   * Returns a uniformly random sample of the graph
   * 
   * @param currentNodeId
   * @return Returns a uniformly random sample of the graph
   */
  private Integer[] getSample(int currentNodeId) {
    int count = config.getUniformRandomSampleSize();
    int rndId;
    int size = entireGraph.size();
    ArrayList<Integer> rndIds = new ArrayList<Integer>();

    while (true) {
      rndId = nodeIds.get(RandNoGenerator.nextInt(size));
      if (rndId != currentNodeId && !rndIds.contains(rndId)) {
        rndIds.add(rndId);
        count--;
      }

      if (count == 0)
        break;
    }

    Integer[] ids = new Integer[rndIds.size()];
    return rndIds.toArray(ids);
  }

  /**
   * Get random neighbors. The number of random neighbors is controlled using
   * -closeByNeighbors command line argument which can be obtained from the config
   * using {@link Config#getRandomNeighborSampleSize()}
   * 
   * @param node
   * @return
   */
  private Integer[] getNeighbors(Node node) {
    ArrayList<Integer> list = node.getNeighbours();
    int count = config.getRandomNeighborSampleSize();
    int rndId;
    int index;
    int size = list.size();
    ArrayList<Integer> rndIds = new ArrayList<Integer>();

    if (size <= count)
      rndIds.addAll(list);
    else {
      while (true) {
        index = RandNoGenerator.nextInt(size);
        rndId = list.get(index);
        if (!rndIds.contains(rndId)) {
          rndIds.add(rndId);
          count--;
        }

        if (count == 0)
          break;
      }
    }

    Integer[] arr = new Integer[rndIds.size()];
    return rndIds.toArray(arr);
  }

  /**
   * Generate a report which is stored in a file in the output dir.
   *
   * @throws IOException
   */
  private void report() throws IOException {
    int grayLinks = 0;
    int migrations = 0; // number of nodes that have changed the initial color
    int size = entireGraph.size();

    for (int i : entireGraph.keySet()) {
      Node node = entireGraph.get(i);
      int nodeColor = node.getColor();
      ArrayList<Integer> nodeNeighbours = node.getNeighbours();

      if (nodeColor != node.getInitColor()) {
        migrations++;
      }

      if (nodeNeighbours != null) {
        for (int n : nodeNeighbours) {
          Node p = entireGraph.get(n);
          int pColor = p.getColor();

          if (nodeColor != pColor)
            grayLinks++;
        }
      }
    }

    int edgeCut = grayLinks / 2;

    if (edgeCut < minEdgeCut) {
      minEdgeCut = edgeCut;
      System.out.println("New min edge cut: " + minEdgeCut + " at round " + round);
    }

    logger.info("round: " + round +
        ", edge cut:" + edgeCut +
        ", swaps: " + numberOfSwaps +
        ", migrations: " + migrations);

    saveToFile(edgeCut, migrations);
  }

  private void saveToFile(int edgeCuts, int migrations) throws IOException {
    String delimiter = "\t\t";
    String outputFilePath;

    // output file name
    File inputFile = new File(config.getGraphFilePath());
    outputFilePath = config.getOutputDir() +
        File.separator +
        inputFile.getName() + "_" +
        "NS" + "_" + config.getNodeSelectionPolicy() + "_" +
        "GICP" + "_" + config.getGraphInitialColorPolicy() + "_" +
        "T" + "_" + config.getTemperature() + "_" +
        "D" + "_" + config.getDelta() + "_" +
        "RNSS" + "_" + config.getRandomNeighborSampleSize() + "_" +
        "URSS" + "_" + config.getUniformRandomSampleSize() + "_" +
        "A" + "_" + config.getAlpha() + "_" +
        "R" + "_" + config.getRounds() + ".txt";

    if (!resultFileCreated) {
      File outputDir = new File(config.getOutputDir());
      if (!outputDir.exists()) {
        if (!outputDir.mkdir()) {
          throw new IOException("Unable to create the output directory");
        }
      }
      // create folder and result file with header
      String header = "# Migration is number of nodes that have changed color.";
      header += "\n\nRound" + delimiter + "Edge-Cut" + delimiter + "Swaps" + delimiter + "Migrations" + delimiter
          + "Skipped" + "\n";
      FileIO.write(header, outputFilePath);
      resultFileCreated = true;
    }

    FileIO.append(round + delimiter + (edgeCuts) + delimiter + numberOfSwaps + delimiter + migrations + "\n",
        outputFilePath);
  }
}
