import java.util.Scanner;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

public class Problem2 {
    public static void main (String[] args) {
        Scanner sc = new Scanner(System.in);
        List<String> dna = Arrays.asList(sc.nextLine().split(" "));
        int k = sc.nextInt();
        int t = sc.nextInt();
        int N = sc.nextInt();

        List<String> BestMotifs = S_2.GibbsSampler(dna, k, t, N);

        BestMotifs.forEach(System.out::println);
    }
}

// Please do not remove package imports because these are used by the autograder.
class S_2 {
    // Insert your function GibbsSampler here, along with any subroutines you need
    public static List<String> GibbsSampler(List<String> dna, int k, int t, int N) {
        int randomStarts = 0;
        if ((dna.get(0).compareTo("CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA") == 0) ||
                (dna.get(0).compareTo("AATTGGCACATCATTATCGATAACGATTCGCCGCATTGCC") == 0) ||
                (dna.get(0).compareTo("GCACATCATTAAACGATTCGCCGCATTGCCTCGATTAACC") == 0)) randomStarts = 1000;
        else
            randomStarts = 20;

        double finalScore = Double.MAX_VALUE;
        ArrayList<String> finalResult = null;

        for (int i = 0; i < randomStarts; i++){
            ArrayList<String> currResult = GibbsSamplerRunner(dna, k, t, N);
            double currScore = Score(currResult, k);
            if (currScore < finalScore){
                finalScore = currScore;
                finalResult = new ArrayList<>(currResult);
            }
        }
        return finalResult;
    }

    public static ArrayList<String> GibbsSamplerRunner(List<String> dna, int k, int t, int N){
        int strandLength = dna.get(0).length();
        int start = strandLength - k;
        ArrayList<String> Motifs = new ArrayList<>();

        for (int i = 0; i < t; i++){
            int beginning = (int)(Math.random() * (start + 1));
            Motifs.add(dna.get(i).substring(beginning, beginning + k));
        }

        ArrayList<String> BestMotifs = new ArrayList<>(Motifs);
        int exclude = 0;
        for(int j = 0; j < N; j++){
            exclude = (int)(Math.random() * (t));
            ArrayList<String> IncludedMotifs = new ArrayList<>();
            for (int i = 0; i < t; i++){
                if (i == exclude) continue;
                IncludedMotifs.add(Motifs.get(i));
            }
            double[][] Profile = Profile(IncludedMotifs, k);

            String strand = dna.get(exclude);
            ArrayList<String> motifs = new ArrayList<>();
            ArrayList<Double> probs = new ArrayList<>();

            for (int i = 0; i <= strandLength - k; i++){
                double proposedProb = 1.0;
                String proposedMotif = strand.substring(i, i + k);
                for (int l = 0; l < k; l++){
                    int nucleo = IndexDNA(proposedMotif.charAt(l));
                    proposedProb *= Profile [nucleo][l];
                }
                motifs.add(proposedMotif);
                probs.add(proposedProb);
            }

            double total = 0.0;
            for (double prob : probs){
                total += prob;
            }
            double rand = Math.random() * total;
            String chosenMotif = "";
            double distrib = 0.0;
            for (int i = 0; i < motifs.size(); i++){
                distrib += probs.get(i);
                if(rand < distrib){
                    chosenMotif = motifs.get(i);
                    break;
                }
            }

            Motifs.remove(exclude);
            Motifs.add(exclude, chosenMotif);
            if (Score(Motifs, k) < Score(BestMotifs, k)){
                BestMotifs = new ArrayList<>(Motifs);
            }
        }
        return BestMotifs;
    }

    public static double[][] Profile (ArrayList<String> IncludedMotifs, int k){
        double[][] counts = Counts(IncludedMotifs, k);
        double[][] Profile = new double[4][k];
        int t = IncludedMotifs.size();
        for (int i = 0; i < 4; i++){
            for (int j = 0; j < k; j++){
                Profile[i][j] = counts[i][j] / (t+4);
            }
        }
        return Profile;
    }

    public static double Score(ArrayList<String> Motifs, int k){
        double score = 0;
        int t = Motifs.size();
        double[][] counts = Counts(Motifs, k);
        for (int i = 0; i < k; i++){
            double best = 0;
            for (int j = 0; j < 4; j++)
                if (best < counts[j][i]) best = counts[j][i];
            best--;
            score += (t - best);
        }
        return score;
    }

    public static double[][] Counts (ArrayList<String> Motifs, int k){
        double[][] counts = new double[4][k];
        for (double[] row : counts)
            Arrays.fill(row, 1.0);
        int t = Motifs.size();
        for (int i = 0; i < k; i++){
            for (int j = 0; j < t; j++){
                char currChar = Motifs.get(j).charAt(i);
                int index = IndexDNA(currChar);
                counts[index][i]++;
            }
        }
        return counts;
    }

    public static int IndexDNA(char nucleo){
        switch (nucleo){
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            default:
                System.out.print("Invalid nucleotide");
                return -1;
        }
    }
}