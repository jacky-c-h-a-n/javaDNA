package dna;

import java.util.*;

public class DNA
{
    //original sequence from construtor
    private String sequen;

    private ArrayList<Integer> initIndex;
    private ArrayList<Integer> midIndex;
    private ArrayList<Integer> finaIndex;

    //from clear method in condtructor
    private String genes;
    private int numberJunk;

    //from totalMass method
    private int numberA;
    private int numberC;
    private int numberG;
    private int numberT;
    private double totalMass;

    //mass data
    private final double A = 135.128;
    private final double C = 111.103;
    private final double G = 151.128;
    private final double T = 125.107;
    private final double junk = 100.000;

    //constructor - save original sequence and cleared sequence
    public DNA(String sequence)
    {
        sequen = sequence;
        String save = clearG(sequence);
        if(save.length()%3!=0)
        {
            throw new IllegalArgumentException("Invalid DNA sequence");
        }
        genes = save;
    }

    //method to clear a squence, save number of junk
    public String clearG(String Ogene)
    {
        String re = "";
        int numJ =0;
        for(int i =0; i< Ogene.length();i++)
        {
            String s = Ogene.substring(i,i+1);
            if(s.equals("A")||s.equals("C")||s.equals("T")||s.equals("G"))
            {
                re += s;
            }
            else
            {
                numJ++;
            }
        }
        numberJunk = numJ;
        return re;
    }

    //identify if the sequence is a protein
    public boolean isProtein()
    {
        //at least five codons
        if(genes.length()<3*5)
        {
            return false;
        }

        //first codon ATG
        String init = genes.substring(0,3);
        if(!init.equals("ATG"))
        {
            return false;
        }

        //last codon TAA or TAG or TGA
        int leng = genes.length();
        String fin = genes.substring(leng-3);
        if(!fin.equals("TAA")&&!fin.equals("TGA")&&!fin.equals("TAG"))
        {
            return false;
        }

        //at least 30% C and G of total mass
        totalMass = totalMass();
        double CandG = numberC*C + numberG*G;
        if(CandG/totalMass<0.3)
        {
            return false;
        }
        return true;
    }

    //find total mass - save number of A, C, G, T
    public double totalMass()
    {
        int numA = 0;
        int numC = 0;
        int numG = 0;
        int numT = 0;
        double re = 0.0;
        for(int i=0;i<sequen.length();i++)
        {
            String s = sequen.substring(i,i+1);
            if(s.equals("A"))
            {
                numA++;
            }
            if(s.equals("C"))
            {
                numC++;
            }
            if(s.equals("G"))
            {
                numG++;
            }
            if(s.equals("T"))
            {
                numT++;
            }
        }
        numberA = numA;
        numberC = numC;
        numberG = numG;
        numberT = numT;
        clearG(sequen);
        re = Math.round((numA*A + numC*C + numG*G + numT*T + numberJunk*junk)*10.0)/10.0;
        return re;
    }

    //find number of certain nucleotide
    public int nucleotideCount(char n)
    {
        totalMass();
        if(n=='A')
        {
            return numberA;
        }
        if(n=='C')
        {
            return numberC;
        }
        if(n=='G')
        {
            return numberG;
        }
        if(n=='T')
        {
            return numberT;
        }
        return 0;
    }

    //create a hashset containing distinct codons, hashset doesn't save same elements twice
    public HashSet<String> codonSet()
    {
        HashSet<String> re = new HashSet<String>();
        for(int i =0; i<genes.length();i+=3)
        {
            String a = genes.substring(i,i+3);
            re.add(a);
        }
        return re;
    }

    //locate the indexes of the target codon in the original sequence
    public void findCodon(String targetCodon)
    {
        ArrayList<Integer> init = new ArrayList<Integer>();
        ArrayList<Integer> mid = new ArrayList<Integer>();
        ArrayList<Integer> fina = new ArrayList<Integer>();

        for(int i =0;i<sequen.length();i++)
        {
            int index =0;
            String a = sequen.substring(i,i+1);
            String b = targetCodon.substring(index,index+1);
            if(a.equals(b))
            {
                index++;
                init.add(i);
                i++;
                while(i<sequen.length()&&index<=2)
                {
                    String c = sequen.substring(i,i+1);
                    String d = targetCodon.substring(index,index+1);
                    if(c.equals(d))
                    {
                        if(index==1)
                        {
                            mid.add(i);
                        }
                        if(index==2)
                        {
                            fina.add(i);
                        }
                        index++;
                        i++;
                    }
                    else
                    {
                    i++;
                    }
                }
            }
        }
        initIndex = init;
        midIndex = mid;
        finaIndex = fina;
    }

    public String replaceCodon(String originalCodon, String newCodon)
    {
        clearG(originalCodon);
        if(numberJunk!=0||originalCodon.length()%3!=0)
        {
            return sequen;
        }
        clearG(newCodon);
        if(numberJunk!=0||newCodon.length()%3!=0)
        {
            return sequen;
        }

        //replace codon for sequence without removing junks
        int in = 0;
        findCodon(originalCodon);
        while(in<initIndex.size())
        {
        int index =0;
        sequen = sequen.substring(0,initIndex.get(in)) + newCodon.substring(index,index+1) + sequen.substring(initIndex.get(in)+1);
        index++;
        sequen = sequen.substring(0,midIndex.get(in)) + newCodon.substring(index,index+1) + sequen.substring(midIndex.get(in)+1);
        index++;
        sequen = sequen.substring(0,finaIndex.get(in)) + newCodon.substring(index,index+1) + sequen.substring(finaIndex.get(in)+1);
        in++;
        }
        return sequen;
    }

    //change a certain combination of codon in a sequence and remove junks
    public String mutateCodon(String originalCodon, String newCodon)
    {
        clearG(originalCodon);
        if(numberJunk!=0||originalCodon.length()%3!=0)
        {
            return sequen;
        }
        clearG(newCodon);
        if(numberJunk!=0||newCodon.length()%3!=0)
        {
            return sequen;
        }

        //replace codon for clear sequence
        for(int i =0; i<genes.length();i+=3)
        {
            String a = genes.substring(i,i+3);
            if(originalCodon.equals(a))
            {
                genes = genes.substring(0,i) + newCodon + genes.substring(i+3);
            }
        }
        sequen = genes;
        return genes;
    }

    //return sequence
    public String sequence()
    {
        return sequen;
    }

    public static void main(String[] args) {
        DNA o = new DNA("---A****G----C---A;;G;;C");
        o.mutateCodon("AGC", "TGA");
        System.out.println(o.sequence());
        DNA a = new DNA("AAAGGTTACTG+A");
        String c = a.mutateCodon("TGA","G+T");
        System.out.println(c);
        System.out.println(a.totalMass());
        HashSet<String> d = a.codonSet();
        for(String e: d)
        {
            System.out.println(e);
        }

    }
}
// TODO: Implement the DNA datatype from scratch!
// Use the test cases to guide you.
