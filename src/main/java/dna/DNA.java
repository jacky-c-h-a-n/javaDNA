package dna;

public class DNA
{
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

    //constructor - save cleared sequence
    public DNA(String sequence)
    {
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
        if(init!="ATG")
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
        for(int i=0;i<genes.length();i++)
        {
            String s = genes.substring(i,i+1);
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
}
// TODO: Implement the DNA datatype from scratch!
// Use the test cases to guide you.
