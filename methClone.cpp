#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <cmath>
#include <fstream>
#include "gzstream.h"
#include "api/BamReader.h" 
#include "OptionParser.h"
#include <algorithm>
#include <string>

using namespace std;

bool oneSample = false;

struct parameters{
    string file1;
    string file2;
    string outfile;
    string sample;
    int m;
    int d;
    int c;
};
int getAllPatterns(std::map <std::string, int> & patterns, std::map <std::string, int> & meth)
{
    for (int i=0; i<2; i++) {
        for (int k=0; k<2; k++) {
            for (int j=0; j<2; j++) {
                for (int n=0; n<2; n++) {
                    std::ostringstream x;
                    x << i << k <<j<<n;
                    patterns[x.str()]=0;
                    meth[x.str()]=i+k+j+n;
                    //std::cout << x.str() << ":" << meth[x.str()]<< std::endl;
                }
            }
        }
    }
    return 0;
}

// take pattern vector and return all 16 patterns' frequency
int methPatterns(std::map <std::string, int> & allp, 
                 std::vector<std::string> p, 
                 std::map <std::string, int> patternMeth, 
                 double & meth,
                 std::string & maxp)
{
    int freq=p.size();
    //count the freqeuncy of different patterns.
    for (int i=0; i<freq; i++) {
        allp[p[i]]++;
    }
    int count=0;
    int f=0;
    int maxf=0;
    for (std::map <std::string, int>::const_iterator it=allp.begin(); it!=allp.end(); it++) {
        f=it->second;
        count += f * (patternMeth[it->first]);
        if (f >= maxf) {
            maxp = it->first;
            maxf = f;
        }
//        std::cout << it -> first << " "  << it->second << "*" << patternMeth[it->first] << " " << freq << std::endl;
    }
    meth=100.0 * double(count)/double(freq)/4.0;
    return 0;
}
int normalize(std::map<std::string, int> allp, std::map<std::string, double> & allpNormed, int freq)
{
    for (std::map <std::string, int>::const_iterator i=allp.begin(); i!=allp.end(); i++) {
        // normalize to 100 first
        allpNormed[i->first]=100.0*double(i->second)/double(freq);
    }
    return 0;
}

// group pattern vector by combining two pattern vectors
int groupPatterns(std::map <std::string, double> p1, std::map <std::string, double> p2, std::map <std::string, double> & p)
{
    for (std::map <std::string, double>::const_iterator i=p1.begin(); i!=p1.end(); i++) {
        p[i->first]=0.5* (i->second +p2.find(i->first)->second);
    }
    return 0;
}

// calculate entropy and normalized pattern for one sample
int entropy(std::map<std::string, double> allp,double & entropy, std::string & patterns)
{
    std::map <std::string, double> p;
    entropy=0;
    for (std::map <std::string, double>::const_iterator i=allp.begin(); i!=allp.end(); i++) {
        // entropy function
        double n=i->second;
        entropy=entropy+lgamma(n+1);
        // get normalized pattern
        std::ostringstream oss;
        oss << n << "\t";
        patterns+=oss.str();
    }
    return 0;
}

// calculate combinatorial entropy and normalized patterns
int combEntropy_one_sample(std::string loci, 
                int freq1, 
                std::vector<std::string> p1, 
                std::map <std::string, int> allpatterns,
                std::map <std::string, int> patternMeth,
                ogzstream & myfile)
{
    double meth1,meth2;
    std::string maxp1, maxp2;
    std::map <std::string, int> allp1=allpatterns;
    methPatterns(allp1, p1, patternMeth, meth1, maxp1);
    std::map <std::string, double> allpn, allpn1, allpn2;
    normalize(allp1, allpn1, freq1);
    groupPatterns(allpn1, allpn2, allpn);
    double e1;
    std::string x,x1,x2;
    entropy(allpn1, e1, x1);
    myfile << loci << "\t" <<  freq1 <<  "\t" << meth1 <<  "\t" << maxp1 << "\t" << x1 <<std::endl;
 
    return 0;
}
// iterate lociMeth1 and lociMeth2, and print loci that shared.
int interLoci_one_sample(std::map<std::string, std::vector<std::string> > lociMeth1, 
              std::map <std::string, int> allpatterns,
              std::map <std::string, int> patternMeth,
              int freq, ogzstream & myfile)
{
    std::map<std::string, std::vector<std::string> >::const_iterator it;
    for (it=lociMeth1.begin(); it!=lociMeth1.end(); it++)
    {
        std::string loci=it->first;
   	std::vector<std::string> p1=it->second;
    	int freq1=p1.size();
    if (freq1>=freq)
	   combEntropy_one_sample(loci, freq1, p1, allpatterns, patternMeth, myfile);
    
        
    }
    return 0;
}

int byread(BamTools::BamAlignment al, 
           int d, 
           std::map<std::string, std::vector<std::string> > & lociMeth, 
           const BamTools::RefVector refs, 
           std::string sample)
{
    if (al.IsMapped() && al.HasTag("XM:Z")){
        std::vector<int> loc;
        std::vector<int> m;
        std::string Meth;
        
        loc.clear(); m.clear();// make sure that the vector is clear of content.
        al.GetTag("XM:Z:", Meth);
        std::string strand="+";
        if (al.IsReverseStrand()) {
            std::reverse(Meth.begin(), Meth.end());
            strand="-";
        }
        for (int i = 0; i < al.QueryBases.size(); i++) {
            if (Meth[i] == 'z')
            {
                loc.push_back(al.Position+i+1);
                m.push_back(0);
            } 
            else if (Meth[i] == 'Z')
            {
                loc.push_back(al.Position+i+1);
                m.push_back(1);
            } 
        } 
        if (loc.size()>=4) {
            for (int i=3 ; i <loc.size(); i++) {
                int dist=loc[i]-loc[i-3];
                if (dist<= d) {
                    std::ostringstream loci;
                    loci <<refs.at(al.RefID).RefName<< "\t" << loc[i-3] << "\t" << loc[i] << "\t" << sample << "\t" << loc[i]-loc[i-3]  << "\t" << strand << "\t" << loc[i-3]<<":"<<loc[i-2]<<":"<<loc[i-1]<<":"<<loc[i];
                    std::ostringstream meth;
                    meth << m[i-3]<<m[i-2] <<m[i-1] << m[i];
                    lociMeth[loci.str()].push_back(meth.str());
                    //std::cout << loci.str() << "\t" << meth.str()<< std::endl;
                }
            }
        }
    }
    
    return 0;
}
int bamCheck(std::string bamFile, BamTools::BamReader & reader)
{
    if (!reader.Open(bamFile)){ //bam file
        cerr << "Could not open input BAM file." << std::endl;
        reader.Close();
        return false;
    } else {
        return 0;
    }
    
}

int readerToMeth(BamTools::BamReader & reader1, 
                 std::map<std::string, std::vector<std::string> > & lociMeth1, 
                 BamTools::RefVector::const_iterator i, 
                 int d, 
                 const BamTools::RefVector refs,
                 std::string sample)
{
    const int r1=reader1.GetReferenceID(i->RefName);
    //const int r2=reader2.GetReferenceID(i->RefName);
    const int rl=i->RefLength;
    if(reader1.SetRegion(r1,0, r1, rl))
    {
        //Rcerr << "Processing " << i->RefName << std::endl;
        BamTools::BamAlignment al;            
        while (reader1.GetNextAlignment(al)){
            byread(al, d, lociMeth1, refs, sample);
        }
    }
    reader1.Rewind();       
    return 0;
}
int header_one_sample(ogzstream & myfile, std::map <std::string, int> allpatterns)
{
    myfile << "chr\tstart\tend\tsample\tdistance\tstrand\tloci\tread1\tmeth1\tpattern1\t";
    
        for (std::map <std::string, int>::const_iterator it=allpatterns.begin(); it!=allpatterns.end(); it++) {
            myfile << "s0" << ":" << it->first << "\t";
        }
    
    myfile << std::endl;
    return 0;
}
int finished()
{
    time_t now = time(0);
    //Rcerr << "********************************\n" << ctime(&now) << "************Finished************" << std::endl;
    return 0;
}


//Functions for original MethCone package 



int combEntropy(std::string loci, 
                int freq1, 
                int freq2, 
                std::vector<std::string> p1, 
                std::vector<std::string> p2, 
                std::map <std::string, int> allpatterns,
                std::map <std::string, int> patternMeth,
                int methdiff,
                //double fnorm1k2,
                ogzstream & myfile)
{
    double meth1,meth2;
    std::string maxp1, maxp2;
    std::map <std::string, int> allp1=allpatterns;
    methPatterns(allp1, p1, patternMeth, meth1, maxp1);
    std::map <std::string, int> allp2=allpatterns;
    methPatterns(allp2, p2, patternMeth, meth2, maxp2);
    std::map <std::string, double> allpn, allpn1, allpn2;
    normalize(allp1, allpn1, freq1);
    normalize(allp2, allpn2, freq2);
    //nestedNormalize(allp1, freq1, allp2, freq2, fnorm1k2, allpn1, allpn2);
    groupPatterns(allpn1, allpn2, allpn);
    double e,e1,e2;
    std::string x,x1,x2;
    entropy(allpn, e, x);
    entropy(allpn1, e1, x1);
    entropy(allpn2, e2, x2);
    double combp= 2.0*e - e1 - e2;
    int methd=meth1-meth2;
    if (methd >= methdiff || methd <= - methdiff) {
        myfile << loci << "\t"<< combp << "\t" <<  freq1 << "\t" << freq2 << "\t" << meth1 << "\t" << meth2 << "\t" << maxp1 << "\t" << maxp2 << "\t" << x1 << x2 <<std::endl;
    }
    return 0;
}

int interLoci(std::map<std::string, std::vector<std::string> > lociMeth1, 
              std::map<std::string, std::vector<std::string> > lociMeth2, 
              std::map <std::string, int> allpatterns,
              std::map <std::string, int> patternMeth,
              int freq, int methdiff, ogzstream & myfile)
{
    std::map<std::string, std::vector<std::string> >::const_iterator it;
    for (it=lociMeth1.begin(); it!=lociMeth1.end(); it++)
    {
        std::string loci=it->first;
        if (lociMeth2.find(loci)!=lociMeth2.end()) {
            std::vector<std::string> p1=it->second;
            int freq1=p1.size();
            std::vector<std::string> p2=lociMeth2.find(loci)->second;
            int freq2=p2.size();
            if(freq1 >= freq & freq2 >= freq)
            {
                combEntropy(loci, freq1, freq2, p1, p2, allpatterns, patternMeth, methdiff, myfile);
            }
        }
    }
    return 0;
}

int header(ogzstream & myfile, std::map <std::string, int> allpatterns)
{
    myfile << "chr\tstart\tend\tsample\tdistance\tstrand\tloci\tentropy\tread1\tread2\tmeth1\tmeth2\tpattern1\tpattern2\t";
    for (int k=0; k<2; k++) {
        for (std::map <std::string, int>::const_iterator it=allpatterns.begin(); it!=allpatterns.end(); it++) {
            myfile << "s" << k << ":" << it->first << "\t";
        }
    }
    myfile << std::endl;
    return 0;
}



int MethClone_one_sample(std::string bamFile1,  
                  const char * outFile, 
                  std::string sample, 
                  int d, int freq)
{
    std::map <std::string, int> allpatterns;
    std::map <std::string, int> patternMeth;
    getAllPatterns(allpatterns, patternMeth);
    
    BamTools::BamReader reader1;
    bamCheck(bamFile1, reader1);

    
    // get reference name from first bam
    const BamTools::RefVector refs = reader1.GetReferenceData();
    reader1.LocateIndex();
    //reader2.LocateIndex();

    if (reader1.HasIndex())
    {

        ogzstream myfile;
        myfile.open (outFile);
        header_one_sample(myfile, allpatterns);
        
        for(BamTools::RefVector::const_iterator i = refs.begin(); i != refs.end(); ++i)
        {
            std::map<std::string, std::vector<std::string> > lociMeth1;
            readerToMeth(reader1, lociMeth1,  i, d, refs, sample);
            interLoci_one_sample(lociMeth1, allpatterns, patternMeth, freq, myfile);
        }
        reader1.Close();
        myfile.close();
        finished();
        return 0;
    } else {
        cerr << "Could not load index data for all input BAM file(s)... Aborting." << std::endl;
        return -1;
    }

}


int MethClone_two_samples(std::string bamFile1, 
                  std::string bamFile2, 
                  const char * outFile, 
                  std::string sample, 
                  int d, int freq, int methdiff)
{
    std::map <std::string, int> allpatterns;
    std::map <std::string, int> patternMeth;
    getAllPatterns(allpatterns, patternMeth);
    
    BamTools::BamReader reader1;
    BamTools::BamReader reader2;
    bamCheck(bamFile1, reader1);
    bamCheck(bamFile2, reader2);
    

    // get reference name from first bam
    const BamTools::RefVector refs = reader1.GetReferenceData();
    reader1.LocateIndex();
    reader2.LocateIndex();

    if (reader1.HasIndex() & reader2.HasIndex())
    {

        ogzstream myfile;
        myfile.open (outFile);
        header(myfile, allpatterns);
        
        for(BamTools::RefVector::const_iterator i = refs.begin(); i != refs.end(); ++i)
        {
            std::map<std::string, std::vector<std::string> > lociMeth1;
            std::map<std::string, std::vector<std::string> > lociMeth2;
            readerToMeth(reader1, lociMeth1,  i, d, refs, sample);
            readerToMeth(reader2, lociMeth2,  i, d, refs, sample);
            interLoci(lociMeth1, lociMeth2, allpatterns, patternMeth, freq, methdiff, myfile);
        }
        reader1.Close();
        reader2.Close();
        myfile.close();
        finished();
        return 0;
    } else {
        cerr << "Could not load index data for all input BAM file(s)... Aborting." << std::endl;
        return -1;
    }

}


parameters parse_options(int argc, const char * argv[]){
    const string usage = "usage: methClone [OPTION]";
    const string desc = "methclone: a program to detect the dynamic evolution of clonal epialleles in DNA methylation sequencing data. ";
    const string version = "methClone V 0.2";

    optparse::OptionParser parser = optparse::OptionParser().usage(usage).version(version).description(desc);

    parser.add_option("--f1").help("First methylation bam file for methclone").set_default("");
    parser.add_option("--f2").help("Second methylation bam file for methclone").set_default("");
    parser.add_option("-o").help("Name of the output file for methclone").set_default("");
    parser.add_option("-s").help("Sample name").set_default("sample");
    parser.add_option("-m").type("int").action("store").set_default(0).help("Cutoff for methylation difference between two samples. Default: 0");
    parser.add_option("-d").type("int").action("store").set_default(72).help("Distance cutoff to consider methylated values as one value. Default: 72");
    parser.add_option("-c").type("int").action("store").set_default(60).help("Minimum read coverage. Default: 60");

    // perform basic validations 
    optparse::Values& options = parser.parse_args(argc, argv);

    string file1 = (string) options.get("f1");
    string file2 = (string) options.get("f2");
    string outfile = (string) options.get("o");

     cout<<file1<<endl;

    if (file1.size() == 0){
        cerr<<"Please provide a methyaltion file to methclone. Exiting..."<<endl;
        exit(1);
    }

    if (file1.size() >0 && file2.size() == 0){
        cerr<<"You have provided only one sample. methClone will work in one sample mode."<<endl;
        oneSample = true;
    }

    if (outfile.size() == 0){
        cerr<<"Please provide a proper output file. Exiting..."<<endl;
        exit(1);
    }

    parameters arguments;
    arguments.file1 = file1;
    arguments.file2 = file2;
    arguments.outfile = outfile;
    arguments.sample = (string) options.get("s");
    arguments.m = (int) options.get("m");
    arguments.d = (int) options.get("d");
    arguments.c = (int) options.get("c");

    return arguments;
}


int main(int argc, const char * argv[]){
    parameters args = parse_options(argc,argv);

    if (oneSample){
        MethClone_one_sample(args.file1,  
                  args.outfile.c_str(), 
                  args.sample, 
                  args.d, args.c);
    }
    else{
         MethClone_two_samples(args.file1, 
                               args.file2, 
                               args.outfile.c_str(), 
                               args.sample, 
                                 args.d,  args.c, args.m);
    }
}