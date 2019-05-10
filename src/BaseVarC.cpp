#include <getopt.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <ctime>

#include "htslib/bgzf.h"
#include "RefReader.h"
#include "BamProcess.h"

static const char* BASEVARC_USAGE_MESSAGE = 
"Program: BaseVarC\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC <command> [options]\n\n"
"Commands:\n"
"           basetype       Variants Caller\n"
"           popmatrix      Create population matrix at specific positions.\n" 
"           merge          Merge popmatrix.\n" 
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* BASETYPE_MESSAGE = 
"Program: BaseVarC basetype\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC basetype [options]\n\n"
"Commands:\n"
"  --input,      -l        BAM/CRAM files list, one file per row.\n"
"  --output,     -o        Output path(default stdout)\n"
"  --reference,  -r        Reference file\n"
"  --region,     -g        Samtools-like region\n"
"  --mapq,       -q <INT>  Mapping quality >= INT. [10]\n"
"  --verbose,    -v        Set verbose output\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* POPMATRIX_MESSAGE = 
"Program: BaseVarC popmatrix\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC popmatrix [options]\n\n"
"Commands:\n"
"  --input,      -l        BAM/CRAM files list, one file per row.\n"
"  --posfile,    -p        Position file <CHRID POS REF ALT>\n"
"  --output,     -o        Output path(default stdout)\n"
"  --mapq,       -q <INT>  Mapping quality >= INT. [10]\n"
"  --verbose,    -v        Set verbose output\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* MERGE_MESSAGE =
"Program: BaseVarC popmatrix\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC merge [options]\n\n"
"Commands:\n"
"  --input,      -l       list of matrix files for merge, one file per row.\n"
"  --output,     -o       output path (will be added suffix .gz at the end)\n"
"\nReport bugs to lizilong@bgi.com \n\n";

void runBaseType(int argc, char **argv);
void runPopMatrix(int argc, char **argv);
void runMerge(int argc, char **argv);
void parseOptions(int argc, char **argv, const char* msg);

namespace opt {
    static bool verbose = false;
    static int mapq;
    static std::string input;
    static std::string reference;
    static std::string posfile;
    static std::string region;
    static std::string output;
}

static const char* shortopts = "hvl:r:p:g:o:q:";

static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "input",                   required_argument, NULL, 'l' },
  { "reference",               required_argument, NULL, 'r' },
  { "posfile",                 required_argument, NULL, 'p' },
  { "region",                  required_argument, NULL, 'g' },
  { "output",                  required_argument, NULL, 'o' },
  { "mapq",                    required_argument, NULL, 'q' },
  { NULL, 0, NULL, 0 }
};

int main(int argc, char** argv)
{
    if (argc <= 1 ) {
        std::cerr << BASEVARC_USAGE_MESSAGE;
        return 0;
    } else {
        std::string command(argv[1]);
        if (command == "help" || command == "--help") {
            std::cerr << BASEVARC_USAGE_MESSAGE;
            return 0;
        } else if (command == "basetype") {
            runBaseType(argc - 1, argv + 1);
        } else if (command == "popmatrix") {
            runPopMatrix(argc - 1, argv + 1);
        } else if (command == "merge") {
            runMerge(argc - 1, argv + 1);
        } else {
            std::cerr << BASEVARC_USAGE_MESSAGE;
            return 0;
        }
    }

    return 0;
}

void runMerge(int argc, char **argv)
{
    parseOptions(argc, argv, MERGE_MESSAGE);
    std::string fm = opt::input;
    std::string fo = opt::output;
    clock_t ctb = clock();
    std::vector<std::string> fm_v;
    std::string path;
    std::ifstream ifm(fm);
    while (ifm >> path) fm_v.push_back(path);
    long int i, k, m, count, n = 0, mt = 0;
    long int nm = fm_v.size();
    BGZF* fp = NULL;
    kstring_t ks;
    ks.s = 0; ks.l = 0; ks.m = 0;
    std::string ss;
    std::vector<std::string> Drg;
    std::vector<std::vector<std::string>> D;
    D.reserve(nm);
    for (k = 0; k < nm ; ++k) {
        std::cout << "reading file : " << fm_v[k] << std::endl;
        fp = bgzf_open(fm_v[k].c_str(), "r");
        if (bgzf_getline(fp, '\n', &ks) >= 0) {
            ss = ks.s;
            std::vector<std::string> tokens;
            std::string token;
            std::istringstream ts(ss);
            while (std::getline(ts, token, '\t')) tokens.push_back(token);
            if (n > 0 && n != std::stol(tokens[0])) {
                std::cerr << "error: the number of samples among the inputs are different!" << std::endl;
                exit(EXIT_FAILURE);
            } else {
                n = std::stol(tokens[0]);
                m = std::stol(tokens[1]);
            }
        }
        mt += m;
        count = 0;
        Drg.reserve(n);
        while (bgzf_getline(fp, '\n', &ks) >= 0) {
            ++count;
            Drg.push_back(ks.s);
        }
        if (count != n) {
            std::cerr << "error: the number of samples dont match the header!" << std::endl;
            exit(EXIT_FAILURE);
        }
        Drg.shrink_to_fit();
        D.push_back(Drg);
        Drg.clear();
        if (bgzf_close(fp) < 0) {
            std::cerr << "warning: file cannot be closed!" << std::endl;
        }
    }
    D.shrink_to_fit();

    fp = bgzf_open(fo.c_str(), "w");
    ss = std::to_string(n) + "\t" + std::to_string(mt) + "\n";
    if (bgzf_write(fp, ss.c_str(), ss.length()) != ss.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    // fast concatenation consider using rope data structure
    // @see https://brianbondy.com/blog/tagged/data-structure
    for (i = 0; i < n; ++i){
        ss = "";
        ss.reserve(mt+1);
        for (k = 0; k < nm; ++k) {
            ss += D[k][i];
        }
        ss += "\n";
        if (bgzf_write(fp, ss.c_str(), ss.length()) != ss.length()) {
            std::cerr << "fail to write" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    if (bgzf_close(fp) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cout << "elapsed secs : " << elapsed_secs << std::endl;
    std::cout << "merge done" << std::endl;
}

void runBaseType(int argc, char **argv)
{
    parseOptions(argc, argv, BASETYPE_MESSAGE);
    std::cerr << "basetype start" << std::endl;
    std::ifstream ibam(opt::input);
    std::vector<std::string> bams(std::istream_iterator<Line>{ibam},
    	                          std::istream_iterator<Line>{});
    RefReader fa;
    fa.GetTargetBase(opt::region, opt::reference);
    std::string chr;
    int32_t rg_s, rg_e;
    std::tie(chr, rg_s, rg_e) = BaseVar::splitrg(opt::region);
    std::vector<int32_t> pv;
    for (int i = 0; i < fa.seq.length(); i++) {
        if (fa.seq[i] == 'N') continue;
        pv.push_back(i + rg_s);       // 1-based
    }
    std::vector<PosAlleleMap> allele_mv;
    const uint32_t N = bams.size();
    int count = 0;
    for (int i = 0; i < N; i++) {
        BamProcess reader;
        if (!(++count % 1000)) std::cerr << "Processing the number " << count / 1000 << "k bam" << std::endl;
        if (!reader.Open(bams[i])) {
            std::cerr << "ERROR: could not open file " << bams[i] << std::endl;
            exit(EXIT_FAILURE);
        }
        reader.FindSnpAtPos(opt::region, pv);
        allele_mv.push_back(reader.allele_m);
        if (!reader.Close()) {
            std::cerr << "Warning: could not close file " << bams[i] << std::endl;
        }
    }
    std::stringstream ss;
    for (auto const& p: pv) {
        std::vector<AlleleInfo> aiv;
    	for (auto& m: allele_mv) {
    	    if (m.count(p) == 0) {
                continue;
    	    } else {
                ss << m[p].base;
                // skip N base
                if (m[p].base != 4) aiv.push_back(m[p]);
    	    }
    	}
    	ss << "\n";
        // skip coverage==0
        int j = aiv.size();
        if (j > 0) {
            // here call BaseType
            std::cout << p << "\t" << j << std::endl;
        }
    }
    std::string out = ss.str();
    BGZF* fp = bgzf_open(opt::output.c_str(), "w");
    if (bgzf_write(fp, out.c_str(), out.length()) != out.length()) {
    	std::cerr << "failed to write" << std::endl;
    	exit(EXIT_FAILURE);
    }
    if (bgzf_close(fp) < 0) std::cerr << "failed to close \n";
    std::cerr << "basetype done" << std::endl;
}

void runPopMatrix (int argc, char **argv)
{
    parseOptions(argc, argv, POPMATRIX_MESSAGE);
    std::cerr << "popmatrix start" << std::endl;
    clock_t ctb = clock();
    std::ifstream ibam(opt::input);
    std::ifstream ipos(opt::posfile);
    if (!ibam.is_open() || !ipos.is_open()) {
        std::cerr << "bamlist or posifle can not be opend" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::vector<std::string> bams(std::istream_iterator<Line>{ibam},
    	                          std::istream_iterator<Line>{});
    PosInfoVector pv;
    for (PosInfo p; ipos >> p;) pv.push_back(p);
    pv.shrink_to_fit();        // request for the excess capacity to be released

    const long int N = bams.size();
    const long int M = pv.size();
    std::string out = std::to_string(N) + "\t" + std::to_string(M) + "\n";
    BGZF* fp = bgzf_open(opt::output.c_str(), "w");
    if (bgzf_write(fp, out.c_str(), out.length()) != out.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    // ready for run
    std::string rg = pv.front() + pv.back();
    long int count = 0;
    for (int i = 0; i < N; i++) {
        BamProcess reader;
        if (!(++count % 1000)) std::cerr << "Processing the number " << count / 1000 << "k bam" << std::endl;
        if (!reader.Open(bams[i])) {
            std::cerr << "ERROR: could not open file " << bams[i] << std::endl;
            exit(EXIT_FAILURE);
        }
        reader.FindSnpAtPos(rg, pv);
        std::string out(reader.snps.begin(), reader.snps.end());
        out += "\n";
        if (bgzf_write(fp, out.c_str(), out.length()) != out.length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!reader.Close()) {
            std::cerr << "Warning: could not close file " << bams[i] << std::endl;
        }
    }
    if (bgzf_close(fp) < 0) std::cerr << "failed to close \n";
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cout << "elapsed secs : " << elapsed_secs << std::endl;
    std::cerr << "popmatrix done" << std::endl;
}

void parseOptions(int argc, char **argv, const char* msg)
{
    bool die = false;
    bool help = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case 'v': opt::verbose = true; break;
        case 'q': arg >> opt::mapq; break;
        case 'l': arg >> opt::input; break;
        case 'r': arg >> opt::reference; break;
        case 'p': arg >> opt::posfile; break;
        case 'g': arg >> opt::region; break;
        case 'o': arg >> opt::output; break;
        default: die = true;
        }
    }
    // todo : need more check
    if (die || help || (opt::input.empty() && opt::output.empty())) {
        std::cerr << msg;
        if (die) exit(EXIT_FAILURE);
        else exit(EXIT_SUCCESS);
    }
}