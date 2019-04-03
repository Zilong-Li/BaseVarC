#include <getopt.h>
#include <string>
#include <iterator>
#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>

#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"
#include "BamProcess.h"

static const char* BASEVARC_USAGE_MESSAGE = 
"Program: BaseVarC\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage: BaseVarC <command> [options]\n\n"
"Commands:\n"
"           basetype       Variants Caller\n"
"           popmatrix      Create population matrix at specific positions.\n" 
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* BASETYPE_MESSAGE = 
"Program: BaseVarC basetype\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage: BaseVarC basetype [options]\n\n";

static const char* POPMATRIX_MESSAGE = 
"Program: BaseVarC popmatrix\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage: BaseVarC popmatrix [options]\n\n"
"Commands:\n"
"  --bamlist,    -l        BAM/CRAM files list, one file per row.\n"
"  --posfile,    -p        Position file <CHRID POS REF ALT>\n"
"  --output,     -o        Output path\n"
"  --mapq,       -q <INT>  Mapping quality >= INT. [10]\n"
"  --verbose,    -v        Set verbose output\n"
"\nReport bugs to lizilong@bgi.com \n\n";

//void runBaseType(int argc, char **argv);
void runPopMatrix(int argc, char **argv);
void subPopMatrix (const std::vector<std::string>& , const PosInfoVector&);
void writeOut (const char* , const PosInfoVector& , const int32_t& , const int32_t&);
void parseOptions(int argc, char **argv, const char* msg);

namespace opt {
    static bool verbose = false;
    static int mapq;
    static std::string bamlst;
    static std::string posfile;
    static std::string output;
}

static const char* shortopts = "hv:q:l:p:o";

static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "mapq",                    no_argument, NULL, 'q' },
  { "bamlst",                 required_argument, NULL, 'l' },
  { "posfile",                 required_argument, NULL, 'p' },
  { "output",                  required_argument, NULL, 'o' },
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
	    std::cerr << BASETYPE_MESSAGE;
	    return 0;
	} else if (command == "popmatrix") {
	    runPopMatrix(argc - 1, argv + 1);
	} else {
	    std::cerr << BASEVARC_USAGE_MESSAGE;
	    return 0;
	}
    }

    return 0;
}

void runPopMatrix (int argc, char **argv) {
    parseOptions(argc, argv, POPMATRIX_MESSAGE);
    std::ifstream ibam(opt::bamlst);
    std::ifstream ipos(opt::posfile);
    if (!ibam.is_open() || !ipos.is_open())
      std::cerr << "file can not be opend";
    std::vector<std::string> bams(std::istream_iterator<Line>{ibam},
    	                           std::istream_iterator<Line>{});
    PosInfoVector pv;
    for (PosInfo p; ipos >> p;) pv.push_back(p);
    pv.shrink_to_fit();        // request for the excess capacity to be released
    
    subPopMatrix(bams, pv);
}

void subPopMatrix (const std::vector<std::string>& bams , const PosInfoVector& pv) {
    const int32_t M = bams.size();
    const int32_t N = pv.size();
    char out[N * M];
    std::string rg = pv.front() + pv.back();

    for (int32_t i = 0; i < M; ++i) {
	BamProcess reader;
	if (!reader.Open(bams[i])) {
	    std::cerr << "ERROR: could not open file " << bams[i] << std::endl;
	    exit(EXIT_FAILURE);
	}
	SeqLib::GenomicRegion gr(rg, reader.Header());
	reader.SetRegion(gr);
	reader.findSnpAtPos(pv);
	for (int32_t j = 0; j < N; ++j) {
	    out[j * M + i] = reader.snps[j];
	}
    }

    writeOut(out, pv, N, M);
    
}

void writeOut (const char* out, const PosInfoVector& pv, const int32_t& N, const int32_t& M) {
    char tmp;
    for (int32_t i = 0 ; i < N ; ++i) {
	std::cout << pv[i].chr << " " << pv[i].pos << " " << pv[i].ref << " " << pv[i].alt << " ";
	for (int32_t j = 0; j < M; ++j) {
	    if (out[i * M + j] == pv[i].ref) {
		tmp = '0';
	    }else if (out[i * M + j] == pv[i].alt) {
		tmp = '1';
	    }else{
		tmp = '.';
	    }
	    std::cout << tmp << " ";
	}
	std::cout << std::endl;
    }
}

void parseOptions(int argc, char **argv, const char* msg) {
    bool die = false;
    bool help = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
	std::istringstream arg(optarg != NULL ? optarg : "");
	switch (c) {
	case 'v': opt::verbose = true; break;
	case 'q': arg >> opt::mapq; break;
	case 'l': arg >> opt::bamlst; break;
	case 'p': arg >> opt::posfile; break;
	case 'o': arg >> opt::output; break;
	default: die = true;
	}
    }
    if (die || help || (opt::bamlst.empty() && opt::posfile.empty())) {
	std::cerr << msg;
	if (die)
	  exit(EXIT_FAILURE);
	else
	  exit(EXIT_SUCCESS);
    }
}