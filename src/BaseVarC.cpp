#include <string>
#include <iterator>
#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>

#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"
#include "BamProcess.h"

void subPopMatrix (const std::vector<std::string>& , const PosInfoVector&);

void printOut (const std::vector<char>& , const PosInfoVector& , const int32_t& , const int32_t&);

int main(int argc, char** argv)
{
    if (argc < 2) {
	exit(EXIT_FAILURE);
    }

    std::string bamlst = argv[1];
    std::ifstream ifs1(bamlst);
    std::vector<std::string> files(std::istream_iterator<Line>{ifs1},
    	                           std::istream_iterator<Line>{});

    std::string posfile = argv[2];
    std::ifstream ifs2(posfile);
    PosInfoVector pv;
    for (PosInfo info; ifs2 >> info;) {
    	pv.push_back(info);
    }
    pv.shrink_to_fit();        // request for the excess capacity to be released

    subPopMatrix(files, pv);
}

void subPopMatrix (const std::vector<std::string>& bams , const PosInfoVector& pv) {
    const int32_t M = bams.size();
    const int32_t N = pv.size();
    std::vector<char> out;
    out.reserve(N*M);
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

    printOut(out, pv, N, M);
    
}

void printOut (const std::vector<char>& out, const PosInfoVector& pv, const int32_t& N, const int32_t& M) {
    for (int32_t i = 0 ; i < N ; ++i) {
	std::cout << pv[i].chr << " " << pv[i].pos << " " << pv[i].ref << " " << pv[i].alt << " ";
	for (int32_t j = 0; j < M; ++j) {
	    std::cout << out[i * M + j] << " ";
	}
	std::cout << std::endl;
    }
}