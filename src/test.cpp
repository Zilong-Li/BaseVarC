#include <string>
#include <iterator>
#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>

#include "SeqLib/BamReader.h"
#include "SeqLib/ReadFilter.h"
#include "BamProcess.h"

void subpopmatrix (const std::vector<std::string>& bams , const PosInfoVector& pv) {
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

    for (int32_t i = 0 ; i < N ; ++i) {
	std::cout << pv[i].chr << "\t" << pv[i].pos << "\t" << pv[i].ref << "\t" << pv[i].alt << "\t";
	for (int32_t j = 0; j < M; ++j) {
	    std::cout << out[i * M + j] << ",";
	}
	std::cout << std::endl;
    }
    
}

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
    subpopmatrix(files, pv);
    // std::string rg = pv.front() + pv.back();

    // const int32_t M = files.size();
    // const int32_t N = pv.size();
    // std::vector<char> out;
    // out.reserve(N*M);
    // for (int32_t i = 0; i < M; ++i) {
    // 	BamProcess reader;
    // 	if (!reader.Open(files[i])) {
    // 	    std::cerr << "ERROR: could not open file " << files[i] << std::endl;
    // 	    exit(EXIT_FAILURE);
    // 	}
    // 	SeqLib::GenomicRegion gr(rg, reader.Header());
    // 	reader.SetRegion(gr);
    // 	reader.findSnpAtPos(pv);
    // 	for (int32_t j = 0; j < N; ++j) {
    // 	    out[j * M + i] = reader.snps[j];
    // 	}
    // }

    // for (int32_t i = 0 ; i < N ; ++i) {
    // 	std::cout << pv[i].chr << "\t" << pv[i].pos << "\t" << pv[i].ref << "\t" << pv[i].alt << "\t";
    // 	for (int32_t j = 0; j < M; ++j) {
    // 	    std::cout << out[i * M + j] << ",";
    // 	}
    // 	std::cout << std::endl;
    // }
    
}
