#ifndef __BASEVARC_FASTA_READER_H__
#define __BASEVARC_FASTA_READER_H__

#include "SeqLib/RefGenome.h"

class FastaReader: public SeqLib::RefGenome
{
 public:
    FastaReader(){}
    ~FastaReader(){}
    // rg is samtools-like region
    void GetTargetBase(std::string rg, const std::string& f);
    std::string seq;

};

void FastaReader::GetTargetBase(std::string rg, const std::string& f)
{
    int32_t rg_s = 0, rg_e = 0;
    std::string chr;
    size_t p;
    if ((p = rg.find(":")) != std::string::npos) {
	chr = rg.substr(0, p);
	rg.erase(0, p + 1);
    }
    if ((p = rg.find("-")) != std::string::npos) {
	rg_s = std::stoi(rg.substr(0,p)) - 1;      // make 0-based
	rg.erase(0, p + 1);
	rg_e = std::stoi(rg) - 1;
    }
	
    if (!LoadIndex(f)) {
	std::cerr << "reference must be index with samtools faidx" << std::endl;
	exit(EXIT_FAILURE);
    }
    seq = QueryRegion(chr, rg_s, rg_e);
    for (auto & i: seq) {
	if (i >= 65 && i <= 90)
	    continue;
	i = i ^ 0x20;
    }
}

#endif