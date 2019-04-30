#ifndef __BASEVARC_REF_READER_H__
#define __BASEVARC_REF_READER_H__

#include "SeqLib/RefGenome.h"
#include "BaseVarUtils.h"

class RefReader: public SeqLib::RefGenome
{
 public:
    RefReader(){}
    ~RefReader(){}
    // rg is samtools-like region
    void GetTargetBase(std::string rg, const std::string& f);
    std::string seq;

};

void RefReader::GetTargetBase(std::string rg, const std::string& f)
{
    std::string chr;
    int32_t rg_s, rg_e;
    std::tie(chr, rg_s, rg_e) = BaseVar::splitrg(rg);
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