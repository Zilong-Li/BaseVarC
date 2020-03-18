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
    std::string GetTargetBase(const std::string& rg, const std::string& f);

};

std::string RefReader::GetTargetBase(const std::string& rg, const std::string& f)
{
    std::string chr;
    std::string seq;
    int32_t rg_s, rg_e;
    std::tie(chr, rg_s, rg_e) = BaseVarC::splitrg(rg);
    if (!LoadIndex(f)) {
        throw std::runtime_error("ERROR: reference must be index with samtools faidx");
    }
    // SeqLib will throw exception if something goes wrong.
    seq = QueryRegion(chr, rg_s - 1, rg_e - 1);    // make 0-based
    for (auto & i: seq) {    // may contain '-' character
        if ((i >= 65 && i <= 90) || i == '-') continue;
        i = i ^ 0x20;
    }
    return seq;
}

#endif