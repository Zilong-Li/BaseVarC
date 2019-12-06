#ifndef __BASEVARC_BAM_PROCESS_H__
#define __BASEVARC_BAM_PROCESS_H__

#include "SeqLib/BamReader.h"
#include "BaseVarUtils.h"

struct PosInfo
{
    std::string chr;
    int pos;
    char ref;
    char alt;
    /* provide an overload of operator + to return samtools-like region */
    std::string operator+(const PosInfo& snp) const {
        assert(this->chr == snp.chr);
        const char *colon = ":";
        const char *hyphen = "-";
        auto rg_s = BaseVar::tostring(this->pos);
        auto rg_e = BaseVar::tostring(snp.pos);
        return snp.chr + colon + rg_s + hyphen + rg_e;
    }
    friend std::istream& operator>>(std::istream& is, PosInfo& info) {
        is >> std::ws >> info.chr >> info.pos >> info.ref >> info.alt;
        return is;
    }
};
typedef std::vector<PosInfo> PosInfoVector;

struct AlleleInfo
{
    unsigned int base:  3;      // 0 : A, 1 : C, 2 : G, 3 : T, 4 : N, 5 : .
    unsigned int mapq:  8;
    unsigned int qual:  8;
    unsigned int rpr:   8;
    unsigned int strand:1;      // 0 : -, 1 : +
    unsigned int is_indel:1;      // 0 : -, 1 : +
    std::string indel;
};
typedef std::vector<AlleleInfo> AlleleInfoVector;
typedef std::unordered_map<int32_t, AlleleInfo> PosAlleleMap;

static const std::unordered_map<char, int8_t> BASEM{ {'A', 0},{'C', 1},{'G', 2},{'T', 3} };

class BamProcess: public SeqLib::BamReader
{
 public:
    BamProcess(uint8_t mapq_): mapq(mapq_){}
    ~BamProcess(){
        // if we want to read the same bams many times , we need to reset m_bams, which is _BamMap type
        m_bams.clear();
    }

    bool FindSnpAtPos(int32_t rg_s, const std::string& refseq, const std::string& rg, const std::vector<int32_t>& pv);

    void FindSnpAtPos(const std::string& rg, const PosInfoVector& pv);

    std::string sm;
    std::vector<char> snps;
    PosAlleleMap allele_m;  // allele_m may be uninitialized.

 private:

    const uint8_t mapq;

    int GetOffset(const SeqLib::BamRecord& r, const uint32_t pos) const;

    char GetSnpCode(const SeqLib::BamRecord& r, const PosInfo& s) const;

    void GetAllele(const SeqLib::BamRecord& r, const uint32_t pos, AlleleInfo& ale) const;

    bool GetBRV(const std::string& rg, SeqLib::BamRecordVector& rv);

};


#endif