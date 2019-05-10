#ifndef __BASEVARC_BAM_PROCESS_H__
#define __BASEVARC_BAM_PROCESS_H__

#include "SeqLib/BamReader.h"
#include "BaseVarUtils.h"

struct Line
{
    std::string data;
    operator std::string const&() const {return data;}
    friend std::istream& operator>>(std::istream& is, Line& line) {
        return std::getline(is, line.data);
    }
};

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
    unsigned int base:  3;      // 0 : A, 1 : C, 2 : G, 3 : T
    unsigned int strand:1;      // 0 : -, 1 : +
    unsigned int qual:  8;
    unsigned int mapq:  8;
    unsigned int rpr:   8;
};
typedef std::unordered_map<uint32_t, AlleleInfo> PosAlleleMap;

class BamProcess: public SeqLib::BamReader
{
 public:
    BamProcess(){}
    ~BamProcess(){
        // if we want to read the same bams many times , we need to reset m_bams, which is _BamMap type
        m_bams.clear();
    }

    void FindSnpAtPos(const std::string& rg, const std::vector<int32_t>& pv);

    void FindSnpAtPos(const std::string& rg, const PosInfoVector& pv);

    uint8_t mapq = 10;
    std::vector<char> snps;
    PosAlleleMap allele_m;

 private:
    char GetSnpCode(const SeqLib::BamRecord& r, const PosInfo& s) const;

    void GetAllele(const SeqLib::BamRecord& r, const uint32_t pos, AlleleInfo& ale) const;

    int GetOffset(const SeqLib::BamRecord& r, const uint32_t pos) const;

    std::unordered_map<char, int> base_m{ {'A', 0},{'C', 1},{'G', 2},{'T', 3},{'N', 4} };
};


#endif