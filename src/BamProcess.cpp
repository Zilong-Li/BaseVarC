#include "BamProcess.h"

bool BamProcess::FindSnpAtPos(const std::string& rg, const std::vector<int32_t>& pv)
{
    // check if the BAM is sorted
    std::string hh = Header().AsString(); //std::string(header()->text)
    bool sorted = hh.find("SO:coord") != std::string::npos;
    if (!sorted) {
    	std::cerr << "ERROR: BAM file does not appear to be sorted (no SO:coordinate) found in header." << std::endl;
    	std::cerr << "       Sorted BAMs are required." << std::endl;
    	exit(EXIT_FAILURE);
    }
    // find sm:samplename
    size_t p;
    if ((p = hh.find("SM:")) != std::string::npos) {
        hh.erase(0, p+3);
        if ((p = hh.find("\t")) != std::string::npos) sm = hh.substr(0, p);
    } else {
        std::cerr << "ERROR: No SM tag can be found. Please make sure there is SM tag in the bam header" << std::endl;
    	exit(EXIT_FAILURE);
    }
    SeqLib::GenomicRegion gr(rg, Header());
    if (!SetRegion(gr)) {
        // SetRegion always be true, which is weird
    	std::cerr << sm << ": region " << rg << " is empty." << std::endl;
        return false;
    }
    gr.Pad(1000);
    SeqLib::BamRecord r;
    // filter reads here
    SeqLib::BamRecordVector rv;
    while (GetNextRecord(r)) {
        if (r.MapQuality() < mapq) continue;
        rv.push_back(r);
    }
    // check again;
    if (rv.empty()) return false;
    bool flag = false;
    uint32_t i = 0, j = 0;
    r = rv[i];
    const std::string SKIP = "DPN";
    AlleleInfo ale;
    for (auto const& pos: pv) {
        /* select the first record covering this position */
        if (pos < r.Position() + 1) { // make 1-based
            continue;
        }else{
            // r.PositionEnd is 1-based
            while (pos > r.PositionEnd()) {
                if (i == rv.size() - 1) {flag = true; break;}
                j = ++i;
                r = rv[i];
            }
            // find the proper read here
            while (true) {
                SeqLib::Cigar c = r.GetCigar();
                int track = r.Position();
                bool fail = false;
                for (auto const& cf: c) {
                    auto t = cf.Type();
                    if (t == 'H' || t == 'S') continue;
                    if (t != 'I') track += cf.Length();
                    // skip read if position locus at a SKIP cigar field;
                    if (track >= pos) {
                        if (SKIP.find(t) != std::string::npos) fail = true;
                        break;
                    }
                }
                if (fail && j < rv.size() - 1) {
                    r = rv[++j];
                    if (pos < r.Position() + 1) break;
                } else {
                    break;
                }
            }
            // now we pick this read.
            if (pos < r.Position() + 1 || flag) {
                ;    //continue is a bug; should do nothing
            } else {
                GetAllele(r, pos, ale);
                allele_m.insert({pos, ale});     // key might change to 'chr:pos'
            }
            r = rv[i];     // all back to index i
            j = i;
        }
    }
    return true;
}

void BamProcess::FindSnpAtPos(const std::string& rg, const PosInfoVector& pv)
{
    // check if the BAM is sorted
    std::string hh = Header().AsString(); //std::string(header()->text)
    bool sorted = hh.find("SO:coord") != std::string::npos;
    if (!sorted) {
    	std::cerr << "ERROR: BAM file does not appear to be sorted (no SO:coordinate) found in header." << std::endl;
    	std::cerr << "       Sorted BAMs are required." << std::endl;
    	exit(EXIT_FAILURE);
    }
    SeqLib::GenomicRegion gr(rg, Header());
    SetRegion(gr);
    gr.Pad(1000);
    SeqLib::BamRecord r;
    // filter reads here
    SeqLib::BamRecordVector rv;
    while (GetNextRecord(r)) {
        if (r.MapQuality() < mapq) continue;
        rv.push_back(r);
    }
    if (rv.empty()) {
        for (size_t i = 0; i < pv.size(); ++i) {
            snps.push_back('.');
        }
    } else {
        bool flag = false;
        uint32_t  i = 0, j = 0;
        r = rv[i];
        snps.reserve(pv.size());       // best practice;
        const std::string SKIP = "DPN";
        for (auto const& s: pv) {
            /* select the first record covering this position */
            if (s.pos < r.Position() + 1) { // make 1-based
                snps.push_back('.');
            }else{
                // r.PositionEnd is 1-based
                while (s.pos > r.PositionEnd()) {
                    if (i == rv.size() - 1) {flag = true; break;}
                    j = ++i;
                    r = rv[i];
                }
                // find the proper read here
                while (true) {
                    SeqLib::Cigar c = r.GetCigar();
                    int track = r.Position();
                    bool fail = false;
                    for (auto const& cf: c) {
                        auto t = cf.Type();
                        if (t == 'H' || t == 'S') continue;
                        if (t != 'I') track += cf.Length();
                        // skip read if position locus at a SKIP cigar field;
                        if (track >= s.pos) {
                            if (SKIP.find(t) != std::string::npos) fail = true;
                            break;
                        }
                    }
                    if (fail && j < rv.size() - 1) {
                        r = rv[++j];
                        if (s.pos < r.Position() + 1) break;
                    } else {
                        break;
                    }
                }
                // now we pick this read.
                if (s.pos < r.Position() + 1 || flag) {
                    snps.push_back('.');
                } else {
                    snps.push_back(GetSnpCode(r, s));
                }
                r = rv[i];     // all back to index i
                j = i;
            }
        }
    }
}

char BamProcess::GetSnpCode(const SeqLib::BamRecord& r, const PosInfo& s) const
{
    int offset = GetOffset(r, s.pos);
    char seq[r.Sequence().length() + 1];
    std::strcpy(seq, r.Sequence().c_str());    // must copy r.r.Sequence();
    char x = seq[offset];
    if (x == s.ref) {
        return '0';
    } else if (x == s.alt) {
        return '1';
    } else {
        return '.';
    }
}

void BamProcess::GetAllele(const SeqLib::BamRecord& r, const uint32_t pos, AlleleInfo& ale) const
{
    int offset = GetOffset(r, pos);
    char seq[r.Sequence().length() + 1];
    char qualities[r.Qualities().length() + 1];
    std::strcpy(seq, r.Sequence().c_str());    // must copy r.r.Sequence();
    std::strcpy(qualities, r.Qualities().c_str());
    // assign allele info
    ale.base = base_m.at(seq[offset]);
    // offset = 33 , so need to substract 33;
    ale.qual = qualities[offset] - 33;
    ale.mapq = r.MapQuality();
    ale.rpr = offset + 1;
    if (r.ReverseFlag()) {
        ale.strand = 0;
    } else {
        ale.strand = 1;
    }
}

int BamProcess::GetOffset(const SeqLib::BamRecord& r, const uint32_t pos) const
{
    SeqLib::Cigar c = r.GetCigar();
    // offset is the 0-based
    int offset = pos - (r.Position() + 1);
    uint32_t track = r.Position();
    for (auto const& cf: c) {
        auto t = cf.Type();
        if (t != 'I' && t != 'S' && t != 'H') track += cf.Length();
        if (track < pos) {
            switch (cf.Type()) {
            case 'I': offset += cf.Length(); break;
            case 'S': offset += cf.Length(); break;
            case 'D': offset -= cf.Length(); break;
            case 'P': offset -= cf.Length(); break;
            case 'N': offset -= cf.Length(); break;
            default : break;
            }
        } else {
            break;
        }
    }
    return offset;
}