#include "BamProcess.h"

bool BamProcess::FindSnpAtPos(int32_t rg_s, const std::string& refseq, const std::string& rg, const std::vector<int32_t>& pv)
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
        if ((p = hh.find("\n")) != std::string::npos) {
            hh.erase(hh.begin() + p, hh.end());
            if ((p = hh.find("\t")) != std::string::npos) {
                sm = hh.substr(0, p);
            } else {
                sm = hh;
            }
        }
    } else {
        std::cerr << "ERROR: No SM tag can be found. Please make sure there is SM tag in the bam header" << std::endl;
    	exit(EXIT_FAILURE);
    }
    SeqLib::GenomicRegion gr(rg, Header());
    gr.Pad(1000);
    if (!SetRegion(gr)) {
        // SetRegion always be true, which is weird
    	std::cerr << sm << ": region " << rg << " is empty." << std::endl;
        return false;
    }
    SeqLib::BamRecord r;
    // filter reads here
    SeqLib::BamRecordVector rv;
    while (GetNextRecord(r)) {
        if (r.MapQuality() < mapq) continue;
        rv.push_back(r);
    }
    // check again;
    if (rv.empty()) return false;
 // sk: the index of the CIGAR operator.
 // sx: the reference coordinate of the start of sk
 // sy: the query coordiante of the start of sk
    size_t i = 0, j = 0, k = 0, nc, sk;
    int sx, sy, indel;
    bool eof = false, next = false, is_indel = false, is_del = false, is_refskip = false;
    std::string indel_str;
    AlleleInfo ale; // = {4, 0, 0, 0, 0, 0, "N"};
    SeqLib::Cigar c;
    r = rv[i];
    for (auto const& pos: pv) {
        // select the first record covering this position
        if (pos < r.Position() + 1) { // make 1-based
            continue;
        } else {
            // r.PositionEnd is 1-based
            while (pos > r.PositionEnd()) {
                if (i == rv.size() - 1) {eof = true; break;}
                r = rv[++i]; j = i;
                if (pos < r.Position() + 1) { next = true; break;}
            }
            if (next || eof) { eof = false; next = false; continue; }
            // find proper read here
            while (true) {
                c = r.GetCigar(); nc = c.size();
                for (k = 0, sx = r.Position(), sy = 0; k < nc; ++k) {
                    char op = c[k].Type();
                    int l = c[k].Length();
                    if (op -= 'M' || op == 'I' || op == 'S' || op == 'X') sy += l;
                    if (op == 'H' || op == 'I') continue;
                    else sx += l;
                    if (pos <= sx) break;
                }
                assert(k < nc);
                sk = k;
                // collect indel information
                is_indel = is_del = is_refskip = false; indel = 0;
                {
                    char op = c[sk].Type();
                    if (sx == pos && sk + 1 < nc) { //peek the next operation
                        char op2 = c[sk+1].Type();
                        int l2 = c[sk+1].Length();
                        if (op2 == 'D') indel = -(int)l2, indel_str = "-" + refseq.substr(pos - rg_s + 1, l2);  // start from the next operation position
                        else if (op2 == 'I') indel = l2, indel_str = "+" + r.Sequence().substr(sy, l2);
                        else if (op2 == 'P' && sk + 2 < nc) { // no working for adjacent padding
                            indel_str = 'N';   // this indicates the indel may be adjacent with a padding operation
                            int l3 = 0;
                            for (k = sk + 2; k < nc; ++k) {
                                op2 = c[sk].Type();
                                l2 = c[sk].Length();
                                if (op2 == 'I') l3 += l2;
                                else if (op2 == 'D' || op2 == 'M' || op2 == 'N' || op2 == 'X') break;
                            }
                            if (l3 > 0) indel = l3;
                        }
                    }
                    if (indel != 0) is_indel = true;
                    if (op == 'D') is_del = true;
                    else if (op == 'N') is_refskip = true;
                }
                // here we go
                if (is_indel) {
                    if (r.ReverseFlag() || r.MateReverseFlag()) ale.strand = 0;
                    else ale.strand = 1;
                    ale.base = 5; ale.qual = r.MapQuality(); ale.rpr = 0; ale.is_indel = 1; ale.indel = indel_str;
                    allele_m.insert({pos, ale});
                    break;
                }
                if (is_del == 0 && is_refskip == 0) {
                    GetAllele(r, pos, ale);
                    allele_m.insert({pos, ale});
                    break;
                } else if (j < rv.size() - 1) {
                    r = rv[++j];
                    if (pos < r.Position() + 1 || pos > r.PositionEnd()) break;
                } else {
                    break;
                }
            }
            r = rv[i]; j = i;     // all back to index i
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
    gr.Pad(1000);
    SetRegion(gr);
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
        size_t  i = 0, j = 0;
        r = rv[i];
        snps.reserve(pv.size());
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
                if (s.pos < r.Position() + 1 || s.pos > r.PositionEnd() || flag) {
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
    std::string seq = r.Sequence();
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
    std::string seq = r.Sequence();
    std::string qualities = r.Qualities();
    if (BASEM.count(seq[offset])) {
        ale.base = BASEM.at(seq[offset]);
    } else {
        ale.base = 4;
    }
    ale.qual = qualities[offset] - 33;
    ale.mapq = r.MapQuality();
    ale.rpr = offset + 1;
    ale.is_indel = 0;
    if (r.ReverseFlag() || r.MateReverseFlag()) ale.strand = 0;
    else ale.strand = 1;
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