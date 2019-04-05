#include "BamProcess.h"

void BamProcess::FindSnpAtPos(const SeqLib::GenomicRegion& gr, const PosInfoVector& pv) {

    SetRegion(gr);
    SeqLib::BamRecord r;
    /* check if the BAM is sorted */
    std::string hh = Header().AsString(); //std::string(header()->text)
    bool sorted = hh.find("SO:coord") != std::string::npos;
    if (!sorted) {
    	std::cerr << "ERROR: BAM file does not appear to be sorted (no SO:coordinate) found in header." << std::endl;
    	std::cerr << "       Sorted BAMs are required." << std::endl;
    	exit(EXIT_FAILURE);
    }
    //std::string sm = hh.find("")
    // filter reads here
    SeqLib::BamRecordVector rv;
    while (GetNextRecord(r)) {
	if (r.MapQuality() < mapq) continue;
	rv.push_back(r);
    } 

    bool flag = false;
    unsigned int i = 0;
    r = rv[i];
    snps.reserve(pv.size());       // best practice;
    const std::string SKIP = "SHPN";
    for (auto const& s: pv) {
	/* select the first record covering this position */
	if (s.pos < r.Position() + 1) { // make 1-based
	    snps.push_back('.');
	}else{
	     // r.PositionEnd is 1-based
	    while (s.pos > r.PositionEnd()) {
		if (i == rv.size() - 1) {flag = true; break;}
		i++;
		r = rv[i];
	    } 
	    // filter reads again
	    while (true) {
		SeqLib::Cigar c = r.GetCigar();
		uint32_t offset = s.pos - (r.Position() + 1);
		uint32_t idxl = 0;
		uint32_t idxr = 0;
		bool fail = false;
		for (auto const& cf: c) {
		    auto t = cf.Type();
		    if (SKIP.find(t) != std::string::npos) { fail = true; break;}
		    idxr += cf.Length();
		    if (t == 'D') {
			if (offset >= idxl && offset <= idxr) {fail = true; break;}
		    }
		    idxl += cf.Length();
		}
		if (fail && i < rv.size()) {
		    i++;
		    r = rv[i];
		} else {
		    break;
		}
	    }
	    if (s.pos < r.Position() + 1 || flag) {
		snps.push_back('.');
	    } else {
	 	snps.push_back(GetSnpCode(r, s));
	    }
	}
    }
}

void BamProcess::PrintOut () const {
    std::string sep = " ";
    std::ostringstream tmp;
    for (auto const& s: snps) {
	tmp << sep << s;
    }
    std::cout << tmp.str() << std::endl;
}

inline char BamProcess::GetSnpCode(const SeqLib::BamRecord& r, const PosInfo& s) const {
    SeqLib::Cigar c = r.GetCigar();
    uint32_t offset = s.pos - (r.Position() + 1);
    uint32_t idx = 0;
    for (auto const& cf: c) {
        switch (cf.Type()) {
        case 'I': if (offset >= idx) offset += cf.Length(); break;
        case 'D': if (offset >= idx) offset -= cf.Length(); break;
        default : break;
        }
        idx += cf.Length();
    }
    char seq[r.Sequence().length() + 1]; // must copy r.r.Sequence();
    std::strcpy(seq, r.Sequence().c_str());
    char x = seq[offset];
    if (x == s.ref) {
	return '0';
    } else if (x == s.alt) {
	return '1';
    } else {
	return '.';
    }
}
