#include "BamProcess.h"

void BamProcess::FindSnpAtPos(const SeqLib::GenomicRegion& gr, const PosInfoVector& pv) {
    SetRegion(gr);
    SeqLib::BamRecord r;
    // check if the BAM is sorted 
    std::string hh = Header().AsString(); //std::string(header()->text)
    bool sorted = hh.find("SO:coord") != std::string::npos;
    if (!sorted) {
    	std::cerr << "ERROR: BAM file does not appear to be sorted (no SO:coordinate) found in header." << std::endl;
    	std::cerr << "       Sorted BAMs are required." << std::endl;
    	exit(EXIT_FAILURE);
    }
    // filter reads here
    SeqLib::BamRecordVector rv;
    while (GetNextRecord(r)) {
	if (r.MapQuality() < mapq) continue;
	rv.push_back(r);
    } 

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

void BamProcess::PrintOut () const {
    std::string sep = "\t";
    std::ostringstream tmp;
    for (auto const& s: snps) {
	tmp << s << sep;
    }
    std::cout << tmp.str() << std::endl;
}

char BamProcess::GetSnpCode(const SeqLib::BamRecord& r, const PosInfo& s) const {
    SeqLib::Cigar c = r.GetCigar();
    uint32_t offset = s.pos - (r.Position() + 1);
    int track = r.Position();
    for (auto const& cf: c) {
	auto t = cf.Type();
	if (t != 'I' && t != 'S' && t != 'H') track += cf.Length();
	if (track < s.pos) {
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
