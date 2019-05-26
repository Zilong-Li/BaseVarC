#include <getopt.h>
#include <iterator>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <ctime>

#include "htslib/bgzf.h"
#include "RefReader.h"
#include "BamProcess.h"
#include "BaseType.h"
#include "ThreadPool.h"

static const char* BASEVARC_USAGE_MESSAGE = 
"Program: BaseVarC -- A c version of BaseVar\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC <command> [options]\n\n"
"Commands:\n"
"           basetype       Variants Caller\n"
"           popmatrix      Create population matrix at specific positions.\n" 
"           concat         Concat popmatrix.\n" 
"           merge          Merge vcf/cvg files.\n" 
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* BASETYPE_MESSAGE = 
"Program: BaseVarC basetype\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC basetype [options]\n\n"
"Commands:\n"
"  --input,      -l        BAM/CRAM files list, one file per row.\n"
"  --output,     -o        Output filename prefix\n"
"  --reference,  -r        Reference file\n"
"  --region,     -g        Samtools-like region\n"
"  --mapq,       -q <INT>  Mapping quality >= INT. [10]\n"
"  --thread,     -t <INT>  Number of thread\n"
"  --batch,      -b <INT>  Number of samples each batch\n"
"  --verbose,    -v        Set verbose output\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* POPMATRIX_MESSAGE = 
"Program: BaseVarC popmatrix\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC popmatrix [options]\n\n"
"Commands:\n"
"  --input,      -l        BAM/CRAM files list, one file per row.\n"
"  --posfile,    -p        Position file <CHRID POS REF ALT>\n"
"  --output,     -o        Output filename prefix(.mat.gz will be added auto)\n"
"  --mapq,       -q <INT>  Mapping quality >= INT. [10]\n"
"  --verbose,    -v        Set verbose output\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* CONCAT_MESSAGE =
"Program: BaseVarC popmatrix\n"
"Contact: Zilong Li [lizilong@bgi.com]\n"
"Usage  : BaseVarC concat [options]\n\n"
"Commands:\n"
"  --input,      -l       List of matrix files for concat, one file per row.\n"
"  --output,     -o       Output filename prefix(.gz will be added auto)\n"
"\nReport bugs to lizilong@bgi.com \n\n";

static const char* CVG_HEADER =
"##fileformat=CVGv1.0\n"
"##Group information is the depth of A:C:G:T\n"
"#CHROM\tPOS\tREF\tDepth\tA\tC\tG\tT\tFS\tSOR\tStrand_Coverage(REF_FWD,REF_REV,ALT_FWD,ALT_REV)\n";

static const char* VCF_HEADER =
"##fileformat=VCFv4.2\n"
"##FILTER=<ID=LowQual,Description=\"Low quality (QUAL < 60)\">\n"
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
"##FORMAT=<ID=AB,Number=1,Type=String,Description=\"Allele Base\">\n"
"##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Strand orientation of the mapping base. Marked as + or -\">\n"
"##FORMAT=<ID=BP,Number=1,Type=String,Description=\"Base Probability which calculate by base quality\">\n"
"##INFO=<ID=CM_AF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of allele frequencies base on LRT algorithm\">\n"
"##INFO=<ID=CM_CAF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of allele frequencies just base on read count\">\n"
"##INFO=<ID=CM_AC,Number=A,Type=Integer,Description=\"An ordered, comma delimited allele depth in CMDB\">\n"
"##INFO=<ID=CM_DP,Number=A,Type=Integer,Description=\"Total Depth\">\n"
"##INFO=<ID=SB_REF,Number=A,Type=Integer,Description=\"Read number support REF: Forward,Reverse\">\n"
"##INFO=<ID=SB_ALT,Number=A,Type=Integer,Description=\"Read number support ALT: Forward,Reverse\">\n"
"##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">\n"
"##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Phred-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">\n"
"##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">\n"
"##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Phred-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">\n"
"##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Phred-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">\n"
"##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence Quality by Depth\">\n";
static const char col = ';';
static const char tab = '\t';

typedef std::string String;
typedef std::vector<String> StringV;
typedef std::vector<int32_t> IntV;
typedef std::vector<PosAlleleMap> PosAlleleMapVec;

struct BtRes
{
    String  cvg;
    String  vcf;
};

void runBaseType(int argc, char **argv);
void runPopMatrix(int argc, char **argv);
void runConcat(int argc, char **argv);
void parseOptions(int argc, char **argv, const char* msg);

void bt_read(StringV bams, const String& region, const IntV& pv, String fout);
BtRes bt_f(int32_t p, AlleleInfoVector aiv, DepM idx, int32_t N, String chr, int32_t rg_s, const String& seq);

namespace opt {
    static bool verbose = false;
    static int mapq;
    static int thread;
    static int batch;
    static std::string input;
    static std::string reference;
    static std::string posfile;
    static std::string region;
    static std::string output;
}

static const char* shortopts = "hvl:r:p:g:o:q:t:b:";

static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "input",                   required_argument, NULL, 'l' },
  { "reference",               required_argument, NULL, 'r' },
  { "posfile",                 required_argument, NULL, 'p' },
  { "region",                  required_argument, NULL, 'g' },
  { "output",                  required_argument, NULL, 'o' },
  { "batch",                   required_argument, NULL, 'b' },
  { "thread",                  required_argument, NULL, 't' },
  { "mapq",                    required_argument, NULL, 'q' },
  { NULL, 0, NULL, 0 }
};


int main(int argc, char** argv)
{
    if (argc <= 1 ) {
        std::cerr << BASEVARC_USAGE_MESSAGE;
        return 0;
    } else {
        std::string command(argv[1]);
        if (command == "help" || command == "--help") {
            std::cerr << BASEVARC_USAGE_MESSAGE;
            return 0;
        } else if (command == "basetype") {
            runBaseType(argc - 1, argv + 1);
        } else if (command == "popmatrix") {
            runPopMatrix(argc - 1, argv + 1);
        } else if (command == "concat") {
            runConcat(argc - 1, argv + 1);
        } else {
            std::cerr << BASEVARC_USAGE_MESSAGE;
            return 0;
        }
    }

    return 0;
}

void runConcat(int argc, char **argv)
{
    parseOptions(argc, argv, CONCAT_MESSAGE);
    String fm = opt::input;
    String fo = opt::output;
    clock_t ctb = clock();
    std::ifstream ifm(fm);
    StringV fm_v(std::istream_iterator<BaseVar::Line>{ifm},
	             std::istream_iterator<BaseVar::Line>{});
    int32_t i, k, m, count, n = 0, mt = 0;
    int32_t nm = fm_v.size();
    BGZF* fp = NULL;
    kstring_t ks;
    ks.s = 0; ks.l = 0; ks.m = 0;
    String ss;
    StringV Drg;
    std::vector<StringV> D;
    D.reserve(nm);
    for (k = 0; k < nm ; ++k) {
        std::cout << "reading file : " << fm_v[k] << std::endl;
        fp = bgzf_open(fm_v[k].c_str(), "r");
        if (bgzf_getline(fp, '\n', &ks) >= 0) {
            ss = ks.s;
            StringV tokens;
            String token;
            std::istringstream ts(ss);
            while (std::getline(ts, token, '\t')) tokens.push_back(token);
            if (n > 0 && n != std::stol(tokens[0])) {
                std::cerr << "error: the number of samples among the inputs are different!" << std::endl;
                exit(EXIT_FAILURE);
            } else {
                n = std::stol(tokens[0]);
                m = std::stol(tokens[1]);
                mt += m;
            }
        }
        count = 0;
        Drg.reserve(n);
        while (bgzf_getline(fp, '\n', &ks) >= 0) {
            ++count;
            Drg.push_back(ks.s);
        }
        if (count != n) {
            std::cerr << "error: the number of samples dont match the header!" << std::endl;
            exit(EXIT_FAILURE);
        }
        D.push_back(Drg);
        Drg.clear();
        if (bgzf_close(fp) < 0) {
            std::cerr << "warning: file cannot be closed!" << std::endl;
        }
    }

    fp = bgzf_open(fo.c_str(), "w");
    ss = std::to_string(n) + "\t" + std::to_string(mt) + "\n";
    if (bgzf_write(fp, ss.c_str(), ss.length()) != ss.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    // fast concatenation consider using rope data structure
    // @see https://brianbondy.com/blog/tagged/data-structure
    for (i = 0; i < n; ++i){
        ss = "";
        ss.reserve(mt+1);
        for (k = 0; k < nm; ++k) {
            ss += D[k][i];
        }
        ss += "\n";
        if (bgzf_write(fp, ss.c_str(), ss.length()) != ss.length()) {
            std::cerr << "fail to write" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    if (bgzf_close(fp) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cout << "elapsed secs : " << elapsed_secs << std::endl;
    std::cout << "concat done" << std::endl;
}

void runBaseType(int argc, char **argv)
{
    parseOptions(argc, argv, BASETYPE_MESSAGE);
    std::cerr << "basetype start" << std::endl;
    clock_t ctb = clock();
    std::ifstream ibam(opt::input);
    StringV bams(std::istream_iterator<BaseVar::Line>{ibam},
    	         std::istream_iterator<BaseVar::Line>{});
    const int32_t N = bams.size();
    // fetch ref bases;
    RefReader fa;
    String seq = fa.GetTargetBase(opt::region, opt::reference);
    String chr;
    int32_t rg_s, rg_e;
    std::tie(chr, rg_s, rg_e) = BaseVar::splitrg(opt::region);
    IntV pv;
    for (size_t i = 0; i < seq.length(); i++) {
        if (seq[i] == 'N') continue;
        pv.push_back(i + rg_s);       // 1-based
    }
    // begin to read bams
    String headcvg = String(CVG_HEADER);
    String headvcf = String(VCF_HEADER) + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    int thread = opt::thread;
    int bc = opt::batch;
    int bt = 1 + (N - 1) / bc;    // ceiling
    ThreadPool pool(thread);
    std::vector<std::future<void>> res;
    String region = opt::region;
    String tmp;
    StringV ftmp_v;
    std::cerr << "begin to read bams and save as tmp file" << std::endl;
    for (int i = 0; i < bt; ++i) {
        tmp = opt::output + ".batch." + std::to_string(i) + ".tmp";
        ftmp_v.push_back(tmp);
        if (i == bt - 1) {
            StringV bams_t(bams.begin() + i * bc, bams.end());
            res.emplace_back(pool.enqueue(bt_read, bams_t, std::cref(region), std::cref(pv), tmp));
        } else {
            StringV bams_t(bams.begin() + i * bc, bams.begin() + (i + 1)*bc);
            res.emplace_back(pool.enqueue(bt_read, bams_t, std::cref(region), std::cref(pv), tmp));
        }
    }
    for (auto && r: res) {
        r.get();
    }
    res.clear();
    // begin to call basetype and output
    String vcfout = opt::output + ".vcf.gz";
    String cvgout = opt::output + ".cvg.gz";
    BGZF* fpv = bgzf_open(vcfout.c_str(), "w");
    BGZF* fpc = bgzf_open(cvgout.c_str(), "w");
    std::vector<std::shared_ptr<std::ifstream> > fpiv;
    for (auto & f: ftmp_v) {
        fpiv.push_back(std::make_shared<std::ifstream>(f));
    }
    for (auto & fp: fpiv) {
        if (std::getline(*fp, tmp)) {
            headvcf += tmp;
        }
    }
    headvcf.pop_back();
    headvcf += "\n";
    if (bgzf_write(fpv, headvcf.c_str(), headvcf.length()) != headvcf.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (bgzf_write(fpc, headcvg.c_str(), headcvg.length()) != headcvg.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    size_t pos;
    String token;
    AlleleInfo ai;
    AlleleInfoVector aiv;
    DepM idx;
    int32_t j = 0, k = 0;
    // int buffer = 10000;
    // bt = 1 + (pv.size() - 1) / buffer;
    IntV pv_t;
    std::vector<std::future<BtRes>> res2;
    res2.reserve(pv.size());
    std::cerr << "begin to load data and run basetype" << std::endl;
    // for (int i = 0; i < bt; ++i) {
    //     if (i == bt - 1) {
    //         IntV t(pv.begin() + i * buffer, pv.end());
    //         pv_t = t;
    //     } else {
    //         IntV t(pv.begin() + i * buffer, pv.begin() + (i + 1)*buffer);
    //         pv_t = t;
    //     }
    // }
    for (auto & p : pv) {
        j = 0; k = 0;
        for (auto & fp: fpiv) {
            if (std::getline(*fp, tmp)) {
                while ((pos = tmp.find(' ')) != std::string::npos) {
                    token = tmp.substr(0, pos);
                    if (token != ".") {
                        std::istringstream iss(token);
                        iss >> ai;
                        aiv.push_back(ai);
                        idx.insert({j, k++});
                    }
                    tmp.erase(0, pos + 1);
                    j++;
                }
            }
        }
        res2.emplace_back(pool.enqueue(bt_f, p, aiv, idx, N, chr, rg_s, std::cref(seq)));
        aiv.clear();
        idx.clear();
    }
    for (auto && r: res2) {
        auto btr = r.get();
        if (bgzf_write(fpv, btr.vcf.c_str(), btr.vcf.length()) != btr.vcf.length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (bgzf_write(fpc, btr.cvg.c_str(), btr.cvg.length()) != btr.cvg.length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    res2.clear();
    if (bgzf_close(fpv) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    if (bgzf_close(fpc) < 0) std::cerr << "warning: file cannot be closed" << std::endl;
    // remove tmp file
    for (auto & f: ftmp_v) {
        std::remove(f.c_str());
    }
    // done
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cerr << "basetype done" << std::endl;
    std::cout << "elapsed secs : " << elapsed_secs << std::endl;
}

void bt_read(StringV bams, const String& region, const IntV& pv, String fout)
{
    char space = ' ';
    char comma = ',';
    size_t n = bams.size();
    PosAlleleMapVec allele_mv;
    String names;
    allele_mv.reserve(n);
    for (auto const& bam: bams) {
        BamProcess reader;
        if (!reader.Open(bam)) {
            std::cerr << "ERROR: could not open file " << bam << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!reader.FindSnpAtPos(region, pv)) {
            std::cerr << "Warning: " << reader.sm << " region " << region << " is empty." << std::endl;
        }
        allele_mv.push_back(reader.allele_m);
        names += reader.sm + '\t';
        if (!reader.Close()) {
            std::cerr << "Warning: could not close file " << bam << std::endl;
        }
    }
    assert(allele_mv.size() == n);
    // names.pop_back();
    names += "\n";
    std::ofstream fo(fout);
    fo << names;
    for (auto & p : pv) {
        for (auto const& m : allele_mv) {
            if (m.count(p)) {
                auto const& a = m.at(p);
                fo << a.base << comma << a.mapq << comma << a.qual << comma << a.rpr << comma << a.strand << space;
            } else {
                fo << ". ";
            }
        }
        fo << "\n";
    }
}

BtRes bt_f(int32_t p, AlleleInfoVector aiv, DepM idx, int32_t N, String chr, int32_t rg_s, const String& seq)
{
    int8_t alt_base, ref_base;
    int32_t j = 0, na, nc, ng, nt, ref_fwd, ref_rev, alt_fwd, alt_rev;
    double fs, sor, left_p, right_p, twoside_p;
    double min_af = 0.001;
    BaseV bases, quals;
    IntV tmp;
    std::ostringstream sout;
    BtRes res;
    if (aiv.size() > 0) {
        // output cvg;
        ref_base = BASE_INT8_TABLE[static_cast<int>(seq[p - rg_s])];
        na = 0, nc = 0, ng = 0, nt = 0;
        for (auto const& a: aiv) {
            switch (a.base) {
            case 0 : na += 1; break;
            case 1 : nc += 1; break;
            case 2 : ng += 1; break;
            case 3 : nt += 1; break;
            }
            bases.push_back(a.base);
            quals.push_back(a.qual);
        }
        tmp = {na, nc, ng, nt};
        std::sort(tmp.begin(), tmp.end(), std::greater<int>());
        if (tmp[0] != ref_base) {
            if (tmp[0] == na) alt_base = 0;
            else if (tmp[0] == nc) alt_base = 1;
            else if (tmp[0] == ng) alt_base = 2;
            else alt_base = 3;
        } else{
            if (tmp[1] == na) alt_base = 0;
            else if (tmp[1] == nc) alt_base = 1;
            else if (tmp[1] == ng) alt_base = 2;
            else alt_base = 3;
        }
        ref_fwd = 0, ref_rev = 0, alt_fwd = 0, alt_rev = 0;
        for (auto const& a: aiv) {
            if (a.strand == 1) {
                if (a.base == ref_base) {
                    ref_fwd += 1;
                } else if (a.base == alt_base) {
                    alt_fwd += 1;
                }
            } else if (a.strand == 0) {
                if (a.base == ref_base) {
                    ref_rev += 1;
                } else if (a.base == alt_base) {
                    alt_rev += 1;
                }
            }
        }
        kt_fisher_exact(ref_fwd, ref_rev, alt_fwd, alt_rev, &left_p, &right_p, &twoside_p);
        fs = -10 * log10(twoside_p);
        if (std::isinf(fs)) fs = 10000.0;
        else if (fs == 0) fs = 0.0;
        if (alt_fwd * ref_rev > 0) {
            sor = static_cast<double>(ref_fwd * alt_rev) / (ref_rev * alt_fwd);
        } else {
            sor = 10000.0;
        }
        sout << chr << tab << p << tab << BASE2CHAR[ref_base] << tab << aiv.size() << tab << na << tab << nc << tab << ng << tab << nt << tab << fs << tab << sor << tab << ref_fwd << col << ref_rev << col << alt_fwd << col << alt_rev << "\n";
        res.cvg = sout.str();
        // basetype caller;
        BaseType bt(bases, quals, ref_base, min_af);
        if (bt.LRT()) {
            String outvcf = bt.WriteVcf(bt, chr, p, ref_base, aiv, idx, N);
            res.vcf = outvcf;
        }
    }

    return res;
}

void runPopMatrix (int argc, char **argv)
{
    parseOptions(argc, argv, POPMATRIX_MESSAGE);
    std::cerr << "popmatrix start" << std::endl;
    clock_t ctb = clock();
    std::ifstream ibam(opt::input);
    std::ifstream ipos(opt::posfile);
    if (!ibam.is_open() || !ipos.is_open()) {
        std::cerr << "bamlist or posifle can not be opend" << std::endl;
        exit(EXIT_FAILURE);
    }
    StringV bams(std::istream_iterator<BaseVar::Line>{ibam},
    	         std::istream_iterator<BaseVar::Line>{});
    PosInfoVector pv;
    for (PosInfo p; ipos >> p;) pv.push_back(p);
    pv.shrink_to_fit();        // request for the excess capacity to be released

    const int32_t N = bams.size();
    const int32_t M = pv.size();
    String out = std::to_string(N) + "\t" + std::to_string(M) + "\n";
    BGZF* fp = bgzf_open(opt::output.c_str(), "w");
    if (bgzf_write(fp, out.c_str(), out.length()) != out.length()) {
        std::cerr << "fail to write - exit" << std::endl;
        exit(EXIT_FAILURE);
    }
    // ready for run
    String rg = pv.front() + pv.back();
    int32_t count = 0;
    for (int32_t i = 0; i < N; i++) {
        BamProcess reader;
        if (!(++count % 1000)) std::cerr << "Processing the number " << count / 1000 << "k bam" << std::endl;
        if (!reader.Open(bams[i])) {
            std::cerr << "ERROR: could not open file " << bams[i] << std::endl;
            exit(EXIT_FAILURE);
        }
        reader.FindSnpAtPos(rg, pv);
        String out(reader.snps.begin(), reader.snps.end());
        out += "\n";
        if (bgzf_write(fp, out.c_str(), out.length()) != out.length()) {
            std::cerr << "fail to write - exit" << std::endl;
            exit(EXIT_FAILURE);
        }
        if (!reader.Close()) {
            std::cerr << "Warning: could not close file " << bams[i] << std::endl;
        }
    }
    if (bgzf_close(fp) < 0) std::cerr << "failed to close \n";
    clock_t cte = clock();
    double elapsed_secs = double(cte - ctb) / CLOCKS_PER_SEC;
    std::cout << "elapsed secs : " << elapsed_secs << std::endl;
    std::cerr << "popmatrix done" << std::endl;
}

void parseOptions(int argc, char **argv, const char* msg)
{
    bool die = false;
    bool help = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case 'v': opt::verbose = true; break;
        case 'q': arg >> opt::mapq; break;
        case 'b': arg >> opt::batch; break;
        case 't': arg >> opt::thread; break;
        case 'l': arg >> opt::input; break;
        case 'r': arg >> opt::reference; break;
        case 'p': arg >> opt::posfile; break;
        case 'g': arg >> opt::region; break;
        case 'o': arg >> opt::output; break;
        default: die = true;
        }
    }
    // todo : need more check
    if (die || help || (opt::input.empty() && opt::output.empty())) {
        std::cerr << msg;
        if (die) exit(EXIT_FAILURE);
        else exit(EXIT_SUCCESS);
    }
}
