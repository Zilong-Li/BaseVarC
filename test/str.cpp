#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using std::cout;
using std::endl;

struct pos_info {
    std::string chr;
    int pos;
    int ref;
    friend std::istream& operator>>(std::istream& is, pos_info& info) {
	is >> std::ws; //skip whitespace at the beginning
	is >> info.chr;
	is >> info.pos;
	is >> info.ref;
	return is;
    }
};

int main(int argc, char** argv)
{
    if (argc < 1) {
	cout << "no input line" << endl;
	exit(1);
    }
    std::string f = argv[1];
    std::ifstream ifs;
    ifs.open(f);
    if (!ifs.is_open()) {
	cout << "no input file" << endl;
    }

    // for (std::string item;ifs >> item;) {
    // 	cout << item << endl;
    // }
    pos_info info;
    while (ifs >> info) {
    	cout << info.pos << endl;
    }


    return 0;
}