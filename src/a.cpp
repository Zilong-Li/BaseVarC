#include <iostream>

struct tinyChar {
    unsigned char x:4;
};

int main(){
    tinyChar a;
    std::cout << sizeof(a) << std::endl;
}
