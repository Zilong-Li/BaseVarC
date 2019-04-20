#ifndef BASEVAR_UTILS_H
#define BASEVAR_UTILS_H

namespace BaseVar
{
  template<typename T> 
  inline std::string tostring(T d) { 
      std::stringstream ss;
      ss << d;
    return ss.str();
  }

  struct SamRg {
      std::string chr;
      uint32_t s;
      uint32_t e;
  };

  inline SamRg torg(const std::string& rg_t) {
      std::string rg = rg_t;
      size_t p;
      SamRg srg;
      if ((p = rg.find(":")) != std::string::npos) {
	  srg.chr = rg.substr(0, p);
	  rg.erase(0, p + 1);
      }
      if ((p = rg.find("-")) != std::string::npos) {
	  srg.s = std::stoi(rg.substr(0,p)) - 1;      // make 0-based
	  rg.erase(0, p + 1);
	  srg.e = std::stoi(rg) - 1;
      }
  } 

}

#endif