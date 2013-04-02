
#ifndef __UTILS_DEBUG_H_
#define __UTILS_DEBUG_H_

#include <string>
#include <iostream>
#include <stdint.h>

namespace utils {

struct DbgHelper {
	static bool flag;
	static uint64_t counter_;

	DbgHelper(){ flag = false; }
	static void print(const std::string& s) {
		if (flag) std::cerr << s << std::endl;
	}
	static uint64_t& counter() { return counter_;}
	static bool is_on() { return flag; }
	static void setflag(bool b) { flag = b; }
};


}//namespace
#endif //__UTILS_DEBUG_H_
