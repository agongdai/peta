#include <vector>
#include <cassert>
#include <fstream>
#include <string>
#include <ctime>
#include <algorithm>
#include <stdint.h>
#include "kmers.h"
#include "bitarray/sdarray.h"
#include "bitarray/utest.h"

using namespace std;
using namespace mscds;

void test_sdarray() {
	SDArrayBuilder bd;
	int N = 10000;
	for (int i = 0; i < N; ++i) {
		bd.add(i);
	}
	SDArrayQuery sda;
	bd.build(&sda);

	printf("SDA length: %d \n", sda.length());
	for (int i = 0; i < 10; ++i) {
		printf("SD %d*%d: %d \n", i, i, sda.prefixsum(i * i));
		printf("SD %d*%d: %d \n", i, i, sda.lookup(i * i));
	}
}
