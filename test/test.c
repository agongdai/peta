#include <stdio.h>
#include <unistd.h>

#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))

int main() {
	int test = 100;
	int some = test;
	printf("%d \n", test);
	kroundup32(test);
	printf("%d \n", test);
	return 0;
}
