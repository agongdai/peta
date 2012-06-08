#include <stdio.h>
#include <unistd.h>

int main() {
	int calendar[12][31];
	int month, day, counter = 0;
	for (month = 0; month < 12; month++) {
		for (day = 0; day < 31; day++) {
			calendar[month][day] = counter++;
		}
	}
	printf("[0][0] = %d \n", calendar[0][0]);
	printf("[9][1] = %d \n", calendar[9][1]);
	printf("[9][1] = %d \n", *(*(calendar + 9) + 1));
	printf("[9][1] = %d \n", **(&(calendar[9]) + 1));
	return 0;
}
