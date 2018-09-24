#include "pnt.hpp"

void pprint(int i, int n, bool stay)
{
	int p = (i+1) * 100 / n;
	if (p != i * 100 / n)
		printf("%d%%\r", p);
	if (stay && p == 100)
		putchar('\n');
}