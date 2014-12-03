#include <array>
#include <cstdio>
#include <cstring>
#include <iostream>

#include <unistd.h>
using namespace std;

const int maxlen = 1000;
char in_filename[maxlen];

#include "graph.h"
#include "summ.h"

void usage(char *binary_name) {
    fprintf(stderr, "USAGE: %s edgefile resfile\n", binary_name);
}

int main(int argc, char *argv[]) {
    setbuf(stdout, NULL);   // no buffering for stdout

    if (argc != 3) {
	usage(argv[0]);
	return 1;
    }
    strcpy(in_filename, argv[1]);
    strcpy(summary_filename, argv[2]);
    read_graph();
    read_summary();
    print_summary_as_full_matrix();

    return 0;
}

