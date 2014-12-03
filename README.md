# Graph Summarization with Quality Guarantees

This repository contains the code for the algorithms presented in

Matteo Riondato, David García-Soriano, Francesco Bonchi, "Graph Summarization
with Quality Guarantees", Proceedings of IEEE International Conference on Data
Mining 2014 (ICDM'14) [PDF](http://francescobonchi.com/icdm14CR.pdf).

## License

Everything in this repository is distributed under the Apache License, version
2.0. See the [LICENSE](LICENSE) file and the [NOTICE](NOTICE) file.

## Contacts

For any question or to point out bugs, please contact Matteo Riondato
<rionda@cs.stanford.edu>.

## Compilation

To compile the code, just run `make`. The code is tested on a GNU/Linux system
using the GNU g++ compiler, version 4.9. It requires the [SDPA
libraries](http://sdpa.sourceforge.net/), and the [MUMPS
libraries](http://mumps.enseeiht.fr/) (sequential version).

## Execution

The "summ" binary is the algorithm for graph summarization described in the
article. See the output of `./summ -h` for help about using it.
A typical usage example would be:

```
./summ -a -k 100 -t 2 -d -1 edge_file.txt
```

The "grass" binary is the algorithm for graph summarization by LeFevre and
Terzi. See the output of `./summ -h` for help about using it. 
A typical usage example would be:

```
./grass -a -k 100 -t 2 -c 0.75 edge_file.txt
```

