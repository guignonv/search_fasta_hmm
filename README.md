# Search FASTA HMM

Scans FASTA files with HMM profiles using HMM search.


# Getting Started

Clone the git repository.
Edit `search_gelp_hmm.pl` and adjust the line with your system as needed:

  our $HMMSEARCH_COMMAND = 'hmmsearch';

and you can run search_gelp_hmm:

  perl search_gelp_hmm.pl -...


# Prerequisites

* [HMMER](http://hmmer.org/) 2.3 and above.
* A list of HMM profile files (ending with ".hmm") in a same directory.
* Some FASTA files to scan.


# Usage

Create a directory four output (optional):

  mkdir gelp_hmm_test

Run the script:

  perl search_gelp_hmm.pl -d examples/hmm -f examples/fasta/example_gelp.fa -o ./gelp_hmm_test -a -e 1e-3 -t 8


# Contributing

Alberto CENCI
Valentin GUIGNON
Mathieu ROUARD


# License

Distributed under the GPLv3 License. See `LICENSE` for more information.


# Contact

Valentin GUIGNON - v.guignon@cgiar.org

Project Link: [https://github.com/guignonv/search_fasta_hmm](https://github.com/guignonv/search_fasta_hmm)


# Acknowledgments

[HMMER](http://hmmer.org/)
