#include <iostream>
#include <cassert>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <list>

using namespace std;

struct Alignment {
	string ref_id;
	int strand : 1;
	int q_from : 31;
	int q_to : 32;
	long long ref_from : 48;
	long long ref_to : 48;

	Alignment(const string &ref_id = "", int strand = 0, int q_from = 0, int q_to = 0, long long ref_from = 0, long long ref_to = 0):
		ref_id(ref_id), strand(strand), q_from(q_from), q_to(q_to), ref_from(ref_from), ref_to(ref_to) {}
};

void parseCigar(const string &cigar, int &q_from, int &q_to, int &algn_len) {
	q_from = q_to = algn_len = 0;
	istringstream is(cigar);
	int num;
	char type;
	bool start = false;

	while (is >> num >> type) {
		switch (type) {
			case 'S':
			case 'H':
			if (!start) {
				q_from += num;
				q_to = q_from;
			}
			break;

			case 'M':
			q_to += num;
			algn_len += num;
			break;

			case 'I':
			q_to += num;
			break;

			case 'D':
			case 'N':
			algn_len += num;
			break;

			default:
			assert(false);
		}

		start = true;
	}
};

bool overlap(int f1, int t1, int f2, int t2) {
	return t1 > f2 && t2 > f1;
}

int main(int argc, char **argv) {
	if (argc < 3) {
		fprintf(stderr, "Usage: %s <background.sam> <2.sam>\n", argv[0]);
		return 1;
	}

	map<string, list<Alignment> > sam1;
	ifstream fs(argv[1]);
	string buf;
	string read_id;
	string ref_id;
	int sam_flag;
	long long ref_from;
	string cigar;

	int num_sam1_alignments = 0;

	while (getline(fs, buf)) {
		if (buf.length() == 0 || buf[0] == '@') { continue; }
		istringstream is(buf);
		is >> read_id >> sam_flag >> ref_id >> ref_from >> buf >> cigar;

		if (sam_flag & 0x4) { continue; }

		int q_from, q_to, algn_len;
		int strand = !!(0x16 & sam_flag);
		parseCigar(cigar, q_from, q_to, algn_len);
		sam1[read_id].push_back(Alignment(ref_id, strand, q_from, q_to, ref_from, ref_from + algn_len));
		num_sam1_alignments++;
	}

	ifstream fs2(argv[2]);
	map<string, int> sam2_aligned_reads; // [] = 1: aligned & match at least 1 background alignment; [] = 0: aligned but do not match

	int num_sam2_alignments = 0;
	int num_sam2_alignments_matched = 0;
	int sam2_aligned_reads_matched = 0;

	while (getline(fs2, buf)) {
		if (buf.length() == 0 || buf[0] == '@') { continue; }
		istringstream is(buf);
		is >> read_id >> sam_flag >> ref_id >> ref_from >> buf >> cigar;

		if (sam_flag & 0x4) { continue; }

		int q_from, q_to, algn_len;
		int strand = !!(0x16 & sam_flag);
		parseCigar(cigar, q_from, q_to, algn_len);
		sam2_aligned_reads[read_id] |= 0;

		auto map_it = sam1.find(read_id);
		if (map_it != sam1.end()) {
			for (auto list_it = map_it->second.begin(); list_it != map_it->second.end(); ++list_it) {
				if (list_it->strand == strand && list_it->ref_id == ref_id &&
					overlap(list_it->ref_from, list_it->ref_to, ref_from, ref_from + algn_len) &&
					overlap(list_it->q_from, list_it->q_to, q_from, q_to)) {

					if (sam2_aligned_reads[read_id] == 0) {
						sam2_aligned_reads[read_id] = 1;
						++sam2_aligned_reads_matched;
					}

					num_sam2_alignments_matched++;
					break;
				}
			}
		}

		num_sam2_alignments++;
	}

	cout << argv[1] << ": " << endl;
	cout << "\t# of aligned reads: " << sam1.size() << endl;
	cout << "\t# of alignments: " << num_sam1_alignments << endl;

	cout << argv[2] << ": " << endl;
	cout << "\t# of aligned reads: " << sam2_aligned_reads.size() << endl;
	cout << "\t# of reads that aligned and match at least 1 alignment of " << argv[1] << ": " << sam2_aligned_reads_matched << endl;
	cout << "\t# of alignments: " << num_sam2_alignments << endl;
	cout << "\t# of alignments that match " << argv[1] << ": " << num_sam2_alignments_matched << endl;

	return 0;
}