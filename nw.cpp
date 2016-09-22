//Klaudia Algiz, 333811
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <climits>
#include <cstring>

using namespace std;

const int SEQ_SIZE = 60;
const int GP = 2;
int MIN_WEIGHT = 10;
const int KEYS_TO_SEARCH = 128;
const int TRIGRAMS_SIZE = 58;
const int TRIGRAMS_AMOUNT = 125;

struct Seq
{
	char seq[SEQ_SIZE];
	long int nr;
	bool comple;
};

struct SeqPair
{
	int key[TRIGRAMS_SIZE];
	Seq value;
};

inline bool comparator (const SeqPair& seq1, const SeqPair& seq2)
{
	for (int i = 0; i < TRIGRAMS_SIZE; i++)
	{
		if (seq1.key[i] < seq2.key[i])
		{
			return true;
		}
		else if (seq1.key[i] > seq2.key[i])
		{
			return false;
		}
	}
	return false;
}

int signToInt(char sign)
{
	switch(sign)
	{
		case 'a':
			return 0;
		case 'c':
			return 1;
		case 'g':
			return 2;
		case 't':
			return 3;
		default:
			return 4;
	}
}

int encode(char trigram[3])
{
	return 25*signToInt(trigram[0]) + 5*signToInt(trigram[1]) + signToInt(trigram[2]);
}

std::vector<Seq> buildComplementary(const std::vector<Seq>& fragments)
{
	std::vector<Seq> result;
	for (std::size_t i = 0; i < fragments.size(); i++)
	{
		Seq seq;
		for (int j = 0; j < SEQ_SIZE; j++)
		{
			switch( fragments[i].seq[j] )
			{
				case 'a':
					seq.seq[SEQ_SIZE-j-1] ='t';
					break;
				case 'c':
					seq.seq[SEQ_SIZE-j-1] = 'g';
					break;
				case 'g':
					seq.seq[SEQ_SIZE-j-1] = 'c';
					break;
				case 't':
					seq.seq[SEQ_SIZE-j-1] = 'a';
					break;
				case 'n':
					seq.seq[SEQ_SIZE-j-1] = 'n';
			}
		}
		seq.nr = fragments[i].nr;
		seq.comple = true;
		result.push_back(seq);
	}
	return result;
}

void sort_trigrams(char *seq, int * res)
{
	int occs[TRIGRAMS_AMOUNT];
	for (int i = 0; i < TRIGRAMS_AMOUNT; i++)
	{
		occs[i] = 0;
	}
	for (int i = 0; i < TRIGRAMS_SIZE; i++)
	{
		occs[encode(seq + i)]++;
	}
	
	//getting only existing trigrams
	int index = 0;
	for (int i = 0; i < TRIGRAMS_AMOUNT; i++)
	{
		for (int j = 0; j < occs[i]; j++)
		{
			res[index] = i;
			index++;
		}
	}
	
	//bubble sort
	for (int i = 0; i < TRIGRAMS_SIZE-1; i++)
	{
		for (int j = i+1; j < TRIGRAMS_SIZE; j++)
		{
			if (occs[res[j-1]] < occs[res[j]])
			{
				int tmp = res[j-1];
				res[j-1] = res[j];
				res[j] = tmp;
			}
			else if ((occs[res[j-1]] == occs[res[j]]) && (res[j-1] > res[j]))
			{
				int tmp = res[j-1];
				res[j-1] = res[j];
				res[j] = tmp;
			}
		}
	}
}

void compareSeqences(Seq seq1, Seq seq2)
{
	int A[SEQ_SIZE+1][SEQ_SIZE+1];
	int B[SEQ_SIZE+1][SEQ_SIZE+1];
// 	initialization
	for (int i = 0; i <= SEQ_SIZE; i++)
	{
		A[i][0] = 0;
		A[0][i] = 0;
		for (int j = 0; j <= SEQ_SIZE; j++)
			B[i][j] = 0;
	}
	
	for (int i = 1; i <= SEQ_SIZE; i++)
	{
		for (int j = 1; j <= SEQ_SIZE; j++)
		{
			if (seq1.seq[i-1] == seq2.seq[j-1]) A[i][j] = 1;
			else A[i][j] = -2;
		}
	}
	
// 	computing matrix
	for (int i = 1; i <= SEQ_SIZE; i++)
	{
		for (int j = 1; j <= SEQ_SIZE; j++)
		{
			int maks = std::max(A[i-1][j-1] + A[i][j], std::max(
				A[i][j-1]-GP, A[i-1][j]-GP) );
			if (maks == A[i-1][j-1] + A[i][j])
				B[i][j] = 1;
			else if (maks == A[i][j-1]-GP)
				B[i][j] = 2;
			else
				B[i][j] = 3;
			A[i][j] = maks;
		}
	}
	
// 	finding max weight
	int weight = INT_MIN;
	std::pair<int, int> maxCoordinates;
	for (int i = 1; i <= SEQ_SIZE; i++)
	{
		if (A[i][SEQ_SIZE] > weight) 
		{
			maxCoordinates.first = i;
			maxCoordinates.second = SEQ_SIZE;
			weight = A[i][SEQ_SIZE];
		}
		if (A[SEQ_SIZE][i] > weight)
		{
			maxCoordinates.first = SEQ_SIZE;
			maxCoordinates.second = i;
			weight = A[SEQ_SIZE][i];
		}
	}
	if (weight >= MIN_WEIGHT)
	{
		std::vector<std::pair<int, int> > alignment;
		int row = maxCoordinates.first;
		int col = maxCoordinates.second;
		
		while (row > 0 && col > 0)
		{
			alignment.push_back(std::make_pair(row, col));
			switch (B[row][col])
			{
				case 1:
					row -= 1;
					col -= 1;
					break;
				case 2:
					col -= 1;
					break;
				case 3:
					row -= 1;
					break;
			}
		}
		std::reverse(alignment.begin(), alignment.end());
		std::cout << seq1.nr << "; " << seq2.nr << "; "
			<< weight << "; "
			<< (seq1.comple ? "1" : "0") << ", "
			<< (seq2.comple ? "1" : "0") << "; " << "{";
		for (std::size_t i = 0; i < alignment.size()-1; i++)
			std::cout << alignment[i].first << "-" << alignment[i].second << ";";
		std::cout << alignment[alignment.size()-1].first << "-" << 
			alignment[alignment.size()-1].second << "}" << std::endl;
	}
}

void countDistances(std::vector<SeqPair>& fragments)
{
	for(std::size_t i = 0; i < fragments.size(); i++)
	{
		for (std::size_t j = i+1; j < std::min(fragments.size(), i+KEYS_TO_SEARCH+1); j++)
			compareSeqences(fragments[i].value, fragments[j].value);
	}
}

std::vector<Seq> readEntrance()
{
	std::vector<Seq> result;
	std::string prev;
	long int prevNumber = 0;
	for (std::string line; std::getline(std::cin, line);) {
		if (line.size() > 0 && line[0] == '>')
		{
			if (prev != "") 
			{
				Seq seq;
				strcpy(seq.seq, prev.c_str());
				seq.nr = prevNumber;
				seq.comple = false;
				result.push_back(seq);
			}
			prev = "";
			prevNumber = atoi(line.substr(3).c_str());
		}
		else 
			prev += line;
	}
	Seq seq;
	strcpy(seq.seq, prev.c_str());
	seq.nr = prevNumber;
	seq.comple = false;
	result.push_back(seq);
	return result;
}

int main (int argc, char* argv[]) {
	if (argc != 1 && argc != 2)
	{
		std::cout << "wrong number of parameters\n";
		return 0;
	}
	if (argc == 2)
		MIN_WEIGHT = atoi(argv[1]);
	
	std::vector<Seq> fragments = readEntrance();

	std::vector<Seq> complFr = buildComplementary(fragments);
	fragments.insert(fragments.end(), complFr.begin(), complFr.end());
	
	std::vector<SeqPair> seqMap;
	for (std::size_t i = 0; i < fragments.size(); i++)
	{
		SeqPair sp;
		sort_trigrams(fragments[i].seq, sp.key);
		sp.value = fragments[i];
		seqMap.push_back(sp);
	}
	std::sort(seqMap.begin(), seqMap.end(), comparator);
	countDistances(seqMap);
	
  return 0;
}