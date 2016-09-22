//Klaudia Algiz, 333811
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <climits>
#include <cstring>
#include <stdio.h>

using namespace std;

const int SEQ_SIZE = 60;
const int GP = 2;
int MIN_WEIGHT = 10;
const int KEYS_TO_SEARCH = 128;
const int TRIGRAMS_SIZE = 58;
const int TRIGRAMS_AMOUNT = 125;
const int BLOCK_DIM = KEYS_TO_SEARCH;
const int MAX_AT_TIME = 40000;

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

struct Edge
{
	long int nr1;
	long int nr2;
	int weight;
	bool compl1;
	bool compl2;
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

__global__ void sort_trigrams(char * seq, int * res, unsigned long fragmentsSize)
{
	int idx = blockIdx.x * BLOCK_DIM + threadIdx.x;
	if (idx >= fragmentsSize)
	{
		return;
	}
	
	int occs[TRIGRAMS_AMOUNT]; //how to use shared memory
	for (int i = 0; i < TRIGRAMS_AMOUNT; i++)
	{
		occs[i] = 0;
	}
	for (int i = 0; i < TRIGRAMS_SIZE; i++)
	{
		int vals[3];
		for (int j = 0; j < 3; j++)
		{
			switch(seq[idx*SEQ_SIZE+i+j])
			{
				case 'a':
					vals[j] = 0;
					break;
				case 'c':
					vals[j] = 1;
					break;
				case 'g':
					vals[j] = 2;
					break;
				case 't':
					vals[j] = 3;
					break;
				default:
					vals[j] = 4;
					break;
			}
		}
		int code = 25*vals[0] + 5*vals[1] + vals[2];
		occs[code]++;
	}
	
	//getting only existing trigrams
	int index = 0;
	for (int i = 0; i < TRIGRAMS_AMOUNT; i++)
	{
		for (int j = 0; j < occs[i]; j++)
		{
			res[idx*TRIGRAMS_SIZE+index] = i;
			index++;
		}
	}
	
	//bubble sort
	for (int i = 0; i < TRIGRAMS_SIZE-1; i++)
	{
		for (int j = i+1; j < TRIGRAMS_SIZE; j++)
		{
			if (occs[res[idx*TRIGRAMS_SIZE+j-1]] < 
				occs[res[idx*TRIGRAMS_SIZE+j]])
			{
				int tmp = res[idx*TRIGRAMS_SIZE+j-1];
				res[idx*TRIGRAMS_SIZE+j-1] = res[idx*TRIGRAMS_SIZE+j];
				res[idx*TRIGRAMS_SIZE+j] = tmp;
			}
			else if ((occs[res[idx*TRIGRAMS_SIZE+j-1]] == 
				occs[res[idx*TRIGRAMS_SIZE+j]]) && 
				(res[idx*TRIGRAMS_SIZE+j-1] > res[idx*TRIGRAMS_SIZE+j]))
			{
				int tmp = res[idx*TRIGRAMS_SIZE+j-1];
				res[idx*TRIGRAMS_SIZE+j-1] = res[idx*TRIGRAMS_SIZE+j];
				res[idx*TRIGRAMS_SIZE+j] = tmp;
			}
		}
	}
	
}

std::vector<SeqPair> sort_all_trigrams(std::vector<Seq> fragments)
{
	std::vector<SeqPair> result;
	char * input;
	char * devInput;
	int * output;
	int * devOutput;
	cudaMalloc(&devInput, fragments.size()*SEQ_SIZE*sizeof(char));
	cudaMalloc(&devOutput, fragments.size()*TRIGRAMS_SIZE*sizeof(int));
	cudaMallocHost(&input, fragments.size()*SEQ_SIZE*sizeof(char));
	cudaMallocHost(&output, fragments.size()*TRIGRAMS_SIZE*sizeof(int));
	for (int i = 0; i < fragments.size(); i++)
	{
		for (int j = 0; j < SEQ_SIZE; j++)
		{
			input[i*SEQ_SIZE+j] = fragments[i].seq[j];
		}
	}
	cudaMemcpy(devInput, input, fragments.size()*SEQ_SIZE*sizeof(char), 
							cudaMemcpyHostToDevice);
	int blocksNr = fragments.size()/BLOCK_DIM + (fragments.size()%BLOCK_DIM == 0 ? 0 : 1);
	sort_trigrams<<< blocksNr, BLOCK_DIM >>>(devInput, devOutput, fragments.size());
	cudaMemcpy(output, devOutput, fragments.size()*TRIGRAMS_SIZE*sizeof(int), 
						cudaMemcpyDeviceToHost);
	for (int i = 0; i < fragments.size(); i++)
	{
		SeqPair sp;
		sp.value = fragments[i];
		for (int j = 0; j < TRIGRAMS_SIZE; j++)
		{
			sp.key[j] = output[i*TRIGRAMS_SIZE+j];
		}
		result.push_back(sp);
	}
	cudaFree(devInput);
	cudaFree(devOutput);
	cudaFreeHost(input);
	cudaFreeHost(output);
	return result;
}

__global__ void countDistances(Seq * fragments, Edge * graph, 
										unsigned long fragmentsSize, int loopNr)
{
	int idx1 = blockIdx.x;
	__shared__ Seq seq1;
	if (threadIdx.x == 0)
	{
		seq1.nr = fragments[idx1].nr;
		seq1.comple = fragments[idx1].comple;
		for (int i = 0; i < SEQ_SIZE; i++)
			seq1.seq[i] = fragments[idx1].seq[i];
	}
	syncthreads();
	int idx2 = idx1 + threadIdx.x+1;
	if (idx1 >= fragmentsSize || idx2 >= fragmentsSize)
	{
		Edge * edge = graph + blockIdx.x*BLOCK_DIM+threadIdx.x;
		edge->weight = -1;
		return;
	}
	Seq *seq2 = fragments+idx2;
	
	int A[SEQ_SIZE+1][SEQ_SIZE+1];
 	// initialization
	for (int i = 0; i <= SEQ_SIZE; i++)
	{
		A[i][0] = 0;
		A[0][i] = 0;
	}
 		
	for (int i = 1; i <= SEQ_SIZE; i++)
		for (int j = 1; j <= SEQ_SIZE; j++)
		{
			if (seq1.seq[i-1] == seq2->seq[j-1]) A[i][j] = 1;
			else A[i][j] = -2;
		}
 		
 	// 	computing matrix
 		for (int i = 1; i <= SEQ_SIZE; i++)
		{
 			for (int j = 1; j <= SEQ_SIZE; j++)
 			{
 				A[i][j] = max(A[i-1][j-1] + A[i][j], max(
 					A[i][j-1]-GP, A[i-1][j]-GP) );
 			}
		}
		
 		int maks = INT_MIN;
		for (int i = 0; i <= SEQ_SIZE; i++)
		{
			maks = max(A[SEQ_SIZE][i], maks);
			maks = max(A[i][SEQ_SIZE], maks);
		}
 		Edge * edge = graph + blockIdx.x*BLOCK_DIM+threadIdx.x;
 		edge->nr1 = seq1.nr;
		edge->nr2 = seq2->nr;
		edge->compl1 = seq1.comple;
		edge->compl2 = seq2->comple;
		edge->weight = maks;
}

void createGraph(std::vector<Seq> fragments)
{
	Seq * input;
	Seq * devInput;
	Edge * output;
	Edge * devOutput;
	int graphEdgeNr = ((fragments.size()+MAX_AT_TIME-1)/MAX_AT_TIME)*MAX_AT_TIME*KEYS_TO_SEARCH;
	cudaMallocHost(&input, fragments.size()*sizeof(Seq));
	cudaMallocHost(&output, graphEdgeNr*sizeof(Edge));
	cudaMalloc(&devInput, (MAX_AT_TIME+KEYS_TO_SEARCH)*sizeof(Seq));
	cudaMalloc(&devOutput, MAX_AT_TIME*KEYS_TO_SEARCH*sizeof(Edge));
	for (int i = 0; i < fragments.size(); i++)
	{
		input[i] = fragments[i];
	}
	
	for (int i = 0; i < (fragments.size()+MAX_AT_TIME-1)/MAX_AT_TIME; i++)
	{
		int seqsToCopy = min(MAX_AT_TIME+KEYS_TO_SEARCH, (int)(fragments.size()-i*MAX_AT_TIME));
		cudaMemcpy(devInput, input+i*MAX_AT_TIME, seqsToCopy*sizeof(Seq),
				cudaMemcpyHostToDevice);
		int blocksNr = MAX_AT_TIME;
		countDistances<<< blocksNr, BLOCK_DIM >>>(devInput, devOutput, fragments.size(), i);
		cudaMemcpy(output+i*MAX_AT_TIME*KEYS_TO_SEARCH, devOutput, MAX_AT_TIME*KEYS_TO_SEARCH*sizeof(Edge), 
				cudaMemcpyDeviceToHost);
	}
	
	for (int i = 0; i < graphEdgeNr; i++)
	{
		if (output[i].weight >= MIN_WEIGHT && output[i].weight != -1)
		{
			cout << output[i].nr1 << "; " << output[i].nr2 << "; "
					<< output[i].weight << "; "
					<< (output[i].compl1 ? "1" : "0") << ", "
					<< (output[i].compl2 ? "1" : "0") << ";\n";
		}
	}

	cudaFree(devInput);
	cudaFree(devOutput);
	cudaFreeHost(input);
	cudaFreeHost(output);
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
	
	std::vector<SeqPair> seqMap = sort_all_trigrams(fragments);
	std::sort(seqMap.begin(), seqMap.end(), comparator);
	fragments.clear();
	for (int i = 0; i < seqMap.size(); i++)
	{
		fragments.push_back(seqMap[i].value);
	}
	createGraph(fragments);
	
  return 0;
}