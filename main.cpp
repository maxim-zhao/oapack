/*
 * Copyright (c) 2020-2022 Eugene Larchenko <el6345@gmail.com>. All rights reserved.
 * See the attached LICENSE.txt file.
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h> // clock_t
#include <new> // std::bad_alloc
#include <algorithm> // min,max

// using namespace std;

// for sanity checks
#define throw_impossible() { throw std::exception(); }

// Do not increase this constant, as it may lead to ushort type overflow!
#define MAX_INPUT_SIZE 0x10000
#define MIN_INPUT_SIZE 1

typedef signed char sbyte;
typedef unsigned char byte;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef signed __int64 int64;

#define ARRAYLEN(a) (sizeof a / sizeof a[0])

int FindLongestMatch(const byte* a, int n);
int FindOptimalSolution();
int EmitCompressed();

void PrintVersion()
{
	printf("\n");
	printf("Optimal aPLib compressor");
	if (sizeof(void*) > 4) {
		printf(" (x64)\n");
	} else {
		printf(" (x86)\n");
	}
	printf("version 2022.06.09\n");
	printf("by Eugene Larchenko (https://gitlab.com/eugene77)\n");
	printf("\n");
}

void PrintUsage(bool verbose)
{
	printf("Usage:\n");
	printf("  oapack.exe [-b] <inputfilename> [<outputfilename>]\n");
	if (verbose) {
		printf("  (-b option enables back-to-front compression)\n");
	}
	printf("\n");
}

bool file_exists(const char* path)
{
	FILE* f = fopen(path, "r");
	if (f != NULL) {
		fclose(f);
	}
	return (f != NULL);
}

void reverse_array(byte* a, int count)
{
	for(int i=0, j=count-1; i < j; i++, j--) {
		auto t = a[i]; a[i] = a[j]; a[j] = t;
	}
}

byte* data; // input data
int size; // data size
byte packedData[MAX_INPUT_SIZE * 9/8 + 100];
int packedSize;
int packedBits; // compressed size in bits

bool reverseMode = false;

clock_t starttime, endtime;

bool inprogress = false;
void ReportPosition(int pos)
{
	if (inprogress) {
		printf("\r"); // move cursor back
	}
	printf("%d ", pos);
	inprogress = true;
}
void ReportDone()
{
	if (inprogress) {
		printf("\r"); // move cursor back
		printf("Finished        \n");
		inprogress = false;
	}
}
void ReportAborted()
{
	if (inprogress) {
		printf("\r\n");
		inprogress = false;
	}
}

int main(int argc, const char* argv[])
{
	int retCode = 10;
	try
	{
		PrintVersion();

		int a = 1;
		if (a < argc && strcmp(argv[a], "-b") == 0) {
			a++;
			reverseMode = true;
		}

		// make input file name
		if (a >= argc)
		{
			PrintUsage(true);
			throw 1;
		}
		const char* inputPath = argv[a++];
	
		// make output file name
		bool outputPathSpecified = (a < argc);
		char outputPath[1000];
		const char* s = (a < argc) ? argv[a] : inputPath;
		size_t sl = strlen(s);
		if (sl + 10 > ARRAYLEN(outputPath))
		{
			printf("Path is too long\n");
			throw 2;
		}
		strcpy(outputPath, s);
		if (!outputPathSpecified)
		{
			strcat(outputPath, ".ap");
		}


		data = (byte*)malloc(MAX_INPUT_SIZE + 1);
		if (!data) {
			throw std::bad_alloc();
		}

		FILE* fIn = fopen(inputPath, "rb");
		if (!fIn)
		{
			printf("Error opening input file %s\n", inputPath);
			throw 5;
		}
		else
		{
			size_t fsize = fread(data, 1, MAX_INPUT_SIZE + 1, fIn);
			fclose(fIn);
			if (fsize > MAX_INPUT_SIZE)
			{
				printf("Input file is too large. Max supported file size is %d bytes.\n", MAX_INPUT_SIZE);
				throw 3;
			}
			else if (fsize < 1)
			{
				printf("Input file is too small. Min supported file size is 1 byte.\n");
				throw 3;
			}
			else
			{
				size = (int)fsize;
				data = (byte*)realloc(data, fsize); // helps catching indexing out of data array bugs
				if (fsize && !data) {
					throw std::bad_alloc();
				}

				if (file_exists(outputPath))
				{
					// we don't want overwriting file. Don't waste time, abort now.
					printf("Error: output file already exists: %s\n", outputPath);
					throw 6;
				}

				printf("Using back-to-front compression: %s\n", reverseMode ? "Yes" : "No");

				printf("Compressing file: %s\n", inputPath);

				if (reverseMode)
				{
					reverse_array(data, size);
				}

				packedBits = FindOptimalSolution();
				//ReportDone();
				packedSize = EmitCompressed();

				if (reverseMode)
				{
					reverse_array(packedData, packedSize);
				}

				double duration = (double)(endtime - starttime) / CLOCKS_PER_SEC;
				printf("time: %.0f sec\n", duration);

				double ratio = (double)packedSize / std::max(size, 1);
				if (ratio > 1) {
					ratio = std::max(ratio, 1.001);
				}
				char* warning = (ratio >= 1) ? " (!)" : "";
				printf("compression: %d / %d = %.3f%s\n", packedSize, size, ratio, warning);

				printf("Writing compressed file: %s\n", outputPath);
				if (file_exists(outputPath))
				{
					// we don't want overwriting file
					printf("Error: output file already exists: %s\n", outputPath);
					throw 6;
				}
				else
				{
					FILE* fOut = fopen(outputPath, "wb");
					if (!fOut)
					{
						printf("Error writing output file\n");
						throw 6;
					}
					else
					{
						size_t written = fwrite(packedData, 1, packedSize, fOut);
						fclose(fOut);
						if (written != packedSize)
						{
							// delete incomplete compressed file
							remove(outputPath);
							printf("Error writing output file\n");
							throw 6;
						}

						printf("All OK\n");
						retCode = 0;
					}
				}
			}
		}
	}
	catch (int code)
	{
		ReportAborted();
		retCode = code;
	}
	catch (std::bad_alloc& e)
	{
		ReportAborted();
		printf("\n");
		printf("Error allocating memory (%s)\n", e.what());
		retCode = 8;
	}
	catch (std::exception& e)
	{
		ReportAborted();
		printf("\n");
		printf("Error: %s\n", e.what());
		retCode = 9;
	}

	//printf("\n");
	return retCode;
};


// computes floor(log2(x))
int log2i(int x)
{
	if (x < 1) throw;
	int r = 30;
	for(int b = 1<<30; (x&b)==0; b>>=1) {
		r--;
	}
	return r;
}

template<class T>
T* makearr(int cnt) {
	T* a = new T[cnt];
	if (!a) throw std::bad_alloc(); // doublecheck
	memset(a, 0, cnt*sizeof(a[0]));
	return a;
}


#define LWM_FALSE 2
#define LWM_TRUE 1

#define Put0 1
#define Put1byte 2
#define Match1byte 3
#define MatchShort 4
#define ReuseOffset 5
#define MatchLong 6

#pragma pack(push,1)
struct Op
{
	ushort Len;
	ushort Ofs;
	byte Type;

	Op() {};

	Op(byte type, int len, int o)
	{
		if (len <= 0) throw_impossible();
		if (o < 0 || o == 0 && type != Put0 && type != Put1byte) throw_impossible();

		if (len < 0 || len > 0xFFFF) throw_impossible();
		if (o < 0 || o > 0xFFFF) throw_impossible();
		Len = (ushort)len;
		Ofs = (ushort)o;
		Type = type;
	}
};
#pragma pack(pop)


int* dp2[MAX_INPUT_SIZE + 1]; // result for LWM=2=false
int* dp1[MAX_INPUT_SIZE + 1]; // result for LWM=1=true
Op* dp_op[2][MAX_INPUT_SIZE + 1]; // solution for lwm=1 and lwm=2

int FindOptimalSolution()
{
	// First, find longest possible match. 
	// This helps reducing memory requirements in most practical cases.
	int longestMatch = FindLongestMatch(data, size);
	longestMatch = std::max(1, longestMatch); // at least one position ahead is required for Put0/Put1Byte ops

	int N = size;
	if (N < 1) throw_impossible();

	sbyte* gammalen = new sbyte[N];
	for (int i = 0; i < N; i++)
		gammalen[i] = sbyte(2 * log2i(i + 2));

	int* matchLen = new int[N];

	memset(dp1, 0, sizeof(dp1));
	memset(dp2, 0, sizeof(dp2));
	memset(dp_op, 0, sizeof(dp_op));

	// Preallocate memory, so we don't waste time if there is no enough

	int64 memReq = 0;
	for (int pos = N - 1; pos >= 1; pos--) {
		if (pos + longestMatch >= N) {
			memReq += (int64)2 * pos * sizeof(int);
		}
		memReq += (int64)2 * pos * sizeof(Op);
	}
//	printf("Allocating %I64d MiB of memory...\n", memReq >> 20);

	for (int pos = N - 1; pos >= 1; pos--) {
		if (pos + longestMatch >= N) {
			dp2[pos] = makearr<int>(pos);
			dp1[pos] = makearr<int>(pos);
		}
		dp_op[0][pos] = makearr<Op>(pos);
		dp_op[1][pos] = makearr<Op>(pos);
	}
	dp2[N] = makearr<int>(N); // result for pos=N is 0
	dp1[N] = makearr<int>(N); // result for pos=N is 0

//	printf("Memory allocated successfully\n");

//	printf("Processing position:\n");
	starttime = endtime = clock();

//	double lastreport = 0;
	for (int pos = N - 1; pos >= 1; pos--)
	{
//		double now = clock() / (double)CLOCKS_PER_SEC;
//		if (now < lastreport || now > lastreport + 0.03) { // limit to 30 reports/sec
//			ReportPosition(reverseMode ? N-pos : pos);
//			lastreport = now;
//		}

		matchLen[pos] = -1;
		for (int hl = pos - 1; hl >= 0; hl--)
		{
			int l = 0;
			while (pos + l < N && data[pos + l] == data[hl + l]) {
				l++;
			}
			matchLen[hl] = l;
		}

		// last_offset will be in [0..pos) range, so need arrays of this size
		if (!dp2[pos]) dp2[pos] = makearr<int>(pos);
		if (!dp1[pos]) dp1[pos] = makearr<int>(pos);
		//dp_op[0][pos] = makearr<Op>(pos); // this was preallocated
		//dp_op[1][pos] = makearr<Op>(pos); // this was preallocated

		// Try short match: len=2..3, o=1..127; LWM <- true
		int shortmatch_best = INT_MAX;
		Op shortmatch_bestOp = Op();
		for (int o = 1; o <= 127 && o <= pos; o++)
			for (int l = 2; l <= 3; l++)
			{
				if (pos + l <= N && matchLen[pos - o] >= l)
				{
					int t = 3 + 8 + dp1[pos + l][o]; // LWM_TRUE
					if (t < shortmatch_best)
					{
						shortmatch_best = t; shortmatch_bestOp = Op(MatchShort, l, o);
					}
				}
			}

		// Try long match for lwm=1 & lwm=2
		int longmatch2_best = INT_MAX;
		Op longmatch2_bestOp = Op();
		int longmatch1_best = INT_MAX;
		Op longmatch1_bestOp = Op();
		{
			for (int o = 1; o <= pos; o++)
			{
				int f = o < 128 ? 2
					: o < 1280 ? 0
					: o < 32000 ? 1
					: 2;
				//if (matchLen[pos - o] > N - pos) throw;
				int maxl = matchLen[pos - o];
				for (int l = 2 + f; l <= maxl; l++)
				{
					int t = 2 + gammalen[(o >> 8) + 1 + 1 - 2] + 8 + gammalen[l - f - 2];
					t += dp1[pos + l][o]; // LWM_TRUE
					if (t < longmatch1_best)
					{
						longmatch1_best = t; longmatch1_bestOp = Op(MatchLong, l, o);
					}

					t += gammalen[(o >> 8) + 2 + 1 - 2] - gammalen[(o >> 8) + 1 + 1 - 2];
					if (t < longmatch2_best)
					{
						longmatch2_best = t; longmatch2_bestOp = Op(MatchLong, l, o);
					}
				}
			}
		}

		for (int last_offset = 0; last_offset < pos; last_offset++)
		{
			// Try put 1 byte
			int best1 = 1 + 8 + dp2[pos + 1][last_offset]; // LWM_FALSE
			Op best1Op = Op(Put1byte, 1, 0);

			// Try short match
			if (shortmatch_best < best1)
			{
				best1 = shortmatch_best; best1Op = shortmatch_bestOp;
			}

			// Try put zero byte
			if (data[pos] == 0)
			{
				int t = 3 + 4 + dp2[pos + 1][last_offset]; // LWM_FALSE
				if (t < best1)
				{
					best1 = t; best1Op = Op(Put0, 1, 0);
				}
			}

			// Try 1-byte match
			for (int o = 1; o <= 15 && o <= pos; o++)
				if (matchLen[pos - o] > 0)
				{
					int t = 3 + 4 + dp2[pos + 1][last_offset]; // LWM_FALSE
					if (t < best1)
					{
						best1 = t; best1Op = Op(Match1byte, 1, o);
						break; // there will be no better solution
					}
				}

			// So far solutions for LWM=1 and LWM=2 are the same
			int best2 = best1;
			Op best2Op = best1Op;

			// Try last_offset assuming LWM=2=false
			if (last_offset != 0)
			{
				int o = last_offset;
				//if (matchLen[pos - o] > N - pos) throw;
				int maxl = matchLen[pos - o];
				for (int l = 2; l <= maxl; l++)
				{
					int t = 2 + 2 + gammalen[l - 2];
					t += dp1[pos + l][o]; // LWM_TRUE
					if (t < best2)
					{
						best2 = t; best2Op = Op(ReuseOffset, l, o);
					}
				}
			}

			// Try long match
			if (longmatch1_best < best1)
			{
				best1 = longmatch1_best; best1Op = longmatch1_bestOp;
			}
			if (longmatch2_best < best2)
			{
				best2 = longmatch2_best; best2Op = longmatch2_bestOp;
			}

			dp2[pos][last_offset] = best2;
			dp1[pos][last_offset] = best1;

			dp_op[2-1][pos][last_offset] = best2Op;
			dp_op[1-1][pos][last_offset] = best1Op;

		} // last_offset

		if (pos + longestMatch <= N)
		{
			// we don't need results for these positions anymore, let's reuse memory
			delete[] dp1[pos + longestMatch]; dp1[pos + longestMatch] = NULL;
			delete[] dp2[pos + longestMatch]; dp2[pos + longestMatch] = NULL;
		}

	} // pos

	endtime = clock();

	// return compressed size in bits
	int reslen = dp2[1][0]; // LWM_FALSE
	reslen += 8; // first byte
	reslen += 3+8; // eos

	// We don't need dp1 and dp2 anymore
	for(int i=0; i<=N; i++)
	{
		if (dp1[i]) { delete[] dp1[i]; dp1[i] = NULL; }
		if (dp2[i]) { delete[] dp2[i]; dp2[i] = NULL; }
	}

	return reslen;
}

// Builds final compressed block using precalculations from dp_op array
int EmitCompressed()
{
	int rpos = 0;

	auto emitByte = [&rpos](byte b) {
		packedData[rpos++] = b;
	};

	int apos = -1;
	int acnt = 8;
	auto emitBit = [&apos, &acnt, &rpos](int bit) {
		if (acnt == 8) {
			apos = rpos;
			packedData[rpos++] = 0;
			acnt = 0;
		}
		packedData[apos] = packedData[apos] * 2 + (bit & 1);
		acnt++;
	};

	auto emitGamma = [emitBit](int x) {
		if (x < 2) throw_impossible();
		int b = 1 << 30;
		while ((x & b) == 0) {
			b >>= 1;
		}
		while ((b >>= 1) != 0) {
			emitBit((x & b) == 0 ? 0 : 1);
			emitBit(b == 1 ? 0 : 1);
		}
	};

	int pos = 0;
	int lwm = LWM_FALSE;
	int last_offset = 0;

	if (size < 1) throw_impossible();
	emitByte(data[pos++]); // first byte is simply copied

	while (pos != size)
	{
		if (pos >= size) throw_impossible(); // something is wrong
	
		Op op = dp_op[lwm-1][pos][last_offset];
		switch (op.Type)
		{
			case Put1byte: {
				if (op.Len != 1) throw_impossible();
				emitBit(0);
				emitByte(data[pos]);
				lwm = LWM_FALSE;
				break;
			}
			case Put0: {
				if (op.Len != 1) throw_impossible();
				if (data[pos] != 0) throw_impossible();
				emitBit(1);
				emitBit(1);
				emitBit(1);
				emitBit(0); emitBit(0); emitBit(0); emitBit(0);
				lwm = LWM_FALSE;
				break;
			}
			case Match1byte: {
				if (op.Len != 1) throw_impossible();
				if (op.Ofs < 1 || op.Ofs > 15) throw_impossible();
				emitBit(1);
				emitBit(1);
				emitBit(1);
				for (int i = 3; i >= 0; i--) emitBit(op.Ofs >> i & 1);
				lwm = LWM_FALSE;
				break;
			}
			case MatchShort: {
				if (op.Len < 2 || op.Len > 3) throw_impossible();
				if (op.Ofs < 1 || op.Ofs > 127) throw_impossible();
				emitBit(1);
				emitBit(1);
				emitBit(0);
				emitByte((op.Ofs * 2 + (op.Len - 2)));
				lwm = LWM_TRUE;
				last_offset = op.Ofs;
				break;
			}
			case ReuseOffset: {
				if (lwm != LWM_FALSE || last_offset <= 0) throw_impossible();
				if (op.Len < 2) throw_impossible();
				if (op.Ofs != last_offset) throw_impossible();
				emitBit(1);
				emitBit(0);
				emitBit(0); emitBit(0);
				emitGamma(op.Len);
				lwm = LWM_TRUE;
				break;
			}
			case MatchLong: {
				int f = op.Ofs < 128 ? 2
					: op.Ofs < 1280 ? 0
					: op.Ofs < 32000 ? 1
					: 2;
				if (op.Len - f < 2) throw_impossible();
				if (op.Ofs < 1) throw_impossible();
				emitBit(1);
				emitBit(0);
				emitGamma((op.Ofs >> 8) + lwm + 1);
				emitByte((byte)op.Ofs);
				emitGamma(op.Len - f);
				lwm = LWM_TRUE;
				last_offset = op.Ofs;
				break;
			}
			default: {
				throw_impossible();
			}
		}
		pos += op.Len;
	}

	if (pos != size) {
		throw_impossible();
	}

	// end of stream marker
	emitBit(1);
	emitBit(1);
	emitBit(0);
	emitByte(0);

	// final check
	if (packedBits != rpos * 8 - (8 - acnt)) {
		throw_impossible();
	}

	// finalize bit stream
	while (acnt != 8) {
		emitBit(0);
	}

	int packedSize = rpos;
	if (packedSize > ARRAYLEN(packedData)) {
		throw_impossible() ; // buffer overflow?
	}

	return packedSize;
}


// Find max L, 0<=L<=n-1, such that there exists i,j such that a{i..i+L-1} == a{j..j+L-1}
int FindLongestMatch(const byte* a, int n)
{
	if (n < 2) {
		return 0;
	}

	// Using binary search to find L.
	// Using sliding hash and a hashtable to find matches of given len.
	// Collisions are possible, in which case the result will be greater than it should, which is acceptable.

	const int M1 = 0x000ffffd; // some small prime
	const int M2 = 0x007ffff1; // some bigger prime
	if (n > M1 / 2) {
		throw exception("file is too large"); // need larger hashtable
	}
	int* hashtable = new int[M1];
	if (!hashtable) {
		throw std::bad_alloc();
	}
	int min = 0, max = n - 1;
	while (min < max)
	{
		int L = (min + max + 1) / 2;

		bool foundMatch = false;
		for (int j = 0; j < M1; j++) hashtable[j] = -1;
		int i;
		int hash1 = 0, c1 = 1;
		int hash2 = 0, c2 = 1;
		for (i = 0; i < L; i++)
		{
			hash1 = (hash1 * 256 + a[i]) % M1;
			c1 = (c1 * 256) % M1;
			hash2 = (hash2 * 256 + a[i]) % M2;
			c2 = (c2 * 256) % M2;
		}
		hashtable[hash1] = hash2;
		for (; i < n; i++)
		{
			hash1 = ((hash1 * 256 - c1 * a[i - L] + a[i]) % M1 + M1) % M1;
			hash2 = ((hash2 * 256 - c2 * a[i - L] + a[i]) % M2 + M2) % M2;

			int j;
			for (j = hash1; hashtable[j] >= 0; j = (j + 1) % M1) {
				if (hashtable[j] == hash2) {
					foundMatch = true;
					goto k0;
				}
			}

			hashtable[j] = hash2;
		}

	k0:
		if (foundMatch) {
			min = L;
		} else {
			max = L - 1;
		}
	}

	delete[] hashtable;

	if (min != max)
		throw_impossible();

	return min;
}
