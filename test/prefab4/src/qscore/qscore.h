// qscore.h

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>

#define	all(t, n)		(t *) allocmem((n)*sizeof(t))
#define	reall(p, t, n)		p = (t *) reallocmem(p, (n)*sizeof(t))
#define zero(p,	t, n)	memset(p, 0, (n)*sizeof(t))
void *allocmem(int bytes);
void freemem(void *p);
void *reallocmem(void *p, int bytes);

inline int iabs(int i)
	{
	return i >= 0 ? i : -i;
	}

class MSA;

const double dInsane = double(0xffffffff);
const unsigned uInsane = 987654321;

unsigned CharToLetter(char c);
char LetterToChar(unsigned Letter);
bool IsGap(char c);

void ComparePair(const MSA &msaTest, unsigned uTestSeqIndexA,
  unsigned uTestSeqIndexB, const MSA &msaRef, unsigned uRefSeqIndexA,
  unsigned uRefSeqIndexB, double *ptrdSP, double *ptrdPS, double *ptrdCS);

double ComparePairSP(const MSA &msaTest, unsigned uTestSeqIndexA,
  unsigned uTestSeqIndexB, const MSA &msaRef, unsigned uRefSeqIndexA,
  unsigned uRefSeqIndexB);

void ComparePairMap(const int iTestMapA[], const int iTestMapB[],
  const int iRefMapA[], const int iRefMapB[], int iLengthA, int iLengthB,
  double *ptrdSP, double *ptrdPS, double *ptrdCS);

double ComparePairMapSP(const int iTestMapA[], const int iTestMapB[],
  const int iRefMapA[], const int iRefMapB[], int iLengthA, int iLengthB);

double SumPairs(const int iMapRef[], const int iMapTest[], unsigned uLength);

double ClineShift(const int iTestMapA[], const int iRefMapA[], unsigned uLengthA,
  const int iTestMapB[], const int iRefMapB[], unsigned uLengthB, double dEpsilon = 0.2);

void MakePairMaps(const MSA &msaTest, unsigned uTestSeqIndexA, unsigned uTestSeqIndexB,
  const MSA &msaRef, unsigned uRefSeqIndexA, unsigned uRefSeqIndexB, int **ptriTestMapAr,
  int **ptriTestMapBr, int **ptriRefMapAr, int **ptriRefMapBr);

void Quit(const char *Format, ...);

FILE *OpenStdioFile(const char *FileName);
int GetFileSize(FILE *f);

void ParseOptions(int argc, char *argv[]);
bool FlagOpt(const char *Name);
const char *ValueOpt(const char *Name);
const char *RequiredValueOpt(const char *Name);
void Usage();

void CompareMSA(const MSA &msaTest, const MSA &msaRef, double *ptrdSP,
  double *ptrdPS, double *ptrdCS);
double TC(MSA &msaTest, MSA &msaRef);

void QScore();

#include "msa.h"
#include "seq.h"

extern const char *g_TestFileName;
extern const char *g_RefFileName;
