#include "qscore.h"

void ComparePair(const MSA &msaTest, unsigned uTestSeqIndexA,
  unsigned uTestSeqIndexB, const MSA &msaRef, unsigned uRefSeqIndexA,
  unsigned uRefSeqIndexB, double *ptrdSP, double *ptrdPS, double *ptrdCS)
	{
	const int iLengthA = (int) msaTest.GetSeqLength(uTestSeqIndexA);
	const int iLengthB = (int) msaTest.GetSeqLength(uTestSeqIndexB);
	const int iLengthAr = (int) msaRef.GetSeqLength(uRefSeqIndexA);
	const int iLengthBr = (int) msaRef.GetSeqLength(uRefSeqIndexB);

	if (iLengthA != iLengthAr || iLengthB != iLengthBr)
		Quit("Lengths differ");

	int *iRefMapA = new int[iLengthA];
	int *iRefMapB = new int[iLengthB];
	int *iTestMapA = new int[iLengthA];
	int *iTestMapB = new int [iLengthB];

	msaTest.GetPairMap(uTestSeqIndexA, uTestSeqIndexB, iTestMapA, iTestMapB);
	msaRef.GetPairMap(uRefSeqIndexA, uRefSeqIndexB, iRefMapA, iRefMapB);

	ComparePairMap(iTestMapA, iTestMapB, iRefMapA, iRefMapB, iLengthA, iLengthB,
	  ptrdSP, ptrdPS, ptrdCS);

	delete[] iRefMapA;
	delete[] iRefMapB;
	delete[] iTestMapA;
	delete[] iTestMapB;
	}
