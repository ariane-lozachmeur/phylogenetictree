#include "qscore.h"

void CompareMSA(const MSA &msaTest, const MSA &msaRef, double *ptrdSP,
  double *ptrdPS, double *ptrdCS)
	{
	const unsigned uRefSeqCount = msaRef.GetSeqCount();

	double dTotalSP = 0.0;
	double dTotalPS = 0.0;
	double dTotalCS = 0.0;
	unsigned uPairCount = 0;

	for (unsigned uRefSeqIndexA = 0; uRefSeqIndexA < uRefSeqCount; ++uRefSeqIndexA)
		{
		const char *pstrSeqNameA = msaRef.GetSeqName(uRefSeqIndexA);
		unsigned uTestSeqIndexA;
		bool bFound = msaTest.GetSeqIndex(pstrSeqNameA, &uTestSeqIndexA);
		if (!bFound)
			{
			Quit("Sequence '%s' not found in test alignment", pstrSeqNameA);
			continue;
			}

		for (unsigned uRefSeqIndexB = uRefSeqIndexA + 1; uRefSeqIndexB < uRefSeqCount;
		  ++uRefSeqIndexB)
			{
			const char *pstrSeqNameB = msaRef.GetSeqName(uRefSeqIndexB);
			unsigned uTestSeqIndexB;
			bool bFound = msaTest.GetSeqIndex(pstrSeqNameB, &uTestSeqIndexB);
			if (!bFound)
				{
				Quit("Sequence '%s' not found in test alignment", pstrSeqNameA);
				continue;
				}

			double dSP = dInsane;
			double dPS = dInsane;
			double dCS = dInsane;
			ComparePair(msaTest, uTestSeqIndexA, uTestSeqIndexB, msaRef, uRefSeqIndexA,
			  uRefSeqIndexB, &dSP, &dPS, &dCS);

			dTotalSP += dSP;
			dTotalPS += dPS;
			dTotalCS += dCS;
			++uPairCount;
			}
		}
	if (0 == uPairCount)
		{
		Quit("No sequence pairs in common between test and reference alignment");
		*ptrdSP = 0;
		*ptrdPS = 0;
		*ptrdCS = 0;
		return;
		}

	*ptrdSP = dTotalSP / uPairCount;
	*ptrdPS = dTotalPS / uPairCount;
	*ptrdCS = dTotalCS / uPairCount;
	}
