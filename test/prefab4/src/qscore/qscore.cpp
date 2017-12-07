#include "qscore.h"

const char *g_TestFileName;
const char *g_RefFileName;

static void ToUpper(MSA &msa)
	{
	const int SeqCount = msa.GetSeqCount();
	const int ColCount = msa.GetColCount();

	for (int SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		for (int ColIndex = 0; ColIndex < ColCount; ++ColIndex)
			{
			char c = msa.GetChar(SeqIndex, ColIndex);
			if (isalpha(c))
				{
				c = toupper(c);
				msa.SetChar(SeqIndex, ColIndex, c);
				}
			}
	}

void QScore()
	{
	g_TestFileName = RequiredValueOpt("test");
	g_RefFileName = RequiredValueOpt("ref");

	MSA msaTest;
	MSA msaRef;

	FILE *fTest = OpenStdioFile(g_TestFileName);
	FILE *fRef = OpenStdioFile(g_RefFileName);

	msaTest.FromFASTAFile(fTest);
	msaRef.FromFASTAFile(fRef);

	fclose(fTest);
	fclose(fRef);

	if (FlagOpt("ignoretestcase"))
		ToUpper(msaTest);

	if (FlagOpt("ignorerefcase"))
		ToUpper(msaRef);

	if (0 == msaTest.GetSeqCount())
		Quit("No seqs in test alignment");

	if (0 == msaRef.GetSeqCount())
		Quit("No seqs in ref alignment");

	double dSP = dInsane;
	double dPS = dInsane;
	double dCS = dInsane;

	CompareMSA(msaTest, msaRef, &dSP, &dPS, &dCS);

	double dTC = TC(msaTest, msaRef);

	printf("Test=%s;Ref=%s;Q=%.3g;Modeler=%.3g;Cline=%.3g;TC=%.3g\n",
	  g_TestFileName, g_RefFileName, dSP, dPS, dCS, dTC);
	}
