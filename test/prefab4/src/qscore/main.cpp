#include "qscore.h"

int main(int argc, char *argv[])
	{
	if (argc < 2)
		{
		Usage();
		exit(0);
		}

	ParseOptions(argc - 1, argv + 1);
	if (FlagOpt("version"))
		{
		fprintf(stderr, "qscore v1.1\n");
		exit(0);
		}

	QScore();
	}
