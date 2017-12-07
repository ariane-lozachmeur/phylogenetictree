#include "qscore.h"

/***
Strip trailing whitespace from FASTA annotation and
truncate at first whitespace (-truncname option).
***/
static void FixLabel(char Label[])
	{
// Truncate trailing whitespace
	size_t n = strlen(Label) - 1;
	while (n >= 0 && isspace(Label[n]))
		--n;
	Label[n+1] = 0;

	if (FlagOpt("truncname"))
		{
		n = strcspn(Label, " \t\r\n");
		Label[n] = 0;
		}
	}

#define APPEND_CHAR(c)							\
	{											\
	if (Pos >= BufferSize)						\
		Quit("ReadMFA: buffer too small");		\
	Buffer[Pos++] = (c);						\
	}

#define APPEND_LABEL(c)							\
	{											\
	if (LabelLength >= LabelBufferLength-1)		\
		{										\
		LabelBufferLength += 128;				\
		char *NewLabel = all(char, LabelBufferLength); \
		memcpy(NewLabel, Label, LabelLength);	\
		freemem(Label);							\
		Label = NewLabel;						\
		}										\
	Label[LabelLength] = c;						\
	++LabelLength;								\
	}

#define APPEND_SEQ(Label, Start, Length)			\
	{												\
	if (0 == m_uSeqCount)							\
		m_uColCount = Length;						\
	else if (Length != m_uColCount)					\
		Quit("Sequence lengths differ %s=%d, %s=%d",\
		  m_szNames[0], m_uColCount, Label, Length);\
	if (m_uSeqCount >= m_uCacheSeqCount)			\
		{											\
		m_uCacheSeqCount += 128;					\
		char **Seqs = all(char *, m_uCacheSeqCount);\
		char **SeqNames = all(char *, m_uCacheSeqCount);\
		if (m_uSeqCount > 0)						\
			{										\
			memcpy(Seqs, m_szSeqs, m_uCacheSeqCount*sizeof(char **));\
			memcpy(SeqNames, m_szNames, m_uCacheSeqCount*sizeof(char **));\
			delete[] m_szSeqs;						\
			delete[] m_szNames;						\
			}										\
		m_szSeqs = Seqs;							\
		m_szNames = SeqNames;						\
		}											\
	m_szSeqs[m_uSeqCount] = m_SeqBuffer + Start;	\
	FixLabel(Label);								\
	m_szNames[m_uSeqCount] = strdup(Label);			\
	++m_uSeqCount;									\
	}

void MSA::FromFASTAFile(FILE *f)
	{
	Clear();

	rewind(f);
	int FileSize = GetFileSize(f);
	int BufferSize = FileSize;
	char *Buffer = all(char, BufferSize);
	m_SeqBuffer = Buffer;

	char prev_c = '\n';
	bool InLabel = false;
	int ContigFrom = 0;
	char *Label = 0;
	int LabelLength = 0;
	int LabelBufferLength = 0;
	int Pos = 0;
	int ContigStart = 0;

	for (;;)
		{
		int c = fgetc(f);
		if (EOF == c)
			{
			if (feof(f))
				break;
			Quit("Stream error");
			}
		if (InLabel)
			{
			if (c == '\r')
				continue;
			if ('\n' == c)
				{
				Label[LabelLength] = 0;
				InLabel = false;
				}
			else
				{
				APPEND_LABEL(c)
				}
			}
		else
			{
			if ('>' == c && '\n' == prev_c)
				{
				int ContigLength = Pos - ContigStart;
				if (ContigLength > 0)
					{
					APPEND_SEQ(Label, ContigStart, ContigLength)
					}

				ContigStart = Pos;
				InLabel = true;
				LabelLength = 0;
				}
			else if (!isspace(c))
				{
				APPEND_CHAR(c)
				}
			}
		prev_c = c;
		}

	int ContigLength = Pos - ContigStart;
	if (ContigLength > 0)
		{
		APPEND_SEQ(Label, ContigStart, ContigLength);
		}
	}
