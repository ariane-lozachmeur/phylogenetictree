#ifndef	MSA_h
#define MSA_h

class TextFile;
class Seq;

enum ALPHABET
	{
	ALPHABET_Amino
	};

class MSA
	{
#ifdef	WIN32
	friend void MSA::CopyReversed();
#endif

public:
	MSA();
	virtual ~MSA();

public:
	void FromFASTAFile(FILE *f);

	void SetSeqCount(unsigned uSeqCount);
	bool IsGap(unsigned uSeqIndex, unsigned uColIndex) const;

	void SetChar(unsigned uSeqIndex, unsigned uColIndex, char c);
	char GetChar(unsigned uSeqIndex, unsigned uIndex) const;

	void SetSeqName(unsigned uSeqIndex, const char szName[]);
	const char *GetSeqName(unsigned uSeqIndex) const;

	void GetSeq(unsigned uSeqIndex, Seq &seq) const;
	bool GetSeqIndex(const char *ptrSeqName, unsigned *ptruSeqIndex) const;

	unsigned GetCharCount(unsigned uSeqIndex, unsigned uColIndex) const;
	const char *GetSeqBuffer(unsigned uSeqIndex) const;
	unsigned GetSeqLength(unsigned uSeqIndex) const;

	void GetPairMap(unsigned uSeqIndex1, unsigned uSeqIndex2, int iMap1[],
	  int iMap2[]) const;
	unsigned GetUngappedColIndex(unsigned uSeqIndex, unsigned uColIndex);
	unsigned GetGappedColIndex(unsigned uSeqIndex, unsigned uUngappedColIndex);

	void Free();
	void Clear()
		{
		Free();
		}
	unsigned GetSeqCount() const
		{
		return m_uSeqCount;
		}
	unsigned GetColCount() const
		{
		return m_uColCount;
		}
	char LetterToChar(unsigned uLetter) const
		{
		return ::LetterToChar(uLetter);
		}
	char CharToLetter(char c) const
		{
		return ::CharToLetter(c);
		}

private:
	unsigned m_uSeqCount;
	unsigned m_uColCount;
	unsigned m_uCacheSeqCount;
	char **m_szSeqs;
	char **m_szNames;
	char *m_SeqBuffer;
	unsigned **m_UngapMap;
	};

#endif	// MSA_h
