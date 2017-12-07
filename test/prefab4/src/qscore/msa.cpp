#include "qscore.h"

bool IsGap(char c)
	{
	return '-' == c || '~' == c || '.' == c;
	}

MSA::MSA()
	{
	m_uSeqCount = 0;
	m_uColCount = 0;
	m_uCacheSeqCount = 0;

	m_szSeqs = 0;
	m_szNames = 0;
	m_SeqBuffer = 0;
	m_UngapMap = 0;
	}

MSA::~MSA()
	{
	Free();
	}

void MSA::Free()
	{
	for (unsigned n = 0; n < m_uSeqCount; ++n)
		freemem(m_szNames[n]);

	if (0 != m_UngapMap)
		for (unsigned n = 0; n < m_uSeqCount; ++n)
			delete[] m_UngapMap[n];

	delete[] m_szSeqs;
	delete[] m_szNames;
	delete[] m_SeqBuffer;
	delete[] m_UngapMap;

	m_uSeqCount = 0;
	m_uColCount = 0;
	m_uCacheSeqCount = 0;

	m_SeqBuffer = 0;
	m_szSeqs = 0;
	m_szNames = 0;
	}

bool MSA::GetSeqIndex(const char *ptrSeqName, unsigned *ptruSeqIndex) const
	{
	for (unsigned uSeqIndex = 0; uSeqIndex < GetSeqCount(); ++uSeqIndex)
		if (0 == strcmp(ptrSeqName, GetSeqName(uSeqIndex)))
			{
			*ptruSeqIndex = uSeqIndex;
			return true;
			}
	return false;
	}

const char *MSA::GetSeqName(unsigned uSeqIndex) const
	{
	if (uSeqIndex >= m_uSeqCount)
		Quit("MSA::GetSeqName(%u), count=%u", uSeqIndex, m_uSeqCount);
	if ('>' == m_szNames[uSeqIndex][0])
		return &m_szNames[uSeqIndex][1];
	return m_szNames[uSeqIndex];
	}

/***
It is sometimes very convenient to represent a pairwise alignment
as a "pair map", which works as follows.

Let iPos1 be the index into ungapped sequence 1, similarly for iPos2.

Then if a pair of letters (iPos1, iPos2) is aligned:

	iMap1[iPos1] = iPos2 and iMap2[iPos2] = iPos1.

If iPos1 is not in an aligned column, or is aligned to a gap, then
iMap1[iPos1] = -1, and similarly for iMap2. This overloads the meaning
of the integer value, so is questionable software engineering practice;
however it's a simple and convenient solution for small applications.
***/

void MSA::GetPairMap(unsigned uSeqIndex1, unsigned uSeqIndex2, int iMap1[],
  int iMap2[]) const
	{
	assert(uSeqIndex1 < GetSeqCount());
	assert(uSeqIndex2 < GetSeqCount());

	int iPos1 = 0;
	int iPos2 = 0;
	const unsigned uColCount = GetColCount();
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		char c1 = GetChar(uSeqIndex1, uColIndex);
		char c2 = GetChar(uSeqIndex2, uColIndex);
		bool bIsGap1 = ::IsGap(c1);
		bool bIsGap2 = ::IsGap(c2);
		if (!bIsGap1 && !bIsGap2)
			{
			if (isupper(c1))
				{
				if (!isupper(c2))
					Quit("Both upper and lower case letters (%c,%c) in ref alignment column %d",
					  c1, c2, uColIndex);
				iMap1[iPos1] = iPos2;
				iMap2[iPos2] = iPos1;
				}
			else
				{
				iMap1[iPos1] = -1;
				iMap2[iPos2] = -1;
				}
			++iPos1;
			++iPos2;
			}
		else if (!bIsGap1 && bIsGap2)
			{
			iMap1[iPos1] = -1;
			++iPos1;
			}
		else if (bIsGap1 && !bIsGap2)
			{
			iMap2[iPos2] = -1;
			++iPos2;
			}
		}

#if	_DEBUG
	{
	int iLength1 = iPos1;
	int iLength2 = iPos2;

	for (int iPos1 = 0; iPos1 < iLength1; ++iPos1)
		{
		int iPos2 = iMap1[iPos1];
		if (-1 == iPos2)
			continue;
		assert(iMap2[iPos2] == iPos1);
		}

	for (int iPos2 = 0; iPos2 < iLength2; ++iPos2)
		{
		int iPos1 = iMap2[iPos2];
		if (-1 == iPos1)
			continue;
		assert(iMap1[iPos1] == iPos2);
		}
	}
#endif
	}

unsigned MSA::GetSeqLength(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < GetSeqCount());

	const unsigned uColCount = GetColCount();
	unsigned uLength = 0;
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		if (!IsGap(uSeqIndex, uColIndex))
			++uLength;
	return uLength;
	}

bool MSA::IsGap(unsigned uSeqIndex, unsigned uColIndex) const
	{
	return ::IsGap(GetChar(uSeqIndex, uColIndex));
	}

char MSA::GetChar(unsigned uSeqIndex, unsigned uIndex) const
	{
// TODO: Performance cost?
	if (uSeqIndex >= m_uSeqCount || uIndex >= m_uColCount)
		Quit("MSA::GetChar(%u/%u,%u/%u)", uSeqIndex, m_uSeqCount, uIndex, m_uColCount);

	char c = m_szSeqs[uSeqIndex][uIndex];
//	assert(IsLegalChar(c));
	return c;
	}

void MSA::SetChar(unsigned uSeqIndex, unsigned uIndex, char c)
	{
// TODO: Performance cost?
	if (uSeqIndex >= m_uSeqCount || uIndex >= m_uColCount)
		Quit("MSA::GetLetter(%u/%u,%u/%u)", uSeqIndex, m_uSeqCount, uIndex, m_uColCount);

	m_szSeqs[uSeqIndex][uIndex] = c;
	}

unsigned MSA::GetUngappedColIndex(unsigned uSeqIndex, unsigned uColIndex)
	{
	const unsigned uSeqCount = GetSeqCount();
	const unsigned uColCount = GetColCount();

	assert(uSeqIndex < uSeqCount);
	assert(uColIndex < uColCount);

	if (0 == m_UngapMap)
		{
		m_UngapMap = new unsigned *[uSeqCount];
		memset(m_UngapMap, 0, uSeqCount*sizeof(unsigned *));
		}

	unsigned *ptrMap = m_UngapMap[uSeqIndex];
	if (0 == ptrMap)
		{
		ptrMap = new unsigned[uColCount];
		memset(ptrMap, 0, uColCount*sizeof(unsigned));
		unsigned uUngappedColIndex = 0;
		for (unsigned uGappedColIndex = 0; uGappedColIndex < uColCount;
		  ++uGappedColIndex)
			if (IsGap(uSeqIndex, uGappedColIndex))
				ptrMap[uGappedColIndex] = uInsane;
			else
				{
				ptrMap[uGappedColIndex] = uUngappedColIndex;
				++uUngappedColIndex;
				}
		m_UngapMap[uSeqIndex] = ptrMap;
		}
	unsigned uUngappedColIndex = ptrMap[uColIndex];
	if (uInsane == uUngappedColIndex)
		Quit("GetUngappedColIndex(%u,%u)", uSeqIndex, uUngappedColIndex);
	return uUngappedColIndex;
	}

unsigned MSA::GetGappedColIndex(unsigned uSeqIndex, unsigned uUngappedColIndex)
	{
	unsigned n = 0;
	for (unsigned uColIndex = 0; uColIndex < GetColCount(); ++uColIndex)
		{
		if (!IsGap(uSeqIndex, uColIndex))
			{
			if (n == uUngappedColIndex)
				{
				assert(GetUngappedColIndex(uSeqIndex, uColIndex) == uUngappedColIndex);
				return uColIndex;
				}
			++n;
			}
		}
	assert(false);
	return 0;
	}
