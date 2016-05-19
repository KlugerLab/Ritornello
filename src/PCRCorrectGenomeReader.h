/*
 * PCACorrectGenomeReader.h
 *
 *  Created on: Apr 25, 2016
 *      Author: Kelly P. Stanton
 */

#ifndef PCRCORRECTGENOMEREADER_H_
#define PCRCORRECTGENOMEREADER_H_

#include "BufferedGenomeReader.h"

class PCRCorrectGenomeReader: public BufferedGenomeReader {
public:
	PCRCorrectGenomeReader(long argWindowSize);
	virtual ~PCRCorrectGenomeReader();
	void init(const char* samFileName);
	int next();
	void incrementBuffer();
	long getPos();
	double* getPstrand();
	double* getMstrand();
	long getPstrandReads();
	long getMstrandReads();
	long correctedPReadsAt(long pos);
	long correctedMReadsAt(long pos);
	void onChrChange();
	static long bandwidth;
private:
	double* correctedPstrand;
	double* correctedMstrand;
	double* correctedPstrandBuffer;
	double* correctedMstrandBuffer;
	long _pstrandReadsCorrected;
	long _mstrandReadsCorrected;
};

#endif /* PCRCORRECTGENOMEREADER_H_ */
