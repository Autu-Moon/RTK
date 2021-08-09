#pragma once
#include"EphStore.h"
#include<fstream>
class RxNavReader {
public:
	/*
	* @Description: main funtion to read rinex nav file
	* @Param: fstream of obs file, EphStore
	* @return: 0 for failed
	*/
	virtual int read(std::ifstream& fin, EphStore& ephStore) = 0;
	/*
	* @Description: read rinex head
	* @Param: fstream of obs file, EphStore
	* @return: bool
	*/
	virtual bool readHead(std::ifstream& fin, EphStore& ephStore) = 0;
	/*
	* @Description: read rinex body for a single epoch
	* @Param: fstream of obs file, EphStore
	* @return: 0 for failed
	*/
	virtual int  readData(std::ifstream& fin, EphStore& ephStore) = 0;
};

class Rx304NavReader :RxNavReader
{
	/*read obs data of 3.04 edition*/
public:
	//------------field-------------
	//---------read function--------
	/*
	* @Description: main funtion to read rinex nav file
	* @Param: fstream of obs file, EphStore
	* @return: 0 for failed
	*/
	int read(std::ifstream& fin, EphStore& ephStore)override;
	/*
	* @Description: read rinex head
	* @Param: fstream of obs file, EphStore
	* @return: bool
	*/
	bool readHead(std::ifstream& fin, EphStore& ephStore)override;
	/*
	* @Description: read rinex body for a single epoch
	* @Param: fstream of obs file, EphStore
	* @return: 0 for failed
	*/
	int  readData(std::ifstream& fin, EphStore& ephStore)override;
	//-----------function-----------
	//----------costructor----------
	Rx304NavReader() = default;
	~Rx304NavReader() = default;
};

class Rx302NavReader :RxNavReader
{
	/*read obs data of 3.04 edition*/
public:
	//------------field-------------
	//---------read function--------
	/*
	* @Description: main funtion to read rinex nav file
	* @Param: fstream of obs file, EphStore
	* @return: 0 for failed
	*/
	int read(std::ifstream& fin, EphStore& ephStore)override;
	/*
	* @Description: read rinex head
	* @Param: fstream of obs file, EphStore
	* @return: bool
	*/
	bool readHead(std::ifstream& fin, EphStore& ephStore)override;
	/*
	* @Description: read rinex body for a single epoch
	* @Param: fstream of obs file, EphStore
	* @return: 0 for failed
	*/
	int  readData(std::ifstream& fin, EphStore& ephStore)override;
	//-----------function-----------
	//----------costructor----------
	Rx302NavReader() = default;
	~Rx302NavReader() = default;
	//---------get function---------

};
