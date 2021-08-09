#pragma once
#include<iostream>
#include"ObsStore.h"
#include"SatID.h"
#include<string>
#include<fstream>
using std::string;
class RxObsReader {
public:
	/*
	* @Description: read rinex head
	* @Param: fstream of obs file, obsStore
	* @return: bool
	*/
	virtual bool readHead(std::ifstream& fin, ObsStore& obs) = 0;
	/*
	* @Description: read rinex body for a single epoch
	* @Param: fstream of obs file, obsStore
	* @return: 0 for failed
	*/
	virtual int  readData(std::ifstream& fin, ObsStore& obs) = 0;
	/*
	* @Description: get expected observation type from user
	* @Param: list of string for observation type
	*/
	virtual void setExpectCode(SYSTEM sys, initializer_list<string> code) = 0;
	/*
	* @Description: store obs to ObsStore according to ExpectType
	*/
	virtual void toStored(ObsStore& obs) = 0;
};

class Rx304ObsReader:RxObsReader
{
	/*read obs data of 3.04 edition*/
public:
	//------------field-------------
	struct Observation {		//to store observation for a a epoch sat type
		double obs = OBSEMPTY;	//unit:m cycle
		short LLI = OBSEMPTY;
		short SSI = OBSEMPTY;
		short valid = 1;

	};
	const string TYPE		=	"SYS / # / OBS TYPES";
	const string Interval	=	"INTERVAL";
	const string EOH		=	"END OF HEADER";
	const string EpochMark	=	">";
	const short  MAXTYPENUM	=	25;							//maxium num of type in a single system 
	using ObsCodeMap		=	map<ObsCode, Observation>;
	using EpochObsMap		=	map<SatID, ObsCodeMap>;
	vector<ObsCode>	CodeVec[SYScount];							//typeID Vec for GPS
	vector<ObsCode>	expCodeVec[SYScount];						//expect typeID Vec for GPS

	/*
	* @Description: the main map to store all of obs data of a epoch 
	* map struct is SatID->ObsCode->Observation
	*/
	EpochObsMap		epobsMap;
	
	//---------read function--------
	/*
	* @Description: read rinex head
	* @Param: fstream of obs file, obsStore
	* @return: bool
	*/
	bool readHead(std::ifstream& fin, ObsStore& obs)override;
	/*
	* @Description: read rinex body for a single epoch
	* @Param: fstream of obs file, obsStore
	* @return: 0 for failed
	*/
	int  readData(std::ifstream& fin, ObsStore& obs)override;
	//---------set function---------
	void setExpectCode(SYSTEM sys,initializer_list<string> code)override;

	//-----------function-----------
	/*
	* @Description: add eph to the main map
	* @Param: eph to add
	*/
	void addObservation(const SatID& sat, const ObsCode& code, const Observation& obs);
	/*
	* @Description: turn obs valid in map to zero
	*/
	void invalid();
	void toStored(ObsStore& obsStore)override;

	void fillblank(string& str);
	//----------costructor----------
	Rx304ObsReader() {
		CodeVec[GPS].reserve(MAXTYPENUM);
		CodeVec[BDS].reserve(MAXTYPENUM);
	};
	~Rx304ObsReader() = default;
	
};

