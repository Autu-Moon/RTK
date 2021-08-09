#pragma once
#include<string>
const double V_light = 2.99792458e8;				//光速 unit:m/s

const double Freq_L1 = 1575.42e06;					//各载波频率 unit:Hz
const double Freq_L2 = 1227.6e06;
const double Freq_L5 = 1176.45e06;
const double Freq_B1I = 1561.098e06;			//BDS2 & BDS3
const double Freq_B1C = 1575.42e06;
const double Freq_B2a = 1176.45e06;
const double Freq_B2b = 1207.14e06;
const double Freq_B3I = 1268.52e06;

#define EPHTGDLENGTH 2
#define OBSERVATIONLENGTH 3
#define OBSOCCUPIED 2
#define OBSEMPTY -1
//#define DEBUG 1
enum SYSTEM
{
	GPS = 0,
	BDS,
	SYScount,
	SYSunknown
};
/*
 * @Description: trans SYSTEM to string
 */
extern const std::string systemString[SYScount];
/*
 * @Description: trans string to SYSTEM
 * @return: SYSTEM
 */
SYSTEM str2SYSTEM(const std::string& str);
enum FreqBand {
	L1 = 0,
	L2,
	L5,
	B1I,
	B1C,
	B2a,
	B2b,
	B3I,
	FreqBandcount,
	FreqBandunknown
};
extern const std::string freqbandString[FreqBandcount];

/*
 * @Description: get Frequency base on FreqBand
 * @return: frequency
 */
double getFreq(FreqBand f);
enum  ObsCode {
	ObsCodeunknown = 0,
	_C, _L, _D, _S,
	C1C, L1C, D1C, S1C,
	C1L, L1L, D1L, S1L,
	C1P, L1P, D1P, S1P,
	C1I, L1I, D1I, S1I,

	C2I, L2I, D2I, S2I,
	C2W, L2W, D2W, S2W,
	C2S, L2S, D2S, S2S,

	C5Q, L5Q, D5Q, S5Q,
	C5P, L5P, D5P, S5P,

	C6I, L6I, D6I, S6I,

	C7I, L7I, D7I, S7I,
	ObsCodecount
};
/*
 * @Description: trans ObsCode to string
 */
extern const std::string obscodeString[ObsCodecount];
/*
 * @Description: trans string to ObsCode
 */
ObsCode str2ObsCode(const std::string& str);
/*
 * @Description: get location in Obs[] based on ObsCode
 * @return: location
 */
size_t getLocation(const ObsCode& code, const SYSTEM& sys);
/*
 * @Description: get location in Obs[] based on ObsCode
 * @return: location
 */
size_t getLocation(const FreqBand& band);
/*
 * @Description: get Frequency in Obs[] based on ObsCode
 * @return: location
 */
double getFreq(const size_t& loc, const SYSTEM& sys);
/*
 * @Description: get Freqband in Obs[] based on ObsCode
 * @return: location
 */
FreqBand getFreqBand(const size_t& loc, const SYSTEM& sys);
