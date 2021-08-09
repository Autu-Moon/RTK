#include"SatNavBasic.h"
const std::string systemString[SYScount] = {
	"GPS","BDS"
};
SYSTEM str2SYSTEM(const std::string& str) {
	for (size_t i = 0; i < SYScount; i++) {
		if (str == systemString[i]) {
			return static_cast<SYSTEM>(i);
		}
	}
	return static_cast<SYSTEM>(0);
}
const std::string freqbandString[FreqBandcount] = {
	"L1",
	"L2",
	"L5",
	"B1I",
	"B1C",
	"B2a",
	"B2b",
	"B3I"
};

const std::string obscodeString[ObsCodecount] = {
	"unknown",
	"C",   "L",   "D",   "S",
	"C1C", "L1C", "D1C", "S1C",
	"C1L", "L1L", "D1L", "S1L",
	"C1P", "L1P", "D1P", "S1P",
	"C1I", "L1I", "D1I", "S1I",
	"C2I", "L2I", "D2I", "S2I",
	"C2W", "L2W", "D2W", "S2W",
	"C2S", "L2S", "D2S", "S2S",

	"C5Q", "L5Q", "D5Q", "S5Q",
	"C5P", "L5P", "D5P", "S5P",

	"C6I", "L6I", "D6I", "S6I",

	"C7I", "L7I", "D7I", "S7I"
};

ObsCode str2ObsCode(const std::string& str) {
	for (size_t i = 0; i < ObsCodecount; i++) {
		if (str == obscodeString[i]) {
			return static_cast<ObsCode>(i);
		}
	}
	return static_cast<ObsCode>(0);
}

size_t getLocation(const ObsCode& code, const SYSTEM& sys) {
	std::string codeS = obscodeString[code];
	size_t fband;
	fband = (size_t)stoi(codeS.substr(1, 1));
	if (sys == GPS) {
		if (fband == 5)			//L5
			fband = 2;
		else					//L1 L2
			fband = fband - 1;
	}
	else if (sys == BDS) {
		if (fband == 2)			//B1I
			fband = 0;
		else if (fband == 6)	//B3I
			fband = 1;
		//else if (fband == 7)	//B2
		//	fband = 2;
		else					//not support other freq
			return -1;
	}
	else						//not support other SYSTEM
		return -1;
	return fband;
}

size_t getLocation(const FreqBand& band) {
	switch (band)
	{
	case L1:
		return 0;
	case L2:
		return 1;
	case L5:
		return 2;
	case B1I:
		return 0;
	case B3I:
		return 1;
	default:
		return -1;
	}
}

FreqBand getFreqBand(const size_t& loc, const SYSTEM& sys) {
	if (sys == GPS) {
		switch (loc) {
		case 0:
			return L1;
		case 1:
			return L2;
		case 2:
			return L5;
		default:
			return FreqBandunknown;
		}
	}
	else if (sys == BDS) {
		switch (loc) {
		case 0:
			return B1I;
		case 1:
			return B3I;
		case 2:
			return B2b;
		default:
			return FreqBandunknown;
		}
	}
	else
		return FreqBandunknown;
}

double getFreq(FreqBand f) {
	switch (f) {
	case L1:
		return Freq_L1;
	case L2:
		return Freq_L2;
	case L5:
		return Freq_L5;
	case B1I:
		return Freq_B1I;
	case B1C:
		return Freq_B1C;
	case B2a:
		return Freq_B2a;
	case B2b:
		return Freq_B2b;
	case B3I:
		return Freq_B3I;
	default:
		return 0;
	}
}

double getFreq(const size_t& loc, const SYSTEM& sys) {
	return getFreq(getFreqBand(loc, sys));
}
