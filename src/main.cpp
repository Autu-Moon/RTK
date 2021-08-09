#include<iostream>
#include<iomanip>
#include"MyMatrix.h"
#include"TimeTrans.h"
#include"CoorTrans.h"
#include"SPP.h"
#include"RxObsReader.h"
#include"RxNavReader.h"
#include"UserObsStore.h"
#include<stdlib.h>
#include<string>
#include"RTK.h"
using namespace std;

int main()
{
	//RTK
	Rx304NavReader rxnav;
	Rx304ObsReader Rrxobs;
	Rx304ObsReader Brxobs;
	EphStore ephStore;
	ObsStore RobsStore;
	ObsStore BobsStore;
	UserObsStore RusrStore;
	UserObsStore BusrStore;

	SPP Rspp;
	RTK rtk;
	/*NovAtel 零基线*/
	//ifstream finnav("E:\\Data\\RTKdesignData\\Zero_OEM7_20210412\\brdm1020.21p");
	//ifstream finRover("E:\\Data\\RTKdesignData\\Zero_OEM7_20210412\\202104122230-7200.21O");
	//ifstream finBase("E:\\Data\\RTKdesignData\\Zero_OEM7_20210412\\202104122230-7190.21O");
	//Triple baseCoor(-2267810.196, 5009356.572, 3221000.818);		//位置真值
	//rtk.setRealRecPos(Triple(-2267810.196, 5009356.572, 3221000.818));

	/*Trimble 零基线*/
	/*ifstream finnav("E:\\Data\\RTKdesignData\\Zero_TrimbleAlloy\\brdm1730.20p");
	ifstream finRover("E:\\Data\\RTKdesignData\\Zero_TrimbleAlloy\\ROVE.20o");
	ifstream finBase("E:\\Data\\RTKdesignData\\Zero_TrimbleAlloy\\BASE.20o");
	Triple baseCoor(-2267826.102, 5009338.682, 3220976.493);		//位置真值
	rtk.setRealRecPos(Triple(-2267823.520, 5009337.993, 3220979.409));*/
	
	/*Trimble 短基线*/
	//ifstream finnav("E:\\Data\\RTKdesignData\\Short_TrimbleAlloy_20200621\\brdm1730.20p");
	//ifstream finRover("E:\\Data\\RTKdesignData\\Short_TrimbleAlloy_20200621\\ROVE.20o");
	//ifstream finBase("E:\\Data\\RTKdesignData\\Short_TrimbleAlloy_20200621\\BASE.20o");
	//Triple baseCoor(-2267826.102, 5009338.682, 3220976.493);		//位置真值
	//rtk.setRealRecPos(Triple(-2267823.520, 5009337.993, 3220979.409));
	
	/*NovAtel 0515短基线*/
	ifstream finnav("E:\\Data\\RTKdesignData\\Short_OEM7_20210515\\brdm1350.21p");
	ifstream finRover("E:\\Data\\RTKdesignData\\Short_OEM7_20210515\\CENT_20210515.21O");
	ifstream finBase("E:\\Data\\RTKdesignData\\Short_OEM7_20210515\\SGG_20210515.21O");
	Triple baseCoor(-2267810.196, 5009356.572, 3221000.818);		//位置真值
	rtk.setRealRecPos(Triple(-2267810.910, 5009435.734, 3220956.013));	//OEM short

	/*NovAtel 0514短基线*/
	//ifstream finnav("E:\\Data\\RTKdesignData\\Short_OEM7_20210514\\brdm1340.21p");
	//ifstream finRover("E:\\Data\\RTKdesignData\\Short_OEM7_20210514\\CENT_20210514.21O");
	//ifstream finBase("E:\\Data\\RTKdesignData\\Short_OEM7_20210514\\SGG_20210514.21O");
	//Triple baseCoor(-2267810.196, 5009356.572, 3221000.818);		//位置真值
	//rtk.setRealRecPos(Triple(-2267810.910, 5009435.734, 3220956.013));	//OEM short
	
	//ofstream fout("result\\RTKresult_OEM0515_Short_fixed.txt", ios::out);
	//ofstream fout("result\\RTKresult_OEM_Zero_fixed_3.txt", ios::out);
	ofstream fout("result\\RTKresult.txt", ios::out);

	//read rinex obs file head
	Rrxobs.readHead(finRover, RobsStore);
	Brxobs.readHead(finBase, BobsStore);
	//set obs code used
	Rrxobs.setExpectCode(GPS,{ "C1C", "C2W", "L1C", "L2W", "D1C", "D2W", "S1C", "S2W" });
	Rrxobs.setExpectCode(BDS,{ "C2I", "C6I", "L2I", "L6I", "D2I", "D6I", "S2I", "S6I" });
	Brxobs.setExpectCode(GPS,{ "C1C", "C2W", "L1C", "L2W", "D1C", "D2W", "S1C", "S2W" });
	Brxobs.setExpectCode(BDS,{ "C2I", "C6I", "L2I", "L6I", "D2I", "D6I", "S2I", "S6I" });
	//read rinex navigation file
	rxnav.read(finnav, ephStore);
	LOOPSTART:
	while (!finRover.eof() && !finBase.eof()) {
		//data read
		Rrxobs.readData(finRover, RobsStore);
		Brxobs.readData(finBase, BobsStore);
		//time match
		while (RobsStore.epoch != BobsStore.epoch) {
			if ((RobsStore.epoch - BobsStore.epoch).toSec() > 1e-8) {
				Brxobs.readData(finBase, BobsStore);
			}
			if ((RobsStore.epoch - BobsStore.epoch).toSec() < -1e-8) {
				Rrxobs.readData(finRover, RobsStore);
			}
			if (finRover.eof() || finBase.eof())
				goto LOOPSTART;
		}
		
		//prepare for obs
		RusrStore.setStore(RobsStore, ephStore);
		RusrStore.setFreq({ L1,L2,B1I,B3I });
		BusrStore.setFreq({ L1,L2,B1I,B3I });
		RusrStore.setIFComb({ L1, L2 });
		RusrStore.setIFComb({ B1I, B3I });
		BusrStore.setStore(BobsStore, ephStore);
		//RusrStore.setBanGPS();
		//RusrStore.setBanBDSGEO();
		//spp
		if (!Rspp.Solve(RusrStore)) {
			continue;
		}

		Triple Rpos = Rspp.getPos();
		RusrStore.setElev(Rpos);
		BusrStore.setElev(baseCoor);
		RusrStore.setCutOffElev(10.);
		BusrStore.setCutOffElev(10.);
		//rtk
		if (!rtk.Solve(BusrStore, RusrStore, cout, true, Rpos, baseCoor)) {
			continue;
		}
#ifdef DEBUG_SPP
		if (!Rspp.Solve(RusrStore, cout)) {
			fout << "ERROR" << endl;
		}
		RobsStore.dump(fout);
		RusrStore.dump(fout);
		Rspp.dumpResult(fout);
#endif //DEBUG_SPP
	}
	finnav.close();
	finBase.close();
	finRover.close();
	fout.close();
	system("pause");
	return 0;
}
