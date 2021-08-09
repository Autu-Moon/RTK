#pragma once
#include"SatID.h"
#include"PVT.h"
#include<map>
#include<iostream>
#include<iomanip>
#include"ObsStore.h"
#include"EphStore.h"
class UserObsStore
{
	/*
	* @Description: to store data needed for user solving
	* now include: pos, vel, clock of sat and Observations
	*/
public:
	//------------field-------------
	struct UserObs {
		SatID sat;
		PVT pvt;
		double Tgd[EPHTGDLENGTH];
		double Psr[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };				//伪距		unit:m
		double Adr[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };				//载波相位	unit:cycle
		double Dop[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };				//多普勒
		double CNR[OBSERVATIONLENGTH] = { OBSEMPTY,OBSEMPTY,OBSEMPTY };				//载波信噪比
		double IFComb = OBSEMPTY;											//无电离层组合观测值
		double elev = 0;											//卫星高度角
		short best = 1;
		
	/*
	 * @Description: 判断是否缺少观测值 
	 * @return: bool
	 */		
		bool isLack() const
		{
			for (size_t i = 0; i < OBSOCCUPIED; i++) {
				if (Psr[i] == OBSEMPTY || Adr[i] == OBSEMPTY) {
					return true;
				}
			}
			return false;
		}
	};
	enum CombineObs {
		IonoFree = 0,
		GeometryFree,
		MW,
		WL,
		count
	};
	using UserObsMap = std::map<SatID, UserObs>;
	UserObsMap userobsMap;
	GPSweek epoch;
	double CutOffElev = 0;
	//---------set function---------
	/*
	* @Description: to store this epoch data in the userobsMap
	* @Param: ObsStore EphStore
	* @return: void
	*/
	void setStore(ObsStore& obsStore,const EphStore& ephStore);
	/*
	* @Description: main function to make combination obs
	* @Notice: used after setStore()
	* @Param: CombineObs
	* @return: void
	*/
	bool setComb(CombineObs Comb);
	/*
	* @Description: function to make combination IF obs
	* @Notice: used after setStore()
	* @return: void
	*/
	void setIFComb(initializer_list<FreqBand> _band);
	/*
	 * @Description: get elevating angle based on position
	 * @Param: position
	 * @return: void
	 */
	void setElev(const Triple& pos);
	/*
	 * @Description: set cut off elevating angle, delete data that below it
	 * @Param: cut off elevating angle in degree
	 * @return: void
	 */
	void setCutOffElev(double cutoff);
	/*
	 * @Description: set cut off elevating angle, delete data that below it
	 * @Param: cut off elevating angle in degree
	 * @return: void
	 */
	void setBanSat(initializer_list<SatID> ban);
	/*
	 * @Description: ban BDS sat
	 * @Param:
	 * @return:
	 */
	void setBanBDSGEO();
	void setBanBDSIGSO();
	void setBanBDS2();
	void setBanBDS3();
	void setBanBDS();
	void setBanGPS();

	/*
	 * @Description: set hoped freqband,erase ones that don't have this band
	 * @Param: FreqBand
	 * @return: void
	 */
	void setFreq(initializer_list<FreqBand> band);
	//-----------function-----------
	/*
	* @Description: add data from a sat to the map
	* @Param: SatID,UserObs
	* @return: void
	*/
	void addUserObs(const SatID& sat, const UserObs& data);
	/*
	* @Description: cout all of the map to Console
	*/
	void dump(std::ostream &os)const;
	/*
	 * @Description: find the highest sat in elev of sys and don't lack of obs
	 * @Param: SYSTEM
	 * @return: SatID
	 */
	SatID getBest(SYSTEM sys)const;
	/*
	 * @Description: 
	 * @return: iterator of userobsMap
	 */
	UserObsMap::iterator begin() {
		return userobsMap.begin();
	}
	UserObsMap::iterator end() {
		return userobsMap.end();
	}
	UserObsMap::const_iterator cbegin()const {
		return userobsMap.cbegin();
	}
	UserObsMap::const_iterator cend()const {
		return userobsMap.cend();
	}
	/*
	* @Description: clear the map
	*/
	void clear();
	int size()const;
	UserObs operator[](const SatID& sat)const {
		return userobsMap.at(sat);
	}
	//----------constructor----------
	UserObsStore()=default;
	~UserObsStore()=default;
};

