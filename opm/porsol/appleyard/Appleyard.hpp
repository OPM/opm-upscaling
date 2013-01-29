/*
  Copyright 2011 IRIS - International Research Institute of Stavanger.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef APPLEYARD_CONTROL_HEADER
#define APPLEYARD_CONTROL_HEADER
#include<vector>
#include<cmath>

namespace Opm 
{	
	class AppleyardControl
	{
	
	public:
		AppleyardControl(const std::vector<double> & sat0, const double timeStart, const double timeEnd)
		:dt_(timeEnd-timeStart)
		,dtMin_((timeEnd-timeStart)/128)
		,time_(timeStart)
		,timeEnd_(timeEnd)
		,saturation_(sat0)
		,saturationIntermediate_(sat0)
		,dSat_(sat0.size(),0.0)
		,dSatMax_(1.0)
		,dSatTolerance_(0.001)
		,isChopped_(false) // TODO: Use this in the diagnostics ...
		,isOk_(true)
		,hasConverged_(false)
		,itCnt_(0)
		,itCntMax_(20) 
		{
			assert(dt_ > 0.0);
		}
		
		// Specifies min and max allowed saturation for each phase: <satMin,satMax>
		//   (Black oil: Sw_co <= Sw <= 1-So_r     0 <= Sg <= 1-Sw_co     So <= 1-Sw_co)
		// Layout of the saturation vector: <begin,pastEnd,stride>
		//   (Two options: ...,s1_i,s2_i,s3_i,...   versus  s1_1,...,s1_N,s2_1,...,s2_N,s3_1,...,s3_N,...) 
		void phaseLayout(double satMin, double satMax, unsigned int iBegin, unsigned int iEnd, unsigned int iStride);
		
		void init(const Opm::parameter::ParameterGroup& param) {} // TODO: What need to be tunable from command line?
		void reset(const Opm::parameter::ParameterGroup& param) {} // TODO:  Needed for repeated use without destruct/construct ...
		
		bool continueUpdate();  // TODO: 
		                        //   Reduced timestep: Try to reach timeEnd via micro-steps or immediately 
		                        //     release control for a possible pressure update? Probably a user choice ...
		                        //   Detect when saturation chopping (and not an optimistic timestep) prevents convergence.
		                        //   Etc, etc ...
		void finalizeUpdate(std::vector<double> & sat);
		
		double dt() {return dt_;}
		double time() {return time_;}
		const std::vector<double> & saturation() {return saturation_;}
		const std::vector<double> & saturationIntermediate() {return saturationIntermediate_;}
		std::vector<double> & dSat() {return dSat_;}
		
		bool & isOk() {return isOk_;}
		bool & hasConverged() {return hasConverged_;}
		bool isChopped() {return isChopped_;}
	
	private:
		void update();
		
		double dt_;
		double dtMin_;
		double time_;
		double timeEnd_;
		std::vector<double> saturation_;
		std::vector<double> saturationIntermediate_;
		std::vector<double> dSat_;
		double dSatMax_;
		double dSatTolerance_;
		double isChopped_;
		bool isOk_;
		bool hasConverged_;
		int itCnt_;
		int itCntMax_;
			
		struct AppleyardUpdate 
		{
			AppleyardUpdate(double satMin, double satMax, unsigned int iBegin, unsigned int iEnd, unsigned int iStride=0)
			:satMin_(satMin)
			,satMax_(satMax)
			,iBegin_(iBegin)
			,iEnd_(iEnd)
			,iStride_(iStride)
			,satMaxChange_(0.2)
			,tol_(1.e-6)
			{}
		
			bool update(std::vector<double> & sat, std::vector<double> & dSat);

			double satMin_;
			double satMax_;
			unsigned int iBegin_;
			unsigned int iEnd_;
			unsigned int iStride_;
			double satMaxChange_;
			double tol_;
		
		}; // AppleyardUpdate
		
		std::vector<AppleyardUpdate> phases_;
		
	}; // AppleyardControl

	
	void AppleyardControl::phaseLayout(double satMin, double satMax, unsigned int iBegin, unsigned int iEnd, unsigned int iStride)
	{
		phases_.push_back(AppleyardUpdate(satMin, satMax, iBegin, iEnd, iStride));
	}

	void AppleyardControl::update()
	{
		std::vector<double>::iterator dSatMax = std::max_element(dSat_.begin(), dSat_.end());
		std::vector<double>::iterator dSatMin = std::min_element(dSat_.begin(), dSat_.end());		
		dSatMax_ = std::max(*dSatMax, -(*dSatMin));
		
		isChopped_ = false;
		for (unsigned int i=0; i<phases_.size(); ++i) {
			if (phases_[i].update(saturation_, dSat_))
				isChopped_ = true;
		}
	}
	
	bool AppleyardControl::continueUpdate()
	{
	
		bool reduceTimeStep = false;
		
		if (isOk_) {
			//Update saturation
			update();
			
			std::cout << "  dSatMax: " << dSatMax_ << std::endl;
			for (unsigned int i=0; i<dSat_.size(); ++i)
				std::cout << " +>  " << saturation_[i] << "    " << dSat_[i] << std::endl;
			
			// Diagnostics
			if (dSatMax_ <= dSatTolerance_ && fabs(time_+dt_-timeEnd_)<1.e-8) {
				hasConverged_ = true;
				time_+= dt_;
				std::cout << "--- AppleyardControl: finished, time=" << time_ << std::endl;
				return false;
			}
			else if (dSatMax_ > dSatTolerance_) {
				if (itCnt_ < itCntMax_) 
					++itCnt_;  // Continue iteration
				else 
					reduceTimeStep = true;  //Reduce time step due to convergence problem
			}
			else { // Completed partial time step
				time_ += dt_;
				itCnt_ = 0;
				copy(saturation_.begin(),saturation_.end(),saturationIntermediate_.begin());				
				std::cout << "--- AppleyardControl: partially completed, time=" << time_ << std::endl;
			}				
		}
		else 
			reduceTimeStep = true;  //Reduce time step due to breakdown of increment computation
		
		
		if (reduceTimeStep) {
			dt_*= 0.5;
			if (dt_ < dtMin_) {
				hasConverged_ = false;
				std::cerr << "### AppleyardControl: dt underflow,  dt=" << dt_ << "  time=" << time_ << "  ###" << std::endl;
				return false;
			}
			std::cout << "--- AppleyardControl: reduced stepsize dt=" << dt_ << std::endl;
			itCnt_ = 0;
			copy(saturationIntermediate_.begin(),saturationIntermediate_.end(),saturation_.begin());
		}
		
		return true;
		
	} // AppleyardControl::continueUpdate


	void AppleyardControl::finalizeUpdate(std::vector<double> & sat)
	{
		assert(hasConverged_ && saturation_.size() == sat.size());	
		copy(saturation_.begin(),saturation_.end(),sat.begin());
	}
	
		
	bool
	AppleyardControl::AppleyardUpdate::update(std::vector<double> & sat, std::vector<double> & dSat)
	{
		assert(sat.size() == dSat.size());
		
		bool isChopped = false;

		for (unsigned int i=iBegin_; i < iEnd_; i+=iStride_) {
			
			//First check that the update is not too big. 
			if (dSat[i] > satMaxChange_) 
				dSat[i] = satMaxChange_; 
			else if (dSat[i] < -satMaxChange_) 
				dSat[i] = -satMaxChange_;
			else
				isChopped = true;

			//Now do some critical checks and then perform the chop: 
			if ((sat[i]+dSat[i]) < satMin_) {
				sat[i] = satMin_;
				isChopped = true;
			}
			else if ((sat[i]+dSat[i]) > satMax_) {
				sat[i] = satMax_;
				isChopped = true;
			} 
			else if (sat[i] <= satMin_) { 
				if ((sat[i]+dSat[i]) > satMin_) 
					sat[i] = satMin_+tol_;
				else
					sat[i] += dSat[i];
				isChopped = true;
			} 
			else if ((sat[i] > satMin_)  && (sat[i]+dSat[i])<= satMin_) {
				sat[i] = satMin_+tol_;
				isChopped = true;
			}
			else 		 
				sat[i] += dSat[i];   //The "normal" update ...
		}
		return isChopped;
	} // AppleyardControl::AppleyardUpdate::update
	
	

	
} // namespace Opm




#endif // APPLEYARD_CONTROL_HEADER
