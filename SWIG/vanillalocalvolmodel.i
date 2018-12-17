/*
 Copyright (C) 2018 Sebastian Schlenkrich

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#ifndef quantlib_vanillalocalvolmodel_i
#define quantlib_vanillalocalvolmodel_i

%include common.i
%include date.i
%include daycounters.i
%include types.i
%include currencies.i
%include observer.i
%include marketelements.i
%include interpolation.i
%include indexes.i
%include optimizers.i
%include options.i

%include volatilities.i


%{
#include <ql/experimental/vanillalocalvolmodel/vanillalocalvolmodel.hpp>
#include <ql/experimental/vanillalocalvolmodel/vanillalocalvolsmilesection.hpp>
#include <ql/experimental/vanillalocalvolmodel/vanillalocalvoltermstructures.hpp>

using QuantLib::VanillaLocalVolModel;
using QuantLib::VanillaLocalVolModelSmileSection;
using QuantLib::VanillaLocalVolSwaptionVTS;
%}

// we need to tell SWIG to export shared_ptr with new classes
%template(VanillaLocalVolModelBase) boost::shared_ptr<VanillaLocalVolModel>;
%template(VanillaLocalVolModelSmileSectionBase) boost::shared_ptr<VanillaLocalVolModelSmileSection>;

// we need to tell C++ that our new pointer-based classes are type names
%{
typedef boost::shared_ptr<VanillaLocalVolModel> VanillaLocalVolModelPtr;
typedef boost::shared_ptr<VanillaLocalVolModelSmileSection> VanillaLocalVolModelSmileSectionPtr;
typedef boost::shared_ptr<VanillaLocalVolSwaptionVTS> VanillaLocalVolSwaptionVTSPtr;
%}

// we use an object adapter pattern to extend base class interface by new underlying class methods 

%rename(VanillaLocalVolModel) VanillaLocalVolModelPtr;
class VanillaLocalVolModelPtr : public boost::shared_ptr<VanillaLocalVolModel> {
  public:
    %extend {
        VanillaLocalVolModelPtr(
			const Time                T,
			const Real                S0,
			const Real                sigmaATM,
			const std::vector<Real>&  Sp,
			const std::vector<Real>&  Sm,
			const std::vector<Real>&  Mp,
			const std::vector<Real>&  Mm,
			// controls for calibration
			const Size                maxCalibrationIters = 5,
			const Size                onlyForwardCalibrationIters = 0,
			const bool                adjustATMFlag = true,
			const bool                enableLogging = false,
			const bool                useInitialMu  = false,
			const Real                initialMu     = 0.0 ) {        
            return new VanillaLocalVolModelPtr(
                new VanillaLocalVolModel(T,S0,sigmaATM,Sp,Sm,Mp,Mm,maxCalibrationIters,onlyForwardCalibrationIters,
                                         adjustATMFlag,enableLogging,useInitialMu,initialMu));
        }

        // wrap C++ object back into SWIG object
        VanillaLocalVolModelPtr( const boost::shared_ptr<VanillaLocalVolModel>&  model ) {
            return new VanillaLocalVolModelPtr( model );
        }
        
		// inspectors
		const std::vector<std::string> logging() { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->logging()     ;}
		const Time timeToExpiry()        { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->timeToExpiry()        ;}
		const Real forward()             { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->forward()             ;}
		const Real sigmaATM()            { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->sigmaATM()            ;}
		const Real alpha()               { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->alpha()               ;}
		const Real mu()                  { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->mu()                  ;}
		const Real nu()                  { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->nu()                  ;}
		const Size maxCalibrationIters() { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->maxCalibrationIters() ;}
		const Size onlyForwardCalibrationIters() { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->onlyForwardCalibrationIters(); }
		const bool adjustATMFlag()       { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->adjustATMFlag()       ;}
		const bool enableLogging()       { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->enableLogging()       ;}
		const bool useInitialMu()        { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->useInitialMu()        ;}
		const Real initialMu()           { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->initialMu()           ;}
        
		// attributes in more convenient single-vector format
		const std::vector<Real> underlyingX()   { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->underlyingX()  ;}
		const std::vector<Real> underlyingS()   { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->underlyingS()  ;}
		const std::vector<Real> localVol()      { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->localVol()     ;}
		const std::vector<Real> localVolSlope() { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->localVolSlope();}
        
		// model function evaluations
		const Real localVol(Real S)    { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->localVol(S)    ;}
		const Real underlyingS(Real x) { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->underlyingS(x) ;}
        
		//calculating expectations - that is the actual purpose of that model
        
		// calculate the forward price of an OTM option
		const Real expectation(bool isRightWing, Real strike) { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->expectation(isRightWing, strike); }
        
		// calculate the forward price of an OTM power option with payoff 1_{S>K}(S-K)^2
		const Real variance(bool isRightWing, Real strike) { return boost::dynamic_pointer_cast<VanillaLocalVolModel>(*self)->variance(isRightWing, strike); }
    }    
};

%rename(VanillaLocalVolModelSmileSection) VanillaLocalVolModelSmileSectionPtr;
class VanillaLocalVolModelSmileSectionPtr : public boost::shared_ptr<SmileSection>,
                                            public boost::shared_ptr<VanillaLocalVolModelSmileSection> {
  public:
    %extend {
		VanillaLocalVolModelSmileSectionPtr(
			const Date&                                   expiryDate,
			const Rate&                                   forward,
			const std::vector<Rate>&                      relativeStrikes,
			const std::vector<Volatility>&                smileVolatilities,
			const Real                                    extrapolationRelativeStrike,
			const Real                                    extrapolationSlope,
			bool                                          vegaWeighted = false,
			const boost::shared_ptr<EndCriteria>&         endCriteria = boost::shared_ptr<EndCriteria>(new EndCriteria(100, 10, 1.0e-6, 1.0e-6, 1.0e-6)),
			const boost::shared_ptr<OptimizationMethod>&  method = boost::shared_ptr<OptimizationMethod>(new LevenbergMarquardt(1.0e-6, 1.0e-6, 1.0e-6)),
			const DayCounter&                             dc = Actual365Fixed(),
			const Date&                                   referenceDate = Date(),
			const VolatilityType                          type = Normal,
			const Rate                                    shift = 0.0,
			const boost::shared_ptr<VanillaLocalVolModel>&  model = boost::shared_ptr<VanillaLocalVolModel>(),
			const Real                                    minSlope = -3.0,   //  lower boundary for m in calibration
			const Real                                    maxSlope =  3.0,   //  upper boundary for m in calibration
			const Real                                    alpha = 1.0e-4) {  //  Tikhonov alpha
            return new VanillaLocalVolModelSmileSectionPtr(
                new VanillaLocalVolModelSmileSection(expiryDate,forward,relativeStrikes,smileVolatilities,
                    extrapolationRelativeStrike,extrapolationSlope,vegaWeighted,endCriteria,method,dc,
                    referenceDate,type,shift,model,minSlope,maxSlope,alpha) );
        }
        
		VanillaLocalVolModelSmileSectionPtr(
			const Date&                                   expiryDate,
			const Rate&                                   forward,
			const Volatility&                             atmVolatility,
			const boost::shared_ptr<VanillaLocalVolModelSmileSection>& smile1,
			const boost::shared_ptr<VanillaLocalVolModelSmileSection>& smile2,
			const Real&                                   rho,
			const bool                                    calcSimple = true,  // use only ATM vol for x-grid calculation
			const DayCounter&                             dc = Actual365Fixed(),
			const Date&                                   referenceDate = Date(),
			const VolatilityType                          type = Normal,
			const Rate                                    shift = 0.0) {
            return new VanillaLocalVolModelSmileSectionPtr(
                new VanillaLocalVolModelSmileSection(expiryDate,forward,atmVolatility,smile1,smile2,rho,
                    calcSimple,dc,referenceDate,type,shift) );            
        }

        const boost::shared_ptr<VanillaLocalVolModel>&  model() const { return boost::dynamic_pointer_cast<VanillaLocalVolModelSmileSection>(*self)->model(); }

    }    
};


namespace std {
    %template(VanillaLocalVolModelSmileSectionVector) vector<boost::shared_ptr<VanillaLocalVolModelSmileSection> >;
    %template(VanillaLocalVolModelSmileSectionVectorVector) vector<vector<boost::shared_ptr<VanillaLocalVolModelSmileSection> > >;
}


%rename(VanillaLocalVolSwaptionVTS) VanillaLocalVolSwaptionVTSPtr;
class VanillaLocalVolSwaptionVTSPtr : public boost::shared_ptr<SwaptionVolatilityStructure> {
  public:
    %extend {
		VanillaLocalVolSwaptionVTSPtr(
			const Handle<SwaptionVolatilityStructure>&                                              atmVolTS,
			const std::vector< std::vector< boost::shared_ptr<VanillaLocalVolModelSmileSection> > >&  smiles,
			const std::vector< Period >&                                                            swapTerms,
			const SwapIndexPtr&                                                                     index) {
            boost::shared_ptr<SwapIndex> swapIdx =  boost::dynamic_pointer_cast<SwapIndex>(index);
            return new VanillaLocalVolSwaptionVTSPtr(
                new VanillaLocalVolSwaptionVTS(atmVolTS,smiles,swapTerms,swapIdx) );
        }    
    }    
};


#endif
