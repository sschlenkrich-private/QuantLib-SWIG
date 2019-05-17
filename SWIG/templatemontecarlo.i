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

#ifndef quantlib_template_montecarlo_i
#define quantlib_template_montecarlo_i

%include common.i
%include date.i
%include daycounters.i
%include types.i
%include currencies.i
%include observer.i
%include marketelements.i
%include interpolation.i
%include indexes.i


%{
#include <ql/experimental/templatemodels/montecarlo/montecarlomodells.hpp>

using QuantLib::RealStochasticProcess;
using QuantLib::RealMCSimulation;
using QuantLib::RealMCPayoff;
using QuantLib::RealMCBase;
using QuantLib::RealMCRates;
using QuantLib::RealAMCPricer;
%}

// we need to tell SWIG to export shared_ptr with new classes
%template(RealStochasticProcessBase) boost::shared_ptr<RealStochasticProcess>;
%template(RealMCSimulationBase) boost::shared_ptr<RealMCSimulation>;
%template(RealMCPricerBase) boost::shared_ptr<RealMCPayoff::Pricer>;


// we need to tell C++ that our new pointer-based classes are type names
%{
typedef boost::shared_ptr<RealMCSimulation> RealMCSimulationPtr;
typedef boost::shared_ptr<RealMCPayoff::Pricer> RealMCPricerPtr;
%}

// we use an object adapter pattern to extend base class interface by new underlying class methods 

%rename(RealMCSimulation) RealMCSimulationPtr;
class RealMCSimulationPtr : public boost::shared_ptr<RealMCSimulation> {
  public:
    %extend {
		RealMCSimulationPtr(
		    const boost::shared_ptr<RealStochasticProcess>  process,
			const std::vector<Real>&                        simTimes,
			const std::vector<Real>&                        obsTimes,
			Size                                            nPaths,
			BigNatural                                      seed = 1234,
			bool                                            richardsonExtrapolation = false,
			bool                                            timeInterpolation = true,
			bool                                            storeBrownians = false ) {
                return new RealMCSimulationPtr(
                    new RealMCSimulation(process,simTimes,obsTimes,nPaths,seed,
                            richardsonExtrapolation,timeInterpolation,storeBrownians));
        }
        void simulate() { boost::dynamic_pointer_cast<RealMCSimulation>(*self)->simulate(); }
        
		// inspectors

		const boost::shared_ptr<RealStochasticProcess> process()  { return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->process()  ;}
		const std::vector<Real>&                  simTimes() { return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->simTimes() ;}
		const std::vector<Real>&                  obsTimes() { return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->obsTimes() ;}
		const Size                                nPaths()   { return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->nPaths()   ;}

		const std::vector< std::vector<Real> >&   observedPath(const Size pathIdx) {
            return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->observedPath(pathIdx);
        }

		// find an actual state of the model for a given observation time
		// this method is used by a path object and the result is passed on to
		// the stochastic process for payoff evaluation
		const std::vector<Real> state(const Size pathIdx, const Time t) {
            return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->state(pathIdx,t);     
        }

		const std::vector< std::vector<Real> >& brownian(const Size pathIdx) {
            return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->brownian(pathIdx);                     
        }
        
        // manage MC adjusters
        
        void calculateNumeraireAdjuster(const std::vector<Real>&  numeraireObservTimes){
            boost::dynamic_pointer_cast<RealMCSimulation>(*self)->calculateNumeraireAdjuster(numeraireObservTimes);                             
        }
		
        void calculateZCBAdjuster( const std::vector<Real>& zcbObservTimes, const std::vector<Real>& zcbOffsetTimes ) {
            boost::dynamic_pointer_cast<RealMCSimulation>(*self)->calculateZCBAdjuster(zcbObservTimes,zcbOffsetTimes);                                     
        }
        
        void calculateAssetAdjuster( const std::vector<Real>& assetObservTimes, const std::vector<std::string>& aliases ) {
            boost::dynamic_pointer_cast<RealMCSimulation>(*self)->calculateAssetAdjuster(assetObservTimes,aliases);                                             
        }
        
        const std::vector<Real>& numeraireAdjuster() {
            return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->numeraireAdjuster();
        }

        const std::vector< std::vector<Real> >& zcbAdjuster() {
            return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->zcbAdjuster();
        }
        
        const std::vector<Real>& assetAdjuster(const std::string& alias) {
            return boost::dynamic_pointer_cast<RealMCSimulation>(*self)->assetAdjuster(alias);
        }
    }    
};

%ignore RealMCPayoff;
class RealMCPayoff {
public:
	Real observationTime();

	// calculate observation times recursively
    std::set<Time> observationTimes();
    
	// return a clone but with changed observation time; this effectively allows considering a payoff as an index
	boost::shared_ptr<RealMCPayoff> at(const Time t);
        
};

// we only use shared pointer on this class
%template(RealMCPayoff) boost::shared_ptr<RealMCPayoff>;

// we also want to use vectors of payoffs 
namespace std {
    %template(RealMCPayoffVector) vector<boost::shared_ptr<RealMCPayoff> >;
}


// generic pricer
%rename(RealMCPricer) RealMCPricerPtr;
class RealMCPricerPtr : public boost::shared_ptr<RealMCPayoff::Pricer> {
public:
    %extend {
		RealMCPricerPtr(const std::vector< boost::shared_ptr<RealMCPayoff> >&   payoffs,
				     const boost::shared_ptr<RealMCSimulation>&              simulation) {
            return new RealMCPricerPtr(
                new RealMCPayoff::Pricer(payoffs,simulation));
        }

		Real NPV() { return boost::dynamic_pointer_cast<RealMCPayoff::Pricer>(*self)->NPV(); }
        
        static std::vector<Real> at(const boost::shared_ptr<RealMCPayoff>&      payoff,
			                        const boost::shared_ptr<RealMCSimulation>&  simulation) {
            return RealMCPayoff::Pricer::at(payoff,simulation);
        }
        
		static std::vector<Real> discountedAt(const boost::shared_ptr<RealMCPayoff>&      payoff,
				                              const boost::shared_ptr<RealMCSimulation>&  simulation) {
            return RealMCPayoff::Pricer::discountedAt(payoff,simulation);                                              
        }
                                              
		static Real NPV(const std::vector< boost::shared_ptr<RealMCPayoff> >&  payoffs,
				        const boost::shared_ptr<RealMCSimulation>&             simulation) {
            return RealMCPayoff::Pricer::NPV(payoffs,simulation);
        }
                        
		static std::vector<Real> NPVs(const std::vector< boost::shared_ptr<RealMCPayoff> >&  payoffs,
				                      const boost::shared_ptr<RealMCSimulation>&             simulation) {
            return RealMCPayoff::Pricer::NPVs(payoffs,simulation);
        }        
        
    }
};        


// basic payoffs and operations

// we need to tell C++ that our new pointer-based classes are type names
%{
typedef boost::shared_ptr<RealMCPayoff> RealMCClonePtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCFixedAmountPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCPayPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCCashPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCAssetPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCAssetBarrierNoHitPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCVanillaOptionPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCCachePtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCAxpyPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCMultPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCDivisionPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCMaxPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCMinPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCLogicalPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCIfThenElsePtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCBasketPtr;
%}

%rename(RealMCClone) RealMCClonePtr;
class RealMCClonePtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
	    RealMCClonePtr(const boost::shared_ptr<RealMCPayoff>&   x,
		               const Time                               observationTime) {
            return new RealMCClonePtr(new RealMCBase::Clone(x,observationTime));
        }
    }
};

// a deterministic flow known in advance (undiscounted)
%rename(RealMCFixedAmount) RealMCFixedAmountPtr;
class RealMCFixedAmountPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
	    RealMCFixedAmountPtr(const Real amount) {
            return new RealMCFixedAmountPtr(new RealMCBase::FixedAmount(amount));
        }
    }
};

// (re-)set paydate of a payoff (for discounting)
%rename(RealMCPay) RealMCPayPtr;
class RealMCPayPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
    	RealMCPayPtr(const boost::shared_ptr<RealMCPayoff>&  x,
			         const Time                              payTime) {
            return new RealMCPayPtr(new RealMCBase::Pay(x,payTime));
        }
    }
};

// simple discounted cash payment
%rename(RealMCCash) RealMCCashPtr;
class RealMCCashPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
		RealMCCashPtr(Time obsTime, Time payTime) {
            return new RealMCPayPtr(new RealMCBase::Cash(obsTime,payTime));
        } 
    }
};
	
// 1 unit of modelled asset
%rename(RealMCAsset) RealMCAssetPtr;
class RealMCAssetPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
		RealMCAssetPtr(Time obsTime, const std::string alias) {
            return new RealMCAssetPtr(new RealMCBase::Asset(obsTime,alias));
        }
    }
};

// return the continuous barrier no-hit probability
%rename(RealMCAssetBarrierNoHit) RealMCAssetBarrierNoHitPtr;
class RealMCAssetBarrierNoHitPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
		RealMCAssetBarrierNoHitPtr(Time tStart, Time tEnd, Real downBarrier, Real upBarrier, Real downOrUpOrBoth, const std::string alias) {
            return new RealMCAssetBarrierNoHitPtr(new RealMCBase::AssetBarrierNoHit(tStart,tEnd,downBarrier,upBarrier,downOrUpOrBoth,alias));
        }
    }
};

// 1 unit call or put exercised and settled at observation time
%rename(RealMCVanillaOption) RealMCVanillaOptionPtr;
class RealMCVanillaOptionPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
		RealMCVanillaOptionPtr(Time obsTime, const std::string alias, Real strike, Real callOrPut) {
            return new RealMCVanillaOptionPtr(new RealMCBase::VanillaOption(obsTime,alias,strike,callOrPut));
        }
    }
};

// cache result in case it is requested repeatedly
%rename(RealMCCache) RealMCCachePtr;
class RealMCCachePtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
		RealMCCachePtr(const boost::shared_ptr<RealMCPayoff>& x) {
            return new RealMCCachePtr(new RealMCBase::Cache(x));
        }
    }
};

// arithmetics and functions applied to payoffs

// a x + y  (undiscounted)
%rename(RealMCAxpy) RealMCAxpyPtr;
class RealMCAxpyPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCAxpyPtr(const Real                             a,
		              const boost::shared_ptr<RealMCPayoff>& x,
		              const boost::shared_ptr<RealMCPayoff>& y) {
            return new RealMCAxpyPtr(new RealMCBase::Axpy(a,x,y));                      
        }
    }    
};

// x * y  (undiscounted)		
%rename(RealMCMult) RealMCMultPtr;
class RealMCMultPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCMultPtr(const boost::shared_ptr<RealMCPayoff>& x,
		              const boost::shared_ptr<RealMCPayoff>& y) {
            return new RealMCMultPtr(new RealMCBase::Mult(x,y));                      
        }
    }    
};

// x / y  (undiscounted)		
%rename(RealMCDivision) RealMCDivisionPtr;
class RealMCDivisionPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCDivisionPtr(const boost::shared_ptr<RealMCPayoff>& x,
		                  const boost::shared_ptr<RealMCPayoff>& y) {
            return new RealMCDivisionPtr(new RealMCBase::Division(x,y));                      
        }
    }    
};

// max{x,y}  (undiscounted)
%rename(RealMCMax) RealMCMaxPtr;
class RealMCMaxPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCMaxPtr(const boost::shared_ptr<RealMCPayoff>& x,
		             const boost::shared_ptr<RealMCPayoff>& y) {
            return new RealMCMaxPtr(new RealMCBase::Max(x,y));                      
        }
    }    
};

// min{x,y}  (undiscounted)
%rename(RealMCMin) RealMCMinPtr;
class RealMCMinPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCMinPtr(const boost::shared_ptr<RealMCPayoff>& x,
		             const boost::shared_ptr<RealMCPayoff>& y) {
            return new RealMCMinPtr(new RealMCBase::Min(x,y));                      
        }
    }    
};

// logical operators
%rename(RealMCLogical) RealMCLogicalPtr;
class RealMCLogicalPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCLogicalPtr(const boost::shared_ptr<RealMCPayoff>& x,
		                 const boost::shared_ptr<RealMCPayoff>& y,
                         const std::string&                     op) {
            return new RealMCLogicalPtr(new RealMCBase::Logical(x,y,op));                      
        }
    }    
};

// if-then-else 
%rename(RealMCIfThenElse) RealMCIfThenElsePtr;
class RealMCIfThenElsePtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCIfThenElsePtr(
		              const boost::shared_ptr<RealMCPayoff>& x,
		              const boost::shared_ptr<RealMCPayoff>& y,
                      const boost::shared_ptr<RealMCPayoff>& z) {
            return new RealMCIfThenElsePtr(new RealMCBase::IfThenElse(x,y,z));                      
        }
    }    
};

// basket of underlyings
%rename(RealMCBasket) RealMCBasketPtr;
class RealMCBasketPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCBasketPtr(
            const std::vector<boost::shared_ptr<RealMCPayoff> >& underlyings,
		    const std::vector<Real>                              weights,
				bool                                             rainbow) {
            return new RealMCBasketPtr(new RealMCBase::Basket(underlyings,weights,rainbow));                      
        }
    }    
};


// interest rate payoffs

// we need to tell C++ that our new pointer-based classes are type names
%{
typedef boost::shared_ptr<RealMCPayoff> RealMCSwaptionPtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCSwapRatePtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCLiborRatePtr;
typedef boost::shared_ptr<RealMCPayoff> RealMCAnnuityPtr;
typedef boost::shared_ptr<RealMCPayoff> RealAMCMaxPtr;
typedef boost::shared_ptr<RealMCPayoff> RealAMCMinPtr;

%}


// general swaption instrument
%rename(RealMCSwaption) RealMCSwaptionPtr;
class RealMCSwaptionPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCSwaptionPtr(
            Time                      obsTime,    // observation equals fixing time
			const std::vector<Time>&  floatTimes,
			const std::vector<Real>&  floatWeights,
			const std::vector<Time>&  fixedTimes,
			const std::vector<Real>&  annuityWeights,
			Real                      strikeRate,
			Real                      payOrRec ) {
            return new RealMCSwaptionPtr(new RealMCRates::GeneralSwaption(
                obsTime,floatTimes,floatWeights,fixedTimes,annuityWeights,strikeRate,payOrRec));                      
        }
        RealMCSwaptionPtr(
            Time                               obsTime,    // observation equals fixing time
			const SwapIndexPtr&                index,
			const Handle<YieldTermStructure>&  discYTSH,
			Real                               strikeRate,
			Real                               payOrRec ) {
            boost::shared_ptr<SwapIndex> swapIdx =  boost::dynamic_pointer_cast<SwapIndex>(index);
            return new RealMCSwaptionPtr(new RealMCRates::GeneralSwaption(
                obsTime,swapIdx,discYTSH,strikeRate,payOrRec));
        }
    }    
};

// general swap (or CMS) rate 
%rename(RealMCSwapRate) RealMCSwapRatePtr;
class RealMCSwapRatePtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCSwapRatePtr(
            Time                      obsTime,    // observation equals fixing time
			const std::vector<Time>&  floatTimes,
			const std::vector<Real>&  floatWeights,
			const std::vector<Time>&  fixedTimes,
			const std::vector<Real>&  annuityWeights ) {
            return new RealMCSwapRatePtr(new RealMCRates::SwapRate(
                obsTime,floatTimes,floatWeights,fixedTimes,annuityWeights));                      
        }
        RealMCSwapRatePtr(
            Time                               obsTime,    // observation equals fixing time
			const SwapIndexPtr&                index,
			const Handle<YieldTermStructure>&  discYTSH ) {
            boost::shared_ptr<SwapIndex> swapIdx =  boost::dynamic_pointer_cast<SwapIndex>(index);
            return new RealMCSwapRatePtr(new RealMCRates::SwapRate(obsTime,swapIdx,discYTSH));
        }
    }    
};

// Libor rate 
%rename(RealMCLiborRate) RealMCLiborRatePtr;
class RealMCLiborRatePtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCLiborRatePtr(
            Time                               obsTime,    // observation equals fixing time
			const IborIndexPtr&                index,
			const Handle<YieldTermStructure>&  discYTSH ) {
            boost::shared_ptr<IborIndex> iborIdx =  boost::dynamic_pointer_cast<IborIndex>(index);
            return new RealMCLiborRatePtr(new RealMCRates::LiborRate(obsTime,iborIdx,discYTSH));
        }
    }    
};

// annuity (of a swap rate)
%rename(RealMCAnnuity) RealMCAnnuityPtr;
class RealMCAnnuityPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealMCAnnuityPtr(
            Time                      obsTime,    // observation equals fixing time
			const std::vector<Time>&  fixedTimes,
			const std::vector<Real>&  annuityWeights ) {
            return new RealMCAnnuityPtr(new RealMCRates::Annuity(obsTime,fixedTimes,annuityWeights));                      
        }
    }    
};

// American Monte Carlo payoffs using regression

// maximum of sum of payoffs discounted to obsTime
%rename(RealAMCMax) RealAMCMaxPtr;
class RealAMCMaxPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealAMCMaxPtr(
            const std::vector<boost::shared_ptr<RealMCPayoff>>& x,
            const std::vector<boost::shared_ptr<RealMCPayoff>>& y,
            const std::vector<boost::shared_ptr<RealMCPayoff>>& z,
            const Time                                          observationTime,
            const boost::shared_ptr<RealMCSimulation>           simulation,
            const Size                                          maxPolynDegree) {
            return new RealAMCMaxPtr(new RealAMCPricer::MinMax(
                x,y,z,observationTime,1.0,simulation,maxPolynDegree)); // +1 = max{}           
        }
    }    
};

// minimum of sum of payoffs discounted to obsTime
%rename(RealAMCMin) RealAMCMinPtr;
class RealAMCMinPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend {
        RealAMCMinPtr(
            const std::vector<boost::shared_ptr<RealMCPayoff>>& x,
            const std::vector<boost::shared_ptr<RealMCPayoff>>& y,
            const std::vector<boost::shared_ptr<RealMCPayoff>>& z,
            const Time                                          observationTime,
            const boost::shared_ptr<RealMCSimulation>           simulation,
            const Size                                          maxPolynDegree) {
            return new RealAMCMinPtr(new RealAMCPricer::MinMax(
                x,y,z,observationTime,-1.0,simulation,maxPolynDegree)); // -1 = min{}           
        }
    }    
};


// template
//%rename(RealMC###) RealMC###Ptr;
//class RealMC###Ptr : public boost::shared_ptr<RealMCPayoff> {
//public:
//    %extend {
//        RealMC###Ptr(const Real                             a,
//		              const boost::shared_ptr<RealMCPayoff>& x,
//		              const boost::shared_ptr<RealMCPayoff>& y) {
//            return new RealMC###Ptr(new RealMCPayoff::###( ));                      
//        }
//    }    
//};




#endif
