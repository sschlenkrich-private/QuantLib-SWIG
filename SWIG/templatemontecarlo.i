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
%}

%rename(RealMCClone) RealMCClonePtr;
class RealMCClonePtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
	    RealMCClonePtr(const boost::shared_ptr<RealMCPayoff>&   x,
		               const Time                               observationTime) {
            return new RealMCClonePtr(new RealMCPayoff::Clone(x,observationTime));
        }
    }
};

// a deterministic flow known in advance (undiscounted)
%rename(RealMCFixedAmount) RealMCFixedAmountPtr;
class RealMCFixedAmountPtr : public boost::shared_ptr<RealMCPayoff> {
public:
    %extend{
	    RealMCFixedAmountPtr(const Real amount) {
            return new RealMCFixedAmountPtr(new RealMCPayoff::FixedAmount(amount));
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
            return new RealMCPayPtr(new RealMCPayoff::Pay(x,payTime));
        }
    }
};





#endif
