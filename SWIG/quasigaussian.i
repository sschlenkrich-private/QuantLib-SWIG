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

#ifndef quantlib_quasigaussian_i
#define quantlib_quasigaussian_i

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
%include templatemontecarlo.i



%{
#include <ql/experimental/templatemodels/qgaussian2/qgcalibrator.hpp>

using QuantLib::QuasiGaussianModel;
using QuantLib::QGSwaprateModel;
using QuantLib::QGAverageSwaprateModel;
using QuantLib::QGCalibrator;

%}

// we need to tell SWIG to export shared_ptr with new classes
//%template(QuasiGaussianModelBase) boost::shared_ptr<QuasiGaussianModel>;
//%template(QGSwaprateModelBase) boost::shared_ptr<QGSwaprateModel>;
//%template(QGAverageSwaprateModelBase) boost::shared_ptr<QGAverageSwaprateModel>;
%template(QGCalibratorBase) boost::shared_ptr<QGCalibrator>;

// we need to tell C++ that our new pointer-based classes are type names
%{
typedef boost::shared_ptr<RealStochasticProcess> QuasiGaussianModelPtr;
typedef boost::shared_ptr<RealStochasticProcess> QGSwaprateModelPtr;
typedef boost::shared_ptr<RealStochasticProcess> QGAverageSwaprateModelPtr;
typedef boost::shared_ptr<QGCalibrator> QGCalibratorPtr;
%}

// we use an object adapter pattern to extend base class interface by new underlying class methods 

%rename(QuasiGaussianModel) QuasiGaussianModelPtr;
class QuasiGaussianModelPtr : public boost::shared_ptr<RealStochasticProcess> {
  public:
    %extend {
		QuasiGaussianModelPtr(
			const Handle<YieldTermStructure>&       termStructure,
		    // number of yield curve factors (excluding stoch. vol)
		    const Size                              d,       // (d+1)-dimensional Brownian motion for [x(t), z(t)]^T
		    // unique grid for time-dependent parameters
		    const std::vector<Time>&                times,   // time-grid of left-constant model parameter values
		    // time-dependent parameters, left-piecewise constant on times_-grid
		    const std::vector< std::vector<Real> >& sigma,   // volatility
		    const std::vector< std::vector<Real> >& slope,   // skew
			const std::vector< std::vector<Real> >& curve,   // smile
			const std::vector<Real>&                eta,     // vol-of-vol
		    // time-homogeneous parameters
		    const std::vector<Real>&                delta,   // maturity of benchmark rates f(t,t+delta_i) 		
		    const std::vector<Real>&                chi,     // mean reversions
		    const std::vector< std::vector<Real> >& Gamma,   // (benchmark rate) correlation matrix
			const Real                              theta) { // mean reversion speed
            return new QuasiGaussianModelPtr(
                new QuasiGaussianModel(termStructure,d,times,sigma,slope,curve,eta,delta,chi,Gamma,theta) );
        }
        
        // wrap C++ object back into SWIG object
        QuasiGaussianModelPtr( const boost::shared_ptr<RealStochasticProcess>&  model ) {
            boost::shared_ptr<QuasiGaussianModel> qgmodel = boost::dynamic_pointer_cast<QuasiGaussianModel>(model);
            QL_REQUIRE(qgmodel,"QuasiGaussianModel required");
            return new QuasiGaussianModelPtr( qgmodel );
        }

        
        // inspectors
		const Handle<YieldTermStructure> termStructure()  { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->termStructure() ;}
		const std::vector<Time>&                times()   { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->times()         ;}
		const std::vector< std::vector<Real> >& sigma()   { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->sigma()         ;}
		const std::vector< std::vector<Real> >& slope()   { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->slope()         ;}
		const std::vector< std::vector<Real> >& curve()   { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->curve()         ;}
		const std::vector<Time>&                eta()     { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->eta()           ;}
		const std::vector< std::vector<Real> >& DfT()     { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->DfT()           ;}
		const std::vector< std::vector<Real> >& HHfInv()  { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->HHfInv()        ;}
		const std::vector<Real>&                delta()   { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->delta()         ;}
		const std::vector<Real>&                chi()     { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->chi()           ;}
		const Real                              theta()   { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->theta()         ;}
		const Real                              z0()      { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->z0()            ;}
		const Real                              condHHf() { return boost::dynamic_pointer_cast<QuasiGaussianModel>(*self)->condHHf()       ;}

    }    
};

%rename(QGSwaprateModel) QGSwaprateModelPtr;
class QGSwaprateModelPtr : public boost::shared_ptr<RealStochasticProcess> {
  public:
    %extend {
		QGSwaprateModelPtr (
			const QuasiGaussianModelPtr& model,
			const std::vector<Real>&    floatTimes,      // T[1], ..., T[M]
			const std::vector<Real>&    floatWeights,    // u[1], ..., u[M]
			const std::vector<Real>&    fixedTimes,      // T[1], ..., T[N]
			const std::vector<Real>&    fixedWeights,    // w[1], ..., w[N-1]
			const std::vector<Real>&    modelTimes,      // time grid for numerical integration
			const bool                  useExpectedXY) { // evaluate E^A [ x(t) ], E^A [ y(t) ] as expansion points
            boost::shared_ptr<QuasiGaussianModel> qgmodel = boost::dynamic_pointer_cast<QuasiGaussianModel>(model);
            QL_REQUIRE(qgmodel,"QuasiGaussianModel required");
            return new QGSwaprateModelPtr(                
                new QGSwaprateModel(qgmodel,floatTimes,floatWeights,fixedTimes,fixedWeights,modelTimes,useExpectedXY) );
        }

        Real              annuity     ( const Real t, const std::vector<Real>& x, const std::vector< std::vector<Real> >& y) { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->annuity     (t,x,y);}
        Real              swapRate    ( const Real t, const std::vector<Real>& x, const std::vector< std::vector<Real> >& y) { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->swapRate    (t,x,y);}           
        std::vector<Real> swapGradient( const Real t, const std::vector<Real>& x, const std::vector< std::vector<Real> >& y) { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->swapGradient(t,x,y);}
        Real sigma  (const Real t) { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->sigma (t) ;}
        Real slope  (const Real t) { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->slope (t) ;}
        Real eta    (const Real t) { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->eta   (t) ;}
        Real theta  ()             { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->theta ()  ;}
        Real z0     ()             { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->z0    ()  ;}
        Real rho    ()             { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->rho   ()  ;}
        Real S0     ()             { return boost::dynamic_pointer_cast<QGSwaprateModel>(*self)->S0    ()  ;}
        
    }    
};

%rename(QGAverageSwaprateModel) QGAverageSwaprateModelPtr;
class QGAverageSwaprateModelPtr : public boost::shared_ptr<RealStochasticProcess> {
  public:
    %extend {
		QGAverageSwaprateModelPtr ( const QGSwaprateModelPtr& model ) {
            boost::shared_ptr<QGSwaprateModel> swapRateModel = boost::dynamic_pointer_cast<QGSwaprateModel>(model);
            QL_REQUIRE(swapRateModel,"QuasiGaussianModel required");            
            return new QGAverageSwaprateModelPtr( new QGAverageSwaprateModel(swapRateModel) );
        }        
        // inspectors
        Real sigma() { return boost::dynamic_pointer_cast<QGAverageSwaprateModel>(*self)->sigma() ;}
		Real slope() { return boost::dynamic_pointer_cast<QGAverageSwaprateModel>(*self)->slope() ;}
		Real eta()   { return boost::dynamic_pointer_cast<QGAverageSwaprateModel>(*self)->eta()   ;}

		// undiscounted expectation of vanilla payoff
		Real vanillaOption(const Real strike, const int callOrPut, const Real accuracy = 1.0-6, const Size maxEvaluations = 1000) {
            return boost::dynamic_pointer_cast<QGAverageSwaprateModel>(*self)->vanillaOption(strike,callOrPut,accuracy,maxEvaluations);
        }        
    }    
};

%rename(QGCalibrator) QGCalibratorPtr;
class QGCalibratorPtr : public boost::shared_ptr<QGCalibrator> {
  public:
    %extend {
        QGCalibratorPtr(
		    const QuasiGaussianModelPtr&                model,
			const Handle<SwaptionVolatilityStructure>&  volTS,
			const std::vector< SwapIndexPtr >&          swapIndices,
			const Real                                  modelTimesStepSize,
            const bool                                  useExpectedXY,
            const Real                                  sigmaMax,
            const Real                                  slopeMax,
            const Real                                  etaMax,
            const Real                                  sigmaWeight,
            const Real                                  slopeWeight,
            const Real                                  etaWeight,
			const Real                                  penaltySigma,
			const Real                                  penaltySlope,
			const EndCriteria&                          endCriteria ) {
            boost::shared_ptr<QuasiGaussianModel> qgmodel = boost::dynamic_pointer_cast<QuasiGaussianModel>(model);
            QL_REQUIRE(qgmodel,"QuasiGaussianModel required");            
            std::vector< boost::shared_ptr<SwapIndex> > indices;
            for (Size k=0; k<swapIndices.size(); ++k)
                indices.push_back( boost::dynamic_pointer_cast<SwapIndex>(swapIndices[k]) );
            return new QGCalibratorPtr(   
                new QGCalibrator(qgmodel,volTS,indices,modelTimesStepSize,useExpectedXY,
                    sigmaMax,slopeMax,etaMax,sigmaWeight,slopeWeight,etaWeight,
                    penaltySigma,penaltySlope, boost::make_shared<EndCriteria>(endCriteria)) );
        }
		// a single optimisation run
		Integer calibrate( const std::vector< std::vector< Real > >&  isInput,
			               const std::vector< std::vector< Real > >&  isOutput,
		                   Real                                       epsfcn = 1.0e-4 ) { // delta for finite differences
            return boost::dynamic_pointer_cast<QGCalibrator>(*self)->calibrate(isInput,isOutput,epsfcn);
        }

		// inspectors
		const boost::shared_ptr<RealStochasticProcess> calibratedModel() {
            return boost::dynamic_pointer_cast<RealStochasticProcess>(
                boost::dynamic_pointer_cast<QGCalibrator>(*self)->calibratedModel());
        }
		const std::vector<std::string>& debugLog() {
            return boost::dynamic_pointer_cast<QGCalibrator>(*self)->debugLog();
        }
		void acceptCalibration() { 
            boost::dynamic_pointer_cast<QGCalibrator>(*self)->acceptCalibration();
        }
    }    
};



#endif
