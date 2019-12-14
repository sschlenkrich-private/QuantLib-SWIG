/*
 Copyright (C) 2019 Sebastian Schlenkrich

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

#ifndef quantlib_hybridmodels_i
#define quantlib_hybridmodels_i

%include common.i
%include date.i
%include daycounters.i
%include types.i
%include currencies.i
%include observer.i
%include marketelements.i
%include interpolation.i
%include indexes.i

%include multiassetmodels.i
%include quasigaussian.i


%{

#include <ql/experimental/templatemodels/hybrid/hybridmodels.hpp>

using QuantLib::AssetModel;
using QuantLib::QuasiGaussianModel;
using QuantLib::HybridModel;
using QuantLib::SpreadModel;
%}


%shared_ptr(AssetModel);
class AssetModel : public RealStochasticProcess {
public:
	AssetModel(const Real X0,
	           const Real sigma);

	virtual Size size();
	virtual Size factors();
	virtual std::vector<Real> initialValues();					  

};

namespace std {
    %template(AssetModelVector) vector<boost::shared_ptr<AssetModel> >;
}


%shared_ptr(HybridModel);
class HybridModel : public RealStochasticProcess {
public:
	HybridModel(
        const std::string                                                domAlias,
        const boost::shared_ptr<RealStochasticProcess>                   domRatesModel,
        const std::vector<std::string>&                                  forAliases,
        const std::vector< boost::shared_ptr< AssetModel > >&            forAssetModels,
        const std::vector< boost::shared_ptr< RealStochasticProcess > >& forRatesModels,
        const std::vector< std::vector<Real> >&                          correlations );	

	// inspectors
	const std::string domAlias();
	const boost::shared_ptr<QuasiGaussianModel>& domRatesModel();
	const std::vector<std::string>&  forAliases();
	const std::vector< boost::shared_ptr<AssetModel> >& forAssetModels();
	const std::vector< boost::shared_ptr<QuasiGaussianModel> >& forRatesModels();
	const std::vector< std::vector<Real> >& correlations();
	const std::vector< std::vector<Real> >& L();
	const std::vector<size_t>& modelsStartIdx();

	virtual Size size();
	virtual Size factors();
	virtual std::vector<Real> initialValues();					  

};

%shared_ptr(SpreadModel);
class SpreadModel : public RealStochasticProcess {
public:
	SpreadModel(
        const boost::shared_ptr<RealStochasticProcess>     baseModel,
        const boost::shared_ptr<RealStochasticProcess>     sprdModel,
        const std::vector< std::vector<Real> >&            correlations );	

	// inspectors
	const boost::shared_ptr<RealStochasticProcess>& baseModel();
	const boost::shared_ptr<RealStochasticProcess>& sprdModel();

	virtual Size size();
	virtual Size factors();
	virtual std::vector<Real> initialValues();					  

};


#endif
