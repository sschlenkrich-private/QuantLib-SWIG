
/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007 StatPro Italia srl
 Copyright (C) 2019, 2020 Matthias Lungwitz

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

#ifndef quantlib_vectors_i
#define quantlib_vectors_i

%include stl.i
%include common.i
%include date.i

#if defined(SWIGCSHARP)
SWIG_STD_VECTOR_ENHANCED( std::pair<Date,double> )
#endif

%{
template <class T, class U>
std::vector<T> to_vector(const std::vector<U>& v) {
    std::vector<T> out(v.size());
    std::copy(v.begin(), v.end(), out.begin());
    return out;
}
%}

namespace std {

    %template(IntVector) vector<int>;
    %template(UnsignedIntVector) vector<unsigned int>;
    %template(SizeVector) vector<size_t>;
    %template(DoubleVector) vector<double>;
    %template(StrVector) vector<std::string>;
    %template(BoolVector) vector<bool>;

    %template(DoubleVectorVector) vector<vector<double>>;

    %template(DoublePair) pair<double,double>;
    %template(DoublePairVector) vector<pair<double,double>>;
    %template(PairDoubleVector) pair<vector<double>,vector<double>>;
    %template(UnsignedIntPair) pair<unsigned int,unsigned int>;
    %template(UnsignedIntPairVector) vector<pair<unsigned int,unsigned int>>;

#if !defined(SWIGR)
    %template(NodePair) pair<Date,double>;
    %template(NodeVector) vector<pair<Date,double> >;
#endif

#if defined(SWIGR)
    swigr_list_converter(IntVector,
    _p_std__vectorTint_std__allocatorTint_t_t,
    integer)

    swigr_list_converter(DoubleVector,
    _p_std__vectorTdouble_std__allocatorTdouble_t_t,
    numeric)

    swigr_list_converter(StrVector,
    _p_std__vectorTstd__string_std__allocatorTstd__string_t_t,
    character)
#endif
}

#endif
