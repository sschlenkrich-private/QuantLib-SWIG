# Copyright (C) 2020 StatPro Italia srl
#
# This file is part of QuantLib, a free-software/open-source library
# for financial quantitative analysts and developers - http://quantlib.org/
#
# QuantLib is free software: you can redistribute it and/or modify it under the
# terms of the QuantLib license.  You should have received a copy of the
# license along with this program; if not, please email
# <quantlib-dev@lists.sf.net>. The license is also available online at
# <http://quantlib.org/license.shtml>.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the license for more details.

import QuantLib as ql

today = ql.Date(26, 9, 2019)
spot = ql.TARGET().advance(today, 2, ql.Days)

ql.Settings.instance().evaluationDate = today

bbgDates = [
    ql.Date(31, 3, 2020),
    ql.Date(30, 4, 2020),
    ql.Date(29, 5, 2020),
    ql.Date(30, 6, 2020),
    ql.Date(31, 7, 2020),
    ql.Date(31, 8, 2020),
    ql.Date(30, 9, 2020),
    ql.Date(30, 10, 2020),
    ql.Date(30, 11, 2020),
    ql.Date(31, 12, 2020),
    ql.Date(29, 1, 2021),
    ql.Date(26, 2, 2021),
    ql.Date(31, 3, 2021),
    ql.Date(30, 9, 2021),
    ql.Date(30, 9, 2022),
    ql.Date(29, 9, 2023),
    ql.Date(30, 9, 2024),
    ql.Date(30, 9, 2025),
    ql.Date(30, 9, 2026),
    ql.Date(30, 9, 2027),
    ql.Date(29, 9, 2028),
    ql.Date(28, 9, 2029),
    ql.Date(30, 9, 2030),
    ql.Date(30, 9, 2031),
    ql.Date(29, 9, 2034),
    ql.Date(30, 9, 2039),
    ql.Date(30, 9, 2044),
    ql.Date(30, 9, 2049),
    ql.Date(30, 9, 2054),
    ql.Date(30, 9, 2059),
    ql.Date(30, 9, 2064),
    ql.Date(30, 9, 2069),
]

bbgMktRates = [
    -0.373,
    -0.388,
    -0.402,
    -0.418,
    -0.431,
    -0.441,
    -0.45,
    -0.457,
    -0.463,
    -0.469,
    -0.461,
    -0.463,
    -0.479,
    -0.4511,
    -0.45418,
    -0.439,
    -0.4124,
    -0.37703,
    -0.3335,
    -0.28168,
    -0.22725,
    -0.1745,
    -0.12425,
    -0.07746,
    0.0385,
    0.1435,
    0.17525,
    0.17275,
    0.1515,
    0.1225,
    0.095,
    0.0644,
]

bbgZeroRates = [
    -0.373,
    -0.38058,
    -0.38718,
    -0.39353,
    -0.407,
    -0.41274,
    -0.41107,
    -0.41542,
    -0.41951,
    -0.42329,
    -0.42658,
    -0.42959,
    -0.43297,
    -0.45103,
    -0.4541,
    -0.43905,
    -0.41266,
    -0.3775,
    -0.33434,
    -0.2828,
    -0.2285,
    -0.17582,
    -0.1254,
    -0.0783,
    0.03913,
    0.14646,
    0.17874,
    0.17556,
    0.1531,
    0.12299,
    0.0948,
    0.06383,
]

index = ql.Euribor6M()

# market instruments

helpers = [
    ql.DepositRateHelper(
        bbgMktRates[0] / 100.0, ql.Period(6, ql.Months), 2, ql.TARGET(), ql.ModifiedFollowing, True, ql.Actual360()
    )
]

helpers += [ql.FraRateHelper(r / 100.0, i + 1, index) for i, r in enumerate(bbgMktRates[1:13])]

swapTenors = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20, 25, 30, 35, 40, 45, 50]
helpers += [
    ql.SwapRateHelper(
        r / 100.0, ql.Period(T, ql.Years), ql.TARGET(), ql.Annual, ql.ModifiedFollowing, ql.Thirty360(), index
    )
    for r, T in zip(bbgMktRates[13:32], swapTenors)
]

# synthetic helpers

additional_helpers = [ql.FraRateHelper(-0.004, 12 + i, index) for i in range(7)]
additional_dates = [ql.TARGET().advance(spot, 1 + i, ql.Months) for i in range(5)]


# global bootstrap over both

curve = ql.GlobalLinearSimpleZeroCurve(
    spot, helpers, ql.Actual365Fixed(), ql.GlobalBootstrap(additional_helpers, additional_dates, 1.0e-12)
)
curve.enableExtrapolation()


# report

print("%12s |%12s |%12s |%12s |%12s" % ("target date", "calc. pillar", "target rate", "calculated", "discrepancy"))

for i, (d, h, R) in enumerate(zip(bbgDates, helpers, bbgZeroRates)):
    pillar = h.pillarDate()

    if i < 13:
        day_counter = ql.Actual360()
        compounding = ql.Simple
    else:
        day_counter = ql.Thirty360()
        compounding = ql.SimpleThenCompounded

    r = curve.zeroRate(d, day_counter, compounding, ql.Annual).rate()

    print("%12s |%12s |%10.4f %% |%10.4f %% |%7.3f bps" % (d.ISO(), pillar.ISO(), R, r * 100, (r * 100 - R) * 100))
