###################################################################
Mixed-Model Regression (`astropack.stats`)
###################################################################

Class |MixedModelFit| as a comprehensive but accessible wrapper for statsmodels'
mixed-model regression class, with its impenetrable documentation.
I've translated it from Statsmodelese to a language spoken on planet earth.

Class |LinearFit| offers multivariate linear regression.

Both classes depend on pandas |DataFrame|s for input and output tables.

******************************************
Getting Started with Class |MixedModelFit|
******************************************

>>> from astropack.stats import MixedModelFit
>>> import pandas as pd
>>> import numpy as np
>>> # First, set up regression data:
>>> #    A, B, and C are independent 'X' variables,
>>> #    Ran is a category/group 'random' variable, and
>>> #    Dep is the independent 'Y' variable, with a bit of random error/noise
>>> points = 80
>>> np.random.seed(1234)
>>> df = pd.DataFrame({'A': np.random.randn(points),
         'B': np.random.randn(points),
         'C': np.random.randn(points),
         'Ran': np.random.randint(0, 3, points),
         'Dep': None})
>>> df['Dep'] = 17 + 1*df.A + 2*df.B + 4*df.C + 5*(df.Ran - 1) + \
                1*np.random.randn(len(df))
>>> categories = ['X', 'Y', 'Z']
>>> df['Ran'] = [categories[r] for r in df['Ran']]
>>> # Run the mixed-model regression:
>>> fit = MixedModelFit(df, dep_var='Dep', fixed_vars=['A', 'B', 'C'], group_var='Ran')
>>> # Inspect results:
>>> fit.converged, fit.nobs, fit.sigma
(True, 80, 0.9875046330402963)
>>> fit.df_fixed_effects  # no-noise values would be [17, 1, 2, 4]
               Value     Stdev     Tvalue         Pvalue       Name
Intercept  16.671046  2.874801   5.799025   6.670147e-09  Intercept
A           0.857537  0.115097   7.450542   9.295788e-14          A
B           1.938975  0.111619  17.371308   1.360888e-67          B
C           4.041939  0.115315  35.051338  3.720338e-269          C
>>> fit.df_fixed_effects.loc['B', 'Value']
1.938974721449711
>>> fit.df_random_effects  # no-noise values would be 5 * [-1, 0, 1] = [-5, 0, +5]
      Group GroupName
X -5.237536         X
Y  0.575538         Y
Z  4.661998         Z

******************************************
Getting Started with Class |LinearFit|
******************************************

>>> from astropack.stats import LinearFit
>>> import pandas as pd
>>> import numpy as np
>>> # First, set up regression data:
>>> #    A, B, and C are independent 'X' variables, and
>>> #    Dep is the independent 'Y' variable, with a bit of random error/noise
>>> points = 80
>>> np.random.seed(1234)
>>> d = {'A': np.random.randn(points),
        'B': np.random.randn(points),
        'C': np.random.randn(points),
        'Dep': 0}
>>> df_model = pd.DataFrame(d)
>>> df_model['Dep'] = 13 + 1.5*df_model.A + 2*df_model.B + 4*df_model.C + \
                      1*np.random.randn(len(df_model))
>>> # Run regression, yielding the fit object:
>>> fit = LinearFit(df_model, dep_var='Dep', indep_vars=['A', 'B', 'C'])
>>> # Inspect results:
>>> fit.nobs, fit.sigma, fit.r2
(80.0, 0.922577047445627, 0.960727350867271)
>>> fit.df_indep_vars  # no-noise values would be [13, 1.5, 2, 4]
               Value     Stdev      Tvalue        PValue       Name
Intercept  13.059410  0.107755  121.195947  1.004265e-88  Intercept
A           1.375857  0.109955   12.512852  3.847126e-20          A
B           1.984497  0.106754   18.589413  5.346495e-30          B
C           4.003729  0.109965   36.409117  7.426987e-50          C

***************
Reference/API
***************

.. automodapi:: astropack.stats
   :no-inheritance-diagram:
