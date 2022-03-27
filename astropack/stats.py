""" Module 'astropack.stats'.
    Numerous statistics and regression classes and functions.
"""

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

# Python core packages:
from math import sqrt

# External packages:
import pandas as pd
import statsmodels.regression.mixed_linear_model as sm_mm  # NB: only statsmodels version >= 0.8
import statsmodels.api as sm_api  # sm version >= 0.8


__all__ = ['MixedModelFit', 'LinearFit', 'weighted_mean']


_____REGRESSION_____________________________________________ = 0


class MixedModelFit:
    """ Represents one mixed-model (statsmodel) fit. Generic; not tied to astronomy.

    Results of the fit are made available as attributes.
    Current implementation uses formula form, i.e.,
    statsmodel::sm_mm.MixedLM.from_formula()

    >>> fit = MixedModelFit(df, 'Y', ['X1', 'X2'], 'a_group_type']
    >>> fit = MixedModelFit(df, 'Y', 'X1', 'a_group_type'] (one indep var)

    Parameters
    ----------
    data : |DataFrame|
        Input data, one variable per column, one data point per row.

    dep_var : str
        Name of column in ``data`` that serves as dependent 'Y' variable.

    fixed_vars : list of str, or str (if only one fixed variable, rare)
        One or more names of columns in ``data`` that serve as
        independent 'X' variable(s).

    group_var : str
        Name of column in ``data`` that serves as random-effect (group) variable.

    Attributes
    ----------
    statsmodels_object : ~statsmodels.regression.mixed_linear_model.MixedLMResults
        The complex results object as returned from statsmodels' Mixed Model fit.

    converged : bool
        True if Mixed Model fit successfully converged, else False.

        For mixed-model fits, failure to converge is not necessarily fatal, especially
        for trial fits. When a mixed-model fit fails to converge, the results are
        usually still useful for model refinement and especially for removal of outlier
        data points.

        The most common causes of failure to converge are:
        (1) a significant 'wild' outlier data point, and
        (2) trying to fit using too many independent variables for the quantity
        and quality of data points.

    nobs : int
        Number of data points (observations) used in the fit.

    likelihood : float
        Likelihood of the fit, as returned by statsmodels.

    dep_var : str
        Column name in ``data`` of the dependent 'Y' variable used in the fit,
        from input parameter ``dep_var``.

    fixed_vars : list of str
        Column name(s) in ``data`` of the independent 'X' variable(s) used in the fit,
        from input parameter ``fixed_vars``.

    group_var : str
        Column name in ``data`` of the 'random' group variable used in the fit,
        from input parameter ``group_var``.

    sigma : float
        RMS residual for the fit, in dependent-variable units.

    df_fixed_effects : |DataFrame|
        Dataframe, one row per independent 'X' variable, of data related to those
        variables, including their best-fit values, standard deviation, t-value, and
        P-value. Dataframe row index is the variable's name.

    df_random_effects : |DataFrame|
        Dataframe, one row per 'random-effects' group, of data related to those groups,
        especially the groups' effect sizes.

    df_observations : |DataFrame|
        Dataframe, one per input data point, of data related to each data point,
        including their fitted values and fit residuals.
    """
    # TODO: Replace internal 'formula' API with more robust column-name API.
    def __init__(self, data, dep_var=None, fixed_vars=None, group_var=None):
        if not isinstance(data, pd.DataFrame):
            raise TypeError('Parameter \'data\' must be a pandas Dataframe'
                            ' of input data.')
        if dep_var is None or fixed_vars is None or group_var is None:
            raise ValueError('Must provide all parameters: dep_var, '
                             'fixed_vars, and group_var.')
        if not isinstance(dep_var, str) or not isinstance(group_var, str):
            raise TypeError('Parameters \'dep_var\' and \'group_var\' '
                            'must both be strings.')
        fixed_vars_valid = False  # default if not validated
        if isinstance(fixed_vars, str):
            fixed_vars = [fixed_vars]
            fixed_vars_valid = True
        if isinstance(fixed_vars, list):
            if len(fixed_vars) >= 1:
                if all([isinstance(var, str) for var in fixed_vars]):
                    fixed_vars_valid = True
        if not fixed_vars_valid:
            raise TypeError('Parameter \'fixed_vars\' must be '
                            'a string or a list of strings.')
        formula = dep_var + ' ~ ' + ' + '.join(fixed_vars)

        model = sm_mm.MixedLM.from_formula(formula, groups=data[group_var], data=data)
        fit = model.fit()

        self.statsmodels_object = fit  # instance of class MixedLMResults (py pkg statsmodels)

        # Scalar and naming attributes:
        self.converged = fit.converged  # bool
        self.nobs = fit.nobs  # number of observations used in fit
        self.likelihood = fit.llf
        self.dep_var = dep_var
        self.fixed_vars = fixed_vars
        self.group_var = group_var
        self.sigma = sqrt(sum(fit.resid**2)/(fit.nobs-len(fixed_vars)-2))

        # Fixed-effects dataframe (use joins, so that we don't count on
        # consistent input ordering):
        df = pd.DataFrame({'Value': fit.fe_params})
        # Join on index (enforce consistency):
        df = df.join(pd.DataFrame({'Stdev': fit.bse_fe}))
        df = df.join(pd.DataFrame({'Tvalue': fit.tvalues}))
        df = df.join(pd.DataFrame({'Pvalue': fit.pvalues}))
        df['Name'] = df.index
        self.df_fixed_effects = df.copy()

        # Random-effect dataframe, index=GroupName, cols=GroupName, GroupValue:
        df = pd.DataFrame(fit.random_effects).transpose()  # DataFrame, 1 row/group
        df = df.rename(columns={'groups': 'Group'})  # was 'GroupValue'
        df['GroupName'] = df.index
        self.df_random_effects = df.copy()

        # Observation dataframe (safe to count on consistent input
        # ordering, so fit construction is easier and safer:
        df = pd.DataFrame({'FittedValue': fit.fittedvalues})
        df['Residual'] = fit.resid
        self.df_observations = df.copy()

    def predict(self, df_predict_input, include_random_effect=True):
        """From new data points, renders predicted dependent-variable values.
        Optionally includes effect of groups (random effects), unlike statsmodels
        itself.

        Parameters
        ----------
        df_predict_input : |DataFrame|
            New input data from which to make predictions, one row per new data point
            from which to predict dependent variable.
            Columns must include all independent variables and the 'random' group
            variable used in the fit.

        include_random_effect : bool, optional
            True if the effect of the 'random' group variable is to be included in the
            predicted values. Default is True, which is almost always correct.

        Returns
        -------
        total_prediction : |Series|
            Predictions of dependent-variable values matching rows of new data.
        """
        # First, get predicted values on fixed effects only
        # (per statsmodels' somewhat weird def. of 'predicted'):
        # One column per fixed-effect (independent, 'X') variable.
        fixed_effect_inputs = df_predict_input[self.fixed_vars]
        predicted_on_fixed_only = \
            self.statsmodels_object.predict(exog=fixed_effect_inputs)

        # If requested, add random-effect (group variable) contibutions
        # (which were not included in statsmodels' MixedModels object 'fit'):
        if include_random_effect:
            df_random_effect_inputs = pd.DataFrame(df_predict_input[self.group_var])
            df_random_effect_values = self.df_random_effects[['Group']]
            predicted_on_random_only = pd.merge(df_random_effect_inputs,
                                                df_random_effect_values,
                                                left_on=self.group_var,
                                                right_index=True, how='left',
                                                sort=False)['Group']
            total_prediction = predicted_on_fixed_only + predicted_on_random_only
        else:
            total_prediction = predicted_on_fixed_only

        return total_prediction


class LinearFit:
    """Represents one multivariate linear (statsmodel) fit.
    Generic; not tied to astronomy.

    Internally uses column-name API to statsmodels OLS.

    >>> fit = LinearFit(df_input, 'Y', ['X1', 'X2']]
    >>> fit = LinearFit(df_input, 'Y', 'X1'] (one indep var)

    Parameters
    ----------

    data : |DataFrame|
        Input data, one variable per column, one data point per row.

    dep_var : str
        Name of column in ``data`` that serves as dependent 'Y' variable.

    indep_vars : list of str, or str (if only one independent variable--rare)
        One or more names of columns in ``data`` that serve as
        independent 'X' variable(s).

    Attributes
    ----------

    statsmodels_object : ~statsmodels.regression.linear_model.RegressionResults
        The complex results object as returned from statsmodels' linear (OLS) fit.

    indep_vars : list of str
        Column name(s) in ``data`` of the independent 'X' variable(s) used in the fit,
        from input parameter ``indep_vars``.

    dep_var :str
        Column name in `data` of the dependent 'Y' variable used in the fit,
        from input parameter ``dep_var``.

    nobs : int
        Number of data points (observations) used in the fit.

    sigma : float
        RMS residual for the fit, in dependent-variable units.

    r2 : float
        R^2 correlation coefficient for the fit.

    df_indep_vars : |DataFrame|
        Dataframe, one row per independent 'X' variable, of data related to those
        variables, including their best-fit values, standard deviation, t-value, and
        P-value. Dataframe row index is the variable's name.

    df_observations : |DataFrame|
        Dataframe, one per input data point, of data related to each data point,
        including their fitted values and fit residuals.
    """
    def __init__(self, data, dep_var=None, indep_vars=None):
        if not isinstance(data, pd.DataFrame):
            print('Parameter \'data\' must be a pandas Dataframe of input data.')
            return
        if dep_var is None or indep_vars is None:
            print('Provide parameters: dep_var and indep_vars.')
            return
        if not isinstance(dep_var, str):
            print('Parameter \'dep_var\' must be a string.')
            return
        indep_vars_valid = False  # default if not validated
        if isinstance(indep_vars, str):
            indep_vars = [indep_vars]
            indep_vars_valid = True
        if isinstance(indep_vars, list):
            if len(indep_vars) >= 1:
                if all([isinstance(var, str) for var in indep_vars]):
                    indep_vars_valid = True
        if not indep_vars_valid:
            print('Parameter \'indep_vars\' must be a string or a list of strings.')
            return

        # Build inputs to regression fn and run it:
        y = data[dep_var]
        x = data[indep_vars].copy()
        x = sm_api.add_constant(x)
        model = sm_api.OLS(endog=y, exog=x)
        fit = model.fit()
        self.statsmodels_object = fit

        # Scalar and naming attributes:
        self.indep_vars = indep_vars  # as passed in
        self.dep_var = dep_var        # as passed in
        self.nobs = fit.nobs
        self.sigma = fit.mse_resid
        self.r2 = fit.rsquared_adj

        # Make solution (indep vars) dataframe:
        df = pd.DataFrame({'Value': fit.params})
        df = df.join(pd.DataFrame({'Stdev': fit.bse}))       # use join to enforce consistency
        df = df.join(pd.DataFrame({'Tvalue': fit.tvalues}))  # "
        df = df.join(pd.DataFrame({'PValue': fit.pvalues}))  # "
        df.index = ['Intercept' if x.lower() == 'const' else x for x in df.index]
        df['Name'] = df.index
        self.df_indep_vars = df.copy()

        # Make observation dataframe (rows in same order as in input dataframe):
        df = pd.DataFrame({'FittedValue': fit.fittedvalues})
        df['Residual'] = fit.resid
        self.df_observations = df.copy()

    def predict(self, df_predict_input):
        """From new data points, renders predicted dependent-variable values.

        Parameters
        ----------
        df_predict_input : |DataFrame|
            New input data from which to make predictions, one row per new data point
            from which to predict dependent variable.
            Columns must include all independent variables used in the fit.

        Returns
        -------
        predicted_y_values : |Series|
            Predictions of dependent-variable values matching rows of new data.
        """
        indep_var_inputs = df_predict_input[self.indep_vars]
        predicted_y_values = self.statsmodels_object.predict(exog=indep_var_inputs)
        return predicted_y_values


_____STATISTICAL_FUNCTIONS__________________________________ = 0


def weighted_mean(values, weights):
    """Returns weighted mean, weighted standard deviation of values, and
    weighted standard deviation of the mean.

    Parameters
    ----------
    values : list of float, or |Series| of float
        Values to be averaged.

    weights : list of float, or |Series| of float
        Weights, corresponding to `values`, used in making the weighted average.

        If a |Series|, the weights correspond element-wise to ``values``,
        as though they were lists of values, that is,
        without reference to the Series *index* of either ``values`` or ``weights``.

    Returns
    -------
    w_mean, w_stdev_pop, w_stdev_w_mean : tuple of 3 float
        Weighted mean, weighted standard deviation of the population of values,
        weighted standard deviation of the mean of values.

    Raises
    ------
    VelueError
        Raised when values and weights are unequal in number, as required.

        Raised when sum of weights is not positive, as required.
    """
    if (len(values) != len(weights)) or (len(values) == 0) or (len(weights) == 0):
        raise ValueError('lengths of values & weights must be equal & non-zero.')
    if sum(weights) <= 0:
        raise ValueError('sum of weights must be positive.')
    # Convert to lists, because py list comprehension can mangle pandas Series indices:
    value_list = list(values)
    weight_list = list(weights)
    norm_weights = [wt/sum(weights) for wt in weight_list]
    w_mean = sum(nwt * val for (nwt, val) in zip(norm_weights, value_list))
    n_nonzero_weights = sum(w != 0 for w in weight_list)

    if n_nonzero_weights == 1:
        w_stdev_pop = 0
        w_stdev_w_mean = 0
    else:
        resid2 = [(val-w_mean)**2 for val in value_list]
        nwt2 = sum(nwt**2 for nwt in norm_weights)
        rel_factor = 1.0 / (1.0 - nwt2)  # reliability factor (better than N'/(N'-1))
        w_stdev_pop = sqrt(rel_factor * sum(nwt * r2 for (nwt, r2) in zip(norm_weights, resid2)))
        w_stdev_w_mean = sqrt(nwt2) * w_stdev_pop
    return w_mean, w_stdev_pop, w_stdev_w_mean
