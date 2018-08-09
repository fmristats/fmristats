import pandas as pd

from pandas import Series, DataFrame

# TODO: give also the option to provide a list of ids!

def filter_protocol(df, cohort=None, j=None, k=None, paradigm=None):
    """
    Filter protocol data frame to specified cohort, within cohort id,
    and paradigm.

    Parameters
    ----------
    cohort : str
        Cohort name
    j : int
        Within cohort id of a subject
    paradigm : str
        Paradigm name of the stimulus

    Notes
    -----
    Return a copy of the specified data frame that only contains rows
    which fit to the given cohort name, within cohort id, and paradigm.
    Additionally only rows which marked as valid remain in the data
    frame.

    Returns
    -------
    The filtert data frame
    """
    if cohort:
        assert cohort in df.cohort.cat.categories, '{} is not a listed cohort'.format(cohort)
        df = df[(df.cohort == cohort)]

    if k:
        assert not (cohort is None), 'if k is specified, cohort must be specified'
        assert not (j is None), 'if k is specified, j must be specified'
        df = df[(j <= df.id) & (df.id <= k)]
    elif j:
        assert not (cohort is None), 'if j is specified, cohort must be specified'
        df = df[(df.id == j)]

    if paradigm:
        assert paradigm in df.paradigm.cat.categories, '{} is not a listed paradigm'.format(paradigm)
        df = df[(df.paradigm == paradigm)]

    return df[df.valid == True].copy()

def layout_dummy(df, key, template,
        strftime='%Y-%m-%d-%H%M'):
        df[key] = Series(
                data = [template.format(
                    r.cohort,
                    r.id,
                    r.paradigm,
                    r.date.strftime(strftime))
                    for r in df.itertuples()],
                index = df.index)

def layout_sdummy(df, key, template, urname, scale_type,
        strftime='%Y-%m-%d-%H%M'):
        df[key] = Series(
                data = [template.format(
                    r.cohort,
                    r.id,
                    r.paradigm,
                    r.date.strftime(strftime),
                    urname,
                    scale_type)
                    for r in df.itertuples()],
                index = df.index)

def layout_fdummy(df, key, template, vb, diffeo, scale_type,
        strftime='%Y-%m-%d-%H%M'):
        df[key] = Series(
                data = [template.format(
                    r.cohort,
                    r.id,
                    r.paradigm,
                    r.date.strftime(strftime),
                    vb,
                    diffeo,
                    scale_type)
                    for r in df.itertuples()],
                index = df.index)
