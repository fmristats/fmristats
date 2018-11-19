# Copyright 2016-2017 Thomas W. D. Möbius
#
# This file is part of fmristats.
#
# fmristats is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# fmristats is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# It is not allowed to remove this copy right statement.

"""

Study Layout

"""

from .name import Identifier

from .load import load, load_block_stimulus, load_session, load_refmaps, \
    load_population_map, load_result

import pandas as pd

from pandas import Series, DataFrame

import os

from os.path import isfile, isdir, join

import pickle

def load_verbose(f, verbose=False, name=None):
    try:
        instance = load(f)
        if verbose:
            print('{}: Read {}'.format(name.name(), f))
        return instance
    except Exception as e:
        if verbose > 1:
            print('{}: Unable to read {}, {}'.format(name.name(), f, e))
        return None

class StudyIterator:
    def __init__(self, df, keys, new=None, verbose=True,
            integer_index=False):
        assert type(df) is DataFrame, 'df must be DataFrame'

        self.df = df
        self.keys = keys
        self.new = new
        self.verbose = verbose
        self.integer_index = integer_index

    def __iter__(self):
        self.it = self.df.itertuples()
        return self

    def __next__(self):
        r = next(self.it)
        name = Identifier(cohort=r.cohort, j=r.id, datetime=r.date, paradigm=r.paradigm)
        if self.new is None:
            if self.integer_index is False:
                return name, \
                    {k : load_verbose(getattr(r, k), self.verbose, name) for k in self.keys}
            else:
                return r.Index, name, \
                    {k : load_verbose(getattr(r, k), self.verbose, name) for k in self.keys}
        else:
            if self.integer_index is False:
                return name, \
                    {k : getattr(r, k) for k in self.new}, \
                    {k : load_verbose(getattr(r, k), self.verbose, name) for k in self.keys}
            else:
                return r.Index, name, \
                    {k : getattr(r, k) for k in self.new}, \
                    {k : load_verbose(getattr(r, k), self.verbose, name) for k in self.keys}

class Study:
    """
    Study Layout
    """

    def __init__(self,
            protocol,
            covariates=None,
            vb=None,
            vb_background=None,
            vb_ati=None,
            diffeomorphism=None,
            scale_type=None,
            file_layout=None,
            root_dir=None,
            strftime=None,
            single_subject=False,
            ):
        """
        Parameters
        ----------
        protocol : DataFrame
            A data frame that contains information on the FMRI sessions
            that a subject has been participated.
        covariates : None or DataFrame
            A data frame that contains potential covariates of the
            subjects.
        vb : None or Image
            Template image in standard space.
        vb_background : None or Image
            Background template in standard space.
        vb_ati : None or Image
            ATI reference field.
        file_layout : None or dict
            A dictionary of the file layout.
        strftime : None or str
            Format of the date and time string.
        root_dir : None or str
            The root directory of the study. All paths will always be
            expanded with respect to this root.
        single_subject : False or True or str
            The default is False. If False, a default file layout for
            multiple subject is used. If True, a default file layout for
            a single subject analysis is used. If string, then the
            string is used.
        diffeomorphism : str or None
            Name or type of the family of diffeomorphisms that will map
            from standard space to the respective subjects.
        scale_type : str or float
            Scale type to use.

        Notes
        -----
        The covariates data frame can also be empty (None). If not None,
        though, it will never be allowed to be empty again.
        """
        # TODO: check if protocol has an epi_code column
        # TODO: check if protocol has valid index
        # TODO: check if covariates have valid index
        self.protocol        = protocol
        self.covariates      = covariates
        self.vb              = vb
        self.vb_background   = vb_background
        self.vb_ati          = vb_ati

        if root_dir is None:
            self.root_dir = ''
        else:
            self.root_dir = root_dir

        if single_subject is True:
            self.file_layout = {
                'stimulus'       : '{0}-{1:04d}-{2}-{3}.stm',
                'session'        : '{0}-{1:04d}-{2}-{3}.ses',
                'reference_maps' : '{0}-{1:04d}-{2}-{3}.ref',
                'population_map' : '{0}-{1:04d}-{2}-{3}-{4}.pop',
                'result'         : '{0}-{1:04d}-{2}-{3}-{4}.fit'}
        elif type(single_subject) is str:
            self.file_layout = {
                'stimulus'       : single_subject + '.stm',
                'session'        : single_subject + '.ses',
                'reference_maps' : single_subject + '.ref',
                'population_map' : single_subject + '.pop',
                'result'         : single_subject + '.fit'}
        else:
            self.file_layout = {
                'stimulus'       : 'data/stm/{2}/{0}-{1:04d}-{2}-{3}.stm',
                'session'        : 'data/ses/{2}/{0}-{1:04d}-{2}-{3}.ses',
                'reference_maps' : 'data/ref/{2}/{0}-{1:04d}-{2}-{3}.ref',
                'population_map' : 'data/pop/{2}/{4}/{5}/{0}-{1:04d}-{2}-{3}-{4}.pop',
                'result'         : 'data/fit/{2}/{4}/{5}/{6}/{0}-{1:04d}-{2}-{3}-{4}.fit'}

        if file_layout is not None:
            self.update_layout(file_layout)

        if strftime is None:
            self.strftime='%Y-%m-%d-%H%M'
        elif strftime == 'short':
            self.strftime = '%Y-%m-%d'
        else:
            self.strftime = strftime

        if scale_type is None:
            self.scale_type = 'max'
        else:
            self.scale_type = scale_type

        if diffeomorphism is None:
            self.diffeomorphism = 'identity'
        else:
            self.diffeomorphism = diffeomorphism

    def update_layout(self, file_layout):
        """
        Update the file layout

        Parameters
        ----------
        file_layout : dict
            A dictionary of the file layout.
        """
        self.file_layout.update( (k,v) for k,v in file_layout.items() if
                v is not None)

    def update_protocol(self, df, verbose=True):
        """
        Update fields in the protocol file

        Parameters
        ----------
        df : DataFrame
            Entries
        """
        old_index_names = self.protocol.index.names
        self.protocol.reset_index(inplace=True)
        self.protocol.set_index(df.index.names, inplace=True)
        self.protocol.update(df)
        self.protocol.reset_index(inplace=True)
        self.protocol.set_index(old_index_names, inplace=True)
        self.protocol.valid = self.protocol.valid.astype(bool)

    def update_covariates(self, df, verbose=True):
        """
        Update fields in the covariates file

        Parameters
        ----------
        df : DataFrame
            Entries
        """
        old_index_names = self.covariates.index.names
        self.covariates.reset_index(inplace=True)
        self.covariates.set_index(df.index.names, inplace=True)
        self.covariates.update(df)
        self.covariates.reset_index(inplace=True)
        self.covariates.set_index(old_index_names, inplace=True)
        self.covariates.valid = self.covariates.valid.astype(bool)

    def iterate(self, *keys, new=None,
            vb_name=None, diffeomorphism_name=None, scale_type=None,
            verbose=True, integer_index=False):
        """
        If covariates in not None, then only subjects in the protocol are
        going to be processed which are also marked as valid in the
        covariates file.

        Returns
        -------
        StudyIterator
        """
        if vb_name is None:
            if self.vb is not None:
                vb_name = self.vb.name

        if diffeomorphism_name is None:
            diffeomorphism_name = self.diffeomorphism

        if scale_type is None:
            if type(self.scale_type) is str:
                scale_type = self.scale_type
            elif type(scale_type) is float:
                scale_type = '{:.2f}'.format(
                        self.scale_type).replace('.', 'd')
            else:
                scale_type = 'undefined'

        df = self.protocol[self.protocol.valid == True].copy()

        if self.covariates is not None:
            covariates = self.covariates
            old_index_names = df.index.names
            df.reset_index(inplace=True)

            df = (df.join(covariates, on=['cohort', 'id'], lsuffix='_')
                    .assign(valid=lambda x: x.valid.fillna(False))
                    .assign(valid=lambda x: x.valid & x.valid_)
                    .drop('valid_', axis=1))

            df.set_index(old_index_names, inplace=True)
            df.valid = df.valid.astype(bool)

        else:
            df = protocol.copy()

        df = df[df.valid == True].copy()
        df.reset_index(inplace=True)

        for key in keys:
                df [key] = Series(
                        data = [join(self.root_dir, self.file_layout[key].format(
                            r.cohort,
                            r.id,
                            r.paradigm,
                            r.date.strftime(self.strftime),
                            vb_name,
                            diffeomorphism_name,
                            scale_type))
                            for r in df.itertuples()],
                        index = df.index)
        if new is not None:
            for n in new:
                if n not in keys:
                    df [n] = Series(
                            data = [join(self.root_dir, self.file_layout[n].format(
                                r.cohort,
                                r.id,
                                r.paradigm,
                                r.date.strftime(self.strftime),
                                vb_name,
                                diffeomorphism_name,
                                scale_type))
                                for r in df.itertuples()],
                            index = df.index)

        return StudyIterator(df, keys, new, verbose, integer_index)

    def filter(self, cohort=None, j=None, paradigm=None, inplace=False):
        """
        Filter protocol and covariates

        Will filter the study to only include protocol and the covariate
        entries which match the specified cohort, id, or paradigm.

        Parameters
        ----------
        cohort : str
            Only keep subject that belong to cohort.
        j : int or tuple or list
            Only keep the subject that has this id (if j is int), that
            lies between the tuple of ids (if j is tuple) or that is
            in the list (if j is list)
        paradigm : str
            Only keep protocol entries that belong to this paradigm.

        Returns
        -------
        study : Study or None
            The filtered study or None if inplace is True.
        """
        if cohort is None:
            cohort = slice(None)

        if paradigm is None:
            paradigm = slice(None)

        if j is None:
            j = slice(None)
        elif len(j) == 1:
            j = slice(j[0], j[0])
        elif len(j) == 2:
            j = slice(j[0], j[1])

        protocol = self.protocol.sort_index()

        if len(protocol) < 1:
            print('No valid entries in the protocol')
            return

        protocol = protocol.loc(axis=0)[(cohort, j, paradigm)]

        if len(protocol) < 1:
            print('No entries left in the protocol')
            return

        if self.covariates is not None:
            covariates = self.covariates.sort_index()

            if len(covariates) < 1:
                print('No valid entries in the covariates data set')
                return

            covariates = covariates.loc(axis=0)[(cohort, j)]

            if len(covariates) < 1:
                print('No entries left in the covariates data set')
                return
        else:
            covariates = None

        if inplace is True:
            self.protocol = protocol
            self.covariates = covariates
        else:
            return Study(protocol, covariates,
            vb            = study.vb,
            vb_background = study.vb_background,
            vb_ati        = study.vb_ati,
            file_layout   = study.file_layout,
            strftime      = study.strftime,
            )

    def save(self, file, **kwargs):
        """
        Save instance to disk

        This will save the current instance to disk for later use.

        Parameters
        ----------
        file : str
            File name.
        """
        with open(file, 'wb') as output:
            pickle.dump(self, output, **kwargs)
