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

from .load import load, load_block_stimulus, load_session, \
        load_refmaps, load_population_map, load_result

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
        if verbose > 2:
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
            file_layout=None,
            root_dir=None,
            strftime=None,
            single_subject=False,
            scale_type=None,
            name=None,
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
        single_subject : None or False or True or str
            The default is False. If False, a default file layout for
            multiple subject is used. If True, a default file layout for
            a single subject analysis is used. If string, then the
            string is used.
        scale_type : None or str or float
            Scale type to use.
        name : None or str
            Name of this study

        Notes
        -----
        The covariates data frame can also be empty (None). If not None,
        though, it will never be allowed to be empty again.
        """
        # TODO: assert if protocol has an epi_code column
        # TODO: assert if protocol has valid index
        # TODO: assert if covariates have valid index
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
                'stimulus'       : '{cohort}-{id:04d}-{paradigm}-{date}.stimulus',
                'session'        : '{cohort}-{id:04d}-{paradigm}-{date}.session',
                'reference_maps' : '{cohort}-{id:04d}-{paradigm}-{date}-{rigids}.rigids',
                'population_map' : '{cohort}-{id:04d}-{paradigm}-{date}-{space}-{psi}.popmap',
                'result'         : '{cohort}-{id:04d}-{paradigm}-{date}-{space}-{psi}-{rigids}.fit'}
        elif type(single_subject) is str:
            self.file_layout = {
                'stimulus'       : single_subject + '.stimulus',
                'session'        : single_subject + '.session',
                'reference_maps' : single_subject + '.rigids',
                'population_map' : single_subject + '.popmap',
                'result'         : single_subject + '.fit'}
        else:
            self.file_layout = {
                'stimulus'       : '{study}/{paradigm}/{cohort}/{id:04d}/{cohort}-{id:04d}-{paradigm}-{date}.stimulus',
                'session'        : '{study}/{paradigm}/{cohort}/{id:04d}/{cohort}-{id:04d}-{paradigm}-{date}.session',
                'reference_maps' : '{study}/{paradigm}/{cohort}/{id:04d}/{cohort}-{id:04d}-{paradigm}-{date}-{rigids}.rigids',
                'population_map' : '{study}/{paradigm}/{cohort}/{id:04d}/{cohort}-{id:04d}-{paradigm}-{date}-{space}-{psi}.popmap',
                'result'         : '{study}/{paradigm}/{cohort}/{id:04d}/{cohort}-{id:04d}-{paradigm}-{date}-{space}-{psi}-{rigids}.fit'}

        if file_layout is not None:
            self.update_layout(file_layout)

        self.set_strftime(strftime)
        self.set_scale_type(scale_type)
        self.set_name(name)

    def set_name(self, name):
        if name is None:
            self.name = 'results'
        else:
            self.name = name

    def set_strftime(self, strftime=None):
        if strftime is None:
            self.strftime='%Y-%m-%d-%H%M'
        elif strftime == 'short':
            self.strftime = '%Y-%m-%d'
        else:
            self.strftime = strftime

    def set_scale_type(self, scale_type=None):
        if scale_type is None:
            self.scale_type = 'max'
        else:
            self.scale_type = scale_type


    def set_rigids(self, rigids_name=None):
        """
        Set the rigids field in protocol

        Parameters
        ----------
        rigids : str
            Name of the method that has been used to fit the head
            movements of the subject.

        Note
        ----
        This will overwrite all entries in the protocol! Use with care.
        """
        if rigids_name is None:
            if 'rigids' not in self.protocol.columns:
                self.protocol['rigids'] = 'undefined'
        else:
            self.protocol['rigids'] = rigids_name

    def set_diffeomorphism(self, diffeomorphism_name=None):
        """
        Set the diffeomorphism field in protocol

        Parameters
        ----------
        diffeomorphism : str
            Name of the method that has been used to fit the
            diffeomorphism from standard space to subject reference
            space.

        Note
        ----
        This will overwrite all entries in the protocol! Use with care.
        """
        if diffeomorphism_name is None:
            if 'diffeomorphism' not in self.protocol.columns:
                self.protocol['diffeomorphism'] = 'undefined'
        else:
            self.protocol['diffeomorphism'] = diffeomorphism_name

    def set_standard_space(self, vb_name=None):
        """
        Set the standard space field in protocol

        Parameters
        ----------
        space : str
            Name of the standard space.

        Note
        ----
        This will overwrite all entries in the protocol! Use with care.
        """
        if vb_name is None:
            if 'standard_space' not in self.protocol.columns:
                self.protocol['standard_space'] = 'undefined'
        else:
            self.protocol['standard_space'] = vb_name

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

    def iterate(self, *keys, new=None, lookup=None, vb_name=None,
            diffeomorphism_name=None, rigids_name=None,
            integer_index=False, verbose=True):
        """
        If covariates in not None, then only subjects in the protocol are
        going to be processed which are also marked as valid in the
        covariates file.

        Returns
        -------
        StudyIterator
        """

        ###############################################################
        # Set up data frames
        ###############################################################

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
            df = self.protocol.copy()

        df = df[df.valid == True].copy()
        df.reset_index(inplace=True)

        ###############################################################
        # Create look up entries
        ###############################################################

        if lookup is not None:
            for key in set(lookup):
                df [key] = Series(
                    data = [join(self.root_dir, self.file_layout[key].format(
                        study = self.name,
                        cohort = r.cohort,
                        id = r.id,
                        paradigm = r.paradigm,
                        date = r.date.strftime(self.strftime),
                        space = r.standard_space,
                        psi = r.diffeomorphism,
                        rigids = r.rigids,
                        )) for r in df.itertuples()],
                    index = df.index)

        ###############################################################
        # Define names
        ###############################################################

        if rigids_name is None:
            if ('rigids' not in df.columns):
                df['rigids'] = 'undefined'
        else:
            df['rigids'] = rigids_name

        if vb_name is None:
            if 'standard_space' not in df.columns:
                df['standard_space'] = 'undefined'
        else:
            df['standard_space'] = vb_name

        if diffeomorphism_name is None:
            if 'diffeomorphism' not in df.columns:
                df['diffeomorphism'] = 'undefined'
        else:
            df['diffeomorphism'] = diffeomorphism_name

        ###############################################################
        # Create entries
        ###############################################################

        if lookup is None:
            keys0 = keys
        else:
            keys0 = set(keys) - set(lookup)

        for key in keys0:
            df [key] = Series(
                data = [join(self.root_dir, self.file_layout[key].format(
                    study = self.name,
                    cohort = r.cohort,
                    id = r.id,
                    paradigm = r.paradigm,
                    date = r.date.strftime(self.strftime),
                    space = r.standard_space,
                    psi = r.diffeomorphism,
                    rigids = r.rigids,
                    )) for r in df.itertuples()],
                index = df.index)

        if new is not None:
            for n in new:
                if n not in keys:
                    df [n] = Series(
                        data = [join(self.root_dir, self.file_layout[n].format(
                            study = self.name,
                            cohort = r.cohort,
                            id = r.id,
                            paradigm = r.paradigm,
                            date = r.date.strftime(self.strftime),
                            space = r.standard_space,
                            psi = r.diffeomorphism,
                            rigids = r.rigids,
                            )) for r in df.itertuples()],
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

        if (self.protocol is None) or (len(self.protocol) < 1):
            print('No entries in the protocol')
            return

        protocol = self.protocol.sort_index()

        try:
            protocol = protocol.loc(axis=0)[(cohort, j, paradigm)]
        except KeyError as e:
            print("""
            Unable to process the query combination:
                cohort   = {},
                id       = {},
                paradigm = {}.
            Failed due to {}""".format( cohort, j, paradigm, e))
            raise

        if len(protocol) < 1:
            print('No entries left in the protocol')
            protocol = None

        if self.covariates is not None:
            if len(self.covariates) < 1:
                print('No entries in the covariates set')
                return

            covariates = self.covariates.sort_index()

            try:
                covariates = covariates.loc(axis=0)[(cohort, j)]
            except KeyError as e:
                print("""
                Unable to process the query combination:
                    cohort   = {},
                    id       = {},
                Failed due to {}""".format( cohort, j, e))
                raise

            covariates = covariates.loc(axis=0)[(cohort, j)]

            if len(covariates) < 1:
                print('No entries left in the covariates set')
                covariates = None
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
