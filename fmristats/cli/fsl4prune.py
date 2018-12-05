# Copyright 2016-2018 Thomas W. D. Möbius
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

Prune statistics field from non-brain areas

"""

########################################################################
#
# Command line program
#
########################################################################

from ..epilog import epilog

from .fmriprune import add_arguments

import argparse

def define_parser():
    parser = argparse.ArgumentParser(
            description=__doc__,
            epilog=epilog)

    parser = add_arguments(parser)

    specific = parser.add_argument_group(
        """Prune the statistics field with BET""")

    specific.add_argument('--cmd-bet',
        default='fsl5.0-bet',
        help="""FSL bet command. Must be in your path.""")

    specific.add_argument('--variante',
        default='R',
        help="""FSL bet command. Must be in your path.""")

    specific.add_argument('--vb-file',
        default='pruning/{cohort}-{id:04d}-{paradigm}-{date}-intercept.nii.gz',
        help="""brain mask in image space.""")

    specific.add_argument('--vb-mask',
        default='pruning/{cohort}-{id:04d}-{paradigm}-{date}-mask.nii.gz',
        help="""brain mask in image space.""")

    return parser

from .fmristudy import add_study_arguments

def cmd():
    parser = define_parser()
    add_study_arguments(parser)
    args = parser.parse_args()
    call(args)

cmd.__doc__ = __doc__

########################################################################
#
# Load libraries
#
########################################################################

import os

from os.path import isfile, isdir, join

from multiprocessing.dummy import Pool as ThreadPool

import numpy as np

from .fmristudy import get_study

from ..lock import Lock

from ..study import Study

from ..fsl import bet

########################################################################

def call(args):

    study = get_study(args)

    if study is None:
        print('No study found. Nothing to do.')
        return

    ####################################################################
    # Options
    ####################################################################

    force         = args.force_mask_overwrite
    verbose       = args.verbose
    threshold_df  = args.threshold
    proportion_df = args.fraction

    cmd_bet = args.cmd_bet
    variante = args.variante

    ####################################################################
    # Create the iterator
    ####################################################################

    study.update_layout({
        'vb_file':args.vb_file,
        'vb_mask':args.vb_mask,
        })

    study_iterator = study.iterate('result',
            new=['result', 'vb_file', 'vb_mask'],
            verbose=verbose)

    df = study_iterator.df.copy()

    ####################################################################
    # Wrapper
    ####################################################################

    def wm (result, vb_file, vb_mask, filename, name):

        if type(result) is Lock:
            print('{}: Result file is locked. Still fitting?'.format(
                name.name()))
            return

        if verbose > 2:
            print("""{}: Description of the fit:
                {}
                {}""".format(name.name(), result.describe(),
                    result.population_map.describe()))

        if hasattr(result.population_map, 'vb_mask') and not force:
            print("""{}:
            VB mask already present, use
            -f/--force-mask-overwrite
            to overwrite existing mask""".format(
                name.name()))
            return

        gf = result.get_field('degrees_of_freedom')

        if proportion_df:
            threshold_df = int(proportion_df * np.nanmax(gf.data))

        if verbose:
            print('{}: Lower df threshold: {:d}'.format(name.name(), threshold_df))

        inside = (gf.data >= threshold_df)

        intercept = result.get_field('intercept', 'point')
        intercept.mask(inside)
        intercept = intercept.round()

        if verbose:
            print('{}: Fit brain mask'.format(name.name()))

        template = bet(
                intercept = intercept,
                intercept_file = vb_file,
                mask_file = vb_mask,
                cmd = cmd_bet,
                variante = variante,
                verbose = verbose)

        if template is None:
            print('{}: Unable to bet'.format(name.name()))
            return

        result.population_map.set_vb_mask(gf.data >= threshold_df)

        result.population_map.vb_mask.name = 'FSL_BETTED_INTERCEPT'

        try:
            if verbose:
                print('{}: Save: {}'.format(name.name(), filename))
            result.save(filename)
        except Exception as e:
            print('{}: Unable to write: {}, {}'.format(name.name(),
                filename, e))

    ###################################################################

    if args.cores == 0:
        args.cores = None

    if len(df) > 1 and ((args.cores is None) or (args.cores > 1)):
        try:
            pool = ThreadPool(args.cores)
            for name, files, instances in study_iterator:
                result = instances['result']
                if result is not None:
                    vb_file  = files['vb_file']
                    vb_mask  = files['vb_mask']
                    filename = files['result']
                    pool.apply_async(wm, args=(result, vb_file, vb_mask, filename, name))

            pool.close()
            pool.join()
        except Exception as e:
            pool.close()
            pool.terminate()
            print('Pool execution has been terminated')
            print(e)
        finally:
            pass
    else:
        try:
            print('Process protocol entries sequentially')
            for name, files, instances in study_iterator:
                result = instances['result']
                if result is not None:
                    vb_file  = files['vb_file']
                    vb_mask  = files['vb_mask']
                    filename = files['result']
                    wm(result, vb_file, vb_mask, filename, name)
        finally:
            pass

    ####################################################################
    # Write study to disk
    ####################################################################

    if args.out is not None:
        if args.verbose:
            print('Save: {}'.format(args.out))

        dfile = os.path.dirname(args.out)
        if dfile and not isdir(dfile):
           os.makedirs(dfile)

        study.save(args.out)
