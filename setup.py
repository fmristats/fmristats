from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='fmristats',
    version='0.0.2',
    description='Rigorous statistical estimation of FMRI models',
    long_description=long_description,
    url='https://fmristats.github.io/',
    author='Thomas W. D. Möbius',
    author_email='moebius@medinfo.uni-kiel.de',
    license='GPLv3+',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='fmri mri statistics meta-analysis meta-regression imaging neuroimaging',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['numpy', 'scipy', 'scikit-image', 'pandas',
        'statsmodels', 'matplotlib', 'seaborn', 'nibabel'],
    entry_points={
        'console_scripts': [
            # general api
            'fmriassessment  = fmristats.cmd.api.fmriassessment:cmd',
            'fmriati         = fmristats.cmd.api.fmriati:cmd',
            'fmrifit         = fmristats.cmd.api.fmrifit:cmd',
            'fmriirritation  = fmristats.cmd.api.fmriirritation:cmd',
            'fmririgid       = fmristats.cmd.api.fmririgid:cmd',
            'fmrisample      = fmristats.cmd.api.fmrisample:cmd',
            'fmripop         = fmristats.cmd.api.fmripop:cmd',
            'fmriprune       = fmristats.cmd.api.fmriprune:cmd',

            # wrapper
            'fsl4prune = fmristats.cmd.wrp.fsl4prune:cmd',
            'fsl4pop   = fmristats.cmd.wrp.fsl4pop:cmd',

            # import from third party formats
            'mat2irr = fmristats.cmd.imp.mat2irr:cmd',
            'nii2img = fmristats.cmd.imp.nii2img:cmd',
            'nii2ses = fmristats.cmd.imp.nii2ses:cmd',
            'csv2protocol = fmristats.cmd.imp.csv2protocol:cmd_pkl',
            'csv2covariate = fmristats.cmd.imp.csv2protocol:cmd_cov',

            # plot
            'fit2plot = fmristats.cmd.plt.fit2plot:cmd',
            'ref2plot = fmristats.cmd.plt.ref2plot:cmd',
            #'ses2plot = fmristats.cmd.plt.ses2plot:cmd',

            # query
            'fmrimap         = fmristats.cmd.qry.fmrimap:cmd',
            'fmriinfo        = fmristats.cmd.qry.fmriinfo:cmd',
            'fmriprotocol    = fmristats.cmd.qry.fmriprotocol:cmd',
            'fmricovariate   = fmristats.cmd.qry.fmricovariate:cmd',

            # export to third party formats
            'fit2nii = fmristats.cmd.exp.fit2nii:cmd',

        ],
    },
)
