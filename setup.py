from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='fmristats',
    version='0.1.0',
    description='Modelling the data and not the images in FMRI',
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
        'Programming Language :: Python :: 3.7',
    ],
    keywords='fmri neuroimaging neuroscience statistics',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['numpy', 'scipy', 'scikit-image', 'pandas',
        'statsmodels', 'matplotlib', 'nibabel', 'nipype'],
    entry_points={
        'console_scripts': [
            'ants4pop      = fmristats.cli.ants4pop:cmd',
            'csv2dataframe = fmristats.cli.csv2dataframe:cmd',
            'fmriati       = fmristats.cli.fmriati:cmd',
            'fmriblock     = fmristats.cli.fmriblock:cmd',
            'fmrifit       = fmristats.cli.fmrifit:cmd',
            'fmriinfo      = fmristats.cli.fmriinfo:cmd',
            'fmrimap       = fmristats.cli.fmrimap:cmd',
            'fmripop       = fmristats.cli.fmripop:cmd',
            'fmriprune     = fmristats.cli.fmriprune:cmd',
            'fmririgids    = fmristats.cli.fmririgids:cmd',
            'fmrisample    = fmristats.cli.fmrisample:cmd',
            'fmristudy     = fmristats.cli.fmristudy:cmd',
            'fsl4pop       = fmristats.cli.fsl4pop:cmd',
            'fsl4prune     = fmristats.cli.fsl4prune:cmd',
            'mat2block     = fmristats.cli.mat2block:cmd',
            'nii2image     = fmristats.cli.nii2image:cmd',
            'nii2session   = fmristats.cli.nii2session:cmd',
            'ref2plot      = fmristats.cli.ref2plot:cmd',
        ],
    },
)
