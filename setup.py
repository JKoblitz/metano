from setuptools import setup
from setuptools.command import sdist

sdist.walk_revctrl = lambda **kwargs: []


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="metano",
    version="1.3.0",
    description="metano - A Toolkit for Metabolic Network Analysis and "
    "Optimization",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Julia Helmecke and Alexander Riemer",
    author_email="julia.helmecke@tu-bs.de",
    url='http://metano.tu-bs.de/',
    packages=["metano"],
    package_dir={"metano": "src"},
    package_data={
        "metano": ["example/*.txt", "example/*.sh", "example/*.ini"]
    },
    entry_points={
        "console_scripts": [
            'deadends.py = metano.deadends:main',
            'fba.py = metano.fba:main',
            'fva.py = metano.fva:main',
            'moma.py = metano.moma:main',
            'knockout.py = metano.knockout:main',
            'sbmlparser.py = metano.sbmlparser:main',
            'sbmlwriter.py = metano.sbmlwriter:main',
            'splitratios.py = metano.splitratios:main',
            'mfba.py = metano.mfba:main',
            'mfm.py = metano.mfm:main',
            'modelassert.py = metano.modelassert:main',
            'moment.py = metano.moment:main',
            'to_check_constraints.py = metano.to_check_constraints:main',
            'to_diff_solutions.py = metano.to_diff_solutions:main',
            'to_evaluate_objective.py = metano.to_evaluate_objective:main',
            'to_fix_names.py = metano.to_fix_names:main',
            'to_list_active_metab.py = metano.to_list_active_metab:main',
            'to_rea2m.py = metano.to_rea2m:main',
            'to_scatterplot.py = metano.to_scatterplot:main',
        ]
    },
    classifiers=[
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python :: 3",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="genome-scale metabolic models flux balance analysis fba",
    license="GPL",
    install_requires=[
        "numpy",
        "swiglpk",
        "pymprog>=1.0",
        "ecos",
        "scipy",
        "openopt",
        "cvxpy",
        "cvxopt",
    ]
)
