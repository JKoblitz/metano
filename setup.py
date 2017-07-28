from setuptools import setup
from setuptools.command import sdist

sdist.walk_revctrl = lambda **kwargs: []

setup(
    name = "metano",
    version = "1.2.0",
    description = "metano - A Toolkit for Metabolic Network Analysis and "
                  "Optimization",
    author = "Alexander Riemer and Julia Helmecke",
    author_email = "julia.helmecke@tu-bs.de",
    url='http://metano.tu-bs.de/',
    packages = ["metano", "cvxmod"],
    package_dir = {"metano" : "src", "cvxmod" : "cvxmod"},
    package_data={
        "metano" : ["example/*.txt", "example/*.sh"]
        },
    entry_points = {
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
            'to_fix_names.py = metano.to_fix_names:main',
            'to_list_active_metab.py = metano.to_list_active_metab:main',
            'to_rea2m.py = metano.to_rea2m:main',
  ]
    },
    classifiers = [
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Programming Language :: Python",
        "Development Status :: 5 - Production/Stable",
# Development Status :: 2 - Pre-Alpha
# Development Status :: 3 - Alpha
# Development Status :: 4 - Beta
# Development Status :: 5 - Production/Stable
# Development Status :: 6 - Mature
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords = "genome-scale metabolic models flux balance analysis",
    license = "GPL",
    install_requires = [
        "numpy",
        "swiglpk",
        "pymprog>=1.0",
        "ecos",
        "scipy",
        "openopt",
        "cvxpy",
        "cvxopt<=1.1.8",
]
    )