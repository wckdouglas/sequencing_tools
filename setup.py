from setuptools import find_packages, setup, Extension

# from distutils.core import setup, Extension
import glob

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    raise ImportError(
        "Requires cython to "
        "be installed before running setup.py (pip install cython)"
    )
try:
    import numpy as np
except ImportError:
    raise ImportError(
        "Requires numpy to " "be installed before running setup.py (pip install numpy)"
    )
try:
    import pysam
except ImportError:
    raise ImportError(
        "Requires pysam to " "be installed before running setup.py (pip install pysam)"
    )

include_path = [np.get_include()]
# include_path.append('/services_data/workspace/douglaswu/src/miniconda3/include')
include_path.append("/usr/include")
include_path.append("/usr/include/x86_64-linux-gnu")
include_path.extend(pysam.get_include())
# include_path.extend(os.environ['INCLUDE_PATH'].split(':'))
# include_path = filter(lambda x: x!= '', include_path)
ext_modules = cythonize(
    [
        Extension(
            "*", ["sequencing_tools/*tools/*.pyx"], include_dirs=list(include_path)
        )
    ],
    language_level="3",
)
for e in ext_modules:
    e.cython_directives = {"embedsignature": True}


setup(
    name="sequencing_tools",
    version="0.1",
    description="Tools for different NGS operations",
    url="",
    author="Douglas C. Wu",
    author_email="wckdouglas@gmail.com",
    license="MIT",
    packages=find_packages(),
    zip_safe=False,
    scripts=glob.glob("bin/seqtools"),
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    #    install_requires = [
    #          'cython',
    #          'numpy',
    #          'pysam>0.12.0',
    #          'matplotlib>=2.0.0',
    #          'seaborn>-0.7.1',
    #          'ujson',
    #          'scipy>=0.19.0',
    #          'networkx>=2.0',
    #          'pytest'
    #      ],
)
