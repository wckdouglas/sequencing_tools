from distutils.core import setup, Extension
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    raise ImportError("Requires cython to "
            "be installed before running setup.py (pip install cython)")
try:
    import numpy as np
except ImportError:
    raise ImportError("Requires numpy to "
            "be installed before running setup.py (pip install numpy)")
try:
    import pysam
except ImportError:
    raise ImportError("Requires pysam to "
            "be installed before running setup.py (pip install pysam)")

include_path = [np.get_include()]
include_path.extend(pysam.get_include())
ext_modules=cythonize([
        Extension('*', ['tgirt_seq_tools/*.pyx'],
            include_dirs = include_path)
])


setup(
    name='tgirt_seq_tools',
    version='0.1',
    description='Tools for TGIRT-seq pipeline',
    url='',
    author='Douglas Wu',
    author_email='wckdouglas@gmail.com',
    license='MIT',
    packages=['tgirt_seq_tools'],
    zip_safe=False,
    scripts = ['bin/reduce_multi_reads.py','bin/split_uniq_bam.py',
        'bin/filterSoftClip.py','bin/bam_to_bed.py',
        'bin/depth_to_bigwig.py'],
    ext_modules = ext_modules,
    install_requires=[
          'cython',
          'numpy',
          'pysam>0.11.0'
      ],
    cmdclass = {'build_ext': build_ext}
)
