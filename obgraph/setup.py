from setuptools import setup
from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy as np

setup(name='obgraph',
      version='0.0.34',
      description='obgraph',
      url='http://github.com/ivargr/obgraph',
      author='Ivar Grytten',
      author_email='',
      license='MIT',
      packages=["obgraph"],
      zip_safe=False,
      install_requires=['cython', 'numpy', 'tqdm', 'pathos', 'shared_memory_wrapper>=0.0.18', 'bionumpy'],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points={
            'console_scripts': ['obgraph=obgraph.command_line_interface:main']
      },
      cmdclass = {"build_ext": build_ext},
      ext_modules = cythonize(["obgraph/cython_traversing.pyx"]),
      include_dirs=np.get_include(),
)

"""

rm -rf dist
python3 setup.py sdist
twine upload --skip-existing dist/*


rm -rf dist
python3 setup.py sdist bdist_wheel
auditwheel repair --plat manylinux_2_17_x86_64 dist/obgraph-*-cp38-cp38-linux_x86_64.whl
rm dist/*.whl
mv wheelhouse/* dist
python3 -m twine upload --repository pypi dist/*

"""