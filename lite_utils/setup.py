from setuptools import setup

setup(name='lite_utils',
      version='0.1',
      description='KAGE-lite utility methods',
      long_description_content_type="text/markdown",
      packages=["lite_utils"],
      zip_safe=False,
      install_requires=['numpy', 'tqdm', 'pyfaidx', 'pathos', 'cython', 'scipy', 'scikit-allel',
                        'obgraph>=0.0.32',
                        'graph_kmer_index>=0.0.22',
                        'kmer_mapper>=0.0.30',
                        'graph_read_simulator>=0.0.7',
                        'shared_memory_wrapper>=0.0.27',
                        'bionumpy>=0.2.15',
                        'npstructures>=0.2.9'
                        ],
      include_dirs=["."],
      classifiers=[
            'Programming Language :: Python :: 3'
      ],
      entry_points={
            'console_scripts': ['lite_utils=lite_utils.command_line_interface:main']
      }

)

""""
rm -rf dist
python3 setup.py sdist
twine upload --skip-existing dist/*

"""