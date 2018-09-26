from setuptools import setup
from distutils.util import convert_path

versionholder = {}
ver_path = convert_path('emm_typing/__init__.py')
with open(ver_path) as init_file:
    exec(init_file.read(), versionholder)

setup(name='emm_typing',
      version=versionholder['__version__'],
      description='Typing of the emm gene of Streptococcus pyogenes (and other Streptococci)',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Medical Science Apps.'
      ],
      keywords=['GAS', 'emm', 'typing', 'streptococcus', 'pyogenes'],
      author='Ola Brynildsrud',
      author_email='ola.brynildsrud@fhi.no',
      license='MIT',
      packages=['emm_typing'],
      install_requires=['biopython'],
      entry_points={
          'console_scripts': ['emm_typing = emm_typing.emm_typing:main']
      },
      package_dir={'emm_typing': 'emm_typing'},
      package_data={'emm_typing': ['data/*']}
      )
