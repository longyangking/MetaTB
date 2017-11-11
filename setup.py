from distutils.core import setup

setup(name='metatb',
      version='0.0.1',
      author='Yang Long',
      author_email='longyang_123@yeah.net',
      url='https://longyangking.github.io/',
      download_url='https://github.com/longyangking/MetaTB/blob/master/LICENSE',
      keywords='tight binding, materials science, super crystal, meta atoms',
      py_modules=['metatb'],
      license="GPL 3.0",
      description="Simple solver for tight binding models for materials science and super crystals.",
      long_description="The tight binding method is an powerful theoretical tool for solving for the wave functions for photonics/phononics systems assuming a basis of localized atomic-like orbitals like electron in solid.",
      platforms=['UNIX','MAC OS X','Windows'],
      install_requires=['numpy','matplotlib'],
      )

