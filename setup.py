from setuptools import setup, find_packages


name = 'omamo'
__version__ = None
with open('{:s}/__init__.py'.format(name), 'rt') as fp:
    for line in fp:
        if line.startswith('__version__'):
            exec(line.rstrip())

requirements = ['numpy', 'pandas', 'tables', 'pyoma', 'tqdm']

desc = 'OMAmo - orthology-based model organism selection'
with open("README.md", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name=name,
    version=__version__,
    author='Alina Nicheperovich, Sina Majidian and Adrian Altenhoff',
    email='contact@omabrowser.org',
    url='https://github.com/DessimozLab/omamo',
    description=desc,
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=requirements,
    python_requires=">=3.6",
    license='LGPLv3',
    scripts=['bin/omamo'],
)
