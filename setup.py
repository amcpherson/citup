from setuptools import setup, find_packages

setup(
    name='citup',
    version='0.1.0',
    packages=find_packages(),
    scripts=[
        'citup/run_citup_iter.py',
        'citup/run_citup_qip.py',
    ],
)
