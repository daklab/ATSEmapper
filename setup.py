from setuptools import setup, find_packages

setup(
    name="atsemapper",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "numpy",
        "networkx",
        "matplotlib",
        "seaborn",
        "pyfaidx",
        "gffutils",
        "tqdm",
    ],
    entry_points={
        'console_scripts': [
            'atsemapper=atsemapper.atsemapper.main:main',
        ],
    },
    author="Karin Isaev",
    author_email="karin.isaev@gmail.com",
    description="Utility modules for mapping and visualizing alternative splicing events in single cells or bulk RNA-seq data.",
    keywords="bioinformatics, single-cell, RNA-seq, alternative splicing",
    url="https://github.com/yourusername/ATSEmapper",
)