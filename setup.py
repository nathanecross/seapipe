from pathlib import Path
import re

from setuptools import setup, find_packages

here = Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

version_file = here / "seapipe" / "version.py"
version_text = version_file.read_text(encoding="utf-8")

match = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', version_text, re.M)
if not match:
    raise RuntimeError("Unable to find __version__ in seapipe/version.py")

version = match.group(1)

setup(
    name="seapipe",
    version=version,
    description="Sleep Events Analysis pipeline of EEG data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/nathanecross/seapipe",
    author="Nathan E. Cross",
    author_email="nathan.cross@sydney.edu.au",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Research",
        "Topic :: Scientific/Engineering",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3 :: Only",
    ],
    keywords=[
        "sleep",
        "electrophysiology",
        "detection",
        "signal processing",
        "neuroscience",
        "analysis",
    ],
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=[
        "fooof",
        "mne",
        "numpy<=1.26.4",
        "openpyxl",
        "pandas",
        "psutil",
        "pyedflib",
        "PyWavelets",
        "pinguoin",
        "safepickle",
        "sleepecg",
        "scipy<1.13.0",
        "scikit-learn",
        "tensorpac",
        "wonambi",
        "yasa",
    ],
    package_data={
        "seapipe": ["VERSION"],
    },
    project_urls={
        "Bug Reports": "https://github.com/nathanecross/seapipe/issues",
        "Documentation": "https://seapipe.readthedocs.io/",
        "Source": "https://github.com/nathanecross/seapipe/",
    },
)