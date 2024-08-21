from setuptools import setup, find_packages
import pathlib

from os import path

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")
# Get the current version number from inside the module
with open(path.join('seapipe', 'version.py')) as version_file:
    exec(version_file.read())

setup(
    name = "seapipe",  
    version = __version__,  
    description = "Sleep Events Analysis pipeline of EEG data",  
    long_description=long_description,  
    long_description_content_type="text/markdown",  
    url="https://github.com/nathanecross/seapipe",  
    author="Nathan E. Cross",  
    author_email="nathan.cross@sydney.edu.au",
    classifiers=[  
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Research",
        "Topic :: Neuroscience :: Sleep",
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
    keywords=["sleep", "electrophysiology", "detection", 
              "signal processing", "neuroscience", "analysis"], 
    packages=find_packages(),  
    python_requires=">=3.7",
    install_requires=["fooof",
                      "mne",
                      "numpy<=1.26.4",
                      "openpyxl",
                      "pandas",
                      "pinguoin",
                      "safepickle",
                      "scipy<1.13.0",
                      "tensorpac",
                      "wonambi",
                      "yasa"],  
    package_data={ 
        "seapipe": ["VERSION"],
    },
  
    project_urls={  # Optional
        "Bug Reports": "https://github.com/nathanecross/seapipe/issues",
        "Documentation": "https://seapipe.readthedocs.io/",
        "Source": "https://github.com/nathanecross/seapipe/",
    },
)
