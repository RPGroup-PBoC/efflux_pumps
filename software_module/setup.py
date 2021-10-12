import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="efflux",
    version="0.0.1",
    author="Maria Carilli, Tom Roeschinger",
    author_email="troeschi@caltech.edu",
    description="This repository contains all active research materials for the study of Efflux pumps in E. coli",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/RPGroup-PBoC/efflux_pumps",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)