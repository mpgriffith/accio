from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="accio",
    version="0.1.0",
    author="Marissa Griffith",
    author_email="mpp19@pitt.edu",
    description="A comprehensive tool for plasmid analysis and assignment in bacterial genomes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mpgriffith/accio",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "accio=accio.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "accio": ["data/*.csv", "data/*.txt"],
    },
)