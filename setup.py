import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="kcalign",
    version="0.5.8",
    author="Nick Keener",
    author_email="nickeener@gmail.com",
    description="Kalgin-based codon-aware aligner for multiple sequences",
    install_requires=['biopython'],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/davebx/kc-align",
    packages=['kcalign'],
    entry_points={'console_scripts': ['kc-align=kcalign.cli:main']},
    classifiers=[
        "Programming Language :: Python :: 3",
        'License :: OSI Approved :: Academic Free License (AFL)',
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

