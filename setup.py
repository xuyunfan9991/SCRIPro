import setuptools
import glob
import os

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

setuptools.setup(
    name="scripro",
    version="0.0.5",
    python_requires=">=3.8",
    keywords="Single-Cell Gene Regulatory Network Inference using ChIP-seq for Multi-omics",
    url="https://github.com/xuyunfan9991/SCRIPro",
    license="GPL-3.0+",
    packages=['scripro'],
    install_requires=requirements,
    entry_points={
            'console_scripts': [
            'SCIPRO=scripro.start:main'
            ]
        },
    include_package_data=True,
    data_files=[("", ["requirements.txt", "hg38.genome", "mm10.genome"])],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ]
    
)