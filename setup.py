import setuptools
import glob
import os


with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()


setuptools.setup(
    name="scripro",  
    version="1.0.13",  
    author="Xu Yunfan",  
    author_email="xuyunfan9991@gmail.com",
    description="Single-cell gene regulation network inference by large-scale data integration Pro",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown", 
    python_requires=">=3.8", 
    keywords="Single-Cell Gene Regulatory Network Inference using ChIP-seq for Multi-omics", 
    url="https://github.com/xuyunfan9991/SCRIPro", 
    license="GPL-3.0+",
    packages=setuptools.find_packages(),
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'scripro=scripro.cli:main' 
        ]
    },
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
 
    package_data={
        'scripro': ['data/*'],
    },
)


# setup.py
