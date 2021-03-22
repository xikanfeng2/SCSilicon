import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SCSilicon",
    version="1.0.4",
    author="Xikang Feng",
    author_email="fxk@nwpu.edu.cn",
    description="SCSilicon: a Python package that simulate single-cell DNA sequencing data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xikanfeng2/SCSilicon",
    project_urls={
        "Bug Tracker": "https://github.com/xikanfeng2/SCSilicon/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
    install_requires=[
        'numpy>=1.16.1',
        'pandas>=0.23.4',
        'tasklogger>=0.4.0',
        'wget>=3.2',
        'seaborn>=0.11.1',
        'matplotlib>=3.0.2',
    ],
)
