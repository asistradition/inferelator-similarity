import os
from setuptools import setup, find_packages

install_requires = ["pandas"]
tests_require = ["coverage", "nose"]
version = "0.1.0"

# Description from README.md
base_dir = os.path.dirname(os.path.abspath(__file__))
long_description = "\n\n".join([open(os.path.join(base_dir, "README.md"), "r").read()])

setup(
    name="inferelator_similarity",
    version=version,
    description="Inferelator-Similarity Gene Grouping Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/flatironinstitute/inferelator-similarity",
    author="Chris Jackson",
    author_email="cj59@nyu.edu",
    maintainer="Chris Jackson",
    maintainer_email="cj59@nyu.edu",
    packages=find_packages(include=["inferelator_similarity", "inferelator_similarity.*"]),
    zip_safe=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta"
    ],
    install_requires=install_requires,
    tests_require=tests_require,
)
