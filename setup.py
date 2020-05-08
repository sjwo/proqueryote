import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="proqueryote", # Replace with your own username
    version="0.0.1",
    author="Stephen Wissow",
    author_email="sjw1000@wildcats.unh.edu",
    description="Get prokaryote sequences from NCBI.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sjwo/proqueryote",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
    ],
    python_requires='>=3.7',
)