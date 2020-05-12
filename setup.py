import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="trajectify",
    version="0.0.6",
    author="Adham Elarabawy",
    author_email="adhamelarabawy@gmail.com",
    description="A python package to create continuous quintic hermite splines out of user-selected points.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/adham-elarabawy/trajectify",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
