import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="biodescriptors",
    version="0.0.1",
    author="Igor Glukhov",
    author_email="gluhov2000@gmail.com",
    description="Biodescriptors package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/GlukhovIgor/Descriptors_package",
    project_urls={
        "Bug Tracker": "https://github.com/GlukhovIgor/Descriptors_package/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)