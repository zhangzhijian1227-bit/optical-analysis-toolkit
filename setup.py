from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="optical-analysis-toolkit",
    version="1.0.0",
    author="ZHANG ZHIJIAN",
    author_email="zhangzhijian1227@gmail.com",
    description="A toolkit for analyzing optical experiments including polarization, diffraction, and birefringence",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/zhangzhijian1227-bit/optical-analysis-toolkit",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Physics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
        "scipy>=1.5.0",
        "matplotlib>=3.3.0",
    ],
    extras_require={
        "dev": [
            "jupyter>=1.0.0",
            "pytest>=6.0",
        ],
    },
    keywords="optics physics polarization diffraction birefringence measurement analysis",
    project_urls={
        "Bug Reports": "https://github.com/zhangzhijian1227-bit/optical-analysis-toolkit/issues",
        "Source": "https://github.com/zhangzhijian1227-bit/optical-analysis-toolkit",
    },
)
