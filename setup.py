from setuptools import setup, find_packages
from shutil import copytree
import os


# Install dependencies
setup(
    name='solidstatesynth',
    version='0.0.1',
    description='Package to parse and characterize text-mined solid-state synthesis reactions.',
    author='Jane F. Schlesinger',
    author_email='schle759@umn.edu',
    python_requires='>=3.6.0',
    url='https://github.com/bartel-group/SolidStateSynth',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
)
