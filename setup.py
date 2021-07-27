# coding=utf-8

from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

setup(
    name = 'duet',
    version = '0.1.1',
    description = 'Fast and Accurate Novel SV Calling and Haplotyping',
    author = 'Yekai Zhou',
    author_email = 'yekai.zhou@outlook.com',
    url = 'https://github.com/yekaizhou/duet',
    license = 'BSD',
    packages = find_packages(where = 'src'),
    package_dir = {'': 'src'},
    data_files = [('', ['LICENSE'])],
    scripts = ['src/duet/duet'],
    long_description = readme,
    long_description_content_type = 'text/markdown',
    zip_safe = False,
)
