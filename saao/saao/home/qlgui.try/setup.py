from subprocess import call
from setuptools import setup, find_packages
from setuptools.command.install import install as _install
from distutils.command.build import build as _build
import glob

def find_data_files(dirname, ext):
    files = []
    for fil in glob.glob("{}/*.{}".format(dirname, ext)):
        files.append(fil)
    return files

setup(
    name='qlgui',
    description='QuickLook Software for Cassegrain Spectrograph',
    version='1.0',
    author='David Gilbank',
    author_email='gilbank@saao.ac.za',
    # TODO: license='',
    packages=find_packages(),
    install_requires=['future',],
    scripts=[
        'qlgui/newMain.py',
    cmdclass={
    },
    data_files=[
        ('/usr/local/lib/qlgui/', find_data_files("libs", "lis")), 
        ('/usr/local/lib/qlgui/', find_data_files("libs", "txt")), 
        ('/usr/local/lib/qlgui/CuAr/', find_data_files("libs/CuAr", "txt")), 
        ('/usr/local/lib/qlgui/CuNe/', find_data_files("libs/CuNe", "txt"))
    ],
    entry_points={
        'gui_scripts': [
            'quicklook = newMain:main',
        ],
    }
    
)
