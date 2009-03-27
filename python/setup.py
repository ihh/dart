#!/usr/bin/env python

"""setup of XGram, the python interface to 
Ian Holmes' xgram.


"""
import glob
from distutils.core import setup

setup(name='XGram',
      version='0.1',
      description='Python interface to xgram',
      author='Andreas Heger',
      author_email='andreas.heger@gmail.com',
      url='http://www.biowiki.org',
      packages=['XGram',
                'XGram.Generator',
                'XGram.Generator.Prebuilt', 
                'XGram.Model', 
                'XGram.Run', 
                'XGram.Utils',
                'XGram.test',],
      package_data = {'XGram': ['data/*.stk', 'data/*.eg' ]},
     )

