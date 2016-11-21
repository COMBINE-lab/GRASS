from setuptools import setup

setup(name='grass',
      version='0.1.1',
      scripts=['bin/grass'],
      description='Annotation of de novo assemblies using semi-supervised learning on graphs',
      url='https://github.com/COMBINE-lab/GRASS',
      author='Laraib Malik, Shravya Thatipally, Nikhil Junneti, Rob Patro',
      author_email='rob.patro@cs.stonybrook.edu',
      license='BSD with attribution',
      packages=['grass'],
      install_requires=[
          'PyYAML',
          'coloredlogs',
          'click',
          'networkx',
          'numpy',
          'pandas',
          'tqdm',
          'rapclust',
      ],
      zip_safe=False)
