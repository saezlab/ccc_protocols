from setuptools import setup, find_packages

setup(name='CCC_Protocols',
      version='0.0.1',
      description='LIANA x Tensor Tutorials',
      author='Hratch Baghdassarian, Daniel Dimitrov, Erick Armingol',
      author_email='daniel.dimitrov@uni-heidelberg.de',
      url='https://github.com/hmbaghdassarian/ccc_protocols',
      packages=find_packages(),
      package_data={'liana': ['resource/omni_resource.csv']},
      long_description='LIANA x Tensor Tutorials',
      long_description_content_type="text/markdown",
      install_requires=["ipykernel",
                        ],
      python_requires=">=3.7",
      classifiers=[
          "Programming Language :: Python :: 3",
          "Operating System :: OS Independent"]
      )