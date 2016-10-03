def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.misc_util import get_info


    config = Configuration( parent_package,
                            top_path)
    config.add_extension('fermidirac',['fermidirac.cpp'],
                         extra_compile_args=['-O3 -std=c++11 -fopenmp -Wall'],
                         extra_link_args=['-lgomp'])

    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration)
