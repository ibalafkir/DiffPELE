from setuptools import setup, find_packages

setup(
    name="diffpele",
    author="Ismael Balafkir, Anna M. Diaz-Rovira, Victor Montal Blancafort",
    author_email="ibalafkir@gmail.com, annadiarov@gmail.com, victor.montal@protonmail.com",
    description=("The objective of this pipeline is to optimize protein-protein or antibody-antigen interactions. This is done through interface backbone diffusion and Monte Carlo rotations/translations using RFdiffusion and PELE."),
    license="unlicense",
    keywords="protein-protein interactions, antibody-antigen interactions, interface diffusion, Monte Carlo, RFdiffusion, PELE",
    url="https://github.com/ibalafkir/diffpele",
    packages=find_packages(exclude=('tests', 'docs')),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    install_requires=['pytest']
)
