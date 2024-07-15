from setuptools import find_packages, setup


install_requires = [
    'matplotlib<3.4.0,>=3.1.0',
    'scikit-learn<1.0',
    'numpy>=1.22',
    'scipy<1.6.0',
    'ase==3.22.0',
    'pyyaml',
    'cffi',
    'psutil',
    'tqdm',
    'braceexpand',
]

setup_requires = [
    'cffi',
]

setup(
    name='spinner',
    version='0.2.0',
    description='SPINNER',
    author='Sungwoo Kang, Seungwoo Hwang',
    python_requires='>=3.6',
    packages=find_packages(include=['spinner', 'spinner*']),
    include_package_data=True,
    package_data={
        '': ['configure_default.yaml', 'INCAR_premelt', 'KPOINTS', '*.cpp', '*.h', 'params_*'],
    },
    entry_points={
        'console_scripts':[
            'spinner_auto_md = spinner.auto_md.spinner_auto_md:main',
            'spinner_nnp_train = spinner.nnp_train.initial_NNP_training:main',
            'configure_csp = spinner.utils.configure_csp:main',
            'spinner_csp = spinner.csp.spinner_csp:main',
            #'spinner_dft_relax = spinner.spinner_dft_relax:main',
        ]
    },
    install_requires=install_requires,
    setup_requires=setup_requires,
    cffi_modules=[
        "spinner/simple_nn/features/symmetry_function/libsymf_builder.py:ffibuilder",
        "spinner/simple_nn/utils/libgdf_builder.py:ffibuilder",
    ],
)
