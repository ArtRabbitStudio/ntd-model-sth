import setuptools

setuptools.setup(
    name='sth_simulation',
    version='0.0.3',
    url='https://www.artrabbit.com/',
    maintainer='ArtRabbit',
    maintainer_email='support@artrabbit.com',
    description='STH simulation model',
    long_description='Individual-based model in Medley 1989 thesis and Anderson & Medley 1985.',
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    install_requires=['numpy', 'pandas', 'joblib', 'flask', 'flask_cors', 'google-cloud-storage', 'gcsfs', 'fsspec'],
    include_package_data=True
)
