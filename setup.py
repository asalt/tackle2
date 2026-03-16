from setuptools import setup, find_packages

setup(
    name="tackle2",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "python.report": ["prompt_specs/*.json"],
    },
    install_requires=[
        "Click",
        "jinja2",
        "ollama",
        "pandas",
        "SQLAlchemy",
    ],
    entry_points={
        "console_scripts": [
            "tackle2=python.cli:main",
        ],
    },
)
