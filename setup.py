from setuptools import setup

setup(
    name="pyrbnachem",
    version="1.0.0",
    author="Adam Faulconbridge",
    author_email="afaulconbridge@googlemail.com",
    packages=["pyrbnachem"],
    description="TODO.",
    long_description=open("README.md").read(),
    install_requires=["pyrbn", "pyachem"],
    extras_require={
        "dev": [
            "pytest-cov",
            "flake8",
            "pylint",
            "pip-tools",
            "pipdeptree",
            "pre-commit",
        ],
    },
)
