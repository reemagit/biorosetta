#!/bin/bash
#conda activate biorosetta # change this to environment that has build and twine packages installed
python -m build --outdir=dist/to_upload
python -m twine upload --repository pypi dist/to_upload/*
mv dist/to_upload/* dist/