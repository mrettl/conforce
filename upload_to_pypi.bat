python -m build --sdist
python -m build --wheel
python -m twine check dist/*
python -m twine upload --config-file .pypirc dist/*
rm -r dist
rm -r build