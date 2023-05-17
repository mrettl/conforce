import os
import zipfile
import datetime

version = datetime.datetime.today().strftime('%Y-%m-%d')

with open("cf_shared/__version__.py", "w") as fh:
    fh.write(f"version = '{version}'\n")

os.makedirs("release", exist_ok=True)

with zipfile.ZipFile(f"release/cf_plugin_{version}.zip", "w") as fh:
    for file in ["cf_abq_plugin.py", "README.md", "LICENSE.txt"]:
        with fh.open(f"cf/{file}", "w") as writer, open(file, "rb") as reader:
            print(file)
            writer.write(reader.read())

    for packages in ["cf_abq", "cf_shared"]:
        for root, dirs, files in os.walk(packages):
            if "__pycache__" in dirs:
                dirs.remove("__pycache__")

            for file in files:
                if file.endswith(".pyc"):
                    continue

                path = os.path.relpath(f"{root}/{file}", ".")
                print(path)

                with fh.open("cf/" + path, "w") as writer, open(path, "rb") as reader:
                    writer.write(reader.read())

