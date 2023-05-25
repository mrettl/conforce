import os
import zipfile
import datetime


from conforce_shared.__version__ import version

with open("conforce_shared/__version__.py", "w") as fh:
    (version_id, *modifiers) = version.split("-")
    (major, minor, patch) = version_id.split(".")
    new_patch = int(patch) + 1

    new_version_id = f"{major}.{minor}.{new_patch}"
    version = "-".join([new_version_id] + modifiers)

    fh.write(f"version = '{version}'\n")
    fh.write(f"build_date = '{datetime.datetime.utcnow().strftime('%Y.%m.%d %H:%M:%S')}'\n")

os.makedirs("release", exist_ok=True)

with zipfile.ZipFile(f"release/conforce_plugin_{version}.zip", "w") as fh:
    for file in ["conforce_abq_plugin.py", "README.md", "LICENSE.txt"]:
        with fh.open(f"conforce/{file}", "w") as writer, open(file, "rb") as reader:
            print(file)
            writer.write(reader.read())

    for packages in ["conforce_abq", "conforce_shared"]:
        for root, dirs, files in os.walk(packages):
            if "__pycache__" in dirs:
                dirs.remove("__pycache__")

            for file in files:
                if file.endswith(".pyc"):
                    continue

                path = os.path.relpath(f"{root}/{file}", ".")
                print(path)

                with fh.open("conforce/" + path, "w") as writer, open(path, "rb") as reader:
                    writer.write(reader.read())

