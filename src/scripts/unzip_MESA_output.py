from zipfile import ZipFile
from showyourwork.paths import user as Paths

paths = Paths()

if __name__ == "__main__":
    # loading the temp.zip and creating a zip object
    with ZipFile(paths.data / "data.zip", 'r') as zObject:
        # Extracting all the members of the zip
        # into a specific location.
        zObject.extractall(path=paths.src)
