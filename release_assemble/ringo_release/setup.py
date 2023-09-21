import platform, glob, shutil, os, json
from setuptools import setup

# OS specifics
CUR_OS = platform.system()
SHAREDOBJ_TEMPLATE = {
    'Windows': "ringo_base.cp{py_ver}-win_amd64.pyd",
    'Linux': "ringo_base.cpython-{py_ver}*-x86_64-linux-gnu.so",
}
SO_FINAL_NAMES = {
    'Linux': 'ringo_base.so',
    'Windows': 'ringo_base.pyd',
}
SO_FINAL_NAME = SO_FINAL_NAMES[CUR_OS]
assert CUR_OS in ['Linux', 'Windows'], "Only Linux and Windows platforms are supported"
if CUR_OS == 'Windows':
    DLLDEPS_JSON = 'win_dlldeps.json'
    DLL_STORAGE_DIR = 'win_dlls'

# Python version specifics
python_version_tuple = platform.python_version_tuple()
py_ver = int(f"{python_version_tuple[0]}{python_version_tuple[1]}")

RINGORELEASE_DIR = 'ringo'
INSTALL_ADDITIONAL_FILES = []

# ===============================
# Actual installation starts here
# ===============================
# Find the appripriate build of ringo shared library
ringo_so_list = glob.glob(SHAREDOBJ_TEMPLATE[CUR_OS].format(py_ver=py_ver))
assert len(ringo_so_list) < 2, "Several pre-built libraries were found for your python version. This shouldn't have happend"
assert len(ringo_so_list) == 1, "There is no pre-built library for your version of Python"\
                                f"(your={python_version_tuple[0]}.{python_version_tuple[1]}), "\
                                "supported={3.7, ..., 3.11} (Linux) / {3.8, ..., 3.11} (Windows)"
ringo_object_name = ringo_so_list[0]

# Remove the library copied earlier if exists
so_final_path = os.path.join(RINGORELEASE_DIR, ringo_object_name)
if os.path.exists(so_final_path):
    os.remove(so_final_path)
# Put the build of ringo shared library into the package
shutil.copy2(ringo_object_name, so_final_path)

# Copy all DLLs required for running in Windows
if CUR_OS == 'Windows':
    assert os.path.isfile(DLLDEPS_JSON), f'Required file "{DLLDEPS_JSON}" not found'
    with open(DLLDEPS_JSON, 'r') as f:
        dlldeps_data = json.load(f)
    assert ringo_object_name in dlldeps_data, f"'{ringo_object_name}' is not accounted in {DLLDEPS_JSON}"
    for file in dlldeps_data[ringo_object_name]:
        shutil.copy2(os.path.join(DLL_STORAGE_DIR, file), RINGORELEASE_DIR)
        INSTALL_ADDITIONAL_FILES.append(file)


setup(
    name='ringo',
    version='0.0.1',
    author='Nikolai Krivoshchapov',
    python_requires=f'=={python_version_tuple[0]}.{python_version_tuple[1]}.*',
    install_requires=[
        'numpy',
        'networkx',
    ],
    packages=['ringo'],
    package_data={'ringo': ['__init__.py', ringo_object_name, 'cpppart/*', 'pyutils/**/*', 'pyutils/*'] + INSTALL_ADDITIONAL_FILES}
)