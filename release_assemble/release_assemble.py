import os, sys, shutil, glob, platform, subprocess, ntpath, json

BUILD_FLAGS = "-rmsd -eweak"
RELEASE_DIR = "ringo_release"
RINGORELEASE_DIR = "ringo"
SHAREDOBJ_MASKS = {
    'Windows': "ringo_base.cp{pyver}*.pyd",
    'Linux': "ringo_base.cpython-{pyver}*.so",
}
SHAREDOBJ_MASK = SHAREDOBJ_MASKS[platform.system()]

if platform.system() == "Windows":
    DLL_STORAGE = os.path.join(RELEASE_DIR, 'win_dlls')
    DLLDEPS_JSON = os.path.join(RELEASE_DIR, 'win_dlldeps.json')

    # Have to do TWO build iterations (1) BUILD_DLL_STORAGE = True, (2) BUILD_DLL_STORAGE = False
    # (1) Compiles and stores required DLLs
    # (2) Compiles, installs with DLLs stored at (1) and runs tests
    BUILD_DLL_STORAGE = '-dll-copy' in sys.argv

def modes_from_flags(current_dir, all_envs, flags=""):
    BUILD_COMMAND = f'python build.py {flags}'
    TESTING_TASKS = ['pip install .','cd ../example_scripts', 'python run_ringo.py']

    modes = {}
    for py_version in all_envs:
        # Quick build within current python environment
        if py_version == 'current':
            pyversion_short = ''
            envload = []
        else: # When build requires activation of conda env
            if platform.system() == "Windows":
                envload = ["source /c/Users/Nikolai/anaconda3/etc/profile.d/conda.sh", "conda activate"]
            else:
                envload = [f'source ~/.bashrc']
            pyversion_short = py_version.replace('py', '')
            envload.append(f'conda activate {py_version}')
        
        sharedobj_name = SHAREDOBJ_MASK.format(pyver=pyversion_short)
        modes[py_version] = {
            'envload': envload,
            'build': [
                'cd ..',
                BUILD_COMMAND,
                f'mv {sharedobj_name} {os.path.join(current_dir, RELEASE_DIR)}',
                f'cd {os.path.join(current_dir, RELEASE_DIR)}',
            ]
        }

        # For Windows, do testing only on the 2nd build iteration
        if platform.system() == "Linux" or \
                platform.system() == "Windows" and not BUILD_DLL_STORAGE:
            modes[py_version]['build'] += TESTING_TASKS
        else:
            modes[py_version]['build'].append('cd ..')


    if platform.system() == "Windows":
        modes = {
            py_version: {
                ops_name: [item.replace('\\', '/') for item in ops_list]
                for ops_name, ops_list in mode_data.items()
            }
            for py_version, mode_data in modes.items()
        }
        if BUILD_DLL_STORAGE:
            for py_version, item in modes.items():
                if py_version == 'current':
                    pyversion_short = ''
                else:
                    pyversion_short = py_version.replace('py', '')
                sharedobj_name = SHAREDOBJ_MASK.format(pyver=pyversion_short)
                dll_store_call = function_call_line('dll_store', __file__, sharedobj_name)
                item['build'].append(dll_store_call)
    return modes

def dll_store(dll_name):
    dll_paths = glob.glob(os.path.join(RELEASE_DIR, dll_name))
    assert len(dll_paths) == 1
    dll_path = dll_paths[0]
    assert os.path.isfile(dll_path)
    
    temp_so_name = 'temp.so'
    shutil.copy2(dll_path, temp_so_name)
    sopath = os.path.join(os.getcwd(), temp_so_name)
    sseq = ['ldd', sopath]
    out = subprocess.run(sseq, stdout=subprocess.PIPE)
    outlines = str(out.stdout).split("\\n")
    deps = []
    for line in outlines:
        if len(line) > 10:
            dep_dll = line.split('=>')[1].split('(')[0].strip()
            # assert os.path.isfile(dep_dll), f"DLL dependency {dep_dll} not found"
            if dep_dll.startswith('/mingw64'):
                deps.append(dep_dll)
    os.remove(temp_so_name)

    dlldeps_data = {}
    if os.path.exists(DLLDEPS_JSON):
        with open(DLLDEPS_JSON, 'r') as f:
            dlldeps_data = json.load(f)
    
    library_basename = ntpath.basename(dll_path)
    assert library_basename not in dlldeps_data
    dlldeps_data[library_basename] = []

    for dep in deps:
        print(f"Checking out dep dll = {dep}")
        dll_stored_name = os.path.join(DLL_STORAGE, ntpath.basename(dep))
        if os.path.isfile(dll_stored_name):
            print(f"{ntpath.basename(dll_stored_name)} was already stored")
            # assert filecmp.cmp(dep, dll_stored_name), f"Files are different: {dep} vs. {dll_stored_name}]"
        else:
            if not os.path.isdir(DLL_STORAGE):
                sseq = ['mkdir', DLL_STORAGE]
                subprocess.run(sseq, stdout=subprocess.PIPE)
            sseq = ['cp', dep, dll_stored_name]
            subprocess.run(sseq, stdout=subprocess.PIPE)
        
        dlldeps_data[library_basename].append(ntpath.basename(dep))
    
    with open(DLLDEPS_JSON, 'w') as f:
        json.dump(dlldeps_data, f)


def function_call_line(function_name, script_path, *args):
    script_name = os.path.basename(script_path).replace('.py', '')
    return "python -c 'from {} import {}; {}(*{})'".format(script_name,
                                                           function_name,
                                                           function_name,
                                                           repr(args).replace("'", '\\"'))


def run_separately(script_parts, jobname='Untitled job'):
    # Assume that all elements of script_parts are lists of str
    script_parts = [element for sublist in script_parts for element in sublist]
    script = '&&'.join(script_parts)

    exec_parts = ['bash', '-c', f'"{script}"']
    print(' '.join(exec_parts))
    exit_code = os.system(' '.join(exec_parts))
    assert exit_code == 0, f'nonzero exit code at task "{jobname}"'


def build_ringo(env_name, all_envs, build_flags, script_path):
    # Primary checks
    assert isinstance(build_flags, str)
    script_dir = os.path.dirname(os.path.abspath(script_path))
    base_pyutils = os.path.join(script_dir, '..', 'pyutils')
    use_pyutils = os.path.join(script_dir, RELEASE_DIR, RINGORELEASE_DIR, 'pyutils')
    assert os.path.exists(base_pyutils), f'{base_pyutils} not found'
    # Update (copy a new) pyutils module for ringo
    if os.path.exists(use_pyutils):
        shutil.rmtree(use_pyutils)
    shutil.copytree(base_pyutils, use_pyutils)

    # Build ringo
    execution_mode = modes_from_flags(script_dir, all_envs, build_flags)[env_name]
    script_parts = [
        execution_mode['envload'],
        execution_mode['build']
    ]
    run_separately(script_parts, jobname=f'Ringo build ({env_name})')


if __name__ == "__main__":
    # Prepare
    assert os.path.isdir(f'./{RELEASE_DIR}')
    for sofile in glob.glob(os.path.join('.', RELEASE_DIR, SHAREDOBJ_MASK.format(pyver=''))):
        os.remove(sofile)
    if platform.system() == "Windows":
        if BUILD_DLL_STORAGE and os.path.isfile(DLLDEPS_JSON):
            os.remove(DLLDEPS_JSON)
    
    if '-full' in sys.argv:
        ENVS = ['py37', 'py38', 'py39', 'py310', 'py311']
        if platform.system() == 'Windows':
            del ENVS[0] # Python 3.7 is not supported on Windows
    else:
        ENVS = ['current']
    
    # Build for each linux env
    for env in ENVS:
        build_ringo(env, ENVS, BUILD_FLAGS, __file__)
