import subprocess, os, shutil, glob, platform
import sys, copy, re

NPROC = 6
SHAREDOBJ_MASKS = {
    'Windows': "ringo_base.cp{pyver}*.pyd",
    'Linux': "ringo_base.cpython-{pyver}*.so",
}
assert platform.system() in ['Windows', 'Linux'], f"Running on an unsupported platform ({platform.system()})"
SHAREDOBJ_MASK = SHAREDOBJ_MASKS[platform.system()]

def is_int(x):
    try:
        int(x)
        return True
    except:
        return False
    

class FlagParser:
    def __init__(self, argv_flags):
        self.parsed_flags = []
        self.unparsed_flags = [flag for flag in argv_flags if flag.startswith('-')]
        self.all_flags = copy.deepcopy(argv_flags)

        self.ints_assignment = {}
        for i, flag in enumerate(argv_flags):
            if not flag.startswith('-'):
                assert argv_flags[i-1].startswith('-') and is_int(flag), f"Problem with flag {flag}"
                self.ints_assignment[argv_flags[i - 1]] = int(flag)

    def __contains__(self, flag_name):
        assert flag_name not in self.parsed_flags
        res = flag_name in self.unparsed_flags
        if res:
            self.parsed_flags.append(flag_name)
            del self.unparsed_flags[self.unparsed_flags.index(flag_name)]
        return res

    def __getitem__(self, flag_name):
        return self.ints_assignment[flag_name]

    def parsed_all_flags(self):
        return len(self.unparsed_flags) == 0

def runcmd(cmd):
    # Use regular expression to split by spaces while preserving quotes
    parts = [x for x in re.findall(r'(?:[^\s"]*(?:"[^"]*")?)+', cmd) if len(x)> 0]
    subprocess.call(parts)


def build_ringo(argv_flags):
    flags = []
    flags.append(f"-DBUILDFLAGS=\"{' '.join(argv_flags)}\"")
    argv_flags = FlagParser(argv_flags)

    if '-l' in argv_flags:
        flags.append("-DUNITTEST=1")
        flags.append("-DDEBUG=1")

    if '-t' in argv_flags:
        flags.append("-DDEBUG=1")
        flags.append("-DTLC_UNITTEST=1")

    if '-e' in argv_flags:
        flags.append("-DENDTOEND=1")
    if '-eweak' in argv_flags:
        flags.append("-DENDTOENDWEAK=1")
    
    if '-overlap_mult' in argv_flags:
        flags.append(f"-DOVERLAP_MULTIPLIER={argv_flags['-overlap_mult']}")

    if '-rmsd' in argv_flags:
        flags.append("-DRMSD_FILTER=1")
    
    if '-disable-overlaps' in argv_flags:
        pass
    else:
        flags.append("-DOVERLAP_DETECTION=1")

    if '-overlaps-final' in argv_flags:
        pass
    else:
        flags.append("-DOVERLAP_DETECTION_FINAL=1")
    
    if platform.system() == 'Windows':
        flags.append("-DOPENBLAS_DYNAMIC=1")
    # elif platform.system() == 'Linux':
    #     flags.append("-DOPENBLAS_STATIC=1")
        
    assert argv_flags.parsed_all_flags()

    runcmd("rm -Rf build/")
    runcmd("mkdir build")

    command = "cmake -B build {flags}".format(flags=" ".join(flags))
    runcmd(command) # -DCMAKE_BUILD_TYPE=Debug

    mainwd = os.getcwd()
    os.chdir(os.path.join(mainwd, "build"))
    if platform.system() == 'Windows':
        runcmd(f"ninja -j{NPROC}")
    elif platform.system() == 'Linux':
        runcmd(f"make -j{NPROC}")
    os.chdir(mainwd)

    sofile = glob.glob(os.path.join('./build', SHAREDOBJ_MASK.format(pyver='')))
    assert len(sofile) == 1
    sofile = sofile[0]

    shutil.copy2(sofile, '.')

if __name__ == "__main__":
    build_ringo(sys.argv[1:])
