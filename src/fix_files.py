import os

def modify_cpp_files(directory):
    for subdir, _, files in os.walk(directory):
        for file in files:
            if file.endswith(".cpp") or file.endswith(".h"):
                filepath = os.path.join(subdir, file)
                with open(filepath, 'r') as f:
                    lines = f.readlines()

                i = 0
                while i < len(lines):
                    line = lines[i].strip()
                    if line.startswith("#ifdef KDMOL_LOG"):
                        indent_size = lines[i].index('#')
                        routine_name = None
                        for j in range(i + 1, len(lines)):
                            inner_line = lines[j].strip()
                            if inner_line.startswith("log_routine"):
                                routine_name = inner_line.split("(")[1].split(")")[0].strip('"')
                            elif inner_line.startswith("#endif"):
                                if routine_name is not None:
                                    lines.insert(i, ''.join([' '] * indent_size) + f'checkpoint("{routine_name}");\n')
                                break
                        i = j
                    i += 1

                with open(filepath, 'w') as f:
                    f.writelines(lines)

# Example usage
modify_cpp_files('.')
