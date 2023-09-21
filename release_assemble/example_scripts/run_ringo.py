import math, random, sys
import ringo

# Outcomes of conformer generations are encoded as follows:
APS_RESULTS = {
    0: "success",
    1: "zero_solutions",
    2: "tlcfail",
    3: "overlap",
    4: "validationfail"
}

if __name__ == "__main__":
    # Initialize molecule
    # inpfile = sys.argv[1]
    inpfile = '../test_systems/pdb_7UPJ.sdf'

    print(f"Parsing {inpfile}")
    mol = ringo.Molecule(sdf=inpfile)

    # Obtain kinematically independend dihedrals and set them at random
    dofs_list, dofs_values = mol.get_ps()
    for i, cur_param in enumerate(dofs_list):
        newvalue = random.uniform(-math.pi, math.pi)
        print(f"Setting {repr(cur_param)} to {newvalue}")
        dofs_values[i] = newvalue

    # Obtain array of solution indices (discrete degrees of freedom)
    ddofs_list, ddofs_values = mol.get_discrete_ps()
    for i, cur_param in enumerate(ddofs_list):
        ddofs_values[i] = -1 # -1 requests for random solution of IK
    
    # Attempt to generate a conformer from the values of degrees of freedom
    result = mol.apply_ps()
    print(f"Conformer generation has terminated with flag '{APS_RESULTS[result]}'")

    # On success, print coordinates and save to file
    if result == 0:
        # Access coordinates and atom symbols
        coords = mol.get_xyz()
        symbols = mol.get_symbols()
        natoms = len(symbols)

        # Save XYZ file
        with open("result.xyz", "w") as xyzfile:
            # Write the number of atoms and empty comment line
            xyzfile.write(f"{natoms}\n\n")
            # Write the atomic symbols and coordinates
            for i in range(natoms):
                symbol = symbols[i]
                x, y, z = coords[i]
                xyzfile.write(f"{symbol} {x:.6f} {y:.6f} {z:.6f}\n")

    # Remove temporary file
    ringo.cleanup()
