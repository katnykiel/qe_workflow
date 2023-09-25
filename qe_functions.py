def read_key():
    """
    Read in new Materials Project API key
    """
    import os, stat
    from IPython.display import clear_output

    # Read in new Materials Project API key, if one exists
    try:
        with open(os.path.expanduser('~/.mpkey.txt'), 'r') as f:
            key = f.readlines()[0]
            return key
    except:
        key = ""

    # Check if API key already exists, skip try-except
    if not key:
        # Prompt user for API key
        try:
            user = str(input())
            clear_output()
            if not user.isalnum():
                raise TypeError('Wrong Key')
            if user == None:
                raise TypeError('Empty')
            with open(os.path.expanduser('~/.mpkey.txt'), 'w') as keyfile:
                keyfile.write(user)
            os.chmod(os.path.expanduser('~/.mpkey.txt'), stat.S_IREAD | stat.S_IWRITE)
            del user

            with open(os.path.expanduser('~/.mpkey.txt'),'r') as f:
                key = f.readlines()[0]
                return key
            print("Success")
        except:
            print("Something seems wrong with your key")
            
    
def get_qe_outputs(file):
    """
    Extract outputs (energies, forces, structures) from qe .stdout files
    
    inputs:
        file: path to the file we want to extract outputs from
    outputs:
        dict: dictionary of the extracted outputs
    """
    
    # TODO: this is very VERY hardcoded. you can do better, fix this (past kat to future kat)
    import numpy as np
    
    output = open(file, "r")
    lines = output.readlines()
    iE = [] # energy at each ionic step, Ry
    eE = [[]] # energy at each electronic step, Ry
    P = [] # pressure, kbar
    F = [] # total force, Ry/au
    stresses = [] # stress tensor, kbar
    structures = [] # pymatgen structure objects, angstrom

    from pymatgen.core import Lattice, Structure

    # Check for certain tags on lines, add variables to lists
    for i,line in enumerate(lines):
        if 'total energy' in line and '!' not in line and 'The' not in line:
            eE[-1].append(float(line.split()[3]))
        elif '!' in line:
            eE.append([])
            iE.append(float(line.split()[4]))
        elif 'P=' in line:
            P.append(float(line.split()[5]))
            stresses.append(np.array([lines[i+1].split()[3:6],lines[i+2].split()[3:6],lines[i+3].split()[3:6]]).astype(float))
        elif "Total force" in line:
            F.append(float(line.split()[3]))
        # TODO: come back and fix this, make it more robust 
        # figure out why QE only sometimes gives cell outputs
        elif 'CELL_PARAMETERS' in line:
            try:
                lattice = np.array([lines[i+1].split(),lines[i+2].split(),lines[i+3].split()]).astype(float)
                sites = []
                atoms = []
                j=6
                while ("End" not in lines[i+j]) and (lines[i+j]!=""):
                    sites.append(np.array(lines[i+j].split()[1:]).astype(float))
                    atoms.append(lines[i+j].split()[0])
                    j=j+1
                lattice_obj = Lattice(lattice)
                test_struct = Structure(lattice,atoms,sites)
                structures.append(test_struct)
            except:
                pass
    eE = eE[:-1]

    # return output dictionary
    return {'ionic_energies':iE,'electronic_energies':eE,'pressure':P,'forces':F,'stresses':stresses,'structures':structures}
    

def get_convergence_plots(step_dict, sim_name = ""):
    """
    Plot both ionic and electronic energy at each SCF step
    
    inputs:
        step_dict: dictionary output from get_qe_outputs()
        sim_name: optional name to add to title on plot
    outputs:
        n/a
    """
    
    import plotly.graph_objects as go
    import numpy as np
    # Extract the energies, lining up electronic and ionic steps
    i_energies = step_dict['ionic_energies']
    e_energies_array = step_dict['electronic_energies']
    e_count = [len(e) for e in e_energies_array]
    i_steps = [sum(e_count[0:n+1]) for n in range(len(e_count))]
    e_energies = [item for sublist in e_energies_array for item in sublist]
    e_steps = np.linspace(1,len(e_energies),len(e_energies))

    template='simple_white'
    # Create and save a plotly figure with the energy at each ionic and electronic step
    fig_energies = go.Figure()
    fig_energies.add_trace(go.Scatter(x = e_steps, y = e_energies, name = 'electronic'))
    fig_energies.add_trace(go.Scatter(x = i_steps, y = i_energies, name = 'ionic'))

    scaling_factor = 1.005
    fig_energies.update_layout(
        title = f'{sim_name} energy convergence',
        xaxis_title = 'electronic steps',
        yaxis_title = 'energy (Ry)',
        yaxis_range = [min(i_energies)*scaling_factor,max(i_energies)/scaling_factor],
        template=template
    )
    
    fig_energies.show()
    
    return True

            
def main():
    # Main loop 
    pass

if __name__ == "__main__":
    main()