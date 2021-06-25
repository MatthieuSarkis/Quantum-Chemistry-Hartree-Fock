def xyz_reader(file_name):
    file = open(file_name, 'r')

    number_of_atoms = 0
    atom_type = []
    atom_coordinates = []

    for idx, line in enumerate(file):
        if idx == 0:
            try:
                number_of_atoms = int(line.split()[0])
            except:
                print("file not in correct format.")
        
        if idx == 1:
            continue

        if idx != 0:
            split = line.split()
            atom = split[0]
            coordinates = [float(split[1]),
                           float(split[2]),
                           float(split[3])]
            atom_type.append(atom)
            atom_coordinates.append(coordinates)

    file.close()

    return number_of_atoms, atom_type, atom_coordinates