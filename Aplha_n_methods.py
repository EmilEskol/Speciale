import openmc
import os
import re
import time
import numpy as np

class Alpha_N_calc:
    '''Class containning methods used in alpha_N claculations'''
    
    symbol_to_element_name = {
        "H": "Hydrogen", "He": "Helium", "Li": "Lithium", "Be": "Beryllium",
        "B": "Boron", "C": "Carbon", "N": "Nitrogen", "O": "Oxygen",
        "F": "Fluorine", "Ne": "Neon", "Na": "Sodium", "Mg": "Magnesium",
        "Al": "Aluminum", "Si": "Silicon", "P": "Phosphorus", "S": "Sulfur",
        "Cl": "Chlorine", "Ar": "Argon", "K": "Potassium", "Ca": "Calcium",
        "Sc": "Scandium", "Ti": "Titanium", "V": "Vanadium", "Cr": "Chromium",
        "Mn": "Manganese", "Fe": "Iron", "Co": "Cobalt", "Ni": "Nickel",
        "Cu": "Copper", "Zn": "Zinc", "Ga": "Gallium", "Ge": "Germanium",
        "As": "Arsenic", "Se": "Selenium", "Br": "Bromine", "Kr": "Krypton",
        "Rb": "Rubidium", "Sr": "Strontium", "Y": "Yttrium", "Zr": "Zirconium",
        "Nb": "Niobium", "Mo": "Molybdenum", "Tc": "Technetium", "Ru": "Ruthenium",
        "Rh": "Rhodium", "Pd": "Palladium", "Ag": "Silver", "Cd": "Cadmium",
        "In": "Indium", "Sn": "Tin", "Sb": "Antimony", "Te": "Tellurium",
        "I": "Iodine", "Xe": "Xenon", "Cs": "Cesium", "Ba": "Barium",
        "La": "Lanthanum", "Ce": "Cerium", "Pr": "Praseodymium", "Nd": "Neodymium",
        "Pm": "Promethium", "Sm": "Samarium", "Eu": "Europium", "Gd": "Gadolinium",
        "Tb": "Terbium", "Dy": "Dysprosium", "Ho": "Holmium", "Er": "Erbium",
        "Tm": "Thulium", "Yb": "Ytterbium", "Lu": "Lutetium", "Hf": "Hafnium",
        "Ta": "Tantalum", "W": "Tungsten", "Re": "Rhenium", "Os": "Osmium",
        "Ir": "Iridium", "Pt": "Platinum", "Au": "Gold", "Hg": "Mercury",
        "Tl": "Thallium", "Pb": "Lead", "Bi": "Bismuth", "Po": "Polonium",
        "At": "Astatine", "Rn": "Radon", "Fr": "Francium", "Ra": "Radium",
        "Ac": "Actinium", "Th": "Thorium", "Pa": "Protactinium", "U": "Uranium",
        "Np": "Neptunium", "Pu": "Plutonium", "Am": "Americium", "Cm": "Curium",
        "Bk": "Berkelium", "Cf": "Californium", "Es": "Einsteinium", "Fm": "Fermium",
        "Md": "Mendelevium", "No": "Nobelium", "Lr": "Lawrencium", "Rf": "Rutherfordium",
        "Db": "Dubnium", "Sg": "Seaborgium", "Bh": "Bohrium", "Hs": "Hassium",
        "Mt": "Meitnerium", "Ds": "Darmstadtium", "Rg": "Roentgenium", "Cn": "Copernicium",
        "Nh": "Nihonium", "Fl": "Flerovium", "Mc": "Moscovium", "Lv": "Livermorium",
        "Ts": "Tennessine", "Og": "Oganesson"
    }

    @staticmethod
    def SR_file_write_IN (input_file,mat,E_min=10,E_max=10000,state=0):
        '''
        Method for making txt for .IN files used by SR module
        Parameters
        ----------
        mat: openmc.Materiael
            The material used
        E_min: float
            minimum energy (keV)
        E_max: float
            maximum energy (keV)
        state: int
            solid=0, gas=1
        
        Returns
        -------
        None
        '''
        empty_string = ''
        with open(input_file, "w",newline="\n") as f:
            f.write('---Stopping/Range Input Data (Number-format: Period = Decimal Point)\r\n')
            f.write('---Output File Name\r\n')
            f.write(f'"{mat.name}"\r\n')
            f.write('---Ion(Z), Ion Mass(u)\r\n')
            f.write('2\t4\r\n')  # Helium-4 ion
            f.write('---Target Data: (Solid=0,Gas=1), Density(g/cm3), Compound Corr.\r\n')
            f.write(f'{state}\t {mat.density}\t 0\r\n')  # Note single space after first column
            f.write('---Number of Target Elements\r\n')
            f.write(f'{len(mat.nuclides)}\r\n')
            f.write('---Target Elements: (Z), Target name, Stoich, Target Mass(u)\r\n')
            
            for nuc in mat.nuclides:
                Stoich = nuc[1]
                Z, A, m = openmc.data.zam(nuc[0])
                symbol = re.match('([A-Za-z]+)', nuc[0])[0]
                element_name = Alpha_N_calc.symbol_to_element_name[symbol]
                # Match the working spacing exactly:
                # - 3 spaces after Z
                # - element name padded to 15 characters
                # - tab between name and stoich
                # - spaces between stoich and mass number
                f.write(f'{Z:<3}  "{Alpha_N_calc.symbol_to_element_name[symbol]}"  {empty_string:<8}\t {Stoich:<25}  {A}\r\n')
            
            f.write('---Output Stopping Units (1-8)\r\n')
            f.write('1\r\n')  # Units eV / Angstrom
            f.write('---Ion Energy : E-Min(keV), E-Max(keV)\r\n')
            print('energy span',E_min,E_max)
            f.write(f'{E_min}\t{E_max}')
            f.write(f'\r\n\r\n\r\n\r\n\r\n\r\n')

    
    @staticmethod
    def SR_file_read(fuel_name,shared_folder):
        '''
        Method for reading an file made from SR module
        Parameters
        ----------
        fuel_name: str
            name of the file, which is the fuel name when using make_SR_in
        shared_folder: str
            path to the folder shared between windows and linux
        
        Returns
        -------
        energies :array[float]
            energy array in eV
        stopping_powers: array[float]
            in eV/Angstrom
        '''
        
        #Import file
        output_file = os.path.join(shared_folder, fuel_name)
        
        energies = [] #in Mev
        stopping_powers = [] #in eV/Angstrom
        # Step 3: Read output
        with open(output_file, "r") as f:
            while True:
                line = f.readline()
                try:
                    energy, unit, val1, val2, a1, unit1, a2, unit2, a3, unit3 = line.split()
                    break
                except:
                    n=0
            while n < 100:
                line = f.readline()
                #print(line)
                try:
                    energy, unit, val1, val2, a1, unit1, a2, unit2, a3, unit3 = line.split()
                except:
                    print('end of file')
                    break
                if unit == 'keV':
                    energy = float(energy.replace(",", "."))
                    energies.append(energy*1000)
                elif unit=="MeV":
                    energy = float(energy.replace(",", "."))
                    energies.append(energy*1e6)
                else:
                    print("Something went wrong in line",results)
                val1 = float(val1.replace(",", "."))
                val2 = float(val2.replace(",", "."))
                stopping_powers.append(val1+val2)
                n+=1
        return energies,stopping_powers
        
    @staticmethod
    def SR_file_write_and_read(mat,E_min=10,E_max=10000,shared_folder="/root/SR_Module",new_file=True):
        '''
        Method for write input and read output from SR module that calculates stopping power
            Waits until SR-module has been run in windows
        
        Parameters
        ----------
        mat: openmc.Material
            material, which needs stopping power calculated
        E_min: float
            minimum energy (keV)
        E_max: float
            maximum energy (keV)
        shared_folder: str
            path to the folder shared between windows and linux
        Returns
        -------
        energies :array[float], stopping_powers: array[float]
        '''
        
        input_file = os.path.join(shared_folder, "SR.IN")
        Alpha_N_calc.SR_file_write_IN(input_file,mat,E_min=E_min,E_max=E_max,state=0)
    
        print("Waiting for SRModule to produce output...")
        output_file = os.path.join(shared_folder, mat.name)
        if new_file:
            os.remove(output_file)
        while not os.path.exists(output_file):
            time.sleep(1)  # check every second
    
        energies,stopping_powers = Alpha_N_calc.SR_file_read(mat.name,shared_folder)
        return energies,stopping_powers
    
    @staticmethod
    def make_AZ_Str(AZ):
        '''
        Method for making sure that number has 3 digits as in endf files and turns the numbers into strings
        Parameters
        ----------
        AZ: int
            either takes atom number (Z) or mass number (A)
        Returns
        -------
        AZstr: str
        '''
        #
        if AZ < 10:
            AZstr = f"00{AZ}"
        elif AZ < 100:
            AZstr = f"0{AZ}"
        else:
            AZstr = f"{AZ}"
        return AZstr
    @staticmethod
    def has_alpha_decay(nuclide):
        '''
        Checking if atoms have alpha decay
        ----------
        nuclide: [str,...]
            name of the nuclide in the form of U234 needs to be first entry in a array
        Returns
        -------
            boolean
        '''
        nuclide_name=nuclide[0]
        dec = Alpha_N_calc.get_decay_data(nuclide_name)
        for mode in dec.modes:
            if mode.modes[0] == "alpha":
                return True
        return False  
    
    @staticmethod
    def get_decay_data(nuclide_name):
        '''
        Getter for decay data from endf files using nuclide name
        ----------
        nuclide_name: str
            name of the nuclide in the form of U234
        Returns
        -------
            dec: openmc.data.Decay
        '''
        #Isolating the atomic symbol
        Z, A, m =openmc.data.zam(nuclide_name)
        match = re.match(r"([A-Za-z]+)(\d+)(_)(m\d+)",nuclide_name)
        m1=0
        if match == None:
            match = re.match(r"([A-Za-z]+)",nuclide_name)
            symbol =match.groups()[0]
        else:
            symbol, _,_ ,m1 =match.groups()
        
        Z = Alpha_N_calc.make_AZ_Str(Z)
        A = Alpha_N_calc.make_AZ_Str(A)
    
        #Getting decay data from ground state or excited (m1)
        try:
            if m1 != 0:
                dec = openmc.data.Decay.from_endf(f"../endf-b-vii.1/decay/dec-{Z}_{symbol}_{A}{m1}.endf")
            else:
                dec = openmc.data.Decay.from_endf(f"../endf-b-vii.1/decay/dec-{Z}_{symbol}_{A}.endf")
        except Keyerror:
            print("ERROR",Keyerror)
            dec = None
        return dec

    @staticmethod
    def alpha_decay_values_from_material(material,show_discarded=False):
        '''
        Method for finding aplha source for a given material. Using volume of the openmc.Material
        Parameters
        ----------
        material: openmc.Material
            material to analyze
        geom: openmc.Geometry
            geometry for the given problem used to calculated the amount of atoms
        cell_name: str
            name of the openmc.Geometry cell the material is in
        show_discarded: boolean
            enables the printing of every material discarded do to no alha emission
        
        Returns
        -------
        result:arrays[str,flaot,float,float,float]
            [nuclide_name,activities,energies,energy_sted_devs, mass in kg]
        material_mass: float
            total mass of the material
        None
        '''
        result = []
        material_mass = 0
        material_cell = openmc.Cell(1,'fuel')
        material_cell.fill = material
    
        #Getting the volume from the material
        try:
            material_cell.volume = material.volume
        except Exception as e:
            print(f'Error no material volume defined {e}')
        
        
        for nuclide_name, nuclide_amount_percent,_ in material.nuclides:
            Z, A, m =openmc.data.zam(nuclide_name)
            nuclide_amount = material_cell.atoms[nuclide_name]
            nuclide_mass = openmc.data.atomic_mass(nuclide_name)
            total_mass = nuclide_mass*nuclide_amount
            material_mass += total_mass
            
            dec = Alpha_N_calc.get_decay_data(nuclide_name)
    
            #Seeing if decay constant is possible to return
            try:
                decay_constant = dec.decay_constant.nominal_value #log(2)/half_life
            except ValueError as e:
                if show_discarded:
                    print(f"Skipping {dec.nuclide['name']} amount: {nuclide_amount:.3g}: {e}")
                    decay_constant = 0
            
            #Finding the spectra and calculating the activity
            if len(dec.spectra)!=0 and Alpha_N_calc.has_alpha_decay([nuclide_name]):
                alpha_data=dec.spectra['alpha']
                energies = [float(item['energy'].nominal_value) for item in alpha_data['discrete']] #can be continuous or discrete
                energy_std_devs = [float(item['energy'].std_dev) for item in alpha_data['discrete']]
                intensities  = [float(item['intensity'].nominal_value) for item in alpha_data['discrete']]
    
                #This is total activity
                activities = [decay_constant*nuclide_amount*intensity for intensity in intensities] #Calculation of activities
                result.append([nuclide_name,activities,energies,energy_std_devs, total_mass*1.6605402e-27])
            else:
                #Prints discarded decays due to no or very small activity or no data
                if show_discarded:
                    activity=decay_constant*nuclide_amount
                    print(dec.nuclide['name'],f"activity: {activity:.3g} has no alpha spectra")
    
        material_mass = material_mass*1.6605402e-27
        return result, material_mass

    @staticmethod
    def gaussian(A,a,b,x):
        '''
        Normalized gaussian function scaled with A
        
        Parameters
        ----------
        A: float
            scaling of the normalized guassian
        a: float
            deviation of the gaussian
        b: float
            center of the peak in the gaussian
        x: np.array([float])
    
        Returns
        -------
        fx: np.array([float])
            values for the gaussian function
        '''
        fx=A/(a*np.sqrt(2*np.pi))*np.exp(-(x-b)**2/(2*a**2))
        np.array(fx)
        return fx
    @staticmethod
    def energy_spectra_gaussian(Spectra_data,x=None,min_lim = 3.5e6,max_lim = 5e6,number_of_points = 10000):
        '''
        Method for calculating a continius spectra given the activity, energy and energy deviation 
        Parameters
        ----------
        Spectra_data: arrays[array[str],array[flaot],array[float],array[float],array[float]]
            data from alpha_decay_values_from_material
        min_lim: float
            limit for the lowest value of the gaussian spectra
        max_lim: float
            limit for the highest value of the gaussian spectra
        number_of_points: int
            numbar of points in the spectra
        Returns
        -------
        fx:np.array([float])
            the normilized gaussian spectra of the data given
        '''
        if x is None:
            x = np.linspace(min_lim,max_lim,int(number_of_points))
        elif sum(x)<100:
            x = np.array(x)*1e6
        else:
            x = np.array(x)*1e6
        energy_spectra = np.zeros(len(x), dtype=float)
        for name, activities, energies,energy_devs,_  in Spectra_data:
            for activity,energy,energy_std_dev in zip(activities, energies,energy_devs):
                #print(activity,energy_std_dev,energy)
                energy_spectra += Alpha_N_calc.gaussian(activity,energy_std_dev,energy,x)
        return energy_spectra,x