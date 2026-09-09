import openmc
import os
import re
import time
import numpy as np

from scipy.interpolate import make_interp_spline
from openmc.data.endf import Evaluation, get_head_record, get_tab1_record
from io import StringIO

from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

def init_worker(indexs, alpha_S, E_as, integrents):
        global G_indexs, G_alpha_S, G_E_as, G_integrents
        G_indexs = indexs
        G_alpha_S = alpha_S
        G_E_as = E_as
        G_integrents = integrents


class Nuclide_integrent:
    def __init__(self, E_as, nuclide_name,spline):
        self.nuclide_name = nuclide_name

        self.E_n_maxs = [Alpha_N_calc.neutron_energy_alpha_n(E_a, 0, nuclide_name) for E_a in E_as]
        self.E_n_mins = [min(
            Alpha_N_calc.neutron_energy_alpha_n(E_a, np.pi / 2, nuclide_name),
            Alpha_N_calc.neutron_energy_alpha_n(E_a, np.pi, nuclide_name),
        ) for E_a in E_as]
        self.values = [spline(E_a) / (E_n_max - E_n_min) for E_a,E_n_max,E_n_min in 
                      zip(E_as,self.E_n_maxs,self.E_n_mins)]
    def __call__(self, E_n,i):
        if self.E_n_mins[i] < E_n < self.E_n_maxs[i]:
            return self.values[i]
        else:
            return 0  
            
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
            #print('energy span',E_min,E_max)
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
                    #print('end of file')
                    break
                if unit == 'keV':
                    energy = float(energy.replace(",", "."))
                    energies.append(energy*1000)
                elif unit=="MeV":
                    energy = float(energy.replace(",", "."))
                    energies.append(energy*1e6)
                else:
                    print("Something went wrong in line",n)
                val1 = float(val1.replace(",", "."))
                val2 = float(val2.replace(",", "."))
                stopping_powers.append(val1+val2)
                n+=1
        return energies,stopping_powers
        
    @staticmethod
    def SR_file_write_and_read(mat,E_min=1e4,E_max=1e7,shared_folder="/root/SR_Module",new_file=True):
        '''
        Method for write input and read output from SR module that calculates stopping power
            Waits until SR-module has been run in windows
        
        Parameters
        ----------
        mat: openmc.Material
            material, which needs stopping power calculated
        E_min: float
            minimum energy (eV)
        E_max: float
            maximum energy (eV)
        shared_folder: str
            path to the folder shared between windows and linux
        Returns
        -------
        energies :array[float]
            energies in eV
        stopping_powers: array[float]
            stopping_power in eV/Å
        '''
        E_min=int(E_min/1000)
        E_max=int(E_max/1000)
        if E_min<10:
            E_min=10
        
        input_file = os.path.join(shared_folder, "SR.IN")
        Alpha_N_calc.SR_file_write_IN(input_file,mat,E_min=E_min,E_max=E_max,state=0)
    
        output_file = os.path.join(shared_folder, mat.name)
        if new_file:
            try:
                os.remove(output_file)
            except:
                print('no file to remove')
        n=0        
        while not os.path.exists(output_file):
            time.sleep(1)  # check every second
            n+=1
            if new_file or n>5:
                print("Waiting for SRModule to produce output...")
    
        energies,stopping_powers = Alpha_N_calc.SR_file_read(mat.name,shared_folder)
        energies=np.array(energies)
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
        Parameters
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
    def isolate_atomic_symbol(nuclide_name):
        '''
        Takes nuclide name and isolates into symbol, Z, A as 'Li', '003' and '006' for Li6
        Parameters
        ----------
        nuclide_name: str
            name of the nuclide in the form of U234
        Returns
        -------
        symbol: str
            atomic symbol of the nuclide
        Z: str
            number of protons
        A: str
            number of nucleons
        '''
        #Isolating the atomic symbol
        Z, A, m =openmc.data.zam(nuclide_name)
    
        match = re.match(r"([A-Za-z]+)",nuclide_name)
        symbol =match.groups()[0]
        
        Z = Alpha_N_calc.make_AZ_Str(Z)
        A = Alpha_N_calc.make_AZ_Str(A)
        return symbol,Z,A
    
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

    def print_deacy_data(nuclide_name):
        '''
        print alpha decay data from nuclide name format U234
        ----------
        nuclide_name: str
            name of the nuclide in the form of U234
        Returns
        -------
        None
        '''
        dec = get_decay_data(nuclide_name)
        mass = openmc.data.atomic_mass(nuclide_name)/6.022136651E+26 #Convertion from unit to kg
        print(dec.modes)
        print('Halflife (s)\t',dec.half_life)
        print(f'Halflife (y)\t{dec.half_life.nominal_value/(3600*24*365.25):.3g}')
        print('Decay_constant\t',dec.decay_constant)
        #print('Calculated decay constant (s) ',np.log(2)/dec.half_life.nominal_value)
        
        print('\nEnergy \t\t\tIntensity\t\t decay')
        total = 0
        for item in dec.spectra['alpha']['discrete']:
            print(item['energy'],'\t',item['intensity'],'\t',item['intensity']*dec.decay_constant/mass)
            total+=item['intensity'].nominal_value
        print(f'Total intensity {total:.3g}')

    @staticmethod
    def alpha_decay_values_from_material(material,material_volume=1000,show_discarded=False):
        '''
        Method for finding aplha source for a given material. Using volume of the openmc.Material
        Parameters
        ----------
        material: openmc.Material
            material to analyze
        material_volume:float
            volume of material in cm^3 
        show_discarded: boolean
            enables the printing of every material discarded do to no alha emission
        
        Returns
        -------
        result:arrays[str,flaot,float,float,float]
            [nuclide_name,activities,energies,energy_sted_devs, mass in kg]
        material_mass: float
            total mass of the material
        nuclide_amount: float?
            total number of atoms in the cell with volume given by 
        None
        '''
        result = []
        material_mass = 0
        material_cell = openmc.Cell(1,'fuel')
        material_cell.fill = material
    
        #Getting the volume from the material
        material_cell.volume = material_volume

        
        
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
        return result, material_mass, nuclide_amount

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
            print('given x is used')

        energy_spectra = np.zeros(len(x), dtype=float)
        for name, activities, energies,energy_devs,_  in Spectra_data:
            for activity,energy,energy_std_dev in zip(activities, energies,energy_devs):
                #print(activity,energy_std_dev,energy)
                energy_spectra += Alpha_N_calc.gaussian(activity,energy_std_dev,energy,x)
        return energy_spectra,x
    
    
    @staticmethod
    def neutron_spectra_from_material(mat,E_min=0,E_max=5e6,nr_points=10000,nr_energies=500,printing=True,new_file=True):
        '''
        Method for calculating a continuous neutron spectrum from the
        alpha decay spectrum of a material and its (alpha,n) reaction
        cross sections.
    
        Parameters
        ----------
        mat : openmc.Material
            Material for which the neutron spectrum is calculated.
    
        E_min : float, optional
            Minimum alpha/neutron energy considered, in eV.
            Default is 0 eV.
    
        E_max : float, optional
            Maximum alpha/neutron energy considered, in eV.
            Default is 6e6 eV.
    
        nr_points : int, optional
            Number of energy points used for calculating and interpolating
            the stopping powers and (alpha,n) reaction cross sections.
            Default is 1000.
    
        nr_energies : int, optional
            Number of neutron-energy points at which the final neutron
            spectrum is evaluated.
            Default is 100.
    
        printing : bool, optional
            If True, print the time required for the setup calculation.
            Default is True.
    
        Returns
        -------
        neutron_spectra : np.ndarray
            Calculated neutron spectrum evaluated at each energy in
            ``neutron_energies``. The spectrum contains the contribution
            from (alpha,n) reactions in all nuclides present in ``mat``.
    
        neutron_energies : np.ndarray
            Neutron-energy grid corresponding to ``neutron_spectra``,
            in eV.
    
        Notes
        -----
        The calculation uses the alpha decay spectrum of the material,
        the stopping power of alpha particles in the material, and the
        (alpha,n) reaction cross sections of the constituent nuclides.
    
        The reaction channels are grouped according to the number of
        neutrons produced:
    
            MT = 4, 22, 23, 28, 29 : one neutron
            MT = 11, 16, 24, 30    : two neutrons
            MT = 17, 25            : three neutrons
    
        The calculation is parallelized over the neutron-energy grid
        using ``ProcessPoolExecutor``.
        '''
        #Getting data
        time_setup_start = time.perf_counter()
        energies = np.linspace(E_min,E_max,nr_points)
        
        mat_dens = mat.get_nuclide_atom_densities() #returns nuclide densities in atom/b-cm
    
        #Stopping power splines
        stopping_power_spline = Alpha_N_calc.get_stopping_power_spline(mat,E_min=E_min,E_max=E_max,new_file=new_file)
        stopping_powers= Alpha_N_calc.get_spline_data(stopping_power_spline,energies, extrapolate=True)
        stopping_powers = np.array(stopping_powers)*1e8 #Converting from eV/Å to eV/cm

        #Alpha spectrum
        alpha_decays, mass,nuclide_amount = Alpha_N_calc.alpha_decay_values_from_material(mat)
        alpha_spectra = []
        alpha_energies = []
        for name, a_activities, a_energies,energy_devs,_  in alpha_decays:
            alpha_spectra.extend(a_activities)
            alpha_energies.extend(a_energies)
        #sorting energies
        alpha_spectra = [a/mass for _,a in sorted(zip(alpha_energies,alpha_spectra))]
        alpha_energies = sorted(alpha_energies)
        
        integrent_splines = []
        
        for nuc in mat.nuclides:
            #Getting crosssections
            cross_splines,tab = Alpha_N_calc.get_crossection_from_isotope(nuc.name,k=5,x_max=E_max,plot=False)
            total_crosssection = np.zeros(nr_points)
            for mt in cross_splines:
                spline_data = Alpha_N_calc.get_spline_data(cross_splines[mt],energies, extrapolate=False,non_negative=True)
                
                if mt in [4,22,23,28,29]: #(a,n+anything)
                    total_crosssection += spline_data
                    
                elif mt in [11,16,24,30]: #(a,2n+anything)
                    total_crosssection += 2*spline_data
                    
                elif mt in [17,25]: #(a,3n+anything) for up to mt=30
                    total_crosssection += 3*spline_data
                else:
                    print(f'reaction with mt={mt} is not know')
    
            #Making integrent
            target_density = mat_dens[nuc.name] #atom/b-cm 
            
            #total_crosssection is in barns
            #stopping powers is in eV/cm
            
            integrent_splines.append(
                make_interp_spline(energies,target_density*total_crosssection/stopping_powers,k=5))

        #Setting energy array with all alpha energies and energies in between up
        energies_with_alpha = np.unique(np.concatenate((energies,alpha_energies)))
        indexs = [np.searchsorted(energies_with_alpha, E_a, side="right") for E_a in alpha_energies]
        
        #Making kernels and integrent 
        numeric_integrents = [Nuclide_integrent(energies_with_alpha,nuc.name,spline) 
                              for nuc,spline in zip(mat.nuclides,integrent_splines)]
        
    
        time_setup = time.perf_counter() - time_setup_start
        if printing:
            print(f'setup done in {time_setup} s')
    
        #Neutron energies
        neutrons_total = 0
        neutron_spectra = []
        neutron_energies =  np.linspace(E_min,E_max,nr_energies)

        init_worker(
            indexs,
            alpha_spectra,
            energies_with_alpha,
            numeric_integrents
        )
    
        with ProcessPoolExecutor(initializer=init_worker,
                                 initargs=(indexs, alpha_spectra, energies_with_alpha, numeric_integrents)) as executor:
            neutron_spectra = list(executor.map(Alpha_N_calc.compute_spectrum, neutron_energies))     
        return neutron_spectra,neutron_energies,neutrons_total    
        
    @staticmethod
    def compute_spectrum(E_n):#,indexs,alpha_spectra,energies_with_alpha,integrents):
        neutrons_total = 0
        values = [np.sum([integrent(E_n,i) for integrent in G_integrents]) for i in range(len(G_E_as))]
        
        for i,alpha in zip(G_indexs,G_alpha_S): #over all different alpha energies
            neutrons_total +=  alpha*np.trapezoid(values[:i], G_E_as[:i]) #neutrons
        return neutrons_total
    @staticmethod
    def neutron_energy_alpha_n(E_a,theta,nuclide_name):
        m_a = 4.001506179127#openmc.data.atomic_mass('He4') #alpha mass
        m_n = 1.008664915904 #neutron mass
        m_T = openmc.data.atomic_mass(nuclide_name)#target nuclide
        
        Z, A, m =openmc.data.zam(nuclide_name)
        new_Z = Z+2
        m_R = openmc.data.atomic_mass(f'{openmc.data.ATOMIC_SYMBOL[new_Z]}{A+3}') #residual recoil nucleus
        
        Q_i =m_a+m_T-m_n-m_R  #From energy/mass conservation in m_T+m_a = m_n-m_R
        
        term0 = m_a*m_n*E_a*np.pow(np.cos(theta),2)/np.pow(m_n+m_R,2)
        term1 = (m_R*Q_i+(m_R-m_a)*E_a)/(m_a+m_R)
        term2 = 2*np.cos(theta)/(m_n+m_R)*np.sqrt(m_a*m_n*E_a*(m_R*Q_i+(m_R-m_a)*E_a)/(m_n+m_R))
        
        return term0 + term1 + term2


        
    @staticmethod
    def get_stopping_power_spline(mat,k=5,E_min=1e4,E_max=1e7,new_file=True):
        energies,stopping_powers = Alpha_N_calc.SR_file_write_and_read(mat,E_min=E_min,E_max=E_max,new_file=new_file)
        while k>0:    
            try:    
                spline = make_interp_spline(energies,stopping_powers,k=k,bc_type='periodic')
                break
            except:
                    k-=1
            if k==1:
                spline = make_interp_spline(energies,stopping_powers,k=k)
        return spline
        
    @staticmethod
    def get_spline_data(spline,xs,extrapolate = False,non_negative=False):
        spline_with_nan = spline(xs, extrapolate=extrapolate)
        spline_data = np.nan_to_num(spline_with_nan,nan=0.00)
        
        if non_negative:
            spline_data = np.clip(spline_data,0,None) #sets all negative values to zero
    
        return spline_data

    @staticmethod
    def get_crossection_from_isotope(nuc_name,k=1,xs=None, x_min=0, x_max=8e6, plot=True,print_no_data=False):
        splines={}
        tabs={}
        failed_mts = []
        if plot:
            fig,ax = plt.subplots()
            ax.set_title(rf'{nuc_name} cross section for ($\alpha$,anything)')
            ax.set_xlabel('Energy [eV]')
            ax.set_ylabel(r'cross section [b]')
            
        for mt in np.array([4,16,17,22,23,24,25,28,29]):
            mt= int(mt)
            spline, tab ,xs = Alpha_N_calc.get_crosssection_spline(nuc_name,mt,k,xs,x_min,x_max)
                
            if spline != None:
                splines[mt]=spline
                tabs[mt]=tab
                if plot:
                    ax.plot(tab.x,tab.y,'x-',label=f'mt={mt}')
                    ax.plot(xs,spline(xs))
            else:
                failed_mts.append(mt)
        if len(failed_mts)>0 and print_no_data:  
            print(f'In {nuc_name} following mts had no data',failed_mts)
        if plot:
            ax.legend()
        return splines,tabs
    @staticmethod
    def get_crosssection_spline(nuc_name, mt,k=1, xs=None, x_min=2e6, x_max=8e6):
        #Getting data
        symbol,Z,A = Alpha_N_calc.isolate_atomic_symbol(nuc_name)
        ev = Evaluation(f"/root/TENDL-a/a-{symbol}{A}.tendl")
        
        mf=3 #Reaction cross sections
        try:
            endf_output = StringIO(ev.section[mf,mt])
        except:
            #print(f'mt={mt} does not exist')
            return None,None, None
        head = get_head_record(endf_output)
        params, tab = get_tab1_record(endf_output)
        #print('energy',tab.x)   # energies in eV
        #print('cross section',tab.y)   # cross sections in barns
    
        #Removing data under x_min, over x_max and duplicates 
        mask = (tab.x >= x_min) & (tab.x <= x_max)
        tab.x = tab.x[mask]
        tab.y = tab.y[mask]
        if len(tab.x) ==0:
            #print(f'no data between {x_min} eV and {x_max} eV')
            return None,None,None
        tab.x, idx = np.unique(tab.x, return_index=True)
        tab.y = tab.y[idx]
    
        if xs is None:
            xs = np.linspace(min(tab.x),max(tab.x),int(1e4))
        while k>0:    
            try:    
                spline = make_interp_spline(tab.x,tab.y,k=k)
                break
            except:
                k-=1
                if k==1:
                    spline = make_interp_spline(tab.x,tab.y,k=k)
        return spline, tab, xs