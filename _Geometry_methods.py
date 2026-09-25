import openmc
import numpy as np

class Geometry_helper:
    '''Class with useful openmc methods'''
    @staticmethod
    def set_mat(name,materials=None):
        fuel_mat = None
        water = None
        heavy_water = None
        blanket_CA = None
        air = None
        
        if materials is None:
            print('redefining materials')
            openmc.reset_auto_ids()
            materials = openmc.Materials.from_xml('materials.xml')

        print('Materials defined: ')
        for mat in materials:
            print(mat.name)
            if mat.name==name:
                fuel_mat = mat
            elif mat.name=="h2o":
                water = mat
            elif mat.name=="d2o":
                heavy_water = mat
            elif mat.name == "Air":
                air = mat
            elif mat.name=="Blanket_CA":
                blanket_CA = mat
            
        return fuel_mat, heavy_water, water, blanket_CA, air, materials
    @staticmethod    
    def get_material_by_name(name,materials):
        found_mat = None
        for mat in materials:
            if mat.name==name:
                found_mat = mat
        if found_mat is None:
            print(f'No material was found with the name {name}')
        return found_mat
    @staticmethod
    def create_materials(fuel_name):
        def norm_array(array):
            return np.array(array)/sum(array)

        def weight_to_atompercent(nuclides,weights):
            #Takes array of weights of different isotopes and returns a atom percent for the nuclide
            result = []
            nr_of_atoms = []
            
            for nuc,weight in zip(nuclides,weights):
                atomic_mass = openmc.data.atomic_mass(nuc)#*1.66053906892e-27
                nr_of_atoms.append(weight/atomic_mass)
        
            total_nr_of_atoms = sum(nr_of_atoms)
            
            for atom in nr_of_atoms:
                result.append(atom/total_nr_of_atoms)
            return np.array(result)
        #Constants
        lbft_to_gcm = 1.6018463374/100

        #Making heavy water and air
        heavy_water = openmc.Material(name="d2o")
        heavy_water.add_nuclide('H2', 2.0)
        heavy_water.add_nuclide('O16', 1.0)
        
        heavy_water.set_density('g/cm3', 1.104) #from CA 2 table 1
        heavy_water.temperature = 26.8 + 273.15 #in Kelvin
        heavy_water.depletable = False
        
        heavy_water.add_s_alpha_beta('c_D_in_D2O')

        #Definition of air
        N2 = openmc.Material(name='N2')
        N2.add_elements_from_formula(formula='N2')
        
        O2 = openmc.Material(name='O2')
        O2.add_elements_from_formula(formula='O2')
        
        Ar = openmc.Material(name='Ar')
        Ar.add_elements_from_formula(formula='Ar')

        air = openmc.Material.mix_materials([N2,O2,Ar],[0.7811,0.2096,0.0093],'vo',name='Air')
        
        #Making salts used in fuels and blanket
        LiF = openmc.Material(name='LiF')
        LiF.add_elements_from_formula(formula='LiF',percent_type='ao',enrichment=99.9926,
                                         enrichment_target='Li7',enrichment_type='ao')
        
        LiF_CA = openmc.Material(name='LiF_CA') #LiF with pure Li7
        LiF_CA.add_elements_from_formula(formula='LiF',percent_type='ao',enrichment=99.999,
                                 enrichment_target='Li7',enrichment_type='ao')
        
        ThF4 = openmc.Material(name='ThF4')
        ThF4.add_elements_from_formula('ThF4') 

        UF4_LEU = openmc.Material(name='UF4_LEU')
        UF4_LEU.add_elements_from_formula('UF4',percent_type='ao',enrichment=4.95,enrichment_type='wo')
    
        #Define Copenhagen Atomic blanket salt
        blanket_CA = openmc.Material.mix_materials([LiF,ThF4],norm_array([0.70*2,0.30*5]),'ao',name="Blanket_CA")
        #from CA-2 table 1
        blanket_CA.set_density('g/cm3',4.626)
        blanket_CA.temperature = 626.8 + 273.15 # in Kelvin
        blanket_CA.depletable = True
        
        if fuel_name == 'Fuel_CA_1':
            fuel = openmc.Material.mix_materials([LiF_CA,UF4_LEU,ThF4],
                                                 norm_array([0.73*2,0.23*5,0.04*5]),
                                                 'ao',name="Fuel_CA_1") #from CA_2 table 1
            fuel.set_density('g/cm3',4.905)
            fuel.temperature = 626.8 + 273.15 # in Kelvin
            fuel.depletable = True
        elif fuel_name == 'Fuel_CA_2':
            fuel = openmc.Material.mix_materials([LiF_CA,UF4_LEU],
                                                 norm_array([0.73*2,0.27*5]),
                                                 'ao',name="Fuel_CA_2") #from CA_2 table 1
            fuel.set_density('g/cm3',4.905)
            fuel.temperature = 626.8 + 273.15 # in Kelvin
            fuel.depletable = True
        elif fuel_name == 'Fuel_A':
            #Making UF4 for fuel
            UF4_A = openmc.Material(name='UF4_A')
            #Using mass in kg from ORNL-TM-0611 p.7 to find atom percent 
            atom_percent_U = weight_to_atompercent(['U234','U235','U236','U238'],[0.3,27,0.3,1.5])
            atom_percent_UF4 = atom_percent_U*0.2
            print("ao for UF4 in fuel A",atom_percent_UF4)
            UF4_A.add_nuclide('U234',atom_percent_UF4[0],percent_type='ao')
            UF4_A.add_nuclide('U235',atom_percent_UF4[1],percent_type='ao')
            UF4_A.add_nuclide('U236',atom_percent_UF4[2],percent_type='ao')
            UF4_A.add_nuclide('U238',atom_percent_UF4[3],percent_type='ao')
            UF4_A.add_element('F',0.8,percent_type='ao')

            #Making fuel
            fuel = openmc.Material.mix_materials([LiF,BeF2,ZrF4,ThF4,UF4_A]
                                                 ,norm_array([0.69987*2,0.237*3,0.05*5,0.01*5,0.00313*5])
                                                 ,'ao',name='Fuel_A')
            fuel.set_density('g/cm3',144.5*lbft_to_gcm)
            fuel.depletable = True

        else:
            print(f'No fuel with the name {fuel_name}')

        materials = openmc.Materials([fuel, blanket_CA, heavy_water,air])
        materials.export_to_xml()
        return materials
    
    @staticmethod
    def make_settings(geom):
        inner_fuel_radius = geom.get_surfaces_by_name('fuel_inner')[0].bounding_box('-')[1][0]
        outer_fuel_radius = geom.get_surfaces_by_name('fuel_outer')[0].bounding_box('-')[1][0]
    
        settings = openmc.Settings()
        settings.run_mode = 'eigenvalue'
        #point = openmc.stats.Point()
        src = source = openmc.IndependentSource(
                                space=openmc.stats.Point((0, 0, 
                                (inner_fuel_radius+outer_fuel_radius)/2)),
                                energy=openmc.stats.Discrete(np.linspace(1e6,14e6,100), np.ones(100) / 100))
    
        settings.source = src
        settings.batches = 100
        settings.inactive = 10
        settings.particles = 10000
        return settings
    @staticmethod
    def set_geom_fuel_shell(sphere_ir,sphere_or,water_or,with_blanket=False,blanket_or=0,inner_moderator_height=-1
                            , materials = None,fuel_name = None):
        def shell_vol_calc(outer_radius,inner_radius):
            volume=4/3*np.pi*(outer_radius**3-inner_radius**3)
            return volume
        def sphere_cap_vol(r,h):
            return 1/3*np.pi*h**2*(3*r-h)

        fuel_material, heavy_water, water, blanket_CA, air, mat = Geometry_helper.set_mat(fuel_name, materials)
        eps=0.0001 # constant to make sure geomitries dont overlap

        fuel_material.volume = shell_vol_calc(sphere_or,sphere_ir)
        fuel_or = openmc.Sphere(r=sphere_or, name='fuel_outer' ) #fuel cylinder with outer radius
        fuel_ir = openmc.Sphere(r=sphere_ir, name='fuel_inner') #fuel cylinder with outer radius
        fuel_region = +fuel_ir & -fuel_or    #inside the fuel fuel
        fuel_cell = openmc.Cell(1, 'fuel')
        fuel_cell.fill = fuel_material
        fuel_cell.region = fuel_region
    
        #Define moderator area
        if inner_moderator_height == -1:
            inner_moderator_height = 2*sphere_ir
        inner_moderator_height_plane = openmc.Plane(0,1,0,inner_moderator_height-sphere_ir)
        
        water_inner_sphere= openmc.Sphere(r=sphere_ir)
        water_middle_sphere= openmc.Sphere(r=sphere_or)
        water_outer_sphere= openmc.Sphere(r=water_or, name='blanket_inner')
        air_region = -water_inner_sphere & +inner_moderator_height_plane
        water_region1 = -water_inner_sphere & -inner_moderator_height_plane
        water_region2 = +water_middle_sphere & -water_outer_sphere
    
        moderator_cell = openmc.Cell(2,'moderator')
        moderator_cell.fill = heavy_water
        moderator_cell.region = water_region1 | water_region2
    
        air_cell = openmc.Cell(5,'air')
        air_cell.fill = air
        air_cell.region = air_region
    
        
        #If blanket is a part of the geometry
        if with_blanket:
            blanket_CA.volume = shell_vol_calc(blanket_or,water_or)
            blanket_outer_sphere = openmc.Sphere(r=blanket_or, name='blanket_outer')
            blanket_region = +water_outer_sphere & -blanket_outer_sphere
            blanket_cell = openmc.Cell(3,'blanket')
            blanket_cell.region = blanket_region
            blanket_cell.fill = blanket_CA
    
            boundary = openmc.Sphere(r=blanket_or+eps,boundary_type='reflective')
            outer_cell = openmc.Cell(4,region=+blanket_outer_sphere & -boundary) 
                #Outer_cell is only to sepeate the boundary a bit from the geometry
            root = openmc.Universe(cells=[fuel_cell, moderator_cell,blanket_cell,outer_cell,air_cell])
        
        else:
            #we define boundary condition
            boundary = openmc.Sphere(r = water_or+eps, boundary_type='reflective', name='outer_sphere')
            outer_cell = openmc.Cell(4,region = +water_outer_sphere & -boundary)
            root = openmc.Universe(cells=[fuel_cell, moderator_cell,outer_cell,air_cell])
        
    
        geom = openmc.Geometry()
        geom.root_universe = root
    
        return geom, mat

        @staticmethod
        def shell_vol_calc(outer_radius,inner_radius):
        
            volume=4/3*np.pi*(outer_radius**3-inner_radius**3)
            return volume
        
        @staticmethod
        def sphere_cap_vol(r,h):
            return 1/3*np.pi*h**2*(3*r-h)