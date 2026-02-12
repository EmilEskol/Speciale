import openmc
def set_mat(name,materials=None):
    fuel_mat = None
    water = None
    heavy_water = None
    blanket_CA = None
    air = None
    
    if materials is None:
        openmc.reset_auto_ids()
        materials = openmc.Materials.from_xml('materials.xml')
    
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
        elif mat.name=="Blanket CA":
            blanket_CA = mat
        
    return fuel_mat, heavy_water, water, blanket_CA, air, materials
    
def get_material_by_name(name,materials):
    found_mat = None
    for mat in materials:
        if mat.name==name:
            found_mat = mat
    if found_mat is None:
        print(f'No material was found with the name {name}')
    return found_mat

def make_settings(geom):
    inner_fuel_radius = geom.get_surfaces_by_name('fuel_inner')[0].bounding_box('-')[1][0]
    outer_fuel_radius = geom.get_surfaces_by_name('fuel_outer')[0].bounding_box('-')[1][0]

    settings = openmc.Settings()
    settings.run_mode = 'eigenvalue'
    point = openmc.stats.Point()
    src = source = openmc.IndependentSource(
    space=openmc.stats.Point((0, 0, 
                              inner_fuel_radius+(outer_fuel_radius-inner_fuel_radius)/2)),
                            energy=openmc.stats.Discrete([14e6], [1.0]))

    settings.source = src
    settings.batches = 100
    settings.inactive = 10
    settings.particles = 10000
    return settings

def set_geom_fuel_shell(sphere_ir,sphere_or,water_or,with_blanket=False,blanket_or=0,inner_moderator_height=-1
                        , materials= None,fuel_name="Fuel CA_1"):
    fuel_material, heavy_water, water, blanket_CA, air, mat = set_mat(fuel_name, materials)
    eps=0.0001 # constant to make sure geomitries dont overlap
    
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
    water_outer_sphere= openmc.Sphere(r=water_or)
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
        blanket_outer_sphere = openmc.Sphere(r=blanket_or, name='outer_sphere')
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