# Import statements
import numpy as np
from scipy import integrate

# Global object constants [m]
# Total object
total_length = 60
total_width = 60
total_height = 18

# Beams
beam_length = 60
beam_width = 12.65
beam_height = 6.25

# Pillars
pillar_length = 10.75
pillar_width = pillar_length
pillar_height = total_height - beam_height

# Struts
strut_length = 2
strut_width = total_width - 2 * beam_width
strut_height = 2

# Other global constants
grav_constant = 9.81  # [m/s^2]

# Object lists
pillar_list = [1, 2, 3, 4]
beam_list = [5, 6]
strut_list = [7, 8]


# Object normals
def object_normals():
    """
    Loads the normals of the surfaces on the object.
    :return: normal parameters of the object.
    """

    # Universal normals
    normals = np.array([[0, 0, 1],
                        [1, 0, 0],
                        [0, 1, 0],
                        [-1, 0, 0],
                        [0, -1, 0],
                        [0, 0, -1]])

    return normals


# Centroid modifiers
def centroid_modifiers(length, width, height):
    """
    Loads the centroid modifiers of the object.
    :param length: given length of the object.
    :param width: given width of the object.
    :param height: given height of the object.
    :return: modifiers_centroids: the modifiers required to calculate the surface centroids.
    """

    # Centroid modifiers
    half_length = length / 2
    half_width = width / 2
    half_height = height / 2

    # Universal centroid modifications
    modifiers_centroids = np.array([[-half_length, -half_width, 0],
                                    [0, -half_width, -half_height],
                                    [-half_length, 0, -half_height],
                                    [-length, -half_width, -half_height],
                                    [-half_length, -width, -half_height],
                                    [-half_length, -half_width, -height]])

    return modifiers_centroids


# Start / end modifiers for integration
def start_end_modifiers(length, width, height):
    """
    Loads the start / end integration modifiers of the object.
    :param length: given length of the object.
    :param width: given width of the object.
    :param height: given height of the object.
    :return: surface_start_modifiers, surface_end_modifiers: modifiers to calculate the start / end points of
                                                             the surfaces for integration.
    """

    # Universal start modifiers
    surface_start_modifiers = np.array([[-length, -width, 0],
                                        [0, -width, -height],
                                        [-length, 0, -height],
                                        [-length, -width, -height],
                                        [-length, -width, -height],
                                        [-length, -width, -height]])

    # Universal end modifiers
    surface_end_modifiers = np.array([[0, 0, 0],
                                      [0, 0, 0],
                                      [0, 0, 0],
                                      [-length, 0, 0],
                                      [0, -width, 0],
                                      [0, 0, -height]])

    return surface_start_modifiers, surface_end_modifiers


# Applies modifiers to base distance
def modifier_applier(base_distance, length, width, height):
    """
    Applies modifiers to a base distance for an object.
    :param base_distance: given base distance of the object (corner in first octant).
    :param length: given length of the object.
    :param width: given width of the object.
    :param height: given height of the object.
    :return: centroids, starts, ends: the centroid, start and end points of surfaces on the object.
    """

    # Modifiers for all six surfaces of the object
    modifiers_centroid = centroid_modifiers(length, width, height)
    modifiers_start, modifiers_end = start_end_modifiers(length, width, height)

    # Empty matrices to store data
    centroids = np.zeros(np.shape(modifiers_centroid))
    starts = np.zeros(np.shape(modifiers_start))
    ends = np.zeros(np.shape(modifiers_end))

    # Loop to apply modifiers
    for n in range(6):
        centroids[n] = modifiers_centroid[n] + base_distance
        starts[n] = modifiers_start[n] + base_distance
        ends[n] = modifiers_end[n] + base_distance

    # Return the centroid, start and end locations for all surfaces of a given object
    return centroids, starts, ends


# Surface parameters calculation
def surface_parameters(object_number):
    """
    Loads the different surface parameters (centroids, starts and ends) for a given object.
    :param object_number: int describing object type. 1-4 pillar, 5-6 beam, 7-8 strut.
    :return: area and normal parameters of the object.
    """

    # Object selection:
    # Pillars
    if object_number in pillar_list:
        print("Selected a pillar, with object number", object_number)
        # Base distance per pillar, to the base point
        short = 30 - pillar_width  # The reduced distance from COG if it includes a pillar width [m]
        pillar_base_distances = np.array([[30, 30, 1],
                                          [-short, 30, 1],
                                          [-short, -short, 1],
                                          [30, -short, 1]])

        # The properties of the pillar
        length = pillar_length
        width = pillar_width
        height = pillar_height

        return modifier_applier(pillar_base_distances[object_number - 1], length, width, height)

    # Beams
    elif object_number in beam_list:
        print("Selected a beam, with object number", object_number)
        # Base distances per beam, to the base point
        half_length = 0.5 * beam_length
        beam_top_depth = 1 - 18 + beam_height
        beam_base_distances = np.array([[half_length, half_length, beam_top_depth],
                                        [half_length, -half_length + beam_width, beam_top_depth]])

        # The properties of the beam
        length = beam_length
        width = beam_width
        height = beam_height

        return modifier_applier(beam_base_distances[object_number - 5], length, width, height)

    # Struts
    elif object_number in strut_list:
        print("Selected a strut, with object number", object_number)
        # Base distances per beam
        inter_strut_distance = 0.5 * (total_length - pillar_width)
        vert_strut_distance = 1 - total_height + 7 + 1
        strut_base_distances = np.array([[inter_strut_distance, 0, vert_strut_distance],
                                         [-inter_strut_distance, 0, vert_strut_distance]])

        # The properties of the strut
        length = strut_length
        width = strut_width
        height = strut_height

        return modifier_applier(strut_base_distances[object_number - 7], length, width, height)

    # Error if int
    elif object_number == int:
        raise Exception("Outside of range. Please choose 1-4 for pillar, 5-6 for beam, or 7-8 for strut.")

    # Error if not int
    else:
        raise Exception("Please choose an integer. Choose 1-4 for pillar, 5-6 for beam, or 7-8 for strut.")


# Integration type decision
def integration_type(starts, ends):
    """
    Decides the integration type for a given surface. Chooses x-y, x-z, or y-z.
    :param starts: array describing start coordinates for an object.
    :param ends: array describing end coordinates for an object.
    :return: int_type: array describing the integration type required for the surfaces on an object.
    """

    int_type = np.zeros(6)

    # Loop through all surfaces
    for n in range(6):
        # Difference in start and end coordinates
        difference = starts[n] - ends[n]

        # Same x-coordinates; y-z
        if difference[0] == 0:
            int_type[n] = 1

        # Same y-coordinates, x-z
        elif difference[1] == 0:
            int_type[n] = 2

        # Same z-coordinates, x-y
        elif difference[2] == 0:
            int_type[n] = 3

    return int_type


# Decides on the correct function for the cosine and sine, and integrates
def froude_kriloff_integrator(int_type, start, end, k, mu):
    """
    Initiates the Froude-Kriloff function and integrates it over a surface.
    :param int_type: int describing integration type. 1 (y-z), 2 (x-z), 3 (x-y).
    :param start: array describing start coordinates for a surface.
    :param end: array describing end coordinates for a surface.
    :param k: float describing the wave number in 1/m.
    :param mu: float describing wave angle in radians.
    :return: integrated_cosine, integrated_sine: the integrated values for a given start and end of a surface.
    """

    # Checks int type
    # X is constant
    if int_type == 1:
        # Sets constant X
        local_x = start[0]

        # Creates cosine and sine functions
        cosine = lambda y, z: -grav_constant * np.exp(k * z) * np.cos(k * local_x * np.cos(mu) + k * y * np.sin(mu))
        sine = lambda y, z: -grav_constant * np.exp(k * z) * np.sin(k * local_x * np.cos(mu) + k * y * np.sin(mu))

        # Integrates these functions over the start / end coordinates of the surface
        integrated_cosine = integrate.dblquad(cosine, start[2], end[2], start[1], end[1])[0]
        integrated_sine = integrate.dblquad(sine, start[2], end[2], start[1], end[1])[0]

        return integrated_cosine, integrated_sine

    # Y is constant
    elif int_type == 2:
        # Sets constant Y
        local_y = start[1]

        # Creates cosine and sine functions
        cosine = lambda x, z: -grav_constant * np.exp(k * z) * np.cos(k * x * np.cos(mu) + k * local_y * np.sin(mu))
        sine = lambda x, z: -grav_constant * np.exp(k * z) * np.sin(k * x * np.cos(mu) + k * local_y * np.sin(mu))

        # Integrates these functions over the start / end coordinates of the surface
        integrated_cosine = integrate.dblquad(cosine, start[2], end[2], start[0], end[0])[0]
        integrated_sine = integrate.dblquad(sine, start[2], end[2], start[0], end[0])[0]

        return integrated_cosine, integrated_sine

    # Z is constant
    elif int_type == 3:
        # Sets constant Z
        local_z = start[2]

        # Creates cosine and sine functions
        cosine = lambda x, y: -grav_constant * np.exp(k * local_z) * np.cos(k * x * np.cos(mu) + k * y * np.sin(mu))
        sine = lambda x, y: -grav_constant * np.exp(k * local_z) * np.sin(k * x * np.cos(mu) + k * y * np.sin(mu))

        # Integrates these functions over the start / end coordinates of the surface
        integrated_cosine = integrate.dblquad(cosine, start[1], end[1], start[0], end[0])[0]
        integrated_sine = integrate.dblquad(sine, start[1], end[1], start[0], end[0])[0]

        return integrated_cosine, integrated_sine

    # Error if int
    elif int_type == int:
        raise Exception("Outside of range. Please choose 1 for constant x, 2 for y, and 3 for z.")

    # Error if not int
    else:
        raise Exception("Please choose an integer. Please choose 1 for constant x, 2 for y, and 3 for z.")


# Froude Kriloff integrator
def object_integrator(object_number, wave_angle, k):
    """
    Initiates the froude-kriloff function and integrates it over the surfaces
    :param object_number: int describing object type. 1-4 pillar, 5-6 beam, 7-8 strut.
    :param wave_angle: int describing wave angle in degrees.
    :param k: float describing the wave number in 1/m.
    :return: integrated_cosine_vector, integrated_sine_vector: vectors containing integrated data for the object.
    """

    # Load surface parameter data for centroids, starts and ends
    centroids, starts, ends = surface_parameters(object_number)

    # Load additional normal data
    normals = object_normals()

    # Compute wave angle in radians
    mu = wave_angle / 180 * 2 * np.pi

    # Check integration type
    int_type = integration_type(starts, ends)

    # Create empty vectors to store data
    integrated_cosine_vector = np.zeros(3)
    integrated_sine_vector = np.zeros(3)

    # Cycle through all surfaces
    for n in range(6):
        # Collects the integrated cosine / sine term from the integrator
        integrated_cosine, integrated_sine = froude_kriloff_integrator(int_type[n], starts[n], ends[n], k, mu)

        # Multiplies the term by the normal vector, and adds it to the total cosine / sine vectors
        integrated_cosine_vector += integrated_cosine * normals[n]
        integrated_sine_vector += integrated_sine * normals[n]

    return integrated_cosine_vector, integrated_sine_vector


print(object_integrator(1, 180, 1))






















