# Import statements
import numpy as np

# Global object constants
# Total object
total_length = 60
total_width = 60
total_height = 18

# Beams
beam_length = 60
beam_width = 12.65
beam_height = 6.25

# Pillars
pillar_width = 10.75
pillar_height = total_height - beam_height

# Struts
strut_length = total_width - 2 * beam_width
strut_width = 2

# Object lists
pillar_list = [1, 2, 3, 4]
beam_list = [5, 6]
strut_list = [7, 8]


def object_parameters(object_number):
    """
    Loads the different parameters of the surfaces on the object.
    :param object_number: int describing object type. 1-4 pillar, 5-6 beam, 7-8 strut.
    :return: area and normal parameters of the object.
    """

    # Object selection:
    # Pillars
    if object_number in pillar_list:
        print(object_number)
        pillar_top_area = pillar_width ** 2
        pillar_side_area = pillar_width * pillar_height

        pillar_areas = np.array([[pillar_top_area, pillar_side_area, pillar_side_area, pillar_side_area,
                                  pillar_side_area]])
        pillar_normals = np.array([[0, 0, 1],
                                   [1, 0, 0],
                                   [0, 1, 0],
                                   [-1, 0, 0],
                                   [0, -1, 0]])

        return pillar_areas, pillar_normals

    # Beams
    elif object_number in beam_list:
        beam_top_area = beam_length * beam_width - (2 * pillar_width ** 2)
        beam_end_area = beam_width * beam_height
        beam_side_area = beam_length * beam_height
        beam_bottom_area = beam_length * beam_width

        beam_areas = np.array([[beam_top_area, beam_end_area, beam_side_area, beam_end_area, beam_side_area,
                                beam_bottom_area]])
        beam_normals = np.array([[0, 0, 1],
                                 [1, 0, 0],
                                 [0, 1, 0],
                                 [-1, 0, 0],
                                 [0, -1, 0],
                                 [0, 0, -1]])

        return beam_areas, beam_normals

    # Struts
    elif object_number in strut_list:
        strut_area = strut_width * strut_length

        strut_areas = np.full([4, 1], strut_area)
        strut_normals = np.array([[0, 0, 1],
                                  [1, 0, 0],
                                  [0, 0, -1],
                                  [-1, 0, 0]])

        return strut_areas, strut_normals

    # Error if int
    elif type(object_number) == int:
        raise Exception("Outside of range. Please choose 1-4 for pillar, 5-6 for beam, or 7-8 for strut.")

    # Error if not int
    else:
        raise Exception("Please choose an integer. Choose 1-4 for pillar, 5-6 for beam, or 7-8 for strut.")


def surface_centroids(object_number):
    """
    Loads the different centroids of the surfaces on the object.
    :param object_number: int describing object type. 1-4 pillar, 5-6 beam, 7-8 strut.
    :return: area and normal parameters of the object.
    """

    # Object selection:
    # Pillars
    if object_number in pillar_list:
        # Base distance per pillar, to the base point
        short = 30 - pillar_width  # The reduced distance from COG if it includes a pillar width [m]
        pillar_base_distances = np.array([[30, 30, 1],
                                          [-short, 30, 1],
                                          [-short, -short, 1],
                                          [30, -short, 1]])

        # Modifiers for every pillar
        half = 0.5 * pillar_width   # The half of the width of every pillar [m]
        half_height = 0.5 * pillar_height
        pillar_modifiers = np.array([[-half, -half, 0],
                                     [0, -half, -half_height],
                                     [-half, 0, -half_height],
                                     [-pillar_width, -half, -half_height],
                                     [-half, -pillar_width, -half_height]])

        # Empty matrix to store data
        pillar_centroids = np.zeros(np.shape(pillar_modifiers))

        # Apply base distance on modifiers
        for n in range(np.shape(pillar_modifiers)[0]):
            pillar_centroids[n] = pillar_modifiers[n] + pillar_base_distances[object_number - 1]

        return pillar_centroids

    # Beams
    elif object_number in beam_list:
        # Base distance per beam, to the base point
        half_length = 0.5 * beam_length
        beam_top_depth = 1 - 18 + beam_height
        beam_base_distances = np.array([[half_length, half_length, beam_top_depth],
                                        [half_length, -half_length + beam_width, beam_top_depth]])

        # Beam top y centroid calculation
        top_small_area = (beam_width - pillar_width) * beam_length
        top_small_centroid = -beam_width + ((beam_width - pillar_width) * 0.5)
        top_large_area = (beam_length - 2 * pillar_width) * pillar_width
        top_large_centroid = -pillar_width / 2

        top_y_centroid = (top_small_centroid * top_small_area + top_large_centroid + top_large_area) / \
                         (top_small_area + top_large_area)

        # y top modifier selection based on object number
        y_top_modifier = 0

        if object_number == 5:
            y_top_modifier = top_y_centroid

        elif object_number == 6:
            y_top_modifier = beam_width - top_y_centroid

        # Modifiers for every beam
        half_width = 0.5 * beam_width
        half_height = 0.5 * beam_height
        beam_modifiers = np.array([[-half_length, y_top_modifier, beam_top_depth],
                                   [0, -half_width, -half_height],
                                   [-half_length, 0, -half_height],
                                   [-beam_length, -half_width, -half_height],
                                   [-half_length, -beam_width, -half_height],
                                   [-half_length, -half_width, -beam_height]])

        # Empty matrix to store data
        beam_centroids = np.zeros(np.shape(beam_modifiers))

        # Apply base distance on modifiers
        for n in range(np.shape(beam_modifiers)[0]):
            beam_centroids[n] = beam_modifiers[n] + beam_base_distances[object_number - 5]

        return beam_centroids

    # Struts
    elif object_number in strut_list:
        # Base distance per strut, to the base point
        inter_strut_distance = 0.5 * (total_length - pillar_width)
        vert_strut_distance = 1 - total_height + 7 + 1
        strut_base_distances = np.array([[inter_strut_distance, 0, vert_strut_distance],
                                         [-inter_strut_distance, 0, vert_strut_distance]])

        # Modifiers for every strut
        half = 0.5 * strut_width   # The half of the width of every strut [m]
        strut_modifiers = np.array([[0, 0, half],
                                     [half, 0, 0],
                                     [0, 0, -half],
                                     [-half, 0, 0]])

        # Empty matrix to store data
        strut_centroids = np.zeros(np.shape(strut_modifiers))

        # Apply base distance on modifiers
        for n in range(np.shape(strut_modifiers)[0]):
            strut_centroids[n] = strut_modifiers[n] + strut_base_distances[object_number - 7]

        return strut_centroids

    # Error if int
    elif type(object_number) == int:
        raise Exception("Outside of range. Please choose 1-4 for pillar, 5-6 for beam, or 7-8 for strut.")

    # Error if not int
    else:
        raise Exception("Please choose an integer. Choose 1-4 for pillar, 5-6 for beam, or 7-8 for strut.")








