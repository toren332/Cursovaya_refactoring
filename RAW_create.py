from stl import mesh
import math
import copy
import numpy as np
from box_triangle import intersects_box
from mpl_toolkits import mplot3d
from matplotlib import pyplot
from array import array
import re
import sys
import matplotlib.pyplot as plt
import pickle


def crete_raw(filename, voxel_size, layers, fill_model=True, progress=False, restore_voxels=False):
    def get_model_border():
        """Возвращает максимумы и минимумы по осям, и массив треугольников"""
        triangles_array = []
        maximum_x = -float('inf')
        maximum_y = -float('inf')
        maximum_z = -float('inf')
        minimum_x = float('inf')
        minimum_y = float('inf')
        minimum_z = float('inf')
        for triangle in stl_mesh:
            _ = [[triangle[0], triangle[2], triangle[1]],
                 [triangle[3], triangle[5], triangle[4]],
                 [triangle[6], triangle[8], triangle[7]]]
            for dot in _:
                if dot[0] > maximum_x:
                    maximum_x = dot[0]
                elif dot[0] < minimum_x:
                    minimum_x = dot[0]
                if dot[2] > maximum_z:
                    maximum_z = dot[2]
                elif dot[2] < minimum_z:
                    minimum_z = dot[2]
                if dot[1] > maximum_y:
                    maximum_y = dot[1]
                elif dot[1] < minimum_y:
                    minimum_y = dot[1]
            triangles_array.append(_)
        response = {
            'x_max': maximum_x,
            'y_max': maximum_y,
            'z_max': maximum_z,
            'x_min': minimum_x,
            'y_min': minimum_y,
            'z_min': minimum_z,
            'triangles_array': triangles_array,
        }
        return response

    def real_cube_mesh():
        """Создаёт кубическую сетку для будущей воксельной модели возвращает центры созданных пустых вокселей"""

        def round_segment(x, y):
            """Расширяет отрезок для того, чтобы воксели уместились"""
            length = abs(x - y)
            y = x - x % voxel_size + length - length % voxel_size + 3 * voxel_size
            x = x - x % voxel_size - 2 * voxel_size

            return x, y

        centers = []
        min_x = border_nodes['x_min']
        min_y = border_nodes['y_min']
        min_z = border_nodes['z_min']
        max_x = border_nodes['x_max']
        max_y = border_nodes['y_max']
        max_z = border_nodes['z_max']

        rounded_min_x, rounded_max_x = round_segment(min_x, max_x)
        rounded_min_y, rounded_max_y = round_segment(min_y, max_y)
        rounded_min_z, rounded_max_z = round_segment(min_z, max_z)

        x_centers = []
        y_centers = []
        z_centers = []
        voxel_x_centers = []
        voxel_y_centers = []
        voxel_z_centers = []
        voxel_centers = []

        for i in range(int(round(rounded_min_x / voxel_size)), int(round(rounded_max_x / voxel_size))):
            x_centers.append(i * voxel_size + voxel_size / 2)
            voxel_x_centers.append(round(i - rounded_min_x / voxel_size))

        for i in range(int(round(rounded_min_y / voxel_size)), int(round(rounded_max_y / voxel_size))):
            y_centers.append(i * voxel_size + voxel_size / 2)
            voxel_y_centers.append(round(i - rounded_min_y / voxel_size))

        for i in range(int(round(rounded_min_z / voxel_size)), int(round(rounded_max_z / voxel_size))):
            z_centers.append(i * voxel_size + voxel_size / 2)
            voxel_z_centers.append(round(i - rounded_min_z / voxel_size))

        for x_center in x_centers:
            for y_center in y_centers:
                for z_center in z_centers:
                    centers.append([x_center, y_center, z_center])

        for x_center in voxel_x_centers:
            for y_center in voxel_y_centers:
                for z_center in voxel_z_centers:
                    voxel_centers.append([x_center, y_center, z_center])

        rounded_borders = {
            "rounded_min_x": rounded_min_x,
            "rounded_min_y": rounded_min_y,
            "rounded_min_z": rounded_min_z,
            "rounded_max_x": rounded_max_x,
            "rounded_max_y": rounded_max_y,
            "rounded_max_z": rounded_max_z,
        }
        return centers, rounded_borders, voxel_centers

    def insertions():
        """Выберает воксели, которые появляются на пересечении триангулированной поверхности и сетки"""
        not_empty_voxels = []
        old_cur = -1
        for center in all_centers:
            for triangle in triangle_array:

                if intersects_box(np.array(triangle), np.array(center), voxel_size):
                    not_empty_voxels.append(voxel_centers[all_centers.index(center)])
                    break
            if progress:

                cur = (((all_centers.index(center)/all_centers.__len__())*100)//0.1)/10
                if cur != old_cur:
                    print(str(cur)+'%')
                old_cur = cur

        return not_empty_voxels

    def find_fill_coord():
        for xx in range(x_length):
            for yy in range(y_length):
                for zz in range(z_length):
                    if voxels[xx][yy][zz] == 1:
                        beam1 = ''
                        beam2 = ''
                        for zz_1 in range(zz, z_length):
                            # верхний луч
                            beam1 += str(int(voxels[xx][yy][zz_1]))
                        for zz_2 in range(0, zz + 1):
                            # нижний луч
                            beam2 += str(int(voxels[xx][yy][zz_2]))

                        nbeam1 = beam1[::-1]
                        nbeam2 = beam2[::-1]
                        match1 = re.search(r'2+1+0+[0-1-2]', nbeam1)
                        # match1 = re.search(r'2+1+0+[0-1-2]*1+2+', nbeam1)
                        match2 = re.search(r'2+1+0+[0-1-2]', nbeam2)

                        if match1:
                            first_1 = 0
                            for nchar in range(nbeam1.__len__()):
                                char = nbeam1[nchar]
                                if char == '1' and first_1 == 0:
                                    first_1 = 1
                                if char == '0' and first_1 == 1:
                                    zzz = zz + nbeam1.__len__() - nchar - 1
                                    fill_coords = [xx, yy, zzz]
                                    return fill_coords

                        elif match2:

                            first_1 = 0
                            for nchar in range(nbeam2.__len__()):
                                char = nbeam2[nchar]
                                if char == '1' and first_1 == 0:
                                    first_1 = 1
                                if char == '0' and first_1 == 1:
                                    zzz = nbeam2.__len__() - 1 - nchar
                                    fill_coords = [xx, yy, zzz]
                                    return fill_coords
                        else:
                            pass

    def fill(coord, old_color, new_color):

        stack = [coord]
        while stack:
            point = stack.pop()
            if point[0] < x_length and point[1] < y_length and point[2] < z_length:
                if voxels[point[0]][point[1]][point[2]] == old_color:
                    voxels[point[0]][point[1]][point[2]] = new_color

                    ncoord1 = copy.deepcopy(point)
                    ncoord1[0] = ncoord1[0] - 1
                    if ncoord1[0] >= 0:
                        if voxels[ncoord1[0]][ncoord1[1]][ncoord1[2]] == old_color:
                            stack.append(ncoord1)

                    ncoord2 = copy.deepcopy(point)
                    ncoord2[0] = ncoord2[0] + 1
                    if ncoord2[0] < x_length:
                        if voxels[ncoord2[0]][ncoord2[1]][ncoord2[2]] == old_color:
                            stack.append(ncoord2)

                    ncoord3 = copy.deepcopy(point)
                    ncoord3[1] = ncoord3[1] - 1
                    if ncoord3[1] >= 0:
                        if voxels[ncoord3[0]][ncoord3[1]][ncoord3[2]] == old_color:
                            stack.append(ncoord3)

                    ncoord4 = copy.deepcopy(point)
                    ncoord4[1] = ncoord4[1] + 1
                    if ncoord4[1] < y_length:
                        if voxels[ncoord4[0]][ncoord4[1]][ncoord4[2]] == old_color:
                            stack.append(ncoord4)

                    ncoord5 = copy.deepcopy(point)
                    ncoord5[2] = ncoord5[2] - 1
                    if ncoord5[2] >= 0:
                        if voxels[ncoord5[0]][ncoord5[1]][ncoord5[2]] == old_color:
                            stack.append(ncoord5)

                    ncoord6 = copy.deepcopy(point)
                    ncoord6[2] = ncoord6[2] + 1
                    if ncoord6[2] < z_length:
                        if voxels[ncoord6[0]][ncoord6[1]][ncoord6[2]] == old_color:
                            stack.append(ncoord6)


        pass

    stl_mesh = mesh.Mesh.from_file('models/' + filename + '.stl')
    border_nodes = get_model_border()
    triangle_array = border_nodes.pop('triangles_array')
    all_centers, rounded_borders, voxel_centers = real_cube_mesh()

    good_voxels = insertions()

    rounded_min_x = rounded_borders['rounded_min_x']
    rounded_min_y = rounded_borders['rounded_min_y']
    rounded_min_z = rounded_borders['rounded_min_z']
    rounded_max_x = rounded_borders['rounded_max_x']
    rounded_max_y = rounded_borders['rounded_max_y']
    rounded_max_z = rounded_borders['rounded_max_z']

    x_length = int((rounded_max_x - rounded_min_x) / voxel_size)
    y_length = int((rounded_max_y - rounded_min_y) / voxel_size)
    z_length = int((rounded_max_z - rounded_min_z) / voxel_size)

    voxels = np.zeros((x_length, y_length, z_length))

    for vox in good_voxels:
        x, y, z = vox
        voxels[int(x)][int(y)][int(z)] = 1

    if fill_model:
        fill([0, 0, 0], 0, 2)

        fill_dot = find_fill_coord()
        print(fill_dot)
        if fill_dot:
            fill(fill_dot, 0, 1)
        fill([0, 0, 0], 2, 0)

    raw_list = []
    for z in range(z_length):
        for y in range(y_length):
            for x in range(x_length):
                if layers != -1:
                    if z > layers - 1:
                        voxels[x][y][z] = 0
                if voxels[x][y][z] == 1:
                    raw_list.append(1)
                else:
                    raw_list.append(0)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.voxels(voxels, edgecolor='k')

    output_file = open('models_output/' + filename + '_' + str(x_length) + 'x' + str(y_length) + 'x' + str(z_length) + ".raw", "wb")
    new_array = array('b', raw_list)
    new_array.tofile(output_file)
    output_file.close()

    plt.show()
    return filename + '_' + str(x_length) + 'x' + str(y_length) + 'x' + str(z_length)
