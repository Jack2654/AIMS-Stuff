import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.animation as animation
import numpy as np
import os


def gen_paths(height):
    ret_str = "BAND = "
    for i in range(21):
        ret_str += f'%s %s %s %s %s %s, ' % (round(-0.5 + i * 0.05, 3), -0.5, round(height, 2),
                                             round(-0.5 + i * 0.05, 3), 0.5, round(height, 2))
    ret_str = ret_str[:-2] + "\n"
    ret_str += "BAND_POINTS = 101 \nFREQUENCY_CONVERSION_FACTOR = 64.654178"
    return ret_str


def plot_bands(file, energy):
    x_data = []
    y_data = []
    with open(file, "r") as f:
        lines = f.readlines()
    total_kpoints = int(lines[0].split()[1])
    total_paths = int(lines[1].split()[1])
    points_per_band = int(total_kpoints / total_paths)
    freq_count = 0
    first = 0
    for line in lines:
        if "band" in line and first == 0:
            first = 1
        elif "band" in line and first > 0:
            first = 2
        elif first == 1 and "freq" in line:
            freq_count += 1
    next_point = 0
    for path in range(total_paths):
        frequencies = [[[] for y in range(freq_count)] for x in range(total_paths)]
        current_point = 0
        points = []
        for index, line in enumerate(lines):
            if "band" in line:
                current_point += 1
                if next_point + points_per_band >= current_point > next_point:
                    points.append([float(x.replace(",", "")) for x in lines[index - 2].split()[3:6]])
                    for i in range(freq_count):
                        frequencies[path][i].append(float(lines[index + 2 * (i + 1)].split()[1]))
        locations = [round(-0.5 + i * 0.01, 2) for i in range(101)]
        for level in frequencies[path]:
            crosses = crossings(level, energy)
            if len(crosses) > 0:
                x_data += [-0.5 + 0.05 * path for elem in crosses]
                y_data += crosses
                # plt.scatter([path for elem in crosses], crosses)
            # plt.plot(locations, level)
            # if len(crosses) > 0:
            #    plt.scatter(crosses, [20 for elem in crosses], color='k')
        # plt.hlines(20, -0.5, 0.5)
        # plt.show()
        next_point += points_per_band
    # plt.show()
    return x_data, y_data


def crossings(data, energy):
    locations = [round(-0.5 + i * 0.01, 2) for i in range(101)]
    below = True if data[0] < energy else False
    last = data[0]
    crosses = []
    for index, point in enumerate(data):
        cur = True if point < energy else False
        if below != cur:
            # print(f'crossing between %s and %s' % (last, point))
            # print(f'location between %s and %s' % (locations[index-1], locations[index]))
            locat = locations[index - 1] + (energy - last) * (locations[index] - locations[index - 1]) / (point - last)
            crosses.append(locat)
        last = point
        below = cur
    return crosses


def plot_isosurface(base, energy, save=False):
    x_data = []
    y_data = []
    z_data = []
    for i in range(21):
        height = round(-0.5 + 0.05 * i, 2)
        folder = base + "height_"
        if height < 0:
            folder += "n"
        folder += str(int(abs(height * 100))) + "/band.yaml"
        x_cur, y_cur = plot_bands(folder, energy)
        x_data += x_cur
        y_data += y_cur
        z_data += [height for elem in x_cur]
    all_data = np.array([[x_data[i], y_data[i], z_data[i]] for i in range(len(x_data))])
    # np.save(f'./%s.npy' % str(energy).rjust(2, '0'), all_data)

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x_data, y_data, z_data, color='k')
    ax.set_xlabel('q_x')
    ax.set_ylabel('q_y')
    ax.set_zlabel('q_z')
    ax.set_xlim([-0.5, 0.5])
    ax.set_ylim([-0.5, 0.5])
    ax.set_zlim([-0.5, 0.5])
    ax.set_title(f'%s meV isosurface' % energy)
    if save:
        plt.savefig(base + str(energy) + "_meV.png", dpi=1000, bbox_inches='tight', format="png")
    else:
        plt.show()


for i in range(1, 37):
    continue
    # plot_isosurface("../../Documents/24-25/Classes/ME490/tmp/phonopy_work/new_bands/", i, save=False)


# def update(index):
#     im.set_array(image_array[index])
#     return im,


# files = [f'../../Documents/24-25/Classes/ME490/tmp/phonopy_work/new_bands/%s.npy' % str(i).rjust(2,'0') for i in range(1,37)]
# for energy in range(1, 37):
#     continue
#     data = np.load(files[energy - 1])
#     with open('../../Documents/24-25/Classes/ME490/tmp/phonopy_work/new_bands/%s.xyz' % str(energy).rjust(2, '0'), 'w') as f:
#         for point in data:
#             f.write(f'{point[0]} {point[1]} {point[2]}\n')


# files = [f'../../Documents/24-25/Classes/ME490/tmp/phonopy_work/new_bands_opt_FC/%s_meV.png' % str(i) for i in range(1, 37)]
# files = [f'../../Documents/24-25/Classes/ME490/tmp/vis_unfold_0.6/Meshes/%s.png' % str(i).rjust(2, '0') for i in range(1, 37)]
# folder_path = "../../Documents/24-25/Classes/ME490/tmp/phonopy_work/new_bands/Meshes/"
# files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]

base = "../../Documents/24-25/Classes/ME490/tmp/new_point_cloud/python_meshes/"
files = [f'{base}{str(i).rjust(2, "0")}meV.png' for i in range(1, 37)]
files = [f'../../Documents/24-25/Classes/ME490/tmp/phonopy_work/qpoint_res/optimized/%smeV.png' % str(i).rjust(2, "0")
         for i in range(1, 37)]


def create_gif(image_paths, output_path, duration=500, loop=0):
    """
    Create a GIF from a list of image file paths.xq

    Args:
        image_paths (list of str): Paths to image files.
        output_path (str): File path to save the resulting GIF.
        duration (int): Duration between frames in milliseconds.
        loop (int): Number of loops (0 = infinite).
    """
    # Open images and convert to RGB to ensure consistency
    images = [Image.open(img_path) for img_path in image_paths]

    # Save the first image and append the rest
    images[0].save(
        output_path,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=loop
    )


base = "../../Documents/24-25/Fall/ME490/tmp/point_cloud/new_density/final_figures/scatter_"
files = [f'{base}{str(x).rjust(2, "0")}meV.png' for x in range(1, 37)]
create_gif(files, '../../Documents/24-25/Fall/ME490/tmp/point_cloud/new_density/final_figures/res.gif', duration=1000)

# image_array = []
# for my_file in files:
#     image = Image.open(my_file)
#     image_array.append(image)
# fig, ax = plt.subplots()
# Set the initial image
# im = ax.imshow(image_array[0], animated=True)
# Create the animation object
# animation_fig = animation.FuncAnimation(fig, update, frames=len(image_array),
#                                         interval=1000, blit=True)

# Show the animation
# plt.show()
# animation_fig.save(f'../../Documents/24-25/Classes/ME490/tmp/phonopy_work/qpoint_res/optimized/iso.gif')
# animation_fig.save(f'{base}iso.gif')
