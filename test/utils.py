import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import cv2


def create_uv_pos_video(uv_pos_res, save_path):
    # Assuming read_uv_obj and plot_mesh are already defined
    _, faces, vertex_labels = read_uv_obj("input/torus/simple_torus_square_uv.obj")

    # Parameters for the video
    output_filename = save_path
    fps = 15  # Frames per second
    width, height = 1200, 1200  # Increased resolution for larger canvas

    # Initialize the video writer
    fourcc = cv2.VideoWriter_fourcc(*"mp4v")  # Codec for MP4
    video_writer = cv2.VideoWriter(output_filename, fourcc, fps, (width, height))

    # Create a reusable figure and axis
    fig, ax = plt.subplots(figsize=(width / 100, height / 100), dpi=100)

    # find the max and min on x and y axes of the uv_pos_res respectively
    max_x = max([np.max(uv_pos_res[i][:, 0]) for i in range(len(uv_pos_res))])
    min_x = min([np.min(uv_pos_res[i][:, 0]) for i in range(len(uv_pos_res))])
    x_diff = max_x - min_x
    max_y = max([np.max(uv_pos_res[i][:, 1]) for i in range(len(uv_pos_res))])
    min_y = min([np.min(uv_pos_res[i][:, 1]) for i in range(len(uv_pos_res))])
    y_diff = max_y - min_y

    # Generate frames
    for i, uv_positions in enumerate(uv_pos_res):
        # Clear the axis to prepare for the new frame
        ax.clear()

        # Plot the mesh
        plot_2D_mesh(
            uv_positions,
            faces,
            vertex_labels=None,
            fig=fig,
            ax=ax,
            show_vertex_label=False,
        )
        ax.set_aspect("equal", "box")
        # Add the unit square with dashed black boundary
        unit_square = np.array(
            [[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]]
        )  # Coordinates of the square (closed loop)
        ax.plot(
            unit_square[:, 0],
            unit_square[:, 1],
            linestyle="--",
            color="black",
            zorder=2,
        )
        ax.set_title(f"Frame: {i}", fontsize=16, pad=20)
        # set the x and y limits
        ax.set_xlim(min_x - x_diff * 0.1, max_x + x_diff * 0.1)
        ax.set_ylim(min_y - y_diff * 0.1, max_y + y_diff * 0.1)

        # plt.tight_layout()

        # Save the plot to a NumPy array
        fig.canvas.draw()
        frame = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
        frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (3,))

        # Write the frame to the video
        video_writer.write(cv2.cvtColor(frame, cv2.COLOR_RGB2BGR))

    # Release the video writer and close the figure
    video_writer.release()
    plt.close(fig)

    print(f"Video saved as {output_filename}")


def plotly_draw_3D_mesh(vertices, faces, fig=None, show_now=False):
    """
    Draw a 3D mesh using Plotly.

    Parameters
    ----------
    vertices : numpy.ndarray of shape (n, 3)
        The vertices of the mesh.
    faces : numpy.ndarray of shape (m, 3)
        The faces of the mesh.
    fig : plotly.graph_objects.Figure, optional
        The figure to add the mesh to. If None, a new figure is created.
    show_now : bool, optional
        If True, the figure is shown immediately.
    """

    mesh_3d_plot = go.Mesh3d(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        # colorbar=dict(title=dict(text='z')),
        # colorscale=[[0, 'gold'],
        #             [0.5, 'mediumturquoise'],
        #             [1, 'magenta']],
        # intensity=np.linspace(0, 1, 8, endpoint=True),
        i=faces[:, 0],
        j=faces[:, 1],
        k=faces[:, 2],
        name="Mesh",
        legendgroup="Mesh",
        showlegend=True,
        # showscale=True
    )

    if not fig:
        # Define the mesh (cube)
        fig = go.Figure(data=[mesh_3d_plot])
    else:
        fig.add_trace(mesh_3d_plot)

    # Create the edges
    # Define the edges as pairs of vertex indices
    edges = set()
    for f in faces:
        for i in range(3):
            edge = [f[i], f[(i + 1) % 3]]
            edges.add(tuple(sorted(edge)))
    edge_x = []
    edge_y = []
    edge_z = []

    for edge in edges:
        for point in edge:
            edge_x.append(vertices[point][0])
            edge_y.append(vertices[point][1])
            edge_z.append(vertices[point][2])
        # Add None to prevent incorrect connections
        edge_x.append(None)
        edge_y.append(None)
        edge_z.append(None)

    edges_trace = go.Scatter3d(
        x=edge_x,
        y=edge_y,
        z=edge_z,
        mode="lines",
        line=dict(color="black", width=3),
        name="Edges",
        legendgroup="Edges",
        showlegend=True,
    )
    fig.add_trace(edges_trace)

    # Find the min and max values for each axis
    x_range = [vertices[:, 0].min(), vertices[:, 0].max()]
    y_range = [vertices[:, 1].min(), vertices[:, 1].max()]
    z_range = [vertices[:, 2].min(), vertices[:, 2].max()]

    # Compute the max range to enforce equal scaling
    max_range = max(
        x_range[1] - x_range[0], y_range[1] - y_range[0], z_range[1] - z_range[0]
    )
    center_x = sum(x_range) / 2
    center_y = sum(y_range) / 2
    center_z = sum(z_range) / 2

    equal_x_range = [center_x - max_range / 2, center_x + max_range / 2]
    equal_y_range = [center_y - max_range / 2, center_y + max_range / 2]
    equal_z_range = [center_z - max_range / 2, center_z + max_range / 2]

    # Set equal axis scaling
    fig.update_layout(
        scene=dict(
            xaxis=dict(range=equal_x_range),
            yaxis=dict(range=equal_y_range),
            zaxis=dict(range=equal_z_range),
            aspectmode="cube",
        )
    )
    if show_now:
        # Show the figure
        fig.show()
    return fig


def add_path_on_plotly_mesh(
    fig, path_points, figshow_now=False, show_markers=True, path_color="red"
):
    # Add the path as a line
    fig.add_trace(
        go.Scatter3d(
            x=path_points[:, 0],
            y=path_points[:, 1],
            z=path_points[:, 2],
            mode="lines+markers",
            line=dict(color=path_color, width=5),
            marker=dict(size=1, color="black", symbol="circle"),
            name="Path",
        )
    )
    if figshow_now:
        fig.show()
    return fig


def bary_coordinates_to_euc_point(face_idx, bary_coords, faces, vertex_postions):
    """
    Compute the Eucledian point from the barycentric coordinates.
    Parameters
    ----------
    face_idx : int
        The index of the face.
    bary_coords : numpy.ndarray of shape (3,)
        The barycentric coordinates.
    faces : numpy.ndarray of shape (m, 3)
        The faces of the mesh.
    vertex_postions : numpy.ndarray of shape (n, 3) or (n, 2)
        The vertex positions of the mesh either 2D or 3D.
    """
    if vertex_postions.shape[1] == 2:
        p_position = np.zeros(2)
    else:
        p_position = np.zeros(3)
    for i in range(3):
        p_position += bary_coords[i] * vertex_postions[faces[face_idx][i]]
    return p_position


# Function to plot 2D mesh with vertex labels
def plot_2D_mesh(
    vertices,
    faces,
    vertex_labels=None,
    fig=None,
    ax=None,
    fil_color="blue",
    show_vertex_label=True,
):
    """
    Plot a 2D mesh with vertex labels.
    Parameters
    ----------
    vertices : numpy.ndarray of shape (n, 2) or (n, 3)
        The vertices of the mesh. The third dimension is optional.
    faces : numpy.ndarray of shape (m, 3)
        The faces of the mesh.
    vertex_labels : dict, optional
        The vertex labels.
    fig : matplotlib.figure.Figure, optional
        The figure to plot the mesh on.
    """

    if fig is None or ax is None:
        fig, ax = plt.subplots()

    # plot traingles
    for face in faces:
        polygon = [vertices[idx] for idx in face]
        polygon = np.array(polygon)
        ax.fill(
            polygon[:, 0],
            polygon[:, 1],
            edgecolor="blue",
            fill=True,
            color=fil_color,
            alpha=0.2,
            zorder=1,
        )  # Cyan fill with opacity

    # Annotate each vertex with its label
    if show_vertex_label:
        for idx in range(len(vertices)):
            x, y = vertices[idx]
            ax.text(
                x,
                y,
                vertex_labels[str(idx)] + "/" + str(idx),
                fontsize=10,
                ha="center",
                color="black",
                zorder=3,
            )

    # plt.show()
    return fig, ax


def read_uv_obj(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    vertices = []
    faces = []
    vertex_label = {}
    for line in lines:
        if line[0] == "v":
            # example line is: vt 0.333333 -0.333333
            # so we split the line by space and take the
            # last two elements and convert them to float
            vertices.append([float(i) for i in line.split()[1:]])

        elif line[0] == "f":
            # example line is: f 3/-1 13/0 1/-1
            # we only need the index in front of the /
            # so we split the line by space and then by /
            # and then take the first element of the split
            # and convert it to int
            for i in line.split()[1:]:
                after_cut_index = i.split("/")[0]
                ref_index = i.split("/")[1]
                if after_cut_index not in vertex_label:
                    if ref_index == "-1":
                        vertex_label[after_cut_index] = after_cut_index
                    else:
                        vertex_label[after_cut_index] = ref_index

            faces.append([int(i.split("/")[0]) for i in line.split()[1:]])

    return np.array(vertices), np.array(faces), vertex_label


# read obj file, note that the faces are 1-indexed
def read_obj(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    vertices = []
    faces = []
    for line in lines:
        if line[0] == "v":
            # example line is: v 0.1 -0.2 0.3
            # so we split the line by space and take the
            # last two elements and convert them to float
            vertices.append([float(i) for i in line.split()[1:]])

        elif line[0] == "f":
            # example line is: f 3 2 1
            faces.append([int(i) for i in line.split()[1:]])

    return np.array(vertices), np.array(faces)


def save_obj(vertices, faces, filename):
    # save the vertices and faces to a .obj file
    with open(filename, "w") as f:
        for v in vertices:  # use :.6f to format the vertices
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for f in faces:
            f.write(f"f {f[0]} {f[1]} {f[2]}\n")
