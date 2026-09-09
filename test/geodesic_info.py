# add shapecomp in the build path here
import os
import time 
import sys
sys.path.append('../build')
import shapecomp as sc
import numpy as np
import torch
import argparse
import tqdm

# Fix NumPy randomness
np.random.seed(42)

# Fix PyTorch randomness
torch.manual_seed(42)
torch.cuda.manual_seed(42)
torch.cuda.manual_seed_all(42)  # For multi-GPU

# Enable deterministic behavior in PyTorch
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False


def barycentric_coordinates_torch(A, B, C, P):
    """
    Compute the barycentric coordinates of a point P with respect to a triangle ABC.

    Args:
        A (torch.Tensor): Vertex A of the triangle, shape (2,).
        B (torch.Tensor): Vertex B of the triangle, shape (2,).
        C (torch.Tensor): Vertex C of the triangle, shape (2,).
        P (torch.Tensor): Point P, shape (2,).

    Returns:
        tuple: Barycentric coordinates (lambdaA, lambdaB, lambdaC).
    """
    # Compute reference determinant for triangle ABC
    denom = (B[0] - A[0]) * (C[1] - A[1]) - (C[0] - A[0]) * (B[1] - A[1])

    if denom == 0:
        raise ValueError("Degenerate triangle: The points A, B, and C are collinear.")

    # Compute sub-triangle determinants
    lambdaA = ((B[0] - P[0]) * (C[1] - P[1]) - (C[0] - P[0]) * (B[1] - P[1])) / denom
    lambdaB = ((C[0] - P[0]) * (A[1] - P[1]) - (A[0] - P[0]) * (C[1] - P[1])) / denom
    lambdaC = ((A[0] - P[0]) * (B[1] - P[1]) - (B[0] - P[0]) * (A[1] - P[1])) / denom

    # Sanity check if their sum is close to 1
    sum_coords = lambdaA + lambdaB + lambdaC
    if torch.abs(sum_coords - 1.0) > 1e-4:
        raise ValueError(f"Barycentric coordinates do not sum to 1: diff = {torch.abs(sum_coords - 1.0).item()}")
    if lambdaA < 0 or lambdaB < 0 or lambdaC < 0:
        raise ValueError(f"Barycentric coordinates are negative: lambdaA = {lambdaA.item()}, lambdaB = {lambdaB.item()}, lambdaC = {lambdaC.item()}")

    return lambdaA, lambdaB, lambdaC


def compute_geodesic_lengths(
        domain_edges, 
        codomain_faces, codomain_vertices,
        cutted_codomain_faces,
        uv_images, codomain_planar_locator):
    """
    Compute the geodesic lengths of each image edge 
    Args:
        domain_edges: torch.tensor, shape=(nE, 3), dtype=torch.int32
        codomain_faces: torch.tensor, shape=(nF_2, 3), dtype=torch.int32
            The faces of the codomain mesh.
        codomain_vertices: torch.tensor, shape=(nV_2, 3), dtype=torch.float64
            The vertices of the codomain mesh.
        cutted_codomain_faces: torch.tensor, shape=(nF_2, 3), dtype=torch.int32
            The faces of the cutted codomain mesh. 
        uv_images: torch.tensor, shape=(nV_2_2, 2), dtype=torch.float64
            The image of UV positions of the domain mesh at the fundametal domain 
            of the image mesh. Note: nV_2_2 > nV_2 due to the cut. 
        codomain_planr_locator: PlanarLocator object
            The planar locator that maps a 2D vector to the face index 
            and barycentric coordinates of the image mesh.
    """

    # get codomain uv (Numpy 2D array)
    codomain_uvs = codomain_planar_locator.get_uv_positions()
    codomain_uvs = torch.tensor(codomain_uvs, dtype=torch.float64)

    # Use lists to accumulate results
    surface_points_pairs_collections = []
    # Iterate over faces in subdivided mesh1
    n_edges = domain_edges.shape[0]
    euc_and_geo_lengths = np.zeros((n_edges, 2))
    time0 = time.time()
    for e_id in tqdm.tqdm(range(n_edges)):
        edge = domain_edges[e_id]
        # Get the image of 3 UV positions (2D) of the face
        edge_uvs = uv_images[edge]  # Shape: (2, 2)

        # Use the planar locator to locate the image of each vertex 
        img_edge_vertices = torch.zeros((2, 3), dtype=torch.float64) 

        # store surface points for geodesic computation
        surface_pts = [] # each element is a tuple (face_idx, bary_coords)

        for i in range(2):  # For each vertex in the domain face
            
            # Convert the uv to [0, 1] by %1
            unitized_edge_uv_i = edge_uvs[i] % 1

            # Find the face index and barycentric coordinates containing the vertex
            f2_idx, bary_coords_val_debug = codomain_planar_locator.find_bary_coords(unitized_edge_uv_i.tolist())

            # Get the 3 UV positions of the face on the fundamental domain
            A, B, C = codomain_uvs[cutted_codomain_faces[f2_idx]]

            # Compute the barycentric coordinates of the vertex
            bary_coords = barycentric_coordinates_torch(A, B, C, unitized_edge_uv_i)
            bary_coords_torch = torch.stack(bary_coords)

            # assert the barycentric coordinates the same as the debug one
            if not np.allclose(bary_coords_torch.detach().numpy(), bary_coords_val_debug, atol=1e-4):
                print(f"Barycentric coordinates are not the same: {bary_coords_torch.detach().numpy()} != {bary_coords_val_debug}")
                raise ValueError("Barycentric coordinates are not the same")
            
            # Get the 3D position on the image mesh
            f2_vertices = codomain_vertices[codomain_faces[f2_idx]]

            # Compute the image vertex position
            img_edge_vertices[i] = torch.sum(bary_coords_torch.unsqueeze(1) * f2_vertices, dim=0)
            
            # store the surface point
            surface_pts.append((f2_idx, bary_coords_val_debug))
            

        # compute the Eucledian length of the edge
        euc_edge_len = torch.norm(img_edge_vertices[0] - img_edge_vertices[1])

        # store the surface points pairs
        surface_points_pairs_collections.append(surface_pts)

        euc_and_geo_lengths[e_id, 0] = euc_edge_len.item()
    time1 = time.time()
    print(f"Time takes: {time1 - time0:.4f} seconds")

    # compute the geodesic lengths in parallel
    print("Computing geodesic lengths in parallel...")
    t0 = time.time()
    geo_paths = sc.geodesic_compute_parallel(surface_points_pairs_collections, 
                                            codomain_vertices, codomain_faces)
    # store the longest geodesic path
    longest_geo_path = []
    longest_geo_len = 0
    # store the geodeisc path with the biggest ratio geo/euc
    biggest_geo_path_ratio = []
    biggest_ratio = 0
    for i in range(len(geo_paths)):
        pts_on_path = geo_paths[i]
        geo_len = 0
        for j in range(1, len(pts_on_path)):
            geo_len += np.linalg.norm(pts_on_path[j] - pts_on_path[j-1])
        if geo_len == 0:
            print(f"Warning: geodesic length is 0 for edge {i}")
            print(f"points on path:")
            print(pts_on_path)
        euc_and_geo_lengths[i, 1] = geo_len 

        # update the longest geodesic path
        if geo_len > longest_geo_len:
            longest_geo_len = geo_len
            longest_geo_path = pts_on_path
            index_of_longest_geo_path = i
        if geo_len / euc_and_geo_lengths[i, 0] > biggest_ratio:
            biggest_ratio = geo_len / euc_and_geo_lengths[i, 0]
            biggest_geo_path_ratio = pts_on_path
            index_of_biggest_geo_path_ratio = i

    t1 = time.time()
    print(f"Time takes: {t1 - t0:.4f} seconds")

    # save the longest geodesic path
    file_name = f"longest_geo_path_edge_idx_{index_of_longest_geo_path}.npy"
    print(f"Saving the longest geodesic path to {file_name}")
    np.save(file_name, longest_geo_path)

    # save the longest geodesic path with the biggest ratio
    file_name = f"ratio_geo_path_edge_idx_{index_of_biggest_geo_path_ratio}.npy"
    print(f"Saving the geodesic path with the biggest ratio to {file_name}")
    np.save(file_name, biggest_geo_path_ratio)

    return euc_and_geo_lengths


if __name__ == "__main__":
    # add args in command line
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_iter", type=int, default=1000, help="Number of iterations") 
    parser.add_argument("--lr", type=float, default=1e-4, help="Learning rate")
    args = parser.parse_args()
    n_iter = args.n_iter
    lr = args.lr

    # print(f"Number of iterations: {n_iter}, lr = {lr}")
    
    # set the path of the two meshes
    simple_torus_path = 'input/torus/simple_torus.obj'
    subdivied_simple_torus_path = 'input/torus/simple_torus_subdivision_igl.obj'
    torus_path = 'input/torus/torus.obj'
    subdivied_torus_path = 'input/torus/torus_subdivision_igl.obj'

    # start compute the energy
    obj_path1 = simple_torus_path
    obj_path2 = subdivied_simple_torus_path

    # load two meshes and compute the fundamental domains
    mesh1 = sc.load_mesh(obj_path1)
    planar1 = sc.PlanarLocator(mesh1)
    subdivided_mesh_faces_1 = planar1.get_subdivided_mesh_faces()
    subdivided_mesh_vertices_1 = planar1.get_subdivided_mesh_vertices()
    subdivided_mesh_edges_1 = planar1.get_subdivided_mesh_edges()

    mesh2 = sc.load_mesh(obj_path2)
    planar2 = sc.PlanarLocator(mesh2)
    subdivided_mesh_faces_2 = planar2.get_subdivided_mesh_faces()
    subdivided_mesh_vertices_2 = planar2.get_subdivided_mesh_vertices()
    subdivided_mesh_edges_2 = planar2.get_subdivided_mesh_edges()
    cutted_mesh_vertices_2 = planar2.get_cutted_mesh_vertices()
    cutted_mesh_faces_2 = planar2.get_cutted_mesh_faces()
    cutted_mesh_edges_2 = planar2.get_cutted_mesh_edges()
    print("Mesh 2 size info:")
    print("subdivided mesh:")
    print("V = ", subdivided_mesh_vertices_2.shape[0], "F = ", subdivided_mesh_faces_2.shape[0], "E = ", subdivided_mesh_edges_2.shape[0])
    print("cutted mesh:")
    print("V = ", cutted_mesh_vertices_2.shape[0], "F = ", cutted_mesh_faces_2.shape[0], "E = ", cutted_mesh_edges_2.shape[0])


    # Get the UV positions of the domain mesh
    img_f_uvs = torch.tensor(planar1.get_uv_positions(), dtype=torch.float64, requires_grad=True)

    # compute the PARs
    euc_geo_lens = compute_geodesic_lengths(subdivided_mesh_edges_1,
                                            subdivided_mesh_faces_2, subdivided_mesh_vertices_2,
                                            cutted_mesh_faces_2,
                                            img_f_uvs, planar2)
    

    save_folder_path = "big_torus_output/"
    # if the folder does not exist, create it
    if not os.path.exists(save_folder_path):
        os.makedirs(save_folder_path)
    # save the results
    np.save(save_folder_path + "euc_geo_lens_simple_torus_vs_subdivided_simple_torus.npy", euc_geo_lens)
