from __future__ import annotations

from pathlib import Path

import numpy as np
import shapecomp as sc

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_import_and_load_tracked_torus() -> None:
    mesh_path = REPO_ROOT / "test" / "input" / "torus" / "torus.obj"
    assert mesh_path.is_file(), f"missing fixture: {mesh_path}"

    mesh = sc.load_mesh(str(mesh_path), False)
    vertices = mesh.get_vertices()
    faces = mesh.get_faces()
    assert vertices.shape[1] == 3
    assert faces.shape[1] == 3


def _check_sphere_locator(mesh_path: Path) -> None:
    mesh = sc.load_mesh(str(mesh_path), False)
    assert mesh.eulerCharacteristic == 2

    locator = sc.SphereLocator(mesh, Tutte_embedding_type="coTan")
    sphere = locator.get_sphere_vertices()
    vertices = mesh.get_vertices()
    assert sphere.shape == vertices.shape
    assert np.all(np.isfinite(sphere))
    norms = np.linalg.norm(sphere, axis=1)
    np.testing.assert_allclose(norms, np.ones_like(norms), atol=1e-6)

    pole = locator.get_pole_index()
    assert 0 <= pole < sphere.shape[0]

    face, bary = locator.find_bary_coords(sphere[0])
    assert face >= 0
    bary = np.asarray(bary, dtype=np.float64)
    assert bary.shape == (3,)
    np.testing.assert_allclose(np.sum(bary), 1.0, atol=1e-6)
    assert np.all(bary >= -1e-8)

    hits = locator.find_bary_coords_parallel(sphere[: min(4, sphere.shape[0])])
    assert len(hits) == min(4, sphere.shape[0])


def test_sphere_locator_icosahedron() -> None:
    mesh_path = REPO_ROOT / "test" / "input" / "sphere" / "icosahedron.obj"
    assert mesh_path.is_file(), f"missing fixture: {mesh_path}"
    _check_sphere_locator(mesh_path)


def test_sphere_locator_tetrahedron() -> None:
    mesh_path = REPO_ROOT / "test" / "input" / "sphere" / "tetrahedron.obj"
    assert mesh_path.is_file(), f"missing fixture: {mesh_path}"
    _check_sphere_locator(mesh_path)
