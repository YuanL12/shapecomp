from __future__ import annotations

from pathlib import Path

import shapecomp as sc


def test_import_and_load_tracked_torus() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    mesh_path = repo_root / "test" / "input" / "torus" / "torus.obj"
    assert mesh_path.is_file(), f"missing fixture: {mesh_path}"

    mesh = sc.load_mesh(str(mesh_path), False)
    vertices = mesh.get_vertices()
    faces = mesh.get_faces()
    assert vertices.shape[1] == 3
    assert faces.shape[1] == 3
