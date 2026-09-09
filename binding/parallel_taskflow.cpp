#include <atomic>
#include <algorithm>
#include <taskflow/algorithm/for_each.hpp>
#include <taskflow/algorithm/sort.hpp>
#include <taskflow/core/flow_builder.hpp>
#include <taskflow/taskflow.hpp>  // Include Taskflow

#include "core.hpp"
#include "extern/geometry_central/geodesic.hpp"

namespace py = pybind11;

namespace geomp {
    typedef std::pair<int64_t, std::vector<double>> SurfacePoint;
}

std::vector<GC::DenseMatrix<double>> geodesic_compute_parallel(
    std::vector<std::pair<geomp::SurfacePoint, geomp::SurfacePoint>> surface_points_pairs,
    GC::DenseMatrix<double> verts, GC::DenseMatrix<int64_t> faces) {
    size_t iteration_count = surface_points_pairs.size();
    std::vector<GC::DenseMatrix<double>> results(iteration_count);

    tf::Executor executor;
    tf::Taskflow taskflow;

    size_t worker_count = std::min(iteration_count, executor.num_workers());
    for (size_t w = 0; w < worker_count; ++w) {
        taskflow.emplace([&, w, worker_count]() {
            auto GM = std::make_unique<geomp::GeodesicsManager>(verts, faces);
            for (size_t i = w; i < iteration_count; i += worker_count) {
                auto surface_pts_pair = surface_points_pairs[i];
                geomp::SurfacePoint start_surface_point = surface_pts_pair.first;
                geomp::SurfacePoint end_surface_point = surface_pts_pair.second;
                results[i] = GM->find_exact_geodesic_path(start_surface_point.first, end_surface_point.first,
                                                    start_surface_point.second, end_surface_point.second);
            }
        });
    }


    executor.run(taskflow).get();  // Execute all tasks
    return results;
}

std::vector<GC::DenseMatrix<double>> flip_geodesic_compute_parallel(
    std::vector<std::pair<geomp::SurfacePoint, geomp::SurfacePoint>> surface_points_pairs,
    GC::DenseMatrix<double> verts, GC::DenseMatrix<int64_t> faces) {
    size_t iteration_count = surface_points_pairs.size();
    std::vector<GC::DenseMatrix<double>> results(iteration_count);

    tf::Executor executor;
    tf::Taskflow taskflow;

    // Add atomic counter for progress tracking
    std::atomic<size_t> counter{0};

    taskflow.for_each_index((size_t)0, iteration_count, (size_t)1, [&](size_t i) {
        auto surface_pts_pair = surface_points_pairs[i];
        // get the start and end surface point
        geomp::SurfacePoint start_surface_point = surface_pts_pair.first;
        geomp::SurfacePoint end_surface_point = surface_pts_pair.second;
        try {
            auto GM = std::make_unique<geomp::GeodesicsManager>(verts, faces);
            results[i] = GM->find_flip_geodesic_path(start_surface_point.first, end_surface_point.first, start_surface_point.second, end_surface_point.second);
        } catch (const std::exception &e) {
            std::cerr << "Error computing flip geodesic path for pair " << i << ": " << e.what() << std::endl;
            // Initialize to empty matrix to avoid uninitialized data
            results[i] = GC::DenseMatrix<double>(0, 3);
            throw; // Re-throw to propagate the error
        }
        // Update progress
        size_t completed = ++counter;
    });

    executor.run(taskflow).get();  // Execute all tasks
    return results;
}


std::vector<std::vector<std::pair<int64_t, std::vector<double>>>>
geodesic_compute_parallel_surface_points(
    std::vector<std::pair<geomp::SurfacePoint, geomp::SurfacePoint>> surface_points_pairs,
    GC::DenseMatrix<double> verts, GC::DenseMatrix<int64_t> faces) {
    size_t iteration_count = surface_points_pairs.size();
    std::vector<std::vector<std::pair<int64_t, std::vector<double>>>> results(iteration_count);

    tf::Executor executor;
    tf::Taskflow taskflow;

    size_t worker_count = std::min(iteration_count, executor.num_workers());
    for (size_t w = 0; w < worker_count; ++w) {
        taskflow.emplace([&, w, worker_count]() {
            auto GM = std::make_unique<geomp::GeodesicsManager>(verts, faces);
            for (size_t i = w; i < iteration_count; i += worker_count) {
                auto surface_pts_pair = surface_points_pairs[i];
                auto start_surface_point = surface_pts_pair.first;
                auto end_surface_point = surface_pts_pair.second;
                results[i] = GM->find_surface_points_on_exact_geodesic_path(
                    start_surface_point.first, end_surface_point.first, start_surface_point.second,
                    end_surface_point.second);
            }
        });
    }


    executor.run(taskflow).get();  // Execute all tasks
    return results;
}

// Pybind11 Module Definition
void init_parallel_taskflow(py::module& m) {
    m.def("geodesic_compute_parallel", &geodesic_compute_parallel,
          "Parallel compute geodesic path using Taskflow");
    m.def("flip_geodesic_compute_parallel", &flip_geodesic_compute_parallel,
          "Parallel compute flip geodesic path using Taskflow");
    m.def("geodesic_compute_parallel_surface_points", &geodesic_compute_parallel_surface_points,
          "Parallel compute geodesic path and return surface points using Taskflow");
}