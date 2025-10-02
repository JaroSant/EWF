#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <deque>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "WrightFisher.h"

namespace py = pybind11;

double dt_default = 0.1;
double bt_default = 0.04;

PYBIND11_MODULE(EWF_pybind, m) {
  py::class_<WrightFisher>(m, "WrightFisher")
      .def(py::init<std::vector<double>, std::vector<std::vector<double>>, bool,
                    std::vector<double>, int, std::vector<double>, int,
                    std::vector<std::vector<double>>>(),
           py::arg("changepts"), py::arg("thetaP"), py::arg("non_neut"),
           py::arg("sigma"), py::arg("selectionSetup"), py::arg("dom"),
           py::arg("SelPolyDeg"), py::arg("selCoefs"))
      .def("DiffusionRunner", &WrightFisher::DiffusionRunner, py::arg("nSim"),
           py::arg("x"), py::arg("startT"), py::arg("endT"),
           py::arg("Absorption"), py::arg("Filename"),
           py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default)
      .def("DiffusionRunnerVector", &WrightFisher::DiffusionRunnerVector,
           py::arg("nSim"), py::arg("x"), py::arg("startT"), py::arg("endT"),
           py::arg("Absorption"), py::arg("Filename"),
           py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default)
      .def("DiffusionTrajectoryVector",
           &WrightFisher::DiffusionTrajectoryVector, py::arg("nSim"),
           py::arg("x"), py::arg("times"), py::arg("Absorption"),
           py::arg("Filename"), py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default)
      .def("DiffusionDensityCalculator",
           &WrightFisher::DiffusionDensityCalculator, py::arg("meshSize"),
           py::arg("x"), py::arg("startT"), py::arg("endT"),
           py::arg("Absorption"), py::arg("Filename"), py::arg("verbose"),
           py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default)
      .def("BridgeDiffusionRunner", &WrightFisher::BridgeDiffusionRunner,
           py::arg("nSim"), py::arg("x"), py::arg("z"), py::arg("startT"),
           py::arg("endT"), py::arg("sampleT"), py::arg("Absorption"),
           py::arg("Filename"), py::arg("verbose"),
           py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default)
      .def("BridgeDiffusionDensityCalculator",
           &WrightFisher::BridgeDiffusionDensityCalculator, py::arg("meshSize"),
           py::arg("x"), py::arg("z"), py::arg("startT"), py::arg("endT"),
           py::arg("sampleT"), py::arg("Absorption"), py::arg("Filename"),
           py::arg("verbose"), py::arg("diffusion_threshold") = dt_default,
           py::arg("bridge_threshold") = bt_default);
}