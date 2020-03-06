#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../stella/ses-pipeline/c/cmpfit-1.2/MyFit.h"

namespace py = pybind11;

PYBIND11_MODULE(myMPFit, m) {
    m.doc() = "pybind11 mpfit plugin"; // optional module docstring

    m.def("MPFitGaussFuncNB", (int (*)(int m, int n, double *p, double *dy, double **dvec, void *vars)) &MPFitGaussFuncNB);
//    m.def("readCSVFile", (CSVData (*)(string const&)) &readCSVFile);
//    m.def("convertStringVectortoDoubleVector", (vector<double> (*)(vector<string> const&)) &convertStringVectortoDoubleVector);
}
