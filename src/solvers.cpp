#include "solvers.h"


std::pair<std::vector<int>, std::unordered_map<int, int>> TLCSolverModel::walk_cycle(py::object& gr) {
    const int nnodes = gr.attr("number_of_nodes")().cast<int>();
    std::vector<int> idx2atom; // This is a walk around the cycle
    std::unordered_map<int, int> atom2idx;
    idx2atom.reserve(nnodes);
    auto node_list = gr.attr("nodes").cast<py::list>();
    int min_node = -1;
    for (const auto& py_node : node_list) {
        auto cpp_node = py_node.cast<int>();
        if ((min_node == -1) || (cpp_node < min_node)) {
            min_node = cpp_node;
        }
    }
    idx2atom.push_back(min_node);
    atom2idx[idx2atom[0]] = 0;

    auto nblist = py::handle(gr.attr("neighbors")(idx2atom[0])).cast<py::list>();
    int vA = nblist[0].cast<int>(), vB = nblist[1].cast<int>();
    auto cur_node = std::min(vA, vB);
    int i = 1;
    while (cur_node != idx2atom[0]) {
        idx2atom.push_back(cur_node);
        atom2idx[cur_node] = i;

        auto nbnodes = gr.attr("neighbors")(cur_node).cast<py::list>();
        if (nbnodes[0].cast<int>() == idx2atom[i - 1])
            cur_node = nbnodes[1].cast<int>();
        else
            cur_node = nbnodes[0].cast<int>();
        i++;
    }

    return std::make_pair(idx2atom, atom2idx);
}
