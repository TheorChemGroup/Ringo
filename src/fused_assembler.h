#ifndef ASSEMBLER_H
#define ASSEMBLER_H

#include "utils.h"


template <class MD, class AD>
class PolycycleAssembler {
    using umatrix4d_t = Utils::umatrix4d_t;
    using frame_t = Utils::frame_t;

    public:
        using mol_data_t = MD;
        using assembly_data_t = AD;
        using poly_t = typename mol_data_t::poly_t;
        using coord_container_t = typename assembly_data_t::coord_container_t;

        PolycycleAssembler() { }
        PolycycleAssembler(frame_t ind_frame, frame_t dep_frame, mol_data_t* mol_data, assembly_data_t& assembly_data);
        void configure();
        template <class IndC> void connect_cycles(coord_container_t& dep_coord, IndC& ind_nodes
            #ifdef KDMOL_LOG
                , const int ord
            #endif
            );

    private:
        frame_t ind_frame, dep_frame;
        poly_t *poly;
        se_matrix_t transit_mat;

    public:
        coord_container_t *coords;
};

template <class MD, class AD>
PolycycleAssembler<MD,AD>::PolycycleAssembler(frame_t ind_frame, frame_t dep_frame,
                                       mol_data_t* mol_data,
                                       assembly_data_t& assembly_data) :
    ind_frame(ind_frame),
    dep_frame(dep_frame),
    poly(&mol_data->polys[std::get<0>(ind_frame)])
{
    #ifdef KDMOL_DEBUG
        if (std::get<0>(ind_frame) != std::get<0>(dep_frame))
            throw std::runtime_error(fmt::format("Mismatch of central atoms. ind_frame={}, dep_frame={}", 
                                                 py::handle(py::repr(py::cast(ind_frame))).cast<std::string>(),
                                                 py::handle(py::repr(py::cast(dep_frame))).cast<std::string>()
                                    ));
    #endif
    configure();
}

template <class MD, class AD>
void PolycycleAssembler<MD,AD>::configure() {
    auto ind_matr = poly->get_localframe_3d({std::get<1>(ind_frame), std::get<2>(ind_frame)});
    auto dep_matr = poly->get_localframe_3d({std::get<1>(dep_frame), std::get<2>(dep_frame)});
    
    so_matrix_t inv_dep_matr, temp;
    Utils::invert_matrix(dep_matr, inv_dep_matr);
    noalias(temp) = prod(inv_dep_matr, ind_matr);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            transit_mat(i, j) = temp(i, j);
            
    for (int i = 0; i < 3; i++) {
        transit_mat(i, 3) = 0.0;
        transit_mat(3, i) = 0.0;
    }
    transit_mat(3, 3) = 1.0;

    checkpoint("assembler_config");
    #ifdef KDMOL_LOG
        log_routine("assembler_config");
        log_in("center", poly->get_center());
        log_in("ind_frame", ind_frame);
        log_in("dep_frame", dep_frame);
        log_check("transit_mat", transit_mat);
        log_finish();
    #endif
}

template <class MD, class AD>
template <class IndC>
void PolycycleAssembler<MD,AD>::connect_cycles(coord_container_t& dep_coord, IndC& ind_nodes
    #ifdef KDMOL_LOG
        , const int ord
    #endif
    )
{
    checkpoint("assembler_start");
    #ifdef KDMOL_LOG
        {
        py::list log_depcoords, log_indcoords;

        for (auto& [idx, xyz] : dep_coord) {
            py::list newitem;
            newitem.append(idx);
            newitem.append(xyz(0));
            newitem.append(xyz(1));
            newitem.append(xyz(2));
            log_depcoords.append(newitem);
        }

        for (auto& [ idx, _ ] : ind_nodes) {
            py::list newitem;
            newitem.append(idx);
            newitem.append((*coords)[idx](0));
            newitem.append((*coords)[idx](1));
            newitem.append((*coords)[idx](2));
            log_indcoords.append(newitem);
        }

        log_routine("assembler_start");
        log_in("call_idx", Logger::aps_count);
        log_in("center", poly->get_center());
        log_in("ind_frame", ind_frame);
        log_in("dep_frame", dep_frame);
        log_in("ord", ord);
        log_check("transit_mat", transit_mat);
        log_check("depcoords_SET1", log_depcoords);
        log_check("indcoords_SET1", log_indcoords);
        log_finish();
        }
    #endif

    auto ind_matr = poly->get_frame(*coords, {std::get<1>(ind_frame), std::get<2>(ind_frame)}, std::get<0>(ind_frame));
    auto dep_matr = poly->get_frame(dep_coord, {std::get<1>(dep_frame), std::get<2>(dep_frame)}, std::get<0>(dep_frame));
    se_matrix_t inv_ind_matr, temp, transition;
    Utils::invert_matrix(ind_matr, inv_ind_matr);
    noalias(temp) = prod(dep_matr, transit_mat);
    noalias(transition) = prod(temp, inv_ind_matr);

    for (auto& [ ind_atom, _ ] : ind_nodes) {
        auto& xyz = (*coords)[ind_atom];
        uvector4d_t xyz_ext, res;
        for (int i = 0; i < 3; ++i)
            xyz_ext(i) = xyz(i);
        xyz_ext(3) = 1.0;
        res = prod(transition, xyz_ext);
        xyz = boost::numeric::ublas::subrange(res, 0, 3);
    }

    checkpoint("assembler_res");
    #ifdef KDMOL_LOG
        {
        py::list log_indcoords;
        for (auto& [ idx, _ ] : ind_nodes) {
            py::list newitem;
            newitem.append(idx);
            newitem.append((*coords)[idx](0));
            newitem.append((*coords)[idx](1));
            newitem.append((*coords)[idx](2));
            log_indcoords.append(newitem);
        }

        log_routine("assembler_res");
        log_in("call_idx", Logger::aps_count);
        log_in("center", poly->get_center());
        log_in("ind_frame", ind_frame);
        log_in("dep_frame", dep_frame);
        log_check("indcoords_SET1", log_indcoords);
        log_finish();
        }
    #endif
}

#endif //ASSEMBLER_H
