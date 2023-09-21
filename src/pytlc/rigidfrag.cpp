#include "rigidfrag.h"


void RigidFrag::initialize(AtomItem left_atom, AtomItem right_atom, AtomGroup& catoms, BondGroup& cbonds) {
    int left_idx = left_atom.get_num(),
        right_idx = right_atom.get_num();
    unsigned int length = catoms.get_distance(left_idx, right_idx);
    iatoms = AtomContainerDynamic(length - 3);
    ibonds = BondContainerDynamic(length - 4);

    eatoms.set(&catoms[left_atom.get_num()], 0);
    eatoms.set(&catoms[right_atom.get_num()], 1);
    ebonds.set(&left_atom.get_right_bond(), 0);
    ebonds.set(&right_atom.get_left_bond(), 1);

    if(length == 1) { // right_atom - 1 == left_atom + 1
        mycase = RFCase::atom;
        satoms.set(&catoms[(left_atom + 1).get_num()], 0);
    } else if (length == 2) { // right_atom == left_atom + 3
        mycase = RFCase::bond;
        AtomItem &left_satom = catoms[(left_atom + 1).get_num()],
                 &right_satom = catoms[(right_atom - 1).get_num()];
        satoms.set(&left_satom, 0);
        satoms.set(&right_satom, 1);
        sbonds.set(&left_satom.get_right_bond(), 0);
    } else if (length > 2){ 
        mycase = RFCase::general;
        AtomItem &left_satom = catoms[(left_atom + 1).get_num()],
                 &right_satom = catoms[(right_atom - 1).get_num()];
        satoms.set(&catoms[(left_atom + 1).get_num()], 0);
        satoms.set(&catoms[(right_atom - 1).get_num()], 1);
        sbonds.set(&left_satom.get_right_bond(), 0);
        sbonds.set(&right_satom.get_left_bond(), 1);

        iatoms += catoms[(left_atom + 2).get_num()];
        BondItem cur_pos = (left_satom + 1).get_right_bond(), finish_pos = right_satom.get_left_bond();
        while(cur_pos != finish_pos) {
            ibonds += cbonds[cur_pos.get_num()];
            iatoms += cur_pos.get_right_atom();
            cur_pos++;
        }
        iatoms_xyz = Matrix(iatoms.size(), 3);
        // For quick lookup in a buffer later:
        eatom_addr = eatoms[0].get_num() * 3;
        satomL_addr = satoms[0].get_num() * 3;
        satomR_addr = satoms[1].get_num() * 3;
        firstiatom_addr = iatoms[0].get_num() * 3;
        do_Bseg = false;
        for (std::size_t i = 1; i < iatoms.size(); ++i)
            if (iatoms[i].get_num() == 0) {
                do_Bseg = true;
                segA_length = i;
                segB_length = (iatoms.size() - i);
                break;
            }
        if (!do_Bseg)
            segA_length = iatoms.size();
                
    } else
        throw std::runtime_error("The fragment is too small. This shouldn't have happened");
    
    #ifdef EXT_LOG
        Logger::log_rstart("rfinit");
        py::list l;
        l.append(left_atom.get_num());
        l.append(right_atom.get_num());
        Logger::log(fmt::format(" IN:: axatoms = {}", py::repr(l).cast<std::string>()));

        Logger::log(fmt::format(" CHECK:: eatoms = {}", eatoms.repr_list()));
        if (mycase == RFCase::atom)
            Logger::log(fmt::format(" CHECK:: satoms = {}", satoms.repr_list(1)));
        else
            Logger::log(fmt::format(" CHECK:: satoms = {}", satoms.repr_list()));
        
        Logger::log(fmt::format(" CHECK:: iatoms = {}", iatoms.repr_list()));
        Logger::log(fmt::format(" CHECK:: ebonds = {}", ebonds.repr_list()));
        if (mycase == RFCase::atom)
            Logger::log(fmt::format(" CHECK:: sbonds = {}", sbonds.repr_list(0)));
        else if (mycase == RFCase::bond)
            Logger::log(fmt::format(" CHECK:: sbonds = {}", sbonds.repr_list(1)));
        else
            Logger::log(fmt::format(" CHECK:: sbonds = {}", sbonds.repr_list()));
        
        Logger::log(fmt::format(" CHECK:: ibonds = {}", ibonds.repr_list()));
        Logger::log_rend("rfinit");
    #endif
}

void RigidFrag::write_tlc_data_work(AtomGroup& catoms, BondGroup& cbonds) noexcept {
    #define idx_length 3 * num + 1
    #define idx_langle 3 * num + 1
    #define idx_rangle 3 * num + 2
    #define idx_dihedral num
    if (mycase == RFCase::atom) {
        cbonds.tlc_lengths[idx_length] = 0.0;
        catoms.tlc_vangles[idx_langle] = 90.0;
        catoms.tlc_vangles[idx_rangle] = 90.0;
        cbonds.tlc_dihedrals[idx_dihedral] = catoms.vangles[this->satoms[0].get_num()] * RAD2DEG;
    } else if (mycase == RFCase::bond) {
        cbonds.tlc_lengths[idx_length] = cbonds.bonds[this->sbonds[0].get_num()];
        catoms.tlc_vangles[idx_langle] = catoms.vangles[this->satoms[0].get_num()] * RAD2DEG;
        catoms.tlc_vangles[idx_rangle] = catoms.vangles[this->satoms[1].get_num()] * RAD2DEG;
        cbonds.tlc_dihedrals[idx_dihedral] = cbonds.get_dihedral(this->sbonds[0].get_num()) * RAD2DEG;
    } else { // mycase == RFCase::general
        SEMatrix rt_m = IMatrix(4), temp_rt = IMatrix(4);
        SEMatrix dist_m = IMatrix(4), vang_m = IMatrix(4), tang_m = IMatrix(4);
        dist_m(0, 3) = cbonds.bonds[this->sbonds[0].get_num()];
        XYZMatrix atom_pos;
        Matrix iatoms_temp(iatoms_xyz.size1(), 3);
        atom_pos(0, 0) = cbonds.bonds[this->ebonds[0].get_num()];
        
        double tau = -cbonds.get_dihedral(this->sbonds[0].get_num()), ct = cos(tau), st = sin(tau),
               alpha = M_PI - catoms.vangles[this->satoms[0].get_num()], ca = cos(alpha), sa = sin(alpha);
        
        tang_m(1,1) = ct;
        tang_m(2,2) = ct;
        tang_m(1,2) = -st;
        tang_m(2,1) = st;

        vang_m(0,0) = ca;
        vang_m(1,1) = ca;
        vang_m(0,1) = -sa;
        vang_m(1,0) = sa;

        noalias(rt_m) = prod(vang_m, dist_m);
        noalias(temp_rt) = prod(rt_m, tang_m);

        #ifdef EXT_LOG
            Logger::log_rstart("effparamA");
            Logger::log(fmt::format(" IN:: vangle = {}", catoms.vangles[this->satoms[0].get_num()]));
            Logger::log(fmt::format(" IN:: length = {}", cbonds.bonds[this->sbonds[0].get_num()]));
            Logger::log(fmt::format(" IN:: dihedral = {}", cbonds.get_dihedral(this->sbonds[0].get_num())));
            Logger::log(fmt::format(" CHECK:: RTmat = {}", TlcUtils::repr_matrix(temp_rt)));
            Logger::log_rend("effparamA");
        #endif

        const int n_iter = iatoms.size() / 2;
        int j = 0;
        bool extra_step = (iatoms.size() % 2 == 1);
        for(int i = 0; i < n_iter; ++i) {
            insert_iatom_xyz(iatoms_temp, temp_rt, j);
            dist_m(0, 3) = cbonds.bonds[this->iatoms[j].get_num()];
            tau = -cbonds.get_dihedral(this->iatoms[j].get_num());
            ct = cos(tau);
            st = sin(tau);
            alpha = M_PI - catoms.vangles[this->iatoms[j].get_num()];
            ca = cos(alpha);
            sa = sin(alpha);
        
            tang_m(1,1) = ct;
            tang_m(2,2) = ct;
            tang_m(1,2) = -st;
            tang_m(2,1) = st;

            vang_m(0,0) = ca;
            vang_m(1,1) = ca;
            vang_m(0,1) = -sa;
            vang_m(1,0) = sa;

            noalias(rt_m) = prod(temp_rt, vang_m);
            noalias(temp_rt) = prod(rt_m, dist_m);
            noalias(rt_m) = prod(temp_rt, tang_m);
            #ifdef EXT_LOG
                Logger::log_rstart("effparamX");
                Logger::log(fmt::format(" IN:: vangle = {}", catoms.vangles[this->satoms[0].get_num()]));
                Logger::log(fmt::format(" IN:: length = {}", cbonds.bonds[this->sbonds[0].get_num()]));
                Logger::log(fmt::format(" IN:: dihedral = {}", cbonds.get_dihedral(this->sbonds[0].get_num())));
                Logger::log(fmt::format(" IN:: i = {}", j));
                Logger::log(fmt::format(" CHECK:: RTmat = {}", TlcUtils::repr_matrix(rt_m)));
                Logger::log_rend("effparamX");
            #endif
            j++;
            insert_iatom_xyz(iatoms_temp, rt_m, j);
            
            dist_m(0, 3) = cbonds.bonds[this->iatoms[j].get_num()];
            tau = -cbonds.get_dihedral(this->iatoms[j].get_num());
            ct = cos(tau);
            st = sin(tau);
            alpha = M_PI - catoms.vangles[this->iatoms[j].get_num()];
            ca = cos(alpha);
            sa = sin(alpha);

            tang_m(1,1) = ct;
            tang_m(2,2) = ct;
            tang_m(1,2) = -st;
            tang_m(2,1) = st;

            vang_m(0,0) = ca;
            vang_m(1,1) = ca;
            vang_m(0,1) = -sa;
            vang_m(1,0) = sa;

            noalias(temp_rt) = prod(rt_m, vang_m);
            noalias(rt_m) = prod(temp_rt, dist_m);
            noalias(temp_rt) = prod(rt_m, tang_m);
            #ifdef EXT_LOG
                Logger::log_rstart("effparamX");
                Logger::log(fmt::format(" IN:: vangle = {}", catoms.vangles[this->satoms[0].get_num()]));
                Logger::log(fmt::format(" IN:: length = {}", cbonds.bonds[this->sbonds[0].get_num()]));
                Logger::log(fmt::format(" IN:: dihedral = {}", cbonds.get_dihedral(this->sbonds[0].get_num())));
                Logger::log(fmt::format(" IN:: i = {}", j));
                Logger::log(fmt::format(" CHECK:: RTmat = {}", TlcUtils::repr_matrix(temp_rt)));
                Logger::log_rend("effparamX");
            #endif
            j++;
        }

        if (extra_step) {
            insert_iatom_xyz(iatoms_temp, temp_rt, j);
            dist_m(0, 3) = cbonds.bonds[this->iatoms[j].get_num()];
            tau = -cbonds.get_dihedral(this->iatoms[j].get_num());
            ct = cos(tau);
            st = sin(tau);
            alpha = M_PI - catoms.vangles[this->iatoms[j].get_num()];
            ca = cos(alpha);
            sa = sin(alpha);

            tang_m(1,1) = ct;
            tang_m(2,2) = ct;
            tang_m(1,2) = -st;
            tang_m(2,1) = st;

            vang_m(0,0) = ca;
            vang_m(1,1) = ca;
            vang_m(0,1) = -sa;
            vang_m(1,0) = sa;

            noalias(rt_m) = prod(temp_rt, vang_m);
            noalias(temp_rt) = prod(rt_m, dist_m);
            noalias(rt_m) = prod(temp_rt, tang_m);
            #ifdef EXT_LOG
                Logger::log_rstart("effparamX");
                Logger::log(fmt::format(" IN:: vangle = {}", catoms.vangles[this->satoms[0].get_num()]));
                Logger::log(fmt::format(" IN:: length = {}", cbonds.bonds[this->sbonds[0].get_num()]));
                Logger::log(fmt::format(" IN:: dihedral = {}", cbonds.get_dihedral(this->sbonds[0].get_num())));
                Logger::log(fmt::format(" IN:: i = {}", j));
                Logger::log(fmt::format(" CHECK:: RTmat = {}", TlcUtils::repr_matrix(rt_m)));
                Logger::log_rend("effparamX");
            #endif
        } else
            noalias(rt_m) = temp_rt;
        
        #ifdef EXT_LOG
            Logger::log_rstart("effparamB");
            Logger::log(fmt::format(" IN:: vangle = {}", catoms.vangles[this->satoms[0].get_num()]));
            Logger::log(fmt::format(" IN:: length = {}", cbonds.bonds[this->sbonds[0].get_num()]));
            Logger::log(fmt::format(" IN:: dihedral = {}", cbonds.get_dihedral(this->sbonds[0].get_num())));
            Logger::log(fmt::format(" CHECK:: RTmat = {}", TlcUtils::repr_matrix(rt_m)));
            Logger::log(fmt::format(" CHECK:: xyzs = {}", TlcUtils::repr_matrix(iatoms_temp)));
            Logger::log_rend("effparamB");
        #endif
        
        alpha = M_PI - catoms.vangles[this->satoms[1].get_num()];
        ca = cos(alpha);
        sa = sin(alpha);
        vang_m(0,0) = ca;
        vang_m(1,1) = ca;
        vang_m(0,1) = -sa;
        vang_m(1,0) = sa;
        noalias(temp_rt) = prod(rt_m, vang_m);

        #ifdef EXT_LOG
            Logger::log_rstart("effparamC");
            Logger::log(fmt::format(" IN:: vangle = {}", catoms.vangles[this->satoms[0].get_num()]));
            Logger::log(fmt::format(" IN:: length = {}", cbonds.bonds[this->sbonds[0].get_num()]));
            Logger::log(fmt::format(" IN:: dihedral = {}", cbonds.get_dihedral(this->sbonds[0].get_num())));
            Logger::log(fmt::format(" CHECK:: RTmat = {}", TlcUtils::repr_matrix(temp_rt)));
            Logger::log_rend("effparamC");
        #endif
        
        atom_pos(2, 0) = temp_rt(0, 3);
        atom_pos(2, 1) = temp_rt(1, 3);
        atom_pos(2, 2) = temp_rt(2, 3);
        dist_m(0, 3) = cbonds.bonds[this->ebonds[1].get_num()];
        noalias(rt_m) = prod(temp_rt, dist_m);
        atom_pos(3, 0) = rt_m(0, 3);
        atom_pos(3, 1) = rt_m(1, 3);
        atom_pos(3, 2) = rt_m(2, 3);
        #ifdef EXT_LOG
            Logger::log_rstart("effparamD");
            Logger::log(fmt::format(" IN:: vangle = {}", catoms.vangles[this->satoms[0].get_num()]));
            Logger::log(fmt::format(" IN:: length = {}", cbonds.bonds[this->sbonds[0].get_num()]));
            Logger::log(fmt::format(" IN:: dihedral = {}", cbonds.get_dihedral(this->sbonds[0].get_num())));
            Logger::log(fmt::format(" CHECK:: RTmat = {}", TlcUtils::repr_matrix(rt_m)));
            Logger::log(fmt::format(" CHECK:: atom_pos(2) = {}", TlcUtils::repr_row(atom_pos, 2)));
            Logger::log(fmt::format(" CHECK:: atom_pos(3) = {}", TlcUtils::repr_row(atom_pos, 3)));
            Logger::log_rend("effparamD");
        #endif

        FrameMatrix frame;
        frame(0, 0) = 1.0;
        frame(0, 1) = -atom_pos(2, 0);
        frame(1, 1) = -atom_pos(2, 1);
        frame(2, 1) = -atom_pos(2, 2);
        build_frame_Lopt(frame);
        noalias(iatoms_xyz) = prod(iatoms_temp, frame);

        double firstside_dir[] = {-1.0, 0.0, 0.0};
        cbonds.tlc_lengths[idx_length] = cblas_dnrm2(3, &atom_pos(2, 0), 1);
        catoms.tlc_vangles[idx_langle] = TlcUtils::vangle_from_dirsA(firstside_dir, &atom_pos(2, 0)) * RAD2DEG;
        catoms.tlc_vangles[idx_rangle] = TlcUtils::vangle_from_dirsB(&rt_m(0, 0), &atom_pos(2, 0)) * RAD2DEG;
        cbonds.tlc_dihedrals[idx_dihedral] = TlcUtils::dihedral_from_dirs(firstside_dir, &atom_pos(2, 0), &rt_m(0, 0)) * RAD2DEG;
        #ifdef EXT_LOG
            Logger::log_rstart("effparamE");
            Logger::log(fmt::format(" IN:: vangle = {}", catoms.vangles[this->satoms[0].get_num()]));
            Logger::log(fmt::format(" IN:: length = {}", cbonds.bonds[this->sbonds[0].get_num()]));
            Logger::log(fmt::format(" IN:: dihedral = {}", cbonds.get_dihedral(this->sbonds[0].get_num())));
            Logger::log(fmt::format(" CHECK:: leng = {}", cbonds.tlc_lengths[idx_length]));
            Logger::log(fmt::format(" CHECK:: vang1 = {}", catoms.tlc_vangles[idx_langle]));
            Logger::log(fmt::format(" CHECK:: vang2 = {}", catoms.tlc_vangles[idx_rangle]));
            Logger::log(fmt::format(" CHECK:: tang = {}", -cbonds.tlc_dihedrals[idx_dihedral]));
            Logger::log(fmt::format(" CHECK:: xyzs = {}", TlcUtils::repr_matrix(iatoms_xyz)));
            Logger::log_rend("effparamE");
        #endif
    }
}

void RigidFrag::reconstruct(double* buffer, AtomGroup& catoms, BondGroup& cbonds) noexcept {
    if (mycase != RFCase::general)
        return;
    
    for (std::size_t i = 0; i < segA_length; ++i)
        cblas_dcopy(3, &buffer[satomL_addr], 1, &buffer[firstiatom_addr + i * 3], 1);
    if (do_Bseg)
        for (std::size_t i = 0; i < segB_length; ++i)
            cblas_dcopy(3, &buffer[satomL_addr], 1, &buffer[i * 3], 1);

    double frame[9];
    cblas_dcopy(3, &buffer[satomL_addr], 1, &frame[0], 1);
    cblas_dcopy(3, &buffer[satomL_addr], 1, &frame[3], 1);
    cblas_daxpy(3, -1.0, &buffer[eatom_addr], 1, &frame[0], 1);
    cblas_daxpy(3, -1.0, &buffer[satomR_addr], 1, &frame[3], 1);
    build_frame(frame);
    auto N = iatoms.size();
    
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, segA_length, 3, 3, 1.0, &iatoms_xyz(0, 0), 3, &frame[0], 3, 1.0, &buffer[firstiatom_addr], 3);
    if (do_Bseg)
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, segB_length, 3, 3, 1.0, &iatoms_xyz(segA_length, 0), 3, &frame[0], 3, 1.0, &buffer[0], 3);

    #ifdef EXT_LOG
        py::list res;
        for (std::size_t i = 0; i < segA_length; ++i) {
            py::list temp;
            for (unsigned int j = 0; j < 3; ++j)
                temp.append(buffer[firstiatom_addr + i * 3 + j]);
            res.append(temp);
        }
        if (do_Bseg)
            for (std::size_t i = 0; i < segB_length; ++i) {
                py::list temp;
                for (unsigned int j = 0; j < 3; ++j)
                    temp.append(buffer[i * 3 + j]);
                res.append(temp);
            }
        
        Logger::log_rstart("setsolFrag");
        Logger::log(fmt::format(" IN:: vangle = {}", catoms.vangles[this->satoms[0].get_num()]));
        Logger::log(fmt::format(" IN:: length = {}", cbonds.bonds[this->sbonds[0].get_num()]));
        Logger::log(fmt::format(" IN:: dihedral = {}", cbonds.get_dihedral(this->sbonds[0].get_num())));
        // Logger::log(fmt::format(" CHECK:: M_start = {}", py::repr(res).cast<std::string>()));
        Logger::log(fmt::format(" CHECK:: XYZ_old = {}", TlcUtils::repr_matrix(iatoms_xyz)));
        Logger::log(fmt::format(" CHECK:: FrameMatr = {}", TlcUtils::repr_matrix_buffer(frame, 3, 3)));

        res = py::list();
        for (std::size_t i = 0; i < segA_length; ++i) {
            py::list temp;
            for (unsigned int j = 0; j < 3; ++j)
                temp.append(buffer[firstiatom_addr + i * 3 + j]);
            res.append(temp);
        }
        if (do_Bseg)
            for (std::size_t i = 0; i < segB_length; ++i) {
                py::list temp;
                for (unsigned int j = 0; j < 3; ++j)
                    temp.append(buffer[i * 3 + j]);
                res.append(temp);
            }
        Logger::log(fmt::format(" CHECK:: M = {}", py::repr(res).cast<std::string>()));
        // Logger::log(fmt::format(" CHECK:: xyzs = {}", TlcUtils::repr_matrix(iatoms_xyz)));
        Logger::log_rend("setsolFrag");
    #endif
}
