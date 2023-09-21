#ifndef GENERIC_FRAGMENT_H
#define GENERIC_FRAGMENT_H


template <class G, class MD, class DF, class DD, class D>
class GenericFragment {
    public:
        using main_graph_t = G;
        using mol_data_t = MD;
        using dofs_t = DF;
        using discrete_dofs_t = DD;
        using dfg_t = D;

        GenericFragment () { }

        GenericFragment(mol_data_t* mol_data, dofs_t* dofs, discrete_dofs_t* discrete_dofs) :
            mol_data(mol_data), dofs(dofs), discrete_dofs(discrete_dofs) { }
        
        virtual void configure_assembly(const main_graph_t& molgr, const dfg_t& dfg, const int i) = 0;

        virtual void configure_assembly(const main_graph_t& molgr, const dfg_t& dfg) = 0;

    protected:
        mol_data_t* mol_data;
        dofs_t* dofs;
        discrete_dofs_t* discrete_dofs;
};

#endif // GENERIC_FRAGMENT_H
