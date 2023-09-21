#include "utils.h"


namespace Utils {
    #ifdef KDMOL_LOG
        namespace LoggerWrapper {
            void log_routine(const std::string & routine_name) { Logger::log_routine(routine_name); }
            void log_finish() { Logger::log_finish(); }
        }
    #endif

    py::module_ *PyModules::nx, *PyModules::pyutils;
    std::vector<std::string> messages;
    py::list get_status_feed() {
        return py::cast(messages);
    }
    
    std::string assemble_message(const std::string& subject, const std::string& message, const std::vector<int>& atoms, const std::string& filename, const int linenumber) {
        nlohmann::json message_data = {
            {"message", message},
            {"subject", subject},
            {"atoms", atoms},
            {"file", filename},
            {"line", linenumber},
        };
        return message_data.dump();
    }

    void add_message_to_feed(const std::string& json_contents) {
        auto it = std::find(messages.cbegin(), messages.cend(), json_contents);
        if (it == messages.end())
            messages.push_back(json_contents);
    }

    void add_message_to_feed(const std::string& subject, const std::string& message, const std::vector<int>& atoms, const std::string& filename, const int linenumber) {
        const auto json_contents = assemble_message(subject, message, atoms, filename, linenumber);

        auto it = std::find(messages.cbegin(), messages.cend(), json_contents);
        if (it == messages.cend())
            messages.push_back(json_contents);
    }

    void clear_status_feed() {
        messages.clear();
    }

    #ifdef KDMOL_LOG
        std::string Logger::m_log;
        std::string Logger::m_cur_routine;
        int Logger::aps_count;
    #endif


    std::vector<std::string> readlines(const std::string& filename) {
        std::ifstream infile(filename);
        std::string str;
        std::vector<std::string> reslines;

        if(!infile)
            throw std::runtime_error("File " + filename + " not found");

        while (std::getline(infile, str)) {
			erase_remove_if(str, InvalidChar());
			reslines.push_back(str);
        }
        return reslines;
    }

    uvector3d_t gs_rand(uvector3d_t& x, uvector3d_t& y) {
        using namespace boost::numeric::ublas;
        uvector3d_t z;
        auto norm = norm_2(x);
        for(int i = 0; i < 3; ++i)
            x[i] /= norm;
        auto dot = inner_prod(x, y);
        noalias(y) -= x*dot;
        norm = norm_2(y);
        for(int i = 0; i < 3; ++i)
            y[i] /= norm;
        calc_cross(&x[0], &y[0], &z[0]);
        return z;
    }

    #ifdef BUILD_OVERLAP_DETECTION
        std::unordered_map<std::string, double> VDW_RADII = {
            {"H", 1.1}, {"He", 1.4}, {"Li", 1.82}, {"Be", 1.53}, {"B", 1.92}, {"C", 1.7}, {"N", 1.55}, {"O", 1.52}, {"F", 1.47}, {"Ne", 1.54}, {"Na", 2.27}, {"Mg", 1.73}, {"Al", 1.84}, {"Si", 2.1}, {"P", 1.8}, {"S", 1.8}, {"Cl", 1.75}, {"Ar", 1.88}, {"K", 2.75}, {"Ca", 2.31}, {"Sc", 2.15}, {"Ti", 2.11}, {"V", 2.07}, {"Cr", 2.06}, {"Mn", 2.05}, {"Fe", 2.04}, {"Co", 2.0}, {"Ni", 1.97}, {"Cu", 1.96}, {"Zn", 2.01}, {"Ga", 1.87}, {"Ge", 2.11}, {"As", 1.85}, {"Se", 1.9}, {"Br", 1.85}, {"Kr", 2.02}, {"Rb", 3.03}, {"Sr", 2.49}, {"Y", 2.32}, {"Zr", 2.23}, {"Nb", 2.18}, {"Mo", 2.17}, {"Tc", 2.16}, {"Ru", 2.13}, {"Rh", 2.1}, {"Pd", 2.1}, {"Ag", 2.11}, {"Cd", 2.18}, {"In", 1.93}, {"Sn", 2.17}, {"Sb", 2.06}, {"Te", 2.06}, {"I", 1.98}, {"Xe", 2.16}, {"Cs", 3.43}, {"Ba", 2.68}, {"La", 2.43}, {"Ce", 2.42}, {"Pr", 2.4}, {"Nd", 2.39}, {"Pm", 2.38}, {"Sm", 2.36}, {"Eu", 2.35}, {"Gd", 2.34}, {"Tb", 2.33}, {"Dy", 2.31}, {"Ho", 2.3}, {"Er", 2.29}, {"Tm", 2.27}, {"Yb", 2.26}, {"Lu", 2.24}, {"Hf", 2.23}, {"Ta", 2.22}, {"W", 2.18}, {"Re", 2.16}, {"Os", 2.16}, {"Ir", 2.13}, {"Pt", 2.13}, {"Au", 2.14}, {"Au", 2.14}, {"Hg", 2.23}, {"Tl", 1.96}, {"Pb", 2.02}, {"Bi", 2.07}, {"Po", 1.97}, {"At", 2.02}, {"Rn", 2.2}, {"Fr", 3.48}, {"Ra", 2.83}, {"Ac", 2.47}, {"Th", 2.45}, {"Pa", 2.43}, {"U", 2.41}, {"Np", 2.39}, {"Pu", 2.43}, {"Am", 2.44}
        };

        double RADIUS_MULT = 0.5;
        void set_radius_multiplier(double rad_mult) {
            RADIUS_MULT = rad_mult;
        }
        py::dict get_vdw_radii() {
            return py::cast(VDW_RADII);
        }
        void set_vdw_radii(py::dict new_radii) {
            VDW_RADII = new_radii.cast<std::unordered_map<std::string, double>>();
        }
    #endif

    std::unordered_map<std::string, std::string> warning_codes = {
        {"IK_NOT_APPLIED", "IK Not Applied"},
        {"SUBOPTIMAL_SOLN_SEQ", "Suboptimal solution sequence"},
        {"UNMET_DOF_REQUEST", "Unfulfilled DOF request"},
        {"NO_CONFORMERS", "No conformers found"},
        {"FREQUENT_TLC_FAILS", "Frequent TLC failures"},
        {"FREQUENT_GEOM_FAILS", "Frequent geometry failures"}
        #ifdef BUILD_OVERLAP_DETECTION
            , {"UNKNOWN_ELEMENT", "Unknown element"}
        #endif
    };
}