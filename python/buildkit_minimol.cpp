#include "buildkit-minimol.h"
#include "type_conversions.h"
namespace py = pybind11;
using namespace clipper;
// namespace bk = buildkit;
//  23 Feb 23
//  Need to wrap Spacegroup, Cell, typecast String
//  then whole of minimol.h? for iterables
// refer to gemmi mol.cpp to wrap the add_child/set_child etc

// all these need to replace .begin() with something else like [0], and
// children ... there is no access to children (private)
template <typename Item>
void delitem_at_index(std::vector<Item> &items, pybind11::ssize_t index)
{
    items.erase(items.begin() + index);
}

template <typename Item>
void delitem_range(std::vector<Item> &items, pybind11::ssize_t start, pybind11::ssize_t end)
{
    items.erase(items.begin() + start, items.begin() + end);
}

template <typename Items>
void delitem_slice(Items &items, const pybind11::slice &slice)
{
    py::ssize_t start, stop, step, slice_len;
    if (!slice.compute((py::ssize_t)items.size(), &start, &stop, &step, &slice_len))
        throw py::error_already_set();
    if (step == 1)
    {
        delitem_range(items, start, start + slice_len);
    }
    else
    {
        for (int i = 0; i < slice_len; ++i)
            delitem_at_index(items, start + (step > 0 ? slice_len - 1 - i : i) * step);
    }
}

template <typename T, typename C>
C &add_item(T &container, C child, int pos)
{
    if ((size_t)pos > container.size()) // true for negative pos
        pos = (int)container.size();
    return *container.insert(container.begin() + pos, std::move(child));
}

template <typename P, typename C>
C &add_child(P &parent, C child, int pos)
{
    return add_item(parent.children(), std::move(child), pos);
}

template <typename P, typename C>
C &get_child(P &parent, int index)
{
    auto &children = parent.children();
    return children[normalise_index(index, children)];
}

template <typename P, typename C>
void set_child(P &parent, int index, C &child)
{
    auto &children = parent.children();
    children[normalise_index(index, children)] = child;
}

template <typename P>
void remove_child(P &parent, int index)
{
    auto &children = parent.children();
    children.erase(children.begin() + normalise_index(index, children));
}

template <typename P>
void remove_children(P &parent, py::slice slice)
{
    delitem_slice(parent.children(), slice);
}

void init_minimol(py::module &m)
{
    // "Forward declaration" of python classes to avoid
    // C++ signatures in docstrings ?
    py::class_<MAtom> pyAtom(m, "MAtom");
    py::class_<MResidue> pyResidue(m, "MResidue");
    py::class_<MChain> pyChain(m, "MChain");
    py::class_<MModel> pyModel(m, "MModel");

    py::enum_<MM::MODE>(m, "MODE")
        .value("UNIQUE", MM::MODE::UNIQUE)
        .value("ANY", MM::MODE::ANY);

    py::enum_<MM::COPY>(m, "COPY")
        .value("COPY_NONE", MM::COPY::COPY_NONE)
        .value("COPY_M", MM::COPY::COPY_M)
        .value("COPY_P", MM::COPY::COPY_P)
        .value("COPY_MP", MM::COPY::COPY_MP)
        .value("COPY_C", MM::COPY::COPY_C)
        .value("COPY_MC", MM::COPY::COPY_MC)
        .value("COPY_PC", MM::COPY::COPY_PC)
        .value("COPY_MPC", MM::COPY::COPY_MPC)
        .value("MEMBERS", MM::COPY::MEMBERS)
        .value("PROPERTIES", MM::COPY::PROPERTIES)
        .value("CHILDREN", MM::COPY::CHILDREN);

    py::enum_<MResidue::TYPE>(m, "TYPE")
        .value("Default", MResidue::TYPE::Default)
        .value("Dunbrack", MResidue::TYPE::Dunbrack)
        .value("Richardson", MResidue::TYPE::Richardson);

    pyAtom
        .def(py::init<>())
        .def_property("id", &MAtom::id, &MAtom::set_id)
        .def_property("name", &MAtom::name, &MAtom::set_name)
        .def_property("element", &MAtom::element, &MAtom::set_element)
        .def_property("occupancy", &MAtom::occupancy, &MAtom::set_occupancy)
        .def_property("u_iso", &MAtom::u_iso, &MAtom::set_u_iso)
        .def_property("u_aniso_orth", &MAtom::u_aniso_orth, &MAtom::set_u_aniso_orth)
        .def(
            "b_iso", [](const MAtom &self)
            { return clipper::Util::u2b(self.u_iso()); },
            py::return_value_policy::reference_internal)
        .def(
            "b_iso", [](MAtom &self, const ftype &value)
            { return self.set_u_iso(clipper::Util::b2u(value)); },
            py::arg("b_iso"))
        .def_property("pos", &MAtom::coord_orth, &MAtom::set_coord_orth)
        //{
        // auto coord = self.coord_orth();
        // return std::array<ftype,3>{{coord.x(), coord.y(), coord.z()}}; },
        // py::return_value_policy::reference_internal)
        .def("x", [](const MAtom &self)
             { return self.coord_orth().x(); })
        .def("y", [](const MAtom &self)
             { return self.coord_orth().y(); })
        .def("z", [](const MAtom &self)
             { return self.coord_orth().z(); })
        // transform()
        .def("__repr__", [](const MAtom &self)
             {
        std::stringstream stream;
        auto coord = self.coord_orth();
        stream << "<clipper.MAtom " << self.name().trim();
        stream << " at (" << coord.x() << ", " << coord.y() << ", " << coord.z();
        stream << ")>";
        return stream.str(); })
        .def("copy_from", &MAtom::copy, py::arg("other"), py::arg("mode") = MM::COPY::COPY_C)
        .def_static("id_match", &MAtom::id_match, py::arg("id1"), py::arg("id2"), py::arg("mode"));

    pyResidue
        .def(py::init<>())
        .def_property("id", &MResidue::id, &MResidue::set_id)
        .def_property("type", &MResidue::type, &MResidue::set_type)
        .def_property("seqnum", &MResidue::seqnum, &MResidue::set_seqnum) // py::arg("s"), py::arg("inscode") = "")
        .def("size", &MResidue::size)
        .def("__len__", &MResidue::size)
        .def("__repr__", [](const MResidue &self)
             {
        std::stringstream stream;
        stream << "<clipper.MResidue ";
        stream << self.id().trim() << "(" << self.type() << ") containing ";
        stream << self.size() << " atom(s).>";
        return stream.str(); })
        .def(
            "__getitem__", [](MResidue &self, const int i) -> const MAtom &
            { return self[normalise_index(i, self)]; },
            py::arg("i"), py::return_value_policy::reference_internal)
        .def(
            "__getitem__", [](MResidue &self, const std::string &n) -> const MAtom &
            { return self.find(n); },
            py::arg("n"))
        .def(
            "find", [](MResidue &self, const std::string &n, const MM::MODE mode) -> const MAtom &
            { return self.find(n, mode); },
            py::arg("n"), py::arg("mode") = MM::MODE::UNIQUE, py::return_value_policy::reference_internal)
        .def("__setitem__", [](MResidue &self, const int i, MAtom atm)
             { return self[normalise_index(i, self)] = atm; })
        .def("__setitem__", [](MResidue &self, const std::string &n, MAtom atm)
             { return self.find(n) = atm; })
        // atom_list()
        //  transform(const RTop_orth rt)
        .def("insert", &MResidue::insert, py::arg("add"), py::arg("pos"))
        .def(py::self & py::self)
        .def(py::self | py::self)
        .def("copy_from", &MResidue::copy, py::arg("other"), py::arg("mode") = MM::COPY::COPY_C)
        .def("clone", [](const MResidue &self)
             { return new MResidue(self); })
        .def_static("id_match", &MResidue::id_match, py::arg("id1"), py::arg("id2"), py::arg("mode"))
        // UTILITY
        .def("build_carbonyl_oxygen", (void(MResidue::*)(const MResidue &)) & MResidue::protein_mainchain_build_carbonyl_oxygen, py::arg("next"))
        .def("build_carbonyl_oxygen", (void(MResidue::*)()) & MResidue::protein_mainchain_build_carbonyl_oxygen)
        .def("number_of_rotamers", (int(MResidue::*)(MResidue::TYPE) const) & MResidue::protein_sidechain_number_of_rotamers, py::arg("t") = MResidue::TYPE::Richardson)
        .def("number_of_rotamers", (int(MResidue::*)() const) & MResidue::protein_sidechain_number_of_rotamers)
        .def("build_sidechain_numbered_rotamer", (ftype(MResidue::*)(const int &, MResidue::TYPE)) & MResidue::protein_sidechain_build_rotamer, py::arg("n"), py::arg("t") = MResidue::TYPE::Richardson)
        .def("build_sidechain_numbered_rotamer", (ftype(MResidue::*)(const int &)) & MResidue::protein_sidechain_build_rotamer, py::arg("n"))
        .def_static("protein_peptide_bond", &MResidue::protein_peptide_bond, py::arg("m1"), py::arg("m2"), py::arg("r") = 1.5)
        .def_static("protein_ramachandran_phi", &MResidue::protein_ramachandran_phi, py::arg("m1"), py::arg("m2"))
        .def_static("protein_ramachandran_psi", &MResidue::protein_ramachandran_psi, py::arg("m1"), py::arg("m2"))
        .def_static("default_type", &MResidue::default_type);
    // id_match
    // id_tidy
    // lookup
    //;

    pyChain
        .def(py::init<>())
        .def_property("id", &MChain::id, &MChain::set_id)
        .def("size", &MChain::size)
        .def("__len__", &MChain::size)
        .def("__repr__", [](const MChain &self)
             {
        std::stringstream stream;
        stream << "<clipper.MChain ";
        stream << self.id() << " containing ";
        stream << self.size() << " residue(s).>";
        return stream.str(); })
        .def(
            "__getitem__", [](MChain &self, const int i) -> const MResidue &
            { return self[normalise_index(i, self)]; },
            py::arg("i"), py::return_value_policy::reference_internal)
        .def(
            "__getitem__", [](MChain &self, const std::string &n) -> const MResidue &
            { return self.find(n); },
            py::arg("n"), py::return_value_policy::reference_internal)
        .def(
            "find", [](const MChain &self, const std::string &n, const MM::MODE mode) -> const MResidue &
            { return self.find(n, mode); },
            py::arg("n"), py::arg("mode") = MM::MODE::UNIQUE, py::return_value_policy::reference_internal)
        .def("__setitem__", [](MChain &self, const int i, MResidue res)
             { return self[normalise_index(i, self)] = res; })
        .def("__setitem__", [](MChain &self, const std::string &n, MResidue res)
             { return self.find(n) = res; })
        // atom_list()
        // transform()
        .def(
            "__iter__", [](MChain &self)
            { return py::make_iterator(&self[0], &self[self.size()]); },
            py::keep_alive<0, 1>())
        .def("insert", &MChain::insert, py::arg("add"), py::arg("pos"))
        .def(py::self & py::self)
        .def(py::self | py::self)
        .def("copy_from", &MChain::copy, py::arg("other"), py::arg("mode") = MM::COPY::COPY_C)
        .def("clone", [](const MChain &self)
             { return new MChain(self); });

    pyModel
        .def(py::init<>())
        // atom_list()
        // transform()
        .def("size", &MModel::size)
        .def("__len__", &MModel::size)
        .def("__repr__", [](const MModel &self)
             { return "<clipper.MModel containing " + std::to_string(self.size()) + " chain(s)>"; })
        .def(
            "__getitem__", [](MModel &self, const int i) -> const MChain &
            { return self[normalise_index(i, self)]; },
            py::arg("i"), py::return_value_policy::reference_internal)
        .def(
            "__getitem__", [](MModel &self, const std::string &n) -> const MChain &
            { return self.find(n); },
            py::arg("n"), py::return_value_policy::reference_internal)
        .def(
            "find", [](const MModel &self, const std::string &n, const MM::MODE mode) -> const MChain &
            { return self.find(n, mode); },
            py::arg("n"),
            py::arg("mode") = MM::MODE::UNIQUE, py::return_value_policy::reference_internal)
        .def(
            "__setitem__", [](MModel &self, const int i, MChain chn)
            { return self[normalise_index(i, self)] = chn; })
        .def("__setitem__", [](MModel &self, const std::string &n, MChain chn)
             { return self.find(n) = chn; })
        .def(
            "__iter__", [](MModel &self)
            { return py::make_iterator(&self[0], &self[self.size()]); },
            py::keep_alive<0, 1>())
        .def("insert", &MModel::insert, py::arg("add"), py::arg("pos"))
        .def(py::self & py::self)
        .def(py::self | py::self)
        .def("copy_from", &MModel::copy, py::arg("other"), py::arg("mode") = MM::COPY::COPY_C)
        .def("clone", [](const MModel &self)
             { return new MModel(self); });

    py::class_<MiniMol> minimol(m, "MiniMol");
    minimol
        .def(py::init<>())
        .def(py::init<const Spacegroup &, const Cell &>(), py::arg("spacegroup"), py::arg("cell"))
        .def("init", &MiniMol::init)
        .def("__len__", [](const MiniMol &mmol)
             { return mmol.model().size(); })
        .def("__repr__", [](const MiniMol &self)
             {
        std::stringstream stream;
        stream << "<clipper.MiniMol containing model with ";
        stream << self.model().size() << " chain(s).>";
        return stream.str(); })
        .def_property_readonly("cell", &MiniMol::cell)
        .def_property_readonly("spacegroup", &MiniMol::spacegroup)
        .def(
            "model", [](const MiniMol &self) -> const MModel &
            { return self.model(); },
            py::return_value_policy::reference_internal)
        .def(
            "model", [](MiniMol &self, MModel mol)
            { return self.model() = mol; },
            py::return_value_policy::reference_internal)
        .def("clone", [](const MiniMol &self)
             { return new MiniMol(self); })
        .def_property_readonly("is_null", &MiniMol::is_null);

    //.def_property_readonly("model", py::overload_cast<MModel>(&MiniMol::model, py::const_));
    //.def("model", py::overload_cast<>(&MiniMol::model));
    //.def("model", (MModel)&MiniMol::model);
    // py::keep_alive<0, 1>());
    //.def("model", (std::vector<MChain>(MModel::*)()) & MiniMol::model, py::keep_alive<0, 1>());

    m.def(
        "read_structure", [](const std::string &path, bool enable_messages)
        {
        if (path == "undefined")
        {
            throw std::invalid_argument("No path provided for input model! Aborting...");
        }
        MiniMol mmol;
        BuccaneerUtil::read_model(mmol, path, enable_messages);
        // MiniMol *pymmol = new MiniMol(mmol);
        return mmol; },
        // need to see how to read in spacegroup/cell
        // maybe should update clipper to exchange with gemmi
        // return std::unique_ptr<MiniMol>(new MiniMol(mmol)); }, // pymmol; },
        py::arg("path"),
        py::arg("enable_user_messages") = true, "Reads a coordinate file into MiniMol");
    m.def(
        "read_structure", [](const std::string &path, MiniMol &mmol, bool enable_messages)
        {
        if (path == "undefined")
        {
            throw std::invalid_argument("No path provided for input model! Aborting...");
        }
        BuccaneerUtil::read_model(mmol, path, enable_messages);
        // MiniMol *pymmol = new MiniMol(mmol);
        return mmol; },
        // need to see how to read in spacegroup/cell
        // maybe should update clipper to exchange with gemmi
        // return std::unique_ptr<MiniMol>(new MiniMol(mmol)); }, // pymmol; },
        py::arg("path"),
        py::arg("minimol"),
        py::arg("enable_user_messages") = true, "Reads a coordinate file into MiniMol");
    m.def(
        "write_pdb", [](const std::string &path, MiniMol &mmol)
        {
        clipper::MMDBfile mmdb;
        mmdb.export_minimol(mmol);
        mmdb.write_file(path); },
        py::arg("path"), py::arg("minimol"));
    // py::class_<bk::PyCMiniMol>(m, "PyCMiniMol")
    //     .def(py::init<std::string &, bool>(), py::arg("filepath_to_structure") = "undefined", py::arg("enable_messages") = true) //;
    //     .def("get_mmol", &bk::PyCMiniMol::get_mmol, py::arg("filepath_to_structure") = "undefined");
}