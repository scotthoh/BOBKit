#include "buildkit-minimol.h"
#include "type_conversions.h"
namespace py = pybind11;
using namespace clipper;
// namespace bk = buildkit;
//  23 Feb 23
//  Need to wrap Spacegroup, Cell, typecast String
//  then whole of minimol.h? for iterables
// refer to gemmi mol.cpp to wrap the add_child/set_child etc

template <typename T>
int normalize_index(int index, const T &container)
{
    if (index < 0)
        index += (int)container.size();
    if ((size_t)index >= container.size())
        throw pybind11::index_error();
    return index;
}
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
    return children[normalize_index(index, children)];
}

template <typename P, typename C>
void set_child(P &parent, int index, C &child)
{
    auto &children = parent.children();
    children[normalize_index(index, children)] = child;
}

template <typename P>
void remove_child(P &parent, int index)
{
    auto &children = parent.children();
    children.erase(children.begin() + normalize_index(index, children));
}

template <typename P>
void remove_children(P &parent, py::slice slice)
{
    delitem_slice(parent.children(), slice);
}

void init_minimol(py::module &m)
{
    py::enum_<MM::MODE>(m, "MODE")
        .value("UNIQUE", MM::MODE::UNIQUE)
        .value("ANY", MM::MODE::ANY)
        .export_values();

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
        .value("CHILDREN", MM::COPY::CHILDREN)
        .export_values();

    py::class_<MChain>(m, "MChain")
        .def(py::init<>())
        // clipper string convert
        .def_property("id", &MChain::id, &MChain::set_id)
        .def("size", &MChain::size)
        .def("__len__", &MChain::size)
        .def("__repr__", [](const MChain &self)
             {
            std::stringstream stream;
            stream << "<clipper.MChain, Chain ";
            stream << self.id() << " containing ";
            stream << self.size() << " residues.>";
            return stream.str(); });

    py::class_<MModel>(m, "MModel")
        .def(py::init<>())
        .def("size", &MModel::size)
        .def("__len__", &MModel::size)
        .def("__repr__", [](const MModel &self)
             { return "<clipper.MModel containing " + std::to_string(self.size()) + " chain(s)>"; })
        .def(
            "__getitem__", [](MModel &self, const int index) -> const MChain &
            { return self[normalize_index(index, self)]; },
            py::arg("index"), py::return_value_policy::reference_internal)
        .def(
            "__getitem__", [](MModel &self, const std::string &name) -> const MChain &
            { return self.find(name); },
            py::arg("name"), py::return_value_policy::reference_internal)
        .def(
            "find", [](const MModel &self, const std::string &id, const MM::MODE mode) -> const MChain &
            { return self.find(id, mode); },
            py::arg("id"),
            py::arg("mode") = MM::MODE::UNIQUE, py::return_value_policy::reference_internal)
        .def(
            "__setitem__", [](MModel &self, const int index, MChain chn)
            { return self[normalize_index(index, self)] = chn; })
        .def("__setitem__", [](MModel &self, const std::string &name, MChain chn)
             { return self.find(name) = chn; })
        .def("__iter__", [](MModel &self)
             { return py::make_iterator(&self[0], &self[self.size()]); })
        .def("insert", &MModel::insert, py::arg("chain"), py::arg("pos"));

    py::class_<MiniMol>(m, "MiniMol")
        .def(py::init<>())
        .def(py::init<const Spacegroup &, const Cell &>(), py::arg("spacegroup"), py::arg("cell"))
        .def("init", &MiniMol::init)
        .def("__len__", [](const MiniMol &mmol)
             { return mmol.model().size(); })
        .def_property_readonly("cell", &MiniMol::cell)
        .def_property_readonly("spacegroup", &MiniMol::spacegroup)
        .def(
            "model", [](const MiniMol &mmol) -> const MModel &
            { return mmol.model(); },
            py::return_value_policy::reference_internal)
        .def(
            "model", [](MiniMol &mmol, MModel model)
            { return mmol.model() = model; },
            py::return_value_policy::reference_internal)
        //.def("__repr__", [](const MiniMol &self)
        //     {
        //        std::stringstream stream;
        //        stream << "<clipper.MiniMol containing model with ";
        //        stream << self.model().size() << " chain(s).>";
        //        return stream.str(); })
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
            String clipper_str_path = path;
            BuccaneerUtil::read_model(mmol, clipper_str_path, enable_messages);
            //MiniMol *pymmol = new MiniMol(mmol);
            return std::unique_ptr<MiniMol>(new MiniMol(mmol)); }, // pymmol; },
        py::arg("path"), py::arg("enable_user_messages") = true, "Reads a coordinate file into MiniMol");

    // py::class_<bk::PyCMiniMol>(m, "PyCMiniMol")
    //     .def(py::init<std::string &, bool>(), py::arg("filepath_to_structure") = "undefined", py::arg("enable_messages") = true) //;
    //     .def("get_mmol", &bk::PyCMiniMol::get_mmol, py::arg("filepath_to_structure") = "undefined");
}