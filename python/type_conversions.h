#pragma once
#include <clipper/clipper.h>
#include <pybind11/pybind11.h>

// refer to
// https://pybind11.readthedocs.io/en/stable/advanced/cast/custom.html#custom-type-casters
namespace PYBIND11_NAMESPACE
{
  namespace detail
  {
    template <>
    struct type_caster<clipper::String>
    {
    public:
      PYBIND11_TYPE_CASTER(clipper::String, _("String"));
      // Python -> C++
      bool load(handle src, bool convert)
      {
        PyObject *source = src.ptr();
        PyObject *tmp = PyUnicode_AsUTF8String(source);
        if (!tmp)
          return false;
        char *buffer;
        ssize_t slength;
        if (PYBIND11_BYTES_AS_STRING_AND_SIZE(tmp, &buffer, &slength))
          return false;
        value = clipper::String(std::string(buffer, (size_t)slength));
        Py_DECREF(tmp);
        return true;
      }
      // C++ -> Python
      static handle cast(const clipper::String &src, return_value_policy policy, handle parent)
      {
        return PyUnicode_FromString(src.trim().c_str());
      }
    };

    // for gemmi::Structure C++<-->Python
    // template <> struct type_caster<gemmi::Structure> {
    // public:
    //  /**
    //   * This macro establishes the name 'gemmi::Structure' in
    //   * function signatures and declares a local variable
    //   * 'value' of type gemmi::Structure
    //   */
    //  PYBIND11_TYPE_CASTER(gemmi::Structure, _("Structure"));
    //
    //  /**
    //   * Conversion part 1 (Python->C++): convert a PyObject into a
    //   *gemmi::Structure. The second argument indicates whether implicit
    //   *conversions should be applied.
    //   */
    //  bool load(handle src, bool) {
    //    // PyObject *py_object = src.ptr();
    //    gemmi::Structure s = src.cast<gemmi::Structure>();
    //    if (s)
    //    // PyObject *tmp = src.cast<gemmi::Structure *>();
    //    //  std::cout << "python to c++ conversion" << std::endl; // debug
    //    code
    //    //
    //    //  gemmi::Structure gstructure;
    //    //  value = gstructure;
    //    value = s;
    //    return true;
    //  }
    //
    //  /**
    //   * Conversion part 2 (C++ -> Python): convert an gemmi::Structure
    //   instance
    //   * into a Python object. The second and third arguments are used to
    //   * indicate the return value policy and parent object (for
    //   * ``return_value_policy::reference_internal``) and are generally
    //   * ignored by implicit casters.
    //   */
    //  static handle cast(gemmi::Structure src, return_value_policy /* policy
    //  */,
    //                     handle /* parent */) {
    //
    //    std::cout << "C++ to python conversion" << std::endl; // debug mode
    //    return pybind11::bool_(true).release();
    //  }
    //};
  }
}