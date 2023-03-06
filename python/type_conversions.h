#pragma once
#include <pybind11/pybind11.h>
#include <clipper/clipper.h>

// refer to https://pybind11.readthedocs.io/en/stable/advanced/cast/custom.html#custom-type-casters
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
        return PyUnicode_FromString(src.c_str());
      }
    };
  }
}