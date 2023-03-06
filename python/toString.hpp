// #ifndef BUILDKIT_TOSTRING_HPP_
// #define BUILDKIT_TOSTRING_HPP_
// #include <sstream>
// #include <string>
//
// template <typename T, typename... Args>
// inline std::ostream &operator<<(std::ostream &stream, T &&value, Args &&...args)
//{
//   std::ostringstream out;
//   out << std::forward<T>(value);
//   os << out.str();
//   return os;
// }
//
// #endif