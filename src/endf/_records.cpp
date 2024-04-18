#include <cstdlib>
#include <cstring> // for strlen

#include <pybind11/pybind11.h>

//! Convert string representation of a floating point number into a double
//
//! This function handles converting floating point numbers from an ENDF 11
//! character field into a double, covering all the corner cases. Floating point
//! numbers are allowed to contain whitespace (which is ignored). Also, in
//! exponential notation, it allows the 'e' to be omitted. A field containing
//! only whitespace is to be interpreted as a zero.
//
//! \param buffer character input from an ENDF file
//! \return Floating point number

double cfloat_endf(const char* buffer)
{
  char arr[13]; // 11 characters plus e and a null terminator
  int j = 0; // current position in arr
  int found_significand = 0;
  int found_exponent = 0;

  // limit n to 11 characters
  int n = std::strlen(buffer);
  if (n > 11) n = 11;

  for (int i = 0; i < n; ++i) {
    char c = buffer[i];

    // Skip whitespace characters
    if (c == ' ') continue;
    if (found_significand) {
      if (!found_exponent) {
        if (c == '+' || c == '-') {
          // In the case that we encounter +/- and we haven't yet encountered
          // e/E, we manually add it
          arr[j++] = 'e';
          found_exponent = 1;

        } else if (c == 'e' || c == 'E' || c == 'd' || c == 'D') {
          arr[j++] = 'e';
          found_exponent = 1;
          continue;
        }
      }
    } else if (c == '.' || (c >= '0' && c <= '9')) {
      found_significand = 1;
    }

    // Copy character
    arr[j++] = c;
  }

  // Done copying. Add null terminator and convert to double
  arr[j] = '\0';
  return std::atof(arr);
}

PYBIND11_MODULE(_records, m) {
  m.doc() = "float_endf";
  m.def("float_endf", &cfloat_endf, "Convert string to float");
}
