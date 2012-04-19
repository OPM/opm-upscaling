// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
// $Id$

#ifndef OPM_EXCEPTIONS__HH
#define OPM_EXCEPTIONS__HH

#include <string>
#include <sstream>

namespace Opm {

/*! \defgroup Exceptions Exception handling
  \ingroup Common
  \{

  The Dune-exceptions are designed to allow a simple derivation of subclasses
  and to accept a text written in the '<<' syntax.

  Example of usage:

\code
#include "exceptions_.hh"

...

class FileNotFoundError : public Opm::IOError {};

...

void fileopen (std::string name) {
  std::ifstream file;

  file.open(name.c_str());

  if (file == 0)
    OPM_THROW(FileNotFoundError, "File " << name << " not found!");

  ...

  file.close();
}

...

int main () {
  try {
    ...
  } catch (Opm::IOError &e) {
    std::cerr << "I/O error: " << e << std::endl;
    return 1;
  } catch (Opm::Exception &e) {
    std::cerr << "Generic Dune error: " << e << std::endl;
    return 2;  
  }
}
\endcode

  \see exceptions.hh for detailed info

*/

/*! \file
  \brief A few common exception classes
  
This file defines a common framework for generating exception
subclasses and to throw them in a simple manner

*/

/* forward declarations */
class Exception;

/*! \class Exception
  \brief Base class for Dune-Exceptions
  
  all Dune exceptions are derived from this class via trivial subclassing:

\code
  class MyException : public Opm::Exception {};
\endcode

  You should not \c throw a Opm::Exception directly but use the macro
  OPM_THROW() instead which fills the message-buffer of the exception
  in a standard way and features a way to pass the result in the
  operator<<-style

  \see OPM_THROW, IOError, MathError

*/
class Exception {
public:
  Exception ();
  void message(const std::string &message); //!< store string in internal message buffer
  const std::string& what() const;          //!< output internal message buffer
private:
  std::string _message;
};

/*
  Implementation of Opm::Exception
 */

inline Exception::Exception ()
{
}


inline void Exception::message(const std::string &message)
{
  _message = message;
}

inline const std::string& Exception::what() const
{
  return _message;
}

inline std::ostream& operator<<(std::ostream &stream, const Exception &e)
{
  return stream << e.what();
}

#ifndef DOXYGEN
// the "format" the exception-type gets printed.  __FILE__ and
// __LINE__ are standard C-defines, the GNU cpp-infofile claims that
// C99 defines __func__ as well. __FUNCTION__ is a GNU-extension
#define THROWSPEC(E) #E << " [" << __func__ << ":" << __FILE__ << ":" << __LINE__ << "]: "
#endif // DOXYGEN

/*! Macro to throw an exception

  \code
#include "exceptions_.hh"
  \endcode

  \param E exception class derived from Opm::Exception
  \param m reason for this exception in ostream-notation

  Example:

  \code
  if (filehandle == 0)
    OPM_THROW(FileError, "Could not open " << filename << " for reading!");
  \endcode

  OPM_THROW automatically adds information about the exception thrown
  to the text.

  \note
  you can add a hook to be called before a Opm::Exception is emitted,
  e.g. to add additional information to the exception,
  or to invoke a debugger during parallel debugging. (see Opm::ExceptionHook)
 
 */
// this is the magic: use the usual do { ... } while (0) trick, create
// the full message via a string stream and throw the created object
#define OPM_THROW(E, m) do { E th__ex; std::ostringstream th__out; \
 th__out << THROWSPEC(E) << m; th__ex.message(th__out.str()); throw th__ex; \
 } while (0)

/*! \brief Default exception class for I/O errors

  This is a superclass for any errors dealing with file/socket I/O problems
  like

   - file not found
   - could not write file
   - could not connect to remote socket
 */
class IOError : public Exception {};

/*! \brief Default exception class for mathematical errors

  This is the superclass for all errors which are caused by
  mathematical problems like

   - matrix not invertible
   - not convergent
 */
class MathError : public Exception {};

/*! \brief Default exception class for range errors

  This is the superclass for all errors which are caused because
  the user tries to access data that was not allocated before.
  These can be problems like

   - accessing array entries behind the last entry
   - adding the fourth non zero entry in a sparse matrix
     with only three non zero entries per row
  
 */
class RangeError : public Exception {};

/*! \brief Default exception for dummy implementations 

  This exception can be used for functions/methods

  - that have to be implemented but should never be called
  - that are missing
 */
class NotImplemented : public Exception {};

/*! \brief Default exception class for OS errors

  This class is thrown when a system-call is used and returns an
  error.
  
 */
class SystemError : public Exception {};

/*! \brief Default exception if memory allocation fails

 */
class OutOfMemoryError : public SystemError {};

/*! \brief Default exception if a function was called while
  the object is not in a valid state for that function.
 */
class InvalidStateException : public Exception {};

/*! \brief Default exception if an error in the parallel
  communication of the programm occured
  \ingroup ParallelCommunication
 */
class ParallelError : public Exception {};

} // end namespace

#endif
