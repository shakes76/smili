/*=========================================================================
  The Software is copyright (c) Commonwealth Scientific and Industrial Research Organisation (CSIRO)
  ABN 41 687 119 230.
  All rights reserved.

  Licensed under the CSIRO BSD 3-Clause License
  You may not use this file except in compliance with the License.
  You may obtain a copy of the License in the file LICENSE.md or at

  https://stash.csiro.au/projects/SMILI/repos/smili/browse/license.txt

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
=========================================================================*/
#ifndef __MILXGLOBAL_H
#define __MILXGLOBAL_H

//STL
#include <limits>
#include <stdexcept> //Include stdexcept to handle conversion exceptions
#include <sstream> //includes string
#include <iomanip>
#include <iostream>

#ifndef VTK_ONLY
  //VNL
  #include <vnl/vnl_vector_fixed.h>
  //ITK
  #include <itkCommand.h>
#endif
#ifndef ITK_ONLY
  //VTK
  #include <vtkType.h>
  #include <vtkVersion.h>
  #include <vtkSmartPointer.h>
  #include <vtkCommand.h>
  #include <vtkObject.h>
#endif

#ifdef WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace milx
{

/**
 * \defgroup Global Global Methods
 * \brief Common global functions and macros
 */
//@{
//Global Constants
const float Version = static_cast<float>(2.00);
static bool VerboseMode = true;

//Global Typedefs
typedef double coordinateType; //compatible with >= VTK 5
#ifndef VTK_ONLY
  typedef unsigned vnl_size_t;
  typedef vnl_vector_fixed<coordinateType,3> coordinate;
#endif
#ifdef VTK_ONLY
  typedef coordinateType* coordinate;
#endif
const unsigned imgDimension = 3;

//C++ Helper Macros
// Swap two values
template <typename T>
inline void Swap(T &a, T&b) {
  T t=a;
  a=b;
  b=t;
}

// Calculate the minimum of two values
template <typename T>
inline T Minimum(const T &a, const T &b) {
  return (a < b) ? a : b;
}

// Calculate the maximum of two values
template <typename T>
inline T Maximum(const T &a, const T &b) {
  return (a > b) ? a : b;
}

//Global Functions
/**
* \class BadStringConversion
* \brief Exception Handler for stringify function. Thrown if the conversion from
* numeric to string is unsuccessful.
*/
class BadStringConversion : public std::runtime_error
{
public:
  BadStringConversion(const std::string& s)
      : runtime_error(s)
  { }
};

/**
* \fn NumberToString(double num, unsigned zeroPad = 0)
* \brief Number to string converter
*/
inline std::string NumberToString(double num, unsigned zeroPad = 0)
{
  std::ostringstream toString;
  if (!(toString << std::setw(zeroPad) << std::setfill('0') << num))
    throw BadStringConversion("NumberToString(double)");
  return toString.str();
}

/**
* \fn NumberOfProcessors()
* \brief Number of processors or cores on the machine.
*/
inline unsigned NumberOfProcessors()
{
  unsigned cores = 1;

#ifdef _WIN32
  SYSTEM_INFO sysinfo;
  GetSystemInfo( &sysinfo );

  cores = sysinfo.dwNumberOfProcessors;
#else
  cores = sysconf( _SC_NPROCESSORS_ONLN );
#endif

  return cores;
}

/**
* \fn PrintProgressBar( const size_t percent )
* \brief Displays a generic terminal-based progress bar based on the percent provided

 Usage for counter i would be:
 \code
 milx::PrintProgressBar(static_cast<float>(i)/numberOfPoints*100);
 \endcode
*/
inline void PrintProgressBar( const size_t percent )
{
  std::string bar;

  for (size_t i = 0; i < 50; i++)
  {
    if ( i < (percent/2))
      bar.replace(i,1,"=");
    else if ( i == (percent/2))
      bar.replace(i,1,">");
    else
      bar.replace(i,1," ");
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
  if (percent == 100)
    std::cout << std::endl;
}

/**
* \fn PrintInfo(const std::string msg)
* \brief Displays a generic msg to standard output with carriage return
*/
inline void PrintInfo(const std::string msg)
{
  std::cout << msg << std::endl;
}

/**
* \fn PrintWarning(const std::string msg)
* \brief Displays a generic msg to standard output with carriage return
*/
inline void PrintWarning(const std::string msg)
{
  std::cout << "Warning: " << msg << std::endl;
}

/**
* \fn PrintDebug(const std::string msg)
* \brief Displays a generic msg to standard error with carriage return if in Debug mode
*/
inline void PrintDebug(const std::string msg)
{
  if(VerboseMode)
    std::cerr << "Debug: " << msg << std::endl;
}

/**
* \fn PrintError(const std::string msg)
* \brief Displays a generic msg to standard error with carriage return
*/
inline void PrintError(const std::string msg)
{
  std::cerr << msg << std::endl;
}

#if defined(WIN32)
#if defined(SMILI_DLL)
/**
Used for DLL generation purposes (Windows Specific) Import/Export.
Templates classes cannot be imported hence its own variable.
The export command is
*/
#if defined(SMILI_MAKEDLL)     // create a SMILI DLL library
#define SMILI_EXPORT  __declspec(dllexport)
#else                        // use a SMILI DLL library
#define SMILI_EXPORT  __declspec(dllimport)
#endif
#endif // SMILI_DLL
#endif // WIN32

#ifndef SMILI_EXPORT
/**
    \def SMILI_EXPORT
    \brief DLL Function Symbol for Windows. It is empty for other OSes.
*/
#define SMILI_EXPORT
#endif
//@}
}

/**
 * \defgroup Observers Event Propagators
 * \brief Event propagators for SMILI classes milx::Image and milx::Model.

 * This allows the progress events to be propagated to Qt or other GUIs.
 */
//@{
#ifndef ITK_ONLY
  #ifndef VTK_ONLY //Requires VTK and ITK
  namespace itk
  {
    /**
      \class ProgressUpdates
      \brief Object for intercepting progress events triggered by ITK based filters and converting them to VTK progress events.

      Since the milx::Image is unaware of milxQtImage, this Command object is added as an observer of ITK filters.
      Upon a progress event, the Execute member of this class is triggered, which in turn invokes a VTK progress event of
      a dummy internal object. If this internal object is linked as an VTK observer, the event is effectively passed on.

      Link the internal oberver:
      \code
      linkProgressEventOf(milx::ProgressUpdates->GetUpdateObject()); //link itk observer propagator
      \endcode
      Connect the internal observer usig the linkProgressEventOf() member.
    */
    class SMILI_EXPORT ProgressUpdates : public Command
    {
    public:
      /**
          \fn ProgressUpdates::New()
          \brief Standard ITK Object Factory
      */
      static ProgressUpdates* New()
      {
        return new ProgressUpdates();
      }
      itkTypeMacro(ProgressUpdates, Command);

      inline vtkObject* GetUpdateObject()
      { return Filter.GetPointer(); }

    protected:
      ProgressUpdates()
      { Filter = vtkSmartPointer<vtkObject>::New(); }
      ~ProgressUpdates() {}

      vtkSmartPointer<vtkObject> Filter; //!< Dummy object for passing on events

    private:
      inline void Execute(Object *caller, const EventObject & event)
      { Execute( (const Object *)caller, event);  }

      inline void Execute(const Object * object, const EventObject & event)
      { Filter->InvokeEvent(vtkCommand::ProgressEvent); } //pass on event using the dummy object
    };
  }

  namespace milx
  {
    static itk::SmartPointer<itk::ProgressUpdates> ProgressUpdates = itk::ProgressUpdates::New(); //!< Observer to pass on progress events
  }
  #endif
#endif

#ifndef ITK_ONLY //Requires VTK
  /**
    \class vtkProgressUpdates
    \brief Object for intercepting progress events triggered by VTK based filters.

    Since the milx::Model is unaware of milxQtModel, this vtkCommand object is added as an observer of a filter.
    Upon a progress event, the Execute member of this class is triggered, which in turn invokes another progress event of
    a dummy internal object. If this internal object is linked as an observer, the event is effectively passed on.

    Link the internal oberver:
    \code
    linkProgressEventOf(model.GetObserverObject()); //link internal observer
    \endcode
    Connect the internal observer usig the linkProgressEventOf() member.
  */
  class SMILI_EXPORT vtkProgressUpdates : public vtkCommand
  {
  public:
    /**
        \fn vtkProgressUpdates::New()
        \brief Standard VTK Object Factory
    */
    static vtkProgressUpdates* New()
    {
      return new vtkProgressUpdates();
    }
    vtkTypeMacro(vtkProgressUpdates, vtkCommand);

    inline vtkObject* GetUpdateObject()
    { return Filter.GetPointer(); }

    inline bool IsBeingObserved(unsigned long event = vtkCommand::ProgressEvent)
    { return Filter->HasObserver(event);  }

  protected:
    vtkProgressUpdates()
    { Filter = vtkSmartPointer<vtkObject>::New(); }
    ~vtkProgressUpdates() {}

    vtkSmartPointer<vtkObject> Filter; //!< Dummy object for passing on events

  private:
    inline void Execute(vtkObject *caller, unsigned long observedType, void* message)
    { Filter->InvokeEvent(vtkCommand::ProgressEvent); } //pass on event using the dummy object
  };

  namespace milx
  {
    static vtkSmartPointer<vtkProgressUpdates> VTKProgressUpdates = vtkSmartPointer<vtkProgressUpdates>::New(); //!< Observer to pass on progress events
  }

  /**
  \class vtkErrorWarning
  \brief Object for intercepting errors thrown by VTK based readers.
  */
  class SMILI_EXPORT vtkErrorWarning : public vtkCommand
  {
  public:
    /**
        \fn vtkErrorWarning::New()
        \brief Standard VTK Object Factory
    */
    static vtkErrorWarning* New()
    { return new vtkErrorWarning(); }

    vtkTypeMacro(vtkErrorWarning, vtkCommand);

    inline const char* GetMessage()
    { return m_Message.c_str(); }

    inline bool HasFailed()
    { return m_ErrorEncountered;  }

    inline bool ReportsFailure()
    { return m_ErrorEncountered;  }

    inline void Initialize()
    { m_ErrorEncountered = m_WarningEncountered = false;  }

  protected:
    vtkErrorWarning()
    { Initialize(); }
    ~vtkErrorWarning() {}

    bool m_ErrorEncountered;
    bool m_WarningEncountered;
    std::string m_Message;

  private:
    void Execute(vtkObject *caller, unsigned long observedType, void* message)
    {
      m_Message = (const char*) message;

      if(observedType == vtkCommand::ErrorEvent)
          m_ErrorEncountered = true;
      if(observedType == vtkCommand::WarningEvent)
          m_WarningEncountered = true;
    }
  };
#endif
//@}

#endif //__MILXGLOBAL_H
