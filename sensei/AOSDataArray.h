/*=========================================================================

  Program:   Visualization Toolkit
  Module:    AOSDataArray.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   AOSDataArray
 * @brief   Array-Of-Structs implementation of
 * svtkGenericDataArray.
 *
 *
 * svtkGenericDataArray specialization that stores data array in the traditional
 * SVTK memory layout where a 3 component is stored in contiguous memory as
 * \c A1A2A3B1B2B3C1C2C3 ... where A,B,C,... are tuples.
 *
 * This replaces svtkDataArrayTemplate.
 *
 * @sa
 * svtkGenericDataArray svtkSOADataArrayTemplate
 */

#ifndef AOSDataArray_h
#define AOSDataArray_h

#include "svtkBuffer.h"           // For storage buffer.
#include "hamr_buffer.h"
#include "svtkCommonCoreModule.h" // For export macro
#include "svtkAOSDataArrayTemplate.h"


namespace sensei
{

struct DataArray
{
    using allocator = hamr::buffer_allocator;
};



// The export macro below makes no sense, but is necessary for older compilers
// when we export instantiations of this class from svtkCommonCore.
template <class T>
class SVTKCOMMONCORE_EXPORT AOSDataArray : public svtkAOSDataArrayTemplate<T>
{
  typedef svtkGenericDataArray<AOSDataArray<T>, T> GenericDataArrayType;
public:

  using allocator = sensei::DataArray::allocator;

  typedef AOSDataArray<T> SelfType;
  //svtkTemplateTypeMacro(SelfType, GenericDataArrayType);
  using Superclass = svtkAOSDataArrayTemplate<T>;
  typedef typename Superclass::ValueType ValueType;

  /// construct an empty buffer that will use the passed allocator type
  static AOSDataArray *New(allocator alloc)
  {
      return new AOSDataArray(alloc);
  }

  /// construct a buffer with n_elem size using the passed allocator type
  static AOSDataArray *New(allocator alloc, size_t n_elem);/*
  {
      return new AOSDataArray(alloc, n_elem);
  }*/

  /** construct a buffer with n_elem size initialized to the passed value
   * using the passed allocator type
   */
  static AOSDataArray *New(allocator alloc, size_t n_elem, const T &val);/*
  {
      return new AOSDataArray(alloc, n_elem, val);
  }*/

  /** construct a buffer with n_elem size initialized to the passed value
   * using the passed allocator type
   */
  static AOSDataArray *New(allocator alloc, size_t n_elem, const T *vals);/*
  {
      return new AOSDataArray(alloc, n_elem, vals);
  }*/

  /** construct by directly providing the buffer contents. This can be used
   * for zero-copy transfer of data.  One must also name the allocator type
   * and device owning the data.  In addition for new allocations the
   * allocator type and owner are used internally to know how to
   * automatically move data during inter technology transfers.
   *
   * @param[in] alloc a ::buffer_allocator indicating the technology
   *                  backing the pointer
   * @param[in] size  the number of elements in the array pointed to by ptr
   * @param[in] owner the device owning the memory, -1 for CPU. if the
   *                  allocator is a GPU allocator and -1 is passed the
   *                  driver API is used to determine the device that
   *                  allocated the memory.
   * @param[in] ptr   a pointer to the array
   * @param[in] df    a function `void df(void*ptr)` used to delete the array
   *                  when this instance is finished.
   */
  //template <typename delete_func_t>
  //static AOSDataArray *New(allocator alloc, size_t size, int owner, T *ptr, delete_func_t df);

  /** construct by directly providing the buffer contents. This can be used
   * for zero-copy transfer of data.  One must also name the allocator type
   * and device owning the data.  In addition for new allocations the
   * allocator type and owner are used internally to know how to
   * automatically move data during inter technology transfers.
   * The pass ::buffer_allocator is used to create the deleter that will be
   * called when this instance is finished with the memeory. Use this
   * constructor to transfer ownership of the array.
   *
   * @param[in] alloc a ::buffer_allocator indicating the technology
   *                  backing the pointer
   * @param[in] size  the number of elements in the array pointed to by ptr
   * @param[in] owner the device owning the memory, -1 for CPU. if the
   *                  allocator is a GPU allocator and -1 is passed the
   *                  driver API is used to determine the device that
   *                  allocated the memory.
   * @param[in] ptr   a pointer to the array
   */
  //static AOSDataArray *New(allocator alloc, size_t size, int owner, T *ptr);

  /** construct by directly providing the buffer contents. This can be used
   * for zero-copy transfer of data.  One must also name the allocator type
   * and device owning the data.  In addition for new allocations the
   * allocator type and owner are used internally to know how to
   * automatically move data during inter technology transfers.
   *
   * @param[in] alloc a ::buffer_allocator indicating the technology
   *                  backing the pointer
   * @param[in] size  the number of elements in the array pointed to by ptr
   * @param[in] owner the device owning the memory, -1 for CPU. if the
   *                  allocator is a GPU allocator and -1 is passed the
   *                  driver API is used to determine the device that
   *                  allocated the memory.
   * @param[in] data  a shared pointer managing the data
   */
  //static AOSDataArray *New(allocator alloc, size_t size, int owner,
  //    const std::shared_ptr<T> &data);

#if !defined(SWIG)
    /** @name get_cpu_accessible
     * Returns a pointer to the contents of the buffer accessible on the CPU.
     * If the buffer is currently accessible by codes running on the CPU then
     * this call is a NOOP.  If the buffer is not currently accessible by codes
     * running on the CPU then a temporary buffer is allocated and the data is
     * moved to the CPU.  The returned shared_ptr deals with deallocation of
     * the temporary if needed.
     */
    ///@{
    /// returns a pointer to the contents of the buffer accessible on the CPU.
    //std::shared_ptr<T> get_cpu_accessible();

    /// returns a pointer to the contents of the buffer accessible on the CPU.
    //std::shared_ptr<const T> get_cpu_accessible() const;
    ///@}
#endif

    /// returns true if the data is accessible from codes running on the CPU
    //int cpu_accessible() const;

#if !defined(SWIG)
    /** @name get_cuda_accessible
     * returns a pointer to the contents of the buffer accessible from the
     * active CUDA device.  If the buffer is currently accessible on the named
     * CUDA device then this call is a NOOP.  If the buffer is not currently
     * accessible on the named CUDA device then a temporary buffer is allocated
     * and the data is moved.  The returned shared_ptr deals with deallocation
     * of the temporary if needed.
     */
    ///@{
    ///  returns a pointer to the contents of the buffer accessible from within CUDA
    //std::shared_ptr<T> get_cuda_accessible();

    ///  returns a pointer to the contents of the buffer accessible from within CUDA
    //std::shared_ptr<const T> get_cuda_accessible() const;
    ///@}
#endif

    /// returns true if the data is accessible from CUDA codes
    //int cuda_accessible() const;

#if !defined(SWIG)
    /** @name get_hip_accessible
     * returns a pointer to the contents of the buffer accessible from the
     * active HIP device.  If the buffer is currently accessible on the named
     * HIP device then this call is a NOOP.  If the buffer is not currently
     * accessible on the named HIP device then a temporary buffer is allocated
     * and the data is moved.  The returned shared_ptr deals with deallocation
     * of the temporary if needed.
     */
    ///@{
    ///  returns a pointer to the contents of the buffer accessible from within HIP
    //std::shared_ptr<T> get_hip_accessible();

    ///  returns a pointer to the contents of the buffer accessible from within HIP
    //std::shared_ptr<const T> get_hip_accessible() const;
    ///@}
#endif

    /// returns true if the data is accessible from HIP codes
    //int hip_accessible() const;

#if !defined(SWIG)
    /** @name get_openmp_accessible returns a pointer to the contents of
     * the buffer accessible from the active OpenMP off load device.  If the
     * buffer is currently accessible on the named OpenMP off load device then
     * this call is a NOOP.  If the buffer is not currently accessible on the
     * named OpenMP off load device then a temporary buffer is allocated and
     * the data is moved.  The returned shared_ptr deals with deallocation of
     * the temporary if needed.
     */
    ///@{
    /** returns a pointer to the contents of the buffer accessible from within
     * OpenMP off load
     */
    //std::shared_ptr<T> get_openmp_accessible();

    /** returns a pointer to the contents of the buffer accessible from within
     * OpenMP off load
     */
    //std::shared_ptr<const T> get_openmp_accessible() const;
    ///@}
#endif

    /// returns true if the data is accessible from OpenMP off load codes
    ///int openmp_accessible() const;

    /** @name data
     * return the raw pointer to the buffer contents. Use this when you know
     * that the buffer contents are accessible by the code operating on them to
     * save the cost of a std::shared_ptr copy construct.
     */
    ///@{
    /// return a pointer to the buffer contents
    T *data() { return this->Buffer->get(); }

    /// return a const pointer to the buffer contents
    const T *data() const { return this->Buffer->get(); }
    ///@}

protected:
  AOSDataArray(allocator alloc);
  ~AOSDataArray() override;

  /**
   * Allocate space for numTuples. Old data is not preserved. If numTuples == 0,
   * all data is freed.
   */
  bool AllocateTuples(svtkIdType numTuples);

  /**
   * Allocate space for numTuples. Old data is preserved. If numTuples == 0,
   * all data is freed.
   */
  bool ReallocateTuples(svtkIdType numTuples);

  hamr::buffer<ValueType> *Buffer;
  //svtkBuffer<ValueType>* Buffer;

private:
  AOSDataArray(const AOSDataArray&) = delete;
  void operator=(const AOSDataArray&) = delete;

  friend class svtkGenericDataArray<AOSDataArray<T>, T>;
};



// --------------------------------------------------------------------------
template <typename T>
AOSDataArray<T> *AOSDataArray<T>::New(allocator alloc)
{
    return new AOSDataArray<T>(alloc);
}

// --------------------------------------------------------------------------
template <typename T>
AOSDataArray<T> *AOSDataArray<T>::New(allocator alloc, size_t n_elem)
{
    return new AOSDataArray<T>(alloc, n_elem);
}

// --------------------------------------------------------------------------
template <typename T>
AOSDataArray<T> *AOSDataArray<T>::New(allocator alloc, size_t n_elem, const T &val)
{
    return new AOSDataArray<T>(alloc, n_elem, val);
}

// --------------------------------------------------------------------------
template <typename T>
AOSDataArray<T> *AOSDataArray<T>::New(allocator alloc, size_t n_elem, const T *vals)
{
    return new AOSDataArray<T>(alloc, n_elem, vals);
}

// --------------------------------------------------------------------------
template <typename delete_func_t>
template <typename T>
AOSDataArray<T> *AOSDataArray<T>::New(allocator alloc, size_t size, int owner, T *ptr, delete_func_t df);

// --------------------------------------------------------------------------
template <typename T>
AOSDataArray<T> *AOSDataArray<T>::New(allocator alloc, size_t size, int owner, T *ptr);

// --------------------------------------------------------------------------
template <typename T>
AOSDataArray<T> *AOSDataArray<T>::New(allocator alloc, size_t size, int owner,
    const std::shared_ptr<T> &data);












}

#endif // header guard

/*
#if !defined(SWIG)
// Declare svtkArrayDownCast implementations for AoS containers:
svtkArrayDownCast_TemplateFastCastMacro(AOSDataArray);

// This macro is used by the subclasses to create dummy
// declarations for these functions such that the wrapper
// can see them. The wrappers ignore AOSDataArray.
#define svtkCreateWrappedArrayInterface(T)                                                          \
  int GetDataType() const override;                                                                \
  void GetTypedTuple(svtkIdType i, T* tuple) SVTK_EXPECTS(0 <= i && i < GetNumberOfTuples());        \
  void SetTypedTuple(svtkIdType i, const T* tuple) SVTK_EXPECTS(0 <= i && i < GetNumberOfTuples());  \
  void InsertTypedTuple(svtkIdType i, const T* tuple) SVTK_EXPECTS(0 <= i);                          \
  svtkIdType InsertNextTypedTuple(const T* tuple);                                                  \
  T GetValue(svtkIdType id) const SVTK_EXPECTS(0 <= id && id < GetNumberOfValues());                 \
  void SetValue(svtkIdType id, T value) SVTK_EXPECTS(0 <= id && id < GetNumberOfValues());           \
  bool SetNumberOfValues(svtkIdType number) override;                                               \
  void InsertValue(svtkIdType id, T f) SVTK_EXPECTS(0 <= id);                                        \
  svtkIdType InsertNextValue(T f);                                                                  \
  T* GetValueRange(int comp) SVTK_SIZEHINT(2);                                                      \
  T* GetValueRange() SVTK_SIZEHINT(2);                                                              \
  T* WritePointer(svtkIdType id, svtkIdType number);                                                 \
  T* GetPointer(svtkIdType id);                                                                     \
  void SetArray(SVTK_ZEROCOPY T* array, svtkIdType size, int save);                                  \
  void SetArray(SVTK_ZEROCOPY T* array, svtkIdType size, int save, int deleteMethod)

#endif

#if !defined(SWIG)
// This portion must be OUTSIDE the include blockers. This is used to tell
// libraries other than svtkCommonCore that instantiations of
// AOSDataArray can be found externally. This prevents each library
// from instantiating these on their own.
#ifdef SVTK_AOS_DATA_ARRAY_TEMPLATE_INSTANTIATING
#define SVTK_AOS_DATA_ARRAY_TEMPLATE_INSTANTIATE(T)                                                 \
  template class SVTKCOMMONCORE_EXPORT AOSDataArray<T>
#elif defined(SVTK_USE_EXTERN_TEMPLATE)
#ifndef SVTK_AOS_DATA_ARRAY_TEMPLATE_EXTERN
#define SVTK_AOS_DATA_ARRAY_TEMPLATE_EXTERN
#ifdef _MSC_VER
#pragma warning(push)
// The following is needed when the AOSDataArray is declared
// dllexport and is used from another class in svtkCommonCore
#pragma warning(disable : 4910) // extern and dllexport incompatible
#endif
svtkExternTemplateMacro(extern template class SVTKCOMMONCORE_EXPORT AOSDataArray);
#ifdef _MSC_VER
#pragma warning(pop)
#endif
#endif // SVTK_AOS_DATA_ARRAY_TEMPLATE_EXTERN

// The following clause is only for MSVC
#elif defined(_MSC_VER) && !defined(SVTK_BUILD_SHARED_LIBS)
#pragma warning(push)

// C4091: 'extern ' : ignored on left of 'int' when no variable is declared
#pragma warning(disable : 4091)

// Compiler-specific extension warning.
#pragma warning(disable : 4231)

// We need to disable warning 4910 and do an extern dllexport
// anyway.  When deriving svtkCharArray and other types from an
// instantiation of this template the compiler does an explicit
// instantiation of the base class.  From outside the svtkCommon
// library we block this using an extern dllimport instantiation.
// For classes inside svtkCommon we should be able to just do an
// extern instantiation, but VS complains about missing
// definitions.  We cannot do an extern dllimport inside svtkCommon
// since the symbols are local to the dll.  An extern dllexport
// seems to be the only way to convince VS to do the right
// thing, so we just disable the warning.
#pragma warning(disable : 4910) // extern and dllexport incompatible

// Use an "extern explicit instantiation" to give the class a DLL
// interface.  This is a compiler-specific extension.
svtkInstantiateTemplateMacro(extern template class SVTKCOMMONCORE_EXPORT AOSDataArray);

#pragma warning(pop)

#endif
#endif
// SVTK-HeaderTest-Exclude: AOSDataArray.h
*/
