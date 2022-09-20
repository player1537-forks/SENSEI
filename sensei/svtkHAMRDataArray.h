#ifndef svtkHAMRDataArray_h
#define svtkHAMRDataArray_h

#include "hamr_buffer.h"

#include "svtkDataArray.h"
#include "svtkAOSDataArrayTemplate.h"
#include "svtkCommonCoreModule.h"          // For export macro

/*
class svtkDoubleArray;
class svtkIdList;
class svtkInformationStringKey;
class svtkInformationDoubleVectorKey;
class svtkLookupTable;
class svtkPoints;
*/

using Allocator = hamr::buffer_allocator;

template <typename T>
class SVTKCOMMONCORE_EXPORT svtkHAMRDataArray : public svtkDataArray
{
public:
  svtkTypeMacro(svtkHAMRDataArray, svtkDataArray);
  void PrintSelf(ostream& os, svtkIndent indent) override;

  /// allocate a new empty array
  static svtkHAMRDataArray *New();
  
  /** zero-copy the passed data. the allocator is used to tell where the data
   * resides. the callee (array instance) takes ownership of the pointer.
   */
  static svtkHAMRDataArray *New(T *data, size_t numTuples, int numComps, Allocator alloc, int owner);

  /** zero-copy the passed data. the allocator is used to tell where the data
   * resides.
   */
  static svtkHAMRDataArray *New(const std::shared_ptr<T> &data, size_t numTuples, int numComps,
    Allocator alloc, int owner);

  /** zero-copy the passed data. the allocator is used to tell where the data
   * resides the deleter will be called as void deleter(void *dataPtr) when the
   * data is no longer needed
   */
  template <typename deleter_t>
  static svtkHAMRDataArray *New(T *dataPtr, size_t numTuples, int numComps,
    Allocator alloc, int owner, deleter_t deleter);

  /** Allocate a new array of the specified size using the specified allocator */
  static svtkHAMRDataArray *New(Allocator alloc, size_t nTuples, int nComps);

  /** Allocate a new array of the specified size using the specified allocator
   * initialized to the specified value */
  static svtkHAMRDataArray *New(Allocator alloc, size_t nTuples, int nComps, const T &initVal);

  /** Sets or changes the allocator used to manage the menory, this may move
   * the data from one device to another
   */
  void SetAllocator(Allocator alloc);

  /// resize the internal buffer
  void Resize(size_t numTuples, int numComps);

  /** copy contents of the passed in data. allocator and owner tells the current
   * location of the passed data. Use ::Resize to allocate space.
   */
  void CopyData(size_t destStart, T *srcData, size_t numTuples,
    Allocator srcAlloc = Allocator::malloc, int owner = -1);

  /** append the contents of the passed in data. allocator and owner tells the current
   * location of the passed data.
   */
  void AppendData(T *srcData, size_t numTuples,
    Allocator srcAlloc = Allocator::malloc, int owner = -1);

  /** Convert via a deep copy to a SVYTK data array sub class that is 
   * safe to use with SVTK
   */
  svtkAOSDataArrayTemplate<T> *AsSvtkDataArray();

/** zero-copy the passed data. the allocator is used to tell where the data
   * resides. the callee (array instance) takes ownership of the pointer.
   */
  void SetData(T *data, size_t numTuples, int numComps, Allocator alloc, int owner);

  /** zero-copy the passed data. the allocator is used to tell where the data
   * resides.
   */
  void SetData(const std::shared_ptr<T> &data, size_t numTuples, int numComps,
    Allocator alloc, int owner);

  /** zero-copy the passed data. the allocator is used to tell where the data
   * resides the deleter will be called as void deleter(void *dataPtr) when the
   * data is no longer needed
   */
  template <typename deleter_t>
  void SetData(T *dataPtr, size_t numTuples, int numComps,
    Allocator alloc, int owner, deleter_t deleter);

  /// returns a pointer to the data that is safe to use on the CPU
  std::shared_ptr<T> GetCPUAccessible() { return this->Data->get_cpu_accessible(); }

  /// returns a pointer to the data that is safe to use with CUDA
  std::shared_ptr<T> GetCUDAAccessible() { return this->Data->get_cuda_accessible(); }

  /// returns a pointer to the data that is safe to use with HIP
  std::shared_ptr<T> GetHIPAccessible() { return this->Data->get_hip_accessible(); }

  /// returns a pointer to the data that is safe to use with OpenMP device off load
  std::shared_ptr<T> GetOpenMPAccessible() { return this->Data->get_openmp_accessible(); }

  /// return true if a pooniter to the data is safe to use on the CPU
  bool CPUAccessible() { return this->Data->cpu_accessible(); }

  /// return true if a pooniter to the data is safe to use with CUDA
  bool CUDAAccessible() { return this->Data->cuda_accessible(); }

  /// returns a pointer to the data that is safe to use with HIP
  bool HIPAccessible() { return this->Data->hip_accessible(); }

  /// return true if a pooniter to the data is safe to use with OpenMP device off load
  bool OpenMPAccessible() { return this->Data->openmp_accessible(); }

  /** fast access to the internally managed memory. Use this only when you know
   * where the data resides and will access it in that location. This method
   * saves the cost of a smart_ptr copy construct and the cost of the logic
   * that determines if a temporary is needed. For all other cases use
   * GetXAccessible to access the data.
   */
  T *GetData() { return this->Data->data(); }

  /** fast access to the internally managed memory. Use this only when you know
   * where the data resides and will access it in that location. This method
   * saves the cost of the logic determining if a temporary is needed. For all
   * other cases use GetXAccessible to access the data.
   */
  std::shared_ptr<T> GetDataPointer() { return this->Data->pointer(); }

//// TODO -- need to make use of MaxId and NumberOfCompoents because GetNumberOfTuples is not virtual

  /** @name not implemented
   * These methods are not impelemented. If called they will abort. If you find
   * yourself needing these then you are likely writing VTK code and should
   * convert the svtkHAMRDataArray into a svtkDataArray subclass.
   */
  ///@{
  svtkTypeBool Allocate(svtkIdType numValues, svtkIdType ext) override;

  void Initialize() override;

  int GetDataType() const override;

  int GetDataTypeSize() const override;

  int GetElementComponentSize() const override;

  void SetNumberOfTuples(svtkIdType numTuples) override;

  bool SetNumberOfValues(svtkIdType numValues) override;

  void SetTuple(svtkIdType dstTupleIdx, svtkIdType srcTupleIdx, svtkAbstractArray* source) override;

  void InsertTuple(svtkIdType dstTupleIdx, svtkIdType srcTupleIdx, svtkAbstractArray* source) override;

  void InsertTuples(svtkIdList* dstIds, svtkIdList* srcIds, svtkAbstractArray* source) override;

  void InsertTuples(svtkIdType dstStart, svtkIdType n, svtkIdType srcStart, svtkAbstractArray* source) override;

  svtkIdType InsertNextTuple(svtkIdType srcTupleIdx, svtkAbstractArray* source) override;

  void GetTuples(svtkIdList* tupleIds, svtkAbstractArray* output) override;

  void GetTuples(svtkIdType p1, svtkIdType p2, svtkAbstractArray* output) override;

  bool HasStandardMemoryLayout() const override;

  void* GetVoidPointer(svtkIdType valueIdx) override;

  void DeepCopy(svtkAbstractArray* da) override;

  void InterpolateTuple(svtkIdType dstTupleIdx, svtkIdList* ptIndices, svtkAbstractArray* source, double* weights) override;

  void InterpolateTuple(svtkIdType dstTupleIdx, svtkIdType srcTupleIdx1, svtkAbstractArray* source1, svtkIdType srcTupleIdx2, svtkAbstractArray* source2, double t) override;

  void Squeeze() override;

  svtkTypeBool Resize(svtkIdType numTuples) override;

  void SetVoidArray(void* svtkNotUsed(array), svtkIdType svtkNotUsed(size), int svtkNotUsed(save)) override;

  void SetArrayFreeFunction(void (*callback)(void*)) override;

  void ExportToVoidPointer(void* out_ptr) override;

  unsigned long GetActualMemorySize() const override;

  int IsNumeric() const override;

  SVTK_NEWINSTANCE svtkArrayIterator* NewIterator() override;

  svtkIdType LookupValue(svtkVariant value) override;
  void LookupValue(svtkVariant value, svtkIdList* valueIds) override;
  svtkVariant GetVariantValue(svtkIdType valueIdx) override;

  void InsertVariantValue(svtkIdType valueIdx, svtkVariant value) override;

  void SetVariantValue(svtkIdType valueIdx, svtkVariant value) override;

  void DataChanged() override;

  void ClearLookup() override;

  void GetProminentComponentValues(int comp, svtkVariantArray* values, double uncertainty = 1.e-6, double minimumProminence = 1.e-3) override;

  int CopyInformation(svtkInformation* infoFrom, int deep = 1) override;

  double* GetTuple(svtkIdType tupleIdx) override;

  void GetTuple(svtkIdType tupleIdx, double* tuple) override;

  void SetTuple(svtkIdType tupleIdx, const float* tuple) override;
  void SetTuple(svtkIdType tupleIdx, const double* tuple) override;

  void InsertTuple(svtkIdType tupleIdx, const float* tuple)  override;
  void InsertTuple(svtkIdType tupleIdx, const double* tuple)  override;

  svtkIdType InsertNextTuple(const float* tuple) override;
  svtkIdType InsertNextTuple(const double* tuple) override;

  void RemoveTuple(svtkIdType tupleIdx) override;
  void RemoveFirstTuple() override { this->RemoveTuple(0); }
  void RemoveLastTuple() override;

  double GetComponent(svtkIdType tupleIdx, int compIdx) override;

  void SetComponent(svtkIdType tupleIdx, int compIdx, double value) override;

  void InsertComponent(svtkIdType tupleIdx, int compIdx, double value) override;

  void GetData(svtkIdType tupleMin, svtkIdType tupleMax, int compMin, int compMax, svtkDoubleArray* data) override;

  void DeepCopy(svtkDataArray* da) override;

  void ShallowCopy(svtkDataArray* other) override;

  void FillComponent(int compIdx, double value) override;

  void Fill(double value) override;

  void CopyComponent(int dstComponent, svtkDataArray* src, int srcComponent) override;

  void* WriteVoidPointer(svtkIdType valueIdx, svtkIdType numValues) override;

  double GetMaxNorm() override;

  int GetArrayType() const override;
  //@}





//////////////////////////// from generic data array

  /** Get a pointer to the data. The data must already be on the CPU. Call
   * CPUAccessible() instead. */
  /*void* GetVoidPointer(svtkIdType valueIdx) override
  {
    if (!this->Data->cpu_accessible())
    {
      svtkErrroMacro("Accessing a device pointer on the CPU."
        " Call CpuAccessible() instead to move the data.")
      return nullptr;
    }
    return this->Data->data() + valueIdx;
  }*/

  /** Get a pointer to the data. The data must already be on the CPU. Call
   * CPUAccessible() instead. */
  T* GetPointer(svtkIdType valueIdx)
  {
    return this->GetVoidPointer(0);
  }

  /** Pass a pointer to CPU data.
   * @param[in] ptr a pointer to CPU backed data
   * @param[in] size the number of elements of type T
   * @param[in] save if 1 the caller will delete the data
   * @param[in] deleteMethod a flag indicating how to delete the pointed to data.
   *                         SVTK_DATA_ARRAY_FREE, SVTK_DATA_ARRAY_DELETE
   */
  void SetVoidArray(void* ptr, svtkIdType size, int save, int deleteMethod) override
  {
    if (this->Data)
    {
      delete this->Data;
      this->Data = nullptr;
    }

    if (save)
    {
      this->Data = new hamr::buffer<T>(hamr::buffer_allocator::malloc, size,
                                       -1, (T*)ptr, [](void *) -> void {});
    }
    else
    {
      if (deleteMethod == SVTK_DATA_ARRAY_FREE)
      {
        this->Data = new hamr::buffer<T>(hamr::buffer_allocator::malloc, size, -1, (T*)ptr);
      }
      else if (deleteMethod == SVTK_DATA_ARRAY_DELETE)
      {
        this->Data = new hamr::buffer<T>(hamr::buffer_allocator::cpp, size, -1, (T*)ptr);
      }
      else
      {
        svtkErrorMacro("Unsupported delete method for CPU based data");
        abort();
      }
    }
  }





protected:
  svtkHAMRDataArray();
  ~svtkHAMRDataArray() override;

private:
  hamr::buffer<T> *Data;

private:
  svtkHAMRDataArray(const svtkHAMRDataArray&) = delete;
  void operator=(const svtkHAMRDataArray&) = delete;
};

// --------------------------------------------------------------------------
template<typename T>
svtkHAMRDataArray<T>::svtkHAMRDataArray() : Data(nullptr)
{
}

// --------------------------------------------------------------------------
template<typename T>
svtkHAMRDataArray<T>::~svtkHAMRDataArray()
{
  delete this->Data;
}

// --------------------------------------------------------------------------
template<typename T>
svtkHAMRDataArray<T> *svtkHAMRDataArray<T>::New()
{
  return new svtkHAMRDataArray<T>;
}

// --------------------------------------------------------------------------
template<typename T>
svtkHAMRDataArray<T> *svtkHAMRDataArray<T>::New(Allocator alloc, size_t numTuples,
  int numComps)
{
  svtkHAMRDataArray<T> *tmp = new svtkHAMRDataArray<T>;
  tmp->Data = new hamr::buffer<T>(alloc, numTuples*numComps);
  tmp->MaxId = numTuples - 1;
  tmp->SetNumberOfComponents(numComps);
  return tmp;
}

// --------------------------------------------------------------------------
template<typename T>
svtkHAMRDataArray<T> *svtkHAMRDataArray<T>::New(Allocator alloc,
  size_t numTuples, int numComps, const T &initVal)
{
  svtkHAMRDataArray<T> *tmp = new svtkHAMRDataArray<T>;
  tmp->Data = new hamr::buffer<T>(alloc, numTuples*numComps, initVal);
  tmp->MaxId = numTuples - 1;
  tmp->SetNumberOfComponents(numComps);
  return tmp;
}

// --------------------------------------------------------------------------
template<typename T>
svtkHAMRDataArray<T> *svtkHAMRDataArray<T>::New(T *data, size_t numTuples,
  int numComps, Allocator alloc, int owner)
{
  svtkHAMRDataArray<T> *tmp = new svtkHAMRDataArray<T>;
  tmp->Data = new hamr::buffer<T>(alloc, numTuples*numComps, owner, data);
  tmp->MaxId = numTuples - 1;
  tmp->SetNumberOfComponents(numComps);
  return tmp;
}

// --------------------------------------------------------------------------
template<typename T>
svtkHAMRDataArray<T> *svtkHAMRDataArray<T>::New(const std::shared_ptr<T> &data,
    size_t numTuples, int numComps, Allocator alloc, int owner)
{
  svtkHAMRDataArray<T> *tmp = new svtkHAMRDataArray<T>;
  tmp->Data = new hamr::buffer<T>(alloc, numTuples*numComps, owner, data);
  tmp->MaxId = numTuples - 1;
  tmp->SetNumberOfComponents(numComps);
  return tmp;
}

// --------------------------------------------------------------------------
template<typename T>
template <typename deleter_t>
svtkHAMRDataArray<T> *svtkHAMRDataArray<T>::New(T *data,
  size_t numTuples, int numComps, Allocator alloc, int owner,
  deleter_t deleter)
{
  svtkHAMRDataArray<T> *tmp = new svtkHAMRDataArray<T>;
  tmp->Data = new hamr::buffer<T>(alloc, numTuples*numComps, owner, data, deleter);
  tmp->MaxId = numTuples - 1;
  tmp->SetNumberOfComponents(numComps);
  return tmp;
}

// --------------------------------------------------------------------------
template<typename T>
void svtkHAMRDataArray<T>::SetData(T *data, size_t numTuples,
  int numComps, Allocator alloc, int owner)
{
  delete this->Data;
  this->Data = new hamr::buffer<T>(alloc, numTuples*numComps, owner, data);
  this->MaxId = numTuples - 1;
  this->SetNumberOfComponents(numComps);
  return this;
}

// --------------------------------------------------------------------------
template<typename T>
void svtkHAMRDataArray<T>::SetData(const std::shared_ptr<T> &data,
    size_t numTuples, int numComps, Allocator alloc, int owner)
{
  delete this->Data;
  this->Data = new hamr::buffer<T>(alloc, numTuples*numComps, owner, data);
  this->MaxId = numTuples - 1;
  this->SetNumberOfComponents(numComps);
  return this;
}

// --------------------------------------------------------------------------
template<typename T>
template <typename deleter_t>
void svtkHAMRDataArray<T>::SetData(T *data,
  size_t numTuples, int numComps, Allocator alloc, int owner,
  deleter_t deleter)
{
  delete this->Data;
  this->Data = new hamr::buffer<T>(alloc, numTuples*numComps, owner, data, deleter);
  this->MaxId = numTuples - 1;
  this->SetNumberOfComponents(numComps);
  return this;
}

#define AsSvtkDataArrayImpl(_cls_t, _cpp_t)                     \
  if (std::is_same<_cls_t, _cpp_t>::value)                      \
  {                                                             \
    int nComps = this->GetNumberOfComponents();                 \
    size_t nTups = this->GetNumberOfTuples();                   \
                                                                \
    svtkAOSDataArrayTemplate<_cpp_t> *tmp =                     \
      svtkAOSDataArrayTemplate<_cpp_t>::New();                  \
                                                                \
    tmp->SetNumberOfComponents(nComps);                         \
    tmp->SetName(this->GetName());                              \
                                                                \
    if  (this->Data->cpu_accessible())                          \
    {                                                           \
      tmp->SetVoidArray(this->Data->data(), nTups, 1, 0);       \
    }                                                           \
    else                                                        \
    {                                                           \
      tmp->SetNumberOfTuples(nTups);                            \
      T *ptr = (T*)tmp->GetVoidPointer(0);                      \
      this->Data->get(0, ptr, 0, nTups*nComps);                 \
    }                                                           \
                                                                \
    return tmp;                                                 \
  }

// --------------------------------------------------------------------------
template<typename T>
svtkAOSDataArrayTemplate<T> *svtkHAMRDataArray<T>::AsSvtkDataArray()
{
  AsSvtkDataArrayImpl(T, double)
  else AsSvtkDataArrayImpl(T, float)
  else AsSvtkDataArrayImpl(T, char)
  else AsSvtkDataArrayImpl(T, short)
  else AsSvtkDataArrayImpl(T, int)
  else AsSvtkDataArrayImpl(T, long)
  else AsSvtkDataArrayImpl(T, long long)
  else AsSvtkDataArrayImpl(T, unsigned char)
  else AsSvtkDataArrayImpl(T, unsigned short)
  else AsSvtkDataArrayImpl(T, unsigned int)
  else AsSvtkDataArrayImpl(T, unsigned long)
  else AsSvtkDataArrayImpl(T, unsigned long long)
  else
  {
    svtkErrorMacro("No SVTK type for T"); // TODO -- check at compile time
    abort();
  }
  return nullptr;
}

// --------------------------------------------------------------------------
template<typename T>
svtkArrayIterator* svtkHAMRDataArray<T>::NewIterator()
{
  svtkErrorMacro("Not implemented");
  abort();
  return nullptr;
}

// --------------------------------------------------------------------------
template<typename T>
int svtkHAMRDataArray<T>::GetDataType() const
{
  if (std::is_same<T, double>::value)
  {
    return SVTK_DOUBLE;
  }
  else if (std::is_same<T, float>::value)
  {
    return SVTK_FLOAT;
  }
  else if (std::is_same<T, long long>::value)
  {
    return SVTK_LONG_LONG;
  }
  else if (std::is_same<T, long>::value)
  {
    return SVTK_LONG;
  }
  else if (std::is_same<T, int>::value)
  {
    return SVTK_INT;
  }
  else if (std::is_same<T, short>::value)
  {
    return SVTK_SHORT;
  }
  else if (std::is_same<T, char>::value)
  {
    return SVTK_CHAR;
  }
  else if (std::is_same<T, unsigned long long>::value)
  {
    return SVTK_UNSIGNED_LONG_LONG;
  }
  else if (std::is_same<T, unsigned long>::value)
  {
    return SVTK_UNSIGNED_LONG;
  }
  else if (std::is_same<T, unsigned int>::value)
  {
    return SVTK_UNSIGNED_INT;
  }
  else if (std::is_same<T, unsigned short>::value)
  {
    return SVTK_UNSIGNED_SHORT;
  }
  else if (std::is_same<T, unsigned char>::value)
  {
    return SVTK_UNSIGNED_CHAR;
  }
  else
  {
    svtkErrorMacro("No SVTK type for T");
    abort();
    return -1;
  }
}

// --------------------------------------------------------------------------
template<typename T>
int svtkHAMRDataArray<T>::GetDataTypeSize() const
{
  return sizeof(T);
}

/*
// --------------------------------------------------------------------------
template<typename T>
bool svtkHAMRDataArray<T>::HasStandardMemoryLayout() const
{
  return false;
}
*/
/*
// --------------------------------------------------------------------------
template<typename T>
svtkTypeBool svtkHAMRDataArray<T>::Allocate(svtkIdType size, svtkIdType ext)
{
  (void) ext;
  if (this->Data)
  {
    delete this->Data;
  }

  this->Data = new hamr::buffer<T>(hamr::buffer_allocator::malloc, size);

  return true;
}
*/
/*
// --------------------------------------------------------------------------
template<typename T>
svtkTypeBool svtkHAMRDataArray<T>::Resize(svtkIdType numTuples)
{
  this->Data->resize(numTuples*this->NumberOfComponents);
}
*/
/*
// --------------------------------------------------------------------------
template<typename T>
void svtkHAMRDataArray<T>::SetNumberOfComponents(int numComps)
{
  this->NumberOfComponents = numComps;
}
*/

/*
// --------------------------------------------------------------------------
template<typename T>
void svtkHAMRDataArray<T>::SetNumberOfTuples(svtkIdType numTuples)
{
  if (!this->Data)
  {
      this->Data = new hamr::buffer<T>(hamr::buffer_allocator::malloc,
                                       numTuples*this->NumberOfComponents);
  }
  else
  {
      this->Data->resize(numTuples);
  }

  this->MaxId = numTuples - 1;
}
*/





// --------------------------------------------------------------------------
template <typename T>
svtkTypeBool svtkHAMRDataArray<T>::Allocate(svtkIdType numValues, svtkIdType ext)
{
  svtkErrorMacro("Not implemented");
  abort();
  return false;
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::Initialize()
{
  svtkErrorMacro("Not implemented");
  abort();
}

/*
// --------------------------------------------------------------------------
template <typename T>
int svtkHAMRDataArray<T>::GetDataType() const
{
  svtkErrorMacro("Not implemented");
  abort();
  return 0;
}


// --------------------------------------------------------------------------
template <typename T>
int svtkHAMRDataArray<T>::GetDataTypeSize() const
{
  svtkErrorMacro("Not implemented");
  abort();
  return 0;
}
*/

// --------------------------------------------------------------------------
template <typename T>
int svtkHAMRDataArray<T>::GetElementComponentSize() const
{
  svtkErrorMacro("Not implemented");
  abort();
  return 0;
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::SetNumberOfTuples(svtkIdType numTuples)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
bool svtkHAMRDataArray<T>::SetNumberOfValues(svtkIdType numValues)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::SetTuple(svtkIdType dstTupleIdx, svtkIdType srcTupleIdx, svtkAbstractArray* source)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::InsertTuple(svtkIdType dstTupleIdx, svtkIdType srcTupleIdx, svtkAbstractArray* source)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::InsertTuples(svtkIdList* dstIds, svtkIdList* srcIds, svtkAbstractArray* source)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::InsertTuples(svtkIdType dstStart, svtkIdType n, svtkIdType srcStart, svtkAbstractArray* source)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
svtkIdType svtkHAMRDataArray<T>::InsertNextTuple(svtkIdType srcTupleIdx, svtkAbstractArray* source)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::GetTuples(svtkIdList* tupleIds, svtkAbstractArray* output)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::GetTuples(svtkIdType p1, svtkIdType p2, svtkAbstractArray* output)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
bool svtkHAMRDataArray<T>::HasStandardMemoryLayout() const
{
  svtkErrorMacro("Not implemented");
  abort();
  return true;
}

// --------------------------------------------------------------------------
template <typename T>
void* svtkHAMRDataArray<T>::GetVoidPointer(svtkIdType valueIdx)
{
  svtkErrorMacro("Not implemented");
  abort();
  return nullptr;
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::DeepCopy(svtkAbstractArray* da)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::InterpolateTuple(svtkIdType dstTupleIdx, svtkIdList* ptIndices, svtkAbstractArray* source, double* weights)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::InterpolateTuple(svtkIdType dstTupleIdx, svtkIdType srcTupleIdx1, svtkAbstractArray* source1, svtkIdType srcTupleIdx2, svtkAbstractArray* source2, double t)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::Squeeze()
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
svtkTypeBool svtkHAMRDataArray<T>::Resize(svtkIdType numTuples)
{
  svtkErrorMacro("Not implemented");
  abort();
  return true;
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::SetVoidArray(void* svtkNotUsed(array), svtkIdType svtkNotUsed(size), int svtkNotUsed(save))
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::SetArrayFreeFunction(void (*callback)(void*))
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::ExportToVoidPointer(void* out_ptr)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
unsigned long svtkHAMRDataArray<T>::GetActualMemorySize() const
{
  svtkErrorMacro("Not implemented");
  abort();
  return 0;
}

// --------------------------------------------------------------------------
template <typename T>
int svtkHAMRDataArray<T>::IsNumeric() const
{
  svtkErrorMacro("Not implemented");
  abort();
  return 0;
}

// --------------------------------------------------------------------------
template <typename T>
svtkIdType svtkHAMRDataArray<T>::LookupValue(svtkVariant value)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::LookupValue(svtkVariant value, svtkIdList* valueIds)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
svtkVariant svtkHAMRDataArray<T>::GetVariantValue(svtkIdType valueIdx)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::InsertVariantValue(svtkIdType valueIdx, svtkVariant value)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::SetVariantValue(svtkIdType valueIdx, svtkVariant value)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::DataChanged()
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::ClearLookup()
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::GetProminentComponentValues(int comp, svtkVariantArray* values, double uncertainty, double minimumProminence)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
int svtkHAMRDataArray<T>::CopyInformation(svtkInformation* infoFrom, int deep)
{
  svtkErrorMacro("Not implemented");
  abort();
  return 0;
}

// --------------------------------------------------------------------------
template <typename T>
double* svtkHAMRDataArray<T>::GetTuple(svtkIdType tupleIdx)
{
  svtkErrorMacro("Not implemented");
  abort();
  return nullptr;
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::GetTuple(svtkIdType tupleIdx, double* tuple)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::SetTuple(svtkIdType tupleIdx, const float* tuple)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::SetTuple(svtkIdType tupleIdx, const double* tuple)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::InsertTuple(svtkIdType tupleIdx, const float* tuple)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::InsertTuple(svtkIdType tupleIdx, const double* tuple)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
svtkIdType svtkHAMRDataArray<T>::InsertNextTuple(const float* tuple)
{
  svtkErrorMacro("Not implemented");
  abort();
  return 0;
}

// --------------------------------------------------------------------------
template <typename T>
svtkIdType svtkHAMRDataArray<T>::InsertNextTuple(const double* tuple)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::RemoveTuple(svtkIdType tupleIdx)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::RemoveLastTuple()
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
double svtkHAMRDataArray<T>::GetComponent(svtkIdType tupleIdx, int compIdx)
{
  svtkErrorMacro("Not implemented");
  abort();
  return 0.;
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::SetComponent(svtkIdType tupleIdx, int compIdx, double value)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::InsertComponent(svtkIdType tupleIdx, int compIdx, double value)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::GetData( svtkIdType tupleMin, svtkIdType tupleMax, int compMin, int compMax, svtkDoubleArray* data)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::DeepCopy(svtkDataArray* da)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::ShallowCopy(svtkDataArray* other)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::FillComponent(int compIdx, double value)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::Fill(double value)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::CopyComponent(int dstComponent, svtkDataArray* src, int srcComponent)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
void* svtkHAMRDataArray<T>::WriteVoidPointer(svtkIdType valueIdx, svtkIdType numValues)
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
double svtkHAMRDataArray<T>::GetMaxNorm()
{
  svtkErrorMacro("Not implemented");
  abort();
}

// --------------------------------------------------------------------------
template <typename T>
int svtkHAMRDataArray<T>::GetArrayType() const
{
  svtkErrorMacro("Not implemented");
  abort();
}

// ----------------------------------------------------------------------------
template <typename T>
void svtkHAMRDataArray<T>::PrintSelf(ostream& os, svtkIndent indent)
{
  (void) os;
  (void) indent;

  std::cerr << "this->Data = ";
  if (this->Data)
  {
    this->Data->print();
  }
  else
  {
    std::cerr << "nullptr";
  }
}

#endif
