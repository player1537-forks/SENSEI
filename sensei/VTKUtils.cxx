#include "VTKUtils.h"
#include "MPIUtils.h"
#include "MeshMetadata.h"
#include "Error.h"

#include <vtkDataObject.h>
#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkCompositeDataSet.h>
#include <vtkCompositeDataIterator.h>
#include <vtkImageData.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkOverlappingAMR.h>
#include <vtkRectilinearGrid.h>
#include <vtkAMRBox.h>
#include <vtkDataSetAttributes.h>
#include <vtkFieldData.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkObjectBase.h>
#include <vtkObject.h>
#include <vtkDataArray.h>
#include <vtkAbstractArray.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkIntArray.h>
#include <vtkAOSDataArrayTemplate.h>
#include <vtkSOADataArrayTemplate.h>

#include <functional>
#include <mpi.h>

using vtkDataObjectPtr = vtkSmartPointer<vtkDataObject>;
using vtkCompositeDataIteratorPtr = vtkSmartPointer<vtkCompositeDataIterator>;

namespace sensei
{
namespace VTKUtils
{

//----------------------------------------------------------------------------
int GetAssociation(std::string assocStr, int &assoc)
{
  unsigned int n = assocStr.size();
  for (unsigned int i = 0; i < n; ++i)
    assocStr[i] = tolower(assocStr[i]);

  if (assocStr == "point")
    {
    assoc = vtkDataObject::POINT;
    return 0;
    }
  else if (assocStr == "cell")
    {
    assoc = vtkDataObject::CELL;
    return 0;
    }
  else if (assocStr == "field")
    {
    assoc = vtkDataObject::FIELD;
    return 0;
    }

  SENSEI_ERROR("Invalid association \"" << assocStr << "\"")
  return -1;
}

//----------------------------------------------------------------------------
const char *GetAttributesName(int association)
{
  switch (association)
    {
    case vtkDataObject::POINT:
      return "point";
      break;
    case vtkDataObject::CELL:
      return "cell";
      break;
    case vtkDataObject::FIELD:
      return "field";
      break;
    }
  SENSEI_ERROR("Invalid data set attributes association")
  return "";
}

//----------------------------------------------------------------------------
vtkFieldData *GetAttributes(vtkDataSet *dobj, int association)
{
  switch (association)
    {
    case vtkDataObject::POINT:
      return static_cast<vtkFieldData*>(dobj->GetPointData());
      break;
    case vtkDataObject::CELL:
      return static_cast<vtkFieldData*>(dobj->GetCellData());
      break;
    case vtkDataObject::FIELD:
      return static_cast<vtkFieldData*>(dobj->GetFieldData());
      break;
    }
  SENSEI_ERROR("Invalid data set attributes association")
  return nullptr;
}

//----------------------------------------------------------------------------
int Apply(vtkCompositeDataSet *cd, vtkCompositeDataSet *cdo,
  BinaryDatasetFunction &func)
{
  vtkCompositeDataIteratorPtr cdit;
  cdit.TakeReference(cd->NewIterator());
  while (!cdit->IsDoneWithTraversal())
    {
    vtkDataObject *obj = cd->GetDataSet(cdit);
    vtkDataObject *objOut = cdo->GetDataSet(cdit);

    // recurse through nested composite datasets
    if (vtkCompositeDataSet *cdn = dynamic_cast<vtkCompositeDataSet*>(obj))
      {
      vtkCompositeDataSet*cdnOut = static_cast<vtkCompositeDataSet*>(objOut);
      int ret = Apply(cdn, cdnOut, func);
      if (ret < 0)
        {
        // stop with error
        SENSEI_ERROR("Failed to apply to composite data set index "
          << cdit->GetCurrentFlatIndex())
        return -1;
        }
      else if (ret > 0)
        {
        // stop without error
        return 1;
        }
      }
    // process data set leaves
    else if(vtkDataSet *ds = dynamic_cast<vtkDataSet*>(obj))
      {
      vtkDataSet *dsOut = static_cast<vtkDataSet*>(objOut);
      int ret = func(ds, dsOut);
      if (ret < 0)
        {
        // stop with error
        SENSEI_ERROR("Function failed in apply at data set index "
          << cdit->GetCurrentFlatIndex())
        return -1;
        }
      else if (ret)
        {
        // stop without error
        return 1;
        }
      }
    else if (obj)
      {
      SENSEI_ERROR("Can't apply to " << obj->GetClassName())
      return -1;
      }
    cdit->GoToNextItem();
    }
  return 0;
}

//----------------------------------------------------------------------------
int Apply(vtkDataObject *dobj, vtkDataObject *dobjo,
  BinaryDatasetFunction &func)
{
  if (vtkCompositeDataSet *cd = dynamic_cast<vtkCompositeDataSet*>(dobj))
    {
    vtkCompositeDataSet *cdo = static_cast<vtkCompositeDataSet*>(dobjo);
    if (Apply(cd, cdo, func) < 0)
      {
      return -1;
      }
    }
  else if (vtkDataSet *ds = dynamic_cast<vtkDataSet*>(dobj))
    {
    vtkDataSet *dso = static_cast<vtkDataSet*>(dobjo);
    if (func(ds, dso) < 0)
      {
      return -1;
      }
    }
  else
    {
    SENSEI_ERROR("Unsupoorted data object type " << dobj->GetClassName())
    return -1;
    }

  return 0;
}

//----------------------------------------------------------------------------
int Apply(vtkCompositeDataSet *cd, DatasetFunction &func)
{
  vtkCompositeDataIteratorPtr cdit;
  cdit.TakeReference(cd->NewIterator());
  while (!cdit->IsDoneWithTraversal())
    {
    vtkDataObject *obj = cd->GetDataSet(cdit);
    // recurse through nested composite datasets
    if (vtkCompositeDataSet *cdn = dynamic_cast<vtkCompositeDataSet*>(obj))
      {
      int ret = Apply(cdn, func);
      if (ret < 0)
        {
        // stop with error
        SENSEI_ERROR("Failed to apply to composite data set index "
          << cdit->GetCurrentFlatIndex())
        return -1;
        }
      else if (ret > 0)
        {
        // stop without error
        return 1;
        }
      }
    // process data set leaves
    else if(vtkDataSet *ds = dynamic_cast<vtkDataSet*>(obj))
      {
      int ret = func(ds);
      if (ret < 0)
        {
        // stop with error
        SENSEI_ERROR("Function failed to apply to composite data set index "
          << cdit->GetCurrentFlatIndex())
        return -1;
        }
      else if (ret)
        {
        // stop without error
        return 1;
        }
      }
    else if (obj)
      {
      SENSEI_ERROR("Can't apply to " << obj->GetClassName())
      return -1;
      }
    cdit->GoToNextItem();
    }
  return 0;
}

//----------------------------------------------------------------------------
int Apply(vtkDataObject *dobj, DatasetFunction &func)
{
  if (vtkCompositeDataSet *cd = dynamic_cast<vtkCompositeDataSet*>(dobj))
    {
    if (Apply(cd, func) < 0)
      {
      return -1;
      }
    }
  else if (vtkDataSet *ds = dynamic_cast<vtkDataSet*>(dobj))
    {
    if (func(ds) < 0)
      {
      return -1;
      }
    }
  else
    {
    SENSEI_ERROR("Unsupoorted data object type " << dobj->GetClassName())
    return -1;
    }
  return 0;
}

//----------------------------------------------------------------------------
int GetGhostLayerMetadata(vtkDataObject *mesh,
  int &nGhostCellLayers, int &nGhostNodeLayers)
{
  // get the ghost layer metadata
  vtkFieldData *fd = mesh->GetFieldData();

  vtkIntArray *glmd =
    dynamic_cast<vtkIntArray*>(fd->GetArray("senseiGhostLayers"));

  if (!glmd)
    return -1;

  nGhostCellLayers = glmd->GetValue(0);
  nGhostNodeLayers = glmd->GetValue(1);

  return 0;
}

//----------------------------------------------------------------------------
int SetGhostLayerMetadata(vtkDataObject *mesh,
  int nGhostCellLayers, int nGhostNodeLayers)
{
  // pass ghost layer metadata in field data.
  vtkIntArray *glmd = vtkIntArray::New();
  glmd->SetName("senseiGhostLayers");
  glmd->SetNumberOfTuples(2);
  glmd->SetValue(0, nGhostCellLayers);
  glmd->SetValue(1, nGhostNodeLayers);

  vtkFieldData *fd = mesh->GetFieldData();
  fd->AddArray(glmd);
  glmd->Delete();

  return 0;
}

// --------------------------------------------------------------------------
int GetArrayMetadata(vtkDataSetAttributes *dsa, int centering,
  std::vector<std::string> &arrayNames, std::vector<int> &arrayCen,
  std::vector<int> &arrayComps, std::vector<int> &arrayType,
  int &hasGhostArray)
{
  int na = dsa->GetNumberOfArrays();
  for (int i = 0; i < na; ++i)
    {
    vtkDataArray *da = dsa->GetArray(i);
    const char *name = da->GetName();
    arrayNames.emplace_back((name ? name : "unkown"));
    arrayCen.emplace_back(centering);
    arrayComps.emplace_back(da->GetNumberOfComponents());
    arrayType.emplace_back(da->GetDataType());

    if (!hasGhostArray && name && !strcmp("vtkGhostType", name))
      hasGhostArray = 1;
    }
  return 0;
}

// --------------------------------------------------------------------------
int GetArrayMetadata(vtkDataSet *ds, MeshMetadataPtr &metadata)
{
  VTKUtils::GetArrayMetadata(ds->GetPointData(), vtkDataObject::POINT,
    metadata->ArrayName, metadata->ArrayCentering, metadata->ArrayComponents,
    metadata->ArrayType, metadata->NumGhostNodes);

  VTKUtils::GetArrayMetadata(ds->GetCellData(), vtkDataObject::CELL,
    metadata->ArrayName, metadata->ArrayCentering, metadata->ArrayComponents,
    metadata->ArrayType, metadata->NumGhostCells);

  metadata->NumArrays = metadata->ArrayName.size();

  return 0;
}

// --------------------------------------------------------------------------
int GetBlockMetadata(int rank, int id, vtkDataSet *ds,
  const MeshMetadataFlags &flags, std::vector<int> &blockOwner,
  std::vector<int> &blockIds, std::vector<long> &blockPoints,
  std::vector<long> &blockCells, std::vector<long> &blockCellArraySize,
  std::vector<std::array<int,6>> &blockExtents,
  std::vector<std::array<double,6>> &blockBounds)
{
  if (!ds)
    return -1;

  if (flags.BlockDecompSet())
    {
    blockOwner.emplace_back(rank);
    blockIds.emplace_back(id);
    }

  if (flags.BlockSizeSet())
    {
    long nPts = ds->GetNumberOfPoints();
    blockPoints.emplace_back(nPts);

    long nCells = ds->GetNumberOfCells();
    blockCells.emplace_back(nCells);

    long cellArraySize = 0;

    if (vtkUnstructuredGrid *ug = dynamic_cast<vtkUnstructuredGrid*>(ds))
      cellArraySize = ug->GetCells()->GetSize();
    else if (vtkPolyData *pd = dynamic_cast<vtkPolyData*>(ds))
      cellArraySize = pd->GetVerts()->GetSize() + pd->GetLines()->GetSize()
        + pd->GetPolys()->GetSize() + pd->GetStrips()->GetSize();

    blockCellArraySize.emplace_back(cellArraySize);
    }

  if (flags.BlockExtentsSet())
    {
    std::array<int,6> ext;
    if (vtkImageData *im = dynamic_cast<vtkImageData*>(ds))
      {
      im->GetExtent(ext.data());
      }
    else if (vtkRectilinearGrid *rg = dynamic_cast<vtkRectilinearGrid*>(ds))
      {
      rg->GetExtent(ext.data());
      }
    else if (vtkStructuredGrid *sg = dynamic_cast<vtkStructuredGrid*>(ds))
      {
      sg->GetExtent(ext.data());
      }
    else
      {
      SENSEI_ERROR("Extents requested on non-Cartesian block encountered")
      return -1;
      }
    blockExtents.emplace_back(std::move(ext));
    }

  if (flags.BlockBoundsSet())
    {
    std::array<double,6> bounds;
    ds->GetBounds(bounds.data());
    blockBounds.emplace_back(std::move(bounds));
    }

  return 0;
}

// --------------------------------------------------------------------------
int GetBlockMetadata(int rank, int id, vtkDataSet *ds, MeshMetadataPtr metadata)
{
    return GetBlockMetadata(rank, id, ds, metadata->Flags,
      metadata->BlockOwner, metadata->BlockIds, metadata->BlockNumPoints,
      metadata->BlockNumCells, metadata->BlockCellArraySize,
      metadata->BlockExtents, metadata->BlockBounds);
}

// --------------------------------------------------------------------------
int GetMetadata(MPI_Comm comm, vtkCompositeDataSet *cd, MeshMetadataPtr metadata)
{
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  vtkOverlappingAMR *amrds = dynamic_cast<vtkOverlappingAMR*>(cd);

  metadata->MeshType = amrds ? VTK_OVERLAPPING_AMR : VTK_MULTIBLOCK_DATA_SET;

  // get array metadata
  vtkCompositeDataIterator *cdit = cd->NewIterator();
  if (!cdit->IsDoneWithTraversal())
    {
    vtkDataObject *bobj = cd->GetDataSet(cdit);

    metadata->BlockType = bobj->GetDataObjectType();

    if (vtkDataSet *ds = dynamic_cast<vtkDataSet*>(bobj))
      VTKUtils::GetArrayMetadata(ds, metadata);

    if (vtkPointSet *ps = dynamic_cast<vtkPointSet*>(bobj))
      metadata->CoordinateType = ps->GetPoints()->GetData()->GetDataType();
    }

  // get block metadata
  int numBlocks = 0;
  while (!cdit->IsDoneWithTraversal())
    {
    vtkDataObject *dobj = cd->GetDataSet(cdit);
    int bid = std::max(0, int(cdit->GetCurrentFlatIndex() - 1));

    if (vtkDataSet *ds = dynamic_cast<vtkDataSet*>(dobj))
      {
      numBlocks += 1;

      if (VTKUtils::GetBlockMetadata(rank, bid, ds, metadata))
        {
        SENSEI_ERROR("Failed to get block metadata for block "
         << cdit->GetCurrentFlatIndex())
        cdit->Delete();
        return -1;
        }
      }
    cdit->GoToNextItem();
    }

  // set block counts
  metadata->NumBlocksLocal = {numBlocks};

  metadata->NumBlocks = 0;
  cdit->SetSkipEmptyNodes(0);
  cdit->InitTraversal();
  while (!cdit->IsDoneWithTraversal())
    {
    metadata->NumBlocks += 1;
    cdit->GoToNextItem();
    }

  cdit->Delete();

  // get global bounds and extents
  if (metadata->Flags.BlockBoundsSet())
    MPIUtils::GlobalBounds(comm, metadata->BlockBounds, metadata->Bounds);

  if (metadata->Flags.BlockExtentsSet())
    MPIUtils::GlobalBounds(comm, metadata->BlockExtents, metadata->Extent);

  if (amrds)
    {
    // global view of block owner is always required
    if (metadata->Flags.BlockDecompSet())
      MPIUtils::GlobalViewV(comm, metadata->BlockOwner);

    // these are all always global views
    metadata->NumLevels = amrds->GetNumberOfLevels();
    metadata->BlockLevel.resize(metadata->NumBlocks);
    metadata->BoxArray.resize(metadata->NumBlocks);
    metadata->RefRatio.resize(metadata->NumLevels);
    metadata->BlocksPerLevel.resize(metadata->NumLevels);

    int q = 0;
    for (int i = 0; i < metadata->NumLevels; ++i)
      {
      int rr = amrds->GetRefinementRatio(i);
      metadata->RefRatio[i] = {rr, rr, rr};

      int nb = amrds->GetNumberOfDataSets(i);
      metadata->BlocksPerLevel[i] = nb;

      std::vector<double> bounds;
      for (int j = 0; j < nb; ++j, ++q)
        {
        metadata->BlockLevel[q] = i;

        metadata->BoxArray[q].resize(6);
        int *pbaq = metadata->BoxArray[q].data();

        const vtkAMRBox &box = amrds->GetAMRBox(i, j);
        box.GetDimensions(pbaq, pbaq+3);
        }
      }

    if (metadata->Flags.BlockChildrenSet())
      {
      if (!amrds->HasChildrenInformation())
        amrds->GenerateParentChildInformation();

      int q = 0;
      for (int i = 0; i < metadata->NumLevels; ++i)
        {
        int nb = amrds->GetNumberOfDataSets(i);

        for (int j = 0; j < nb; ++j, ++q)
          {
          // parent child info is always local only
          if (metadata->BlockOwner[q] != rank)
            continue;

          unsigned int nc = 0;
          unsigned int *pch = amrds->GetChildren(i, j, nc);
          if (nc)
            {
            std::vector<int> children(pch, pch+nc);
            metadata->BlockChildren.emplace_back(children);
            }
          }
        }
      }
    if (metadata->Flags.BlockParentsSet())
      {
      if (!amrds->HasChildrenInformation())
        amrds->GenerateParentChildInformation();

      int q = 0;
      for (int i = 0; i < metadata->NumLevels; ++i)
        {
        int nb = amrds->GetNumberOfDataSets(i);

        for (int j = 0; j < nb; ++j, ++q)
          {
          // parent child info is always local only
          if (metadata->BlockOwner[q] != rank)
            continue;

          unsigned int nc = 0;
          unsigned int *pp = amrds->GetParents(i, j, nc);
          if (nc)
            {
            std::vector<int> parent(pp, pp+nc);
            metadata->BlockParents.emplace_back(parent);
            }
          }
        }
      }
    }

  return 0;
}

// --------------------------------------------------------------------------
// note: not intended for use on the blocks of a multiblock
int GetMetadata(MPI_Comm comm, vtkDataSet *ds, MeshMetadataPtr metadata)
{
  int rank = 0;
  int nRanks = 1;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nRanks);

  metadata->MeshType = ds->GetDataObjectType();
  metadata->BlockType = ds->GetDataObjectType();

  VTKUtils::GetArrayMetadata(ds, metadata);

  if (VTKUtils::GetBlockMetadata(rank, 0, ds, metadata))
    {
    SENSEI_ERROR("Failed to get block metadata for block " << rank)
    return -1;
    }

  metadata->NumBlocks = nRanks;
  metadata->NumBlocksLocal = {1};

  // get global bounds and extents
  if (metadata->Flags.BlockBoundsSet())
    MPIUtils::GlobalBounds(comm, metadata->BlockBounds, metadata->Bounds);

  if (metadata->Flags.BlockExtentsSet())
    MPIUtils::GlobalBounds(comm, metadata->BlockExtents, metadata->Extent);

  return 0;
}

// --------------------------------------------------------------------------
vtkCompositeDataSetPtr AsCompositeData(MPI_Comm comm,
  vtkDataObject *dobj, bool take)
{
  // make sure we have composite dataset if not create one
  vtkCompositeDataSetPtr cd;
  vtkCompositeDataSet *tmp = nullptr;
  if ((tmp = dynamic_cast<vtkCompositeDataSet*>(dobj)))
    {
    if (take)
      cd.TakeReference(tmp);
    else
      cd = tmp;
    }
  else
    {
    int rank = 0;
    int nRanks = 1;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nRanks);

    vtkMultiBlockDataSet *mb = vtkMultiBlockDataSet::New();
    mb->SetNumberOfBlocks(nRanks);
    mb->SetBlock(rank, dobj);
    if (take)
      dobj->Delete();
    cd.TakeReference(mb);
    }

  return cd;
}

/*
int arrayCpy(void *&wptr, vtkDataArray *da)
{
  unsigned long nt = da->GetNumberOfTuples();
  unsigned int nc = da->GetNumberOfComponents();

  switch (da->GetDataType())
    {
    vtkTemplateMacro(
      if (vtkAOSDataArrayTemplate<VTK_TT> *aosda =
        dynamic_cast<vtkAOSDataArrayTemplate<VTK_TT>*>(da))
        {
        unsigned long long nb = nt*nc*sizeof(VTK_TT);
        VTK_TT *pda = aosda->GetPointer(0);
        memcpy(wptr, pda, nb);
        ((char*)wptr) += nb;
        }
      else if (vtkSOADataArrayTemplate<VTK_TT> *soada =
        dynamic_cast<vtkSOADataArrayTemplate<VTK_TT>*>(da))
        {
        unsigned long long nb = nt*sizeof(VTK_TT);
        for (unsigned int j = 0; j < nc; ++j)
          {
          VTK_TT *pda = soada->GetComponentArrayPointer(j);
          memcpy(wptr, pda, nb);
          ((char*)wptr) += nb;
          }
        }
      else
        {
        SENSEI_ERROR("Invalid data array type " << da->GetClassName())
        return -1;
        }
    )
    }

  return 0;
}


// --------------------------------------------------------------------------
int DataArraySerializer::operator()(vtkDataSet *ds)
{
  vtkDataArray *da = m_centering == vtkDataObject::POINT ?
    ds->GetPointData()->GetName(m_name) : ds->GetCellData()->GetName(m_name);

  if (!da || arrayCpy(m_write_ptr, da))
    {
    SENSEI_ERROR("Failed to serialize "
      << GetAttributesName(m_centering) << " data array \""
      << m_name << "\"")
    return -1;
    }

  return 0;
}

// --------------------------------------------------------------------------
int PointsSerializer::operator()(vtkDataSet *ds)
{
  vtkPointSet *ps = dynamic_cast<vtkPointSet*>(ps);
  if (!ps)
    {
    SENSEI_ERROR("Invalid dataset type " << ds->GetClassName())
    return -1;
    }

  vtkDataArray *da = ps->GetPoints()->GetData();
  if (!da || arrayCpy(m_write_ptr, da))
    {
    SENSEI_ERROR("Failed to serialize points")
    return -1;
    }

  return 0;
}

// --------------------------------------------------------------------------
int CellTypesSerializer::operator()(vtkDataSet *ds)
{
  vtkDataArray *da = nullptr;
  if (vtkUnstructuredGrid *ug = dynamic_cast<vtkUnstructuredGrid*>(ds))
    {
    da = ug->GetCellTypesArray();
    if (!da || arrayCpy(m_write_ptr, da))
      {
      SENSEI_ERROR("Failed to serialize cell types")
      return -1;
      }
    }
  else if (vtkPolyData *pd = dynamic_cast<vtkPolyData*>(ds))
    {
    vtkIdType nv = pd->GetNumberOfVerts();
    memset(m_write_ptr, nv, VTK_VERTEX);
    m_write_ptr += nv;

    vtkIdType nl = pd->GetNumberOfLines();
    memset(m_write_ptr, nl, VTK_LINE);
    m_write_ptr += nl;

    vtkIdType np = pd->GetNumberOfPolys();
    memset(m_write_ptr, np, VTK_POLYGON);
    m_write_ptr += np;

    vtkIdType ns = pd->GetNumberOfStrips();
    memset(m_write_ptr, ns, VTK_TRIANGLE_STRIP)
    m_write_ptr += ns;
    }
  else
    {
    SENSEI_ERROR("Invalid dataset type " << ds->GetClassName())
    return -1;
    }

  return 0;
}

// --------------------------------------------------------------------------
int CellArraySerializer::operator()(vtkDataSet *ds)
{
  vtkDataArray *da = nullptr;
  if (vtkUnstructuredGrid *ug = dynamic_cast<vtkUnstructuredGrid*>(ds))
    {
    da = ug->GetCells()->GetData();
    if (!da || arrayCpy(m_write_ptr, da))
      {
      SENSEI_ERROR("Failed to serialize cells")
      return -1;
      }
    }
  else if (vtkPolyData *pd = dynamic_cast<vtkPolyData*>(ds))
    {
    da = pd->GetVerts()->GetData();
    if (!da || arrayCpy(m_write_ptr, da))
      {
      SENSEI_ERROR("Failed to serialize verts")
      return -1;
      }

    da = pd->GetLines()->GetData();
    if (!da || arrayCpy(m_write_ptr, da))
      {
      SENSEI_ERROR("Failed to serialize lines")
      return -1;
      }

    da = pd->GetLines()->GetPolys();
    if (!da || arrayCpy(m_write_ptr, da))
      {
      SENSEI_ERROR("Failed to serialize polys")
      return -1;
      }

    da = pd->GetStrips()->GetData();
    if (!da || arrayCpy(m_write_ptr, da))
      {
      SENSEI_ERROR("Failed to serialize strips")
      return -1;
      }
    }
  else
    {
    SENSEI_ERROR("Invalid dataset type " << ds->GetClassName())
    return -1;
    }

  return 0;
}

// --------------------------------------------------------------------------
int GetSizesAndOffsets(MPI_Comm comm,
  const sensei::MeshMetadataPtr &md,
  unsigned long long &num_points_total,
  unsigned long long &num_points_local,
  unsigned long long &point_offset_local,
  unsigned long long &num_cells_total,
  unsigned long long &num_cells_local,
  unsigned long long &cell_offset_local,
  unsigned long long &cell_array_size_total,
  unsigned long long &cell_array_size_local,
  unsigned long long &cell_array_offset_local)
{
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  num_points_total = 0;
  num_points_local = 0;
  point_offset_local = 0;
  num_cells_total = 0;
  num_cells_local = 0;
  cell_offset_local = 0;
  cell_array_size_total = 0;
  cell_array_size_local = 0;
  cell_array_offset_local = 0;

  // calculate the number of points and cells total
  // and the local offset to each type of data
  unsigned int num_blocks = md->NumBlocks;
  for (unsigned int i = 0; i < num_blocks; ++i)
    {
    num_points_total += md->BlockNumPoints[i];
    num_cells_total += md->BlockNumCells[i];

    if ((md->BlockType == VTK_POLYDATA) || (md->MeshType == VTK_POLYDATA) ||
     (md->BlockType == VTK_UNSTRUCTURED_GRID) || (md->MeshType == VTK_UNSTRUCTURED_GRID))
     {
     cell_array_size_total += md->BlockCellArraySize[i];
     }

    if (md->BlockOwner[i] < rank)
      {
      point_offset_local += md->BlockNumPoints[i];
      cell_offset_local += md->BlockNumCells[i];
      if ((md->BlockType == VTK_POLYDATA) || (md->MeshType == VTK_POLYDATA) ||
        (md->BlockType == VTK_UNSTRUCTURED_GRID) || (md->MeshType == VTK_UNSTRUCTURED_GRID))
        {
        cell_array_offset_local += md->BlockCellArraySize[i];
        }
      }
    else if (md->BlockOwner[i] == rank)
      {
      num_points_local += md->BlockNumPoints[i]
      num_cells_local += md->BlockNumCells[i]
      if ((md->BlockType == VTK_POLYDATA) || (md->MeshType == VTK_POLYDATA) ||
        (md->BlockType == VTK_UNSTRUCTURED_GRID) || (md->MeshType == VTK_UNSTRUCTURED_GRID))
        {
        cell_array_size_local += md->BlockCellArraySize[i];
        }
      }
    }

  return 0;
}

// --------------------------------------------------------------------------
int GetLocalGeometrySizes(MPI_Comm comm,
  const sensei::MeshMetadataPtr &md,
  unsigned long long &num_points_local,
  unsigned long long &num_cells_local,
  unsigned long long &cell_array_size_local)
{
  int rank = 0;
  MPI_Comm_rank(comm, &rank);

  unsigned int num_blocks = md->NumBlocks;
  for (unsigned int i = 0; i < num_blocks; ++i)
    {
    if (md->BlockOwner[i] == rank)
      {
      num_points_local += md->BlockNumPoints[i]
      num_cells_local += md->BlockNumCells[i]
      if ((md->BlockType == VTK_POLYDATA) || (md->MeshType == VTK_POLYDATA) ||
        (md->BlockType == VTK_UNSTRUCTURED_GRID) || (md->MeshType == VTK_UNSTRUCTURED_GRID))
        {
        cell_array_size_local += md->BlockCellArraySize[i];
        }
      }
    }

  return 0;
}
*/

}
}
