/********************************************************************************
*  This file is part of CinoLib                                                 *
*  Copyright(C) 2016: Marco Livesu                                              *
*                                                                               *
*  The MIT License                                                              *
*                                                                               *
*  Permission is hereby granted, free of charge, to any person obtaining a      *
*  copy of this software and associated documentation files (the "Software"),   *
*  to deal in the Software without restriction, including without limitation    *
*  the rights to use, copy, modify, merge, publish, distribute, sublicense,     *
*  and/or sell copies of the Software, and to permit persons to whom the        *
*  Software is furnished to do so, subject to the following conditions:         *
*                                                                               *
*  The above copyright notice and this permission notice shall be included in   *
*  all copies or substantial portions of the Software.                          *
*                                                                               *
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR   *
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
*  FITNESS FOR A PARTICULAR PURPOSE AND NON INFRINGEMENT. IN NO EVENT SHALL THE *
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING      *
*  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS *
*  IN THE SOFTWARE.                                                             *
*                                                                               *
*  Author(s):                                                                   *
*                                                                               *
*     Marco Livesu (marco.livesu@gmail.com)                                     *
*     http://pers.ge.imati.cnr.it/livesu/                                       *
*                                                                               *
*     Italian National Research Council (CNR)                                   *
*     Institute for Applied Mathematics and Information Technologies (IMATI)    *
*     Via de Marini, 6                                                          *
*     16149 Genoa,                                                              *
*     Italy                                                                     *
*********************************************************************************/
#include <cinolib/io/read_VTU.h>

#ifdef CINOLIB_USES_VTK
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#endif

#include <unordered_map>
#include <algorithm>


namespace cinolib
{

CINO_INLINE
void hash_combine(size_t &seed, size_t hash)
{
    hash += 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hash;
}

CINO_INLINE
std::size_t cinoHash::operator()(const cinolib::vec3d& v) const
{
    std::size_t h = 0;
    hash_combine(h, std::hash<double>()(v.x()));
    hash_combine(h, std::hash<double>()(v.y()));
    hash_combine(h, std::hash<double>()(v.z()));
    return h;
}

CINO_INLINE
bool faceEqualTo::operator()(const std::vector<uint>& lhs, const std::vector<uint>& rhs) const
{
    if(lhs.size() != rhs.size()) return false;
    std::vector<uint> tmp_lhs(lhs);
    std::vector<uint> tmp_rhs(rhs);
    std::sort(tmp_lhs.begin(), tmp_lhs.end());
    std::sort(tmp_rhs.begin(), tmp_rhs.end());
    return std::equal(tmp_lhs.cbegin(), tmp_lhs.cend(), tmp_rhs.cbegin());
}

CINO_INLINE
std::size_t faceHash::operator()(const std::vector<uint>& face) const
{
    assert(face.size()>2);
    std::vector<uint> tmp(face);
    std::sort(tmp.begin(), tmp.end());
    std::size_t h = 0;
    hash_combine(h, std::hash<uint>()(tmp[0]));
    hash_combine(h, std::hash<uint>()(tmp[1]));
    hash_combine(h, std::hash<uint>()(tmp[2]));
    return h;
}

#ifdef CINOLIB_USES_VTK

CINO_INLINE
void read_VTU(const char                      * filename,
               std::vector<vec3d>             & verts,
               std::vector<std::vector<uint>> & faces,
               std::vector<std::vector<uint>> & polys,
               std::vector<std::vector<bool>> & windings)
{
    setlocale(LC_NUMERIC, "en_US.UTF-8"); // makes sure "." is the decimal separator

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> grid(reader->GetOutput());

    verts.clear();
    faces.clear();
    polys.clear();
    windings.clear();

    std::unordered_map<uint, uint> old2new;
    std::unordered_map<vec3d, uint, cinoHash> unique_verts;
    std::unordered_map< std::vector<uint>, uint, faceHash, faceEqualTo > unique_faces;

    for(uint i=0; i<grid->GetNumberOfPoints(); ++i)
    {
        double pnt[3];
        grid->GetPoint(i, pnt);
        vec3d vert(pnt[0],pnt[1],pnt[2]);
        if(unique_verts.count(vert) == 0)
        {
            unique_verts[vert] = i;
            old2new[i] = i;
            verts.emplace_back(vert);
        }
        else
        {
            old2new[i] = unique_verts[vert];
        }
    }

    for(uint i=0; i<grid->GetNumberOfCells(); ++i)
    {
        std::vector<uint> poly;
        std::vector<bool> winding;
        vtkCell *c = grid->GetCell(i);
        if(c == nullptr) continue;

        for(uint j=0; j<c->GetNumberOfFaces(); ++j)
        {
            std::vector<uint> face;
            vtkCell *f = c->GetFace(j);
            if(f == nullptr) continue;
            for(uint k=0; k<f->GetNumberOfPoints(); ++k)
            {
                uint vid = old2new[static_cast<uint>(f->GetPointId(k))];
                face.emplace_back(vid);
            }
            if(face.empty()) continue;

            auto num_duplicated_faces = unique_faces.count(face);

            if(num_duplicated_faces == 0)
            {
                uint fid = static_cast<uint>(faces.size());
                unique_faces[face] = fid;
                poly.emplace_back(fid);
                winding.emplace_back(true);
                faces.emplace_back(face);
            }
            else if(num_duplicated_faces == 1)
            {
                poly.emplace_back(unique_faces[face]);
                winding.emplace_back(false);
            }
            else
            {
                printf("duplicated face!\n");
                assert(false);
            }
        }

        if(!poly.empty())
        {
            polys.emplace_back(poly);
            windings.emplace_back(winding);
        }
    }
}

CINO_INLINE
void read_VTU(const char                      * filename,
               std::vector<vec3d>             & verts,
               std::vector<std::vector<uint>> & poly)
{
    setlocale(LC_NUMERIC, "en_US.UTF-8"); // makes sure "." is the decimal separator

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> grid(reader->GetOutput());

    for(uint i=0; i<grid->GetNumberOfPoints(); ++i)
    {
        double pnt[3];
        grid->GetPoint(i, pnt);

        verts.push_back(vec3d(pnt[0],pnt[1],pnt[2]));
    }

    for(uint i=0; i<grid->GetNumberOfCells(); ++i)
    {
        vtkCell *c = grid->GetCell(i);

        std::vector<uint> polyhedron;
        switch (c->GetCellType())
        {
            case VTK_TETRA:      for(uint j=0; j<4; ++j) polyhedron.push_back(c->GetPointId(j)); break;
            case VTK_HEXAHEDRON: for(uint j=0; j<8; ++j) polyhedron.push_back(c->GetPointId(j)); break;
        }

        if(!polyhedron.empty()) poly.push_back(polyhedron);
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
void read_VTU(const char                      * filename,
               std::vector<double>            & xyz,
               std::vector<std::vector<uint>> & poly)
{
    setlocale(LC_NUMERIC, "en_US.UTF-8"); // makes sure "." is the decimal separator

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> grid(reader->GetOutput());

    for(uint i=0; i<grid->GetNumberOfPoints(); ++i)
    {
        double pnt[3];
        grid->GetPoint(i, pnt);

        xyz.push_back(pnt[0]);
        xyz.push_back(pnt[1]);
        xyz.push_back(pnt[2]);
    }

    for(uint i=0; i<grid->GetNumberOfCells(); ++i)
    {
        vtkCell *c = grid->GetCell(i);

        std::vector<uint> polyhedron;
        switch (c->GetCellType())
        {
            case VTK_TETRA:      for(uint j=0; j<4; ++j) polyhedron.push_back(c->GetPointId(j)); break;
            case VTK_HEXAHEDRON: for(uint j=0; j<8; ++j) polyhedron.push_back(c->GetPointId(j)); break;
        }

        if(!polyhedron.empty()) poly.push_back(polyhedron);
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
void read_VTU(const char           * filename,
               std::vector<double> & xyz,
               std::vector<uint>   & tets,
               std::vector<uint>   & hexa)
{

    setlocale(LC_NUMERIC, "en_US.UTF-8"); // makes sure "." is the decimal separator

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> grid(reader->GetOutput());

    for(uint i=0; i<grid->GetNumberOfPoints(); ++i)
    {
        double pnt[3];
        grid->GetPoint(i, pnt);

        xyz.push_back(pnt[0]);
        xyz.push_back(pnt[1]);
        xyz.push_back(pnt[2]);
    }

    for(uint i=0; i<grid->GetNumberOfCells(); ++i)
    {
        vtkCell *c = grid->GetCell(i);

        switch (c->GetCellType())
        {
            case VTK_TETRA:      for(uint j=0; j<4; ++j) tets.push_back(c->GetPointId(j)); break;
            case VTK_HEXAHEDRON: for(uint j=0; j<8; ++j) hexa.push_back(c->GetPointId(j)); break;
        }
    }
}

#else

CINO_INLINE
void read_VTU(const char                      * filename,
               std::vector<vec3d>             & verts,
               std::vector<std::vector<uint>> & faces,
               std::vector<std::vector<uint>> & polys,
               std::vector<std::vector<bool>> & windings)
{
    std::cerr << "ERROR : VTK missing. Install VTK and recompile defining symbol CINOLIB_USES_VTK" << std::endl;
    exit(-1);
}

CINO_INLINE
void read_VTU(const char          *,
               std::vector<double> &,
               std::vector<uint>   &,
               std::vector<uint>   &)
{
    std::cerr << "ERROR : VTK missing. Install VTK and recompile defining symbol CINOLIB_USES_VTK" << std::endl;
    exit(-1);
}

CINO_INLINE
void read_VTU(const char                      *,
               std::vector<double>            &,
               std::vector<std::vector<uint>> &)
{
    std::cerr << "ERROR : VTK missing. Install VTK and recompile defining symbol CINOLIB_USES_VTK" << std::endl;
    exit(-1);
}

CINO_INLINE
void read_VTU(const char                      *,
               std::vector<vec3d>             &,
               std::vector<std::vector<uint>> &)
{
    std::cerr << "ERROR : VTK missing. Install VTK and recompile defining symbol CINOLIB_USES_VTK" << std::endl;
    exit(-1);
}

#endif

}
