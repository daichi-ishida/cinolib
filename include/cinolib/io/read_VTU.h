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
#ifndef CINO_READ_VTU_H
#define CINO_READ_VTU_H

#include <sys/types.h>
#include <vector>
#include <cinolib/cino_inline.h>
#include <cinolib/geometry/vec_mat.h>


namespace cinolib
{

struct cinoHash
{
    CINO_INLINE
    std::size_t operator()(const cinolib::vec3d& v) const;
};

struct faceEqualTo
{
    CINO_INLINE
    bool operator()(const std::vector<uint>& lhs, const std::vector<uint>& rhs) const;
};

struct faceHash
{
    CINO_INLINE
    std::size_t operator()(const std::vector<uint>& face) const;
};

CINO_INLINE
void read_VTU(const char          * filename,
               std::vector<double> & xyz,
               std::vector<uint>   & tet,
               std::vector<uint>   & hexa);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
void read_VTU(const char                      * filename,
               std::vector<double>            & xyz,
               std::vector<std::vector<uint>> & poly);


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
void read_VTU(const char                      * filename,
               std::vector<vec3d>             & verts,
               std::vector<std::vector<uint>> & poly);


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

CINO_INLINE
void read_VTU(const char                      * filename,
               std::vector<vec3d>             & verts,
               std::vector<std::vector<uint>> & faces,
               std::vector<std::vector<uint>> & polys,
               std::vector<std::vector<bool>> & windings);
}

#ifndef  CINO_STATIC_LIB
#include "read_VTU.cpp"
#endif

#endif // CINO_READ_VTU
