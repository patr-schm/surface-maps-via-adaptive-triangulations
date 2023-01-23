/*
 * This file is part of
 * Surface Maps via Adaptive Triangulations
 * (https://github.com/patr-schm/surface-maps-via-adaptive-triangulations)
 * and is released under the MIT license.
 *
 * Authors: Patrick Schmidt, Janis Born
 */
#pragma once

#include <vector>
#include <OpenMesh/Core/Mesh/Handles.hh>
#include <Eigen/Dense>
#include "Out.hh"

#define EXTERNAL_PROPERTY_RANGE_CHECKS_ENABLED 1

#if EXTERNAL_PROPERTY_RANGE_CHECKS_ENABLED == 1
#define EXTERNAL_PROPERTY_RANGE_CHECK(handle) \
{ \
    ISM_ASSERT_GEQ((int)handle.idx(), 0);  \
    ISM_ASSERT_L((int)handle.idx(), (int)container_.size());  \
}
#endif

namespace SurfaceMaps
{

namespace
{
using VH = OpenMesh::VertexHandle;
using HEH = OpenMesh::HalfedgeHandle;
using EH = OpenMesh::EdgeHandle;
using FH = OpenMesh::FaceHandle;

// Because partial function template specialization is not possible
// we wrap all functions with a specialized
// HandleType but free ValueType in this struct.
template <typename HandleType>
struct TemplateWorkaround
{
    template <typename Mesh, typename ValueType>
    static void toProperty(Mesh& _mesh, OpenMesh::BasePropHandleT<ValueType>& _ph, const std::vector<ValueType>& _container);

    template <typename Mesh>
    static size_t nElements(const Mesh& _mesh);

    template<typename Mesh, typename ValueType>
    static void init_container(const Mesh& _mesh, const OpenMesh::BasePropHandleT<ValueType>& _ph, std::vector<ValueType>& _container);

    template <typename Mesh, typename Function>
    static void foreach(const Mesh& _mesh, Function _f);
};

template <>
struct TemplateWorkaround<VH>
{
    template <typename Mesh, typename ValueType>
    static void toProperty(Mesh& _mesh, OpenMesh::BasePropHandleT<ValueType>& _ph, const std::vector<ValueType>& _container)
    {
        auto ph = static_cast<OpenMesh::VPropHandleT<ValueType>>(_ph);
        if (!ph.is_valid())
            _mesh.add_property(ph);
        for (const auto vh : _mesh.vertices())
            _mesh.property(ph, vh) = _container[vh.idx()];

        _ph = ph;
    }

    template <typename Mesh>
    static size_t nElements(const Mesh& _mesh) { return _mesh.n_vertices(); }

    template<typename Mesh, typename ValueType>
    static void init_container(const Mesh& _mesh, const OpenMesh::BasePropHandleT<ValueType>& _ph, std::vector<ValueType>& _container)
    {
        auto ph = static_cast<OpenMesh::VPropHandleT<ValueType>>(_ph);
        ISM_ASSERT(ph.is_valid());

        _container.clear();
        _container.resize(nElements(_mesh));
        for (const auto vh : _mesh.vertices())
            _container[vh.idx()] = _mesh.property(ph, vh);
    }

    template <typename Mesh, typename Function>
    static void foreach(const Mesh& _mesh, Function _f)
    {
        for (auto v : _mesh.vertices())
            _f(v);
    }
};

template <>
struct TemplateWorkaround<HEH>
{
    template <typename Mesh, typename ValueType>
    static void toProperty(Mesh& _mesh, OpenMesh::BasePropHandleT<ValueType>& _ph, const std::vector<ValueType>& _container)
    {
        auto ph = static_cast<OpenMesh::HPropHandleT<ValueType>>(_ph);
        if (!ph.is_valid())
            _mesh.add_property(ph);
        for (const auto heh : _mesh.halfedges())
            _mesh.property(ph, heh) = _container[heh.idx()];

        _ph = ph;
    }

    template <typename Mesh>
    static size_t nElements(const Mesh& _mesh) { return _mesh.n_halfedges(); }

    template<typename Mesh, typename ValueType>
    static void init_container(const Mesh& _mesh, const OpenMesh::BasePropHandleT<ValueType>& _ph, std::vector<ValueType>& _container)
    {
        auto ph = static_cast<OpenMesh::HPropHandleT<ValueType>>(_ph);
        ISM_ASSERT(ph.is_valid());

        _container.clear();
        _container.resize(nElements(_mesh));
        for (const auto heh : _mesh.halfedges())
            _container[heh.idx()] = _mesh.property(ph, heh);
    }

    template <typename Mesh, typename Function>
    static void foreach(const Mesh& _mesh, Function _f)
    {
        for (auto h : _mesh.halfedges())
            _f(h);
    }
};

template <>
struct TemplateWorkaround<EH>
{
    template <typename Mesh, typename ValueType>
    static void toProperty(Mesh& _mesh, OpenMesh::BasePropHandleT<ValueType>& _ph, const std::vector<ValueType>& _container)
    {
        auto ph = static_cast<OpenMesh::EPropHandleT<ValueType>>(_ph);
        if (!ph.is_valid())
            _mesh.add_property(ph);
        for (const auto eh : _mesh.edges())
            _mesh.property(ph, eh) = _container[eh.idx()];

        _ph = ph;
    }

    template <typename Mesh>
    static size_t nElements(const Mesh& _mesh) { return _mesh.n_edges(); }

    template<typename Mesh, typename ValueType>
    static void init_container(const Mesh& _mesh, const OpenMesh::BasePropHandleT<ValueType>& _ph, std::vector<ValueType>& _container)
    {
        auto ph = static_cast<OpenMesh::EPropHandleT<ValueType>>(_ph);
        ISM_ASSERT(ph.is_valid());

        _container.clear();
        _container.resize(nElements(_mesh));
        for (const auto eh : _mesh.edges())
            _container[eh.idx()] = _mesh.property(ph, eh);
    }

    template <typename Mesh, typename Function>
    static void foreach(const Mesh& _mesh, Function _f)
    {
        for (auto e : _mesh.edges())
            _f(e);
    }
};

template <>
struct TemplateWorkaround<FH>
{
    template <typename Mesh, typename ValueType>
    static void toProperty(Mesh& _mesh, OpenMesh::BasePropHandleT<ValueType>& _ph, const std::vector<ValueType>& _container)
    {
        auto ph = static_cast<OpenMesh::FPropHandleT<ValueType>>(_ph);
        if (!ph.is_valid())
            _mesh.add_property(ph);
        for (const auto fh : _mesh.faces())
            _mesh.property(ph, fh) = _container[fh.idx()];

        _ph = ph;
    }

    template <typename Mesh>
    static size_t nElements(const Mesh& _mesh) { return _mesh.n_faces(); }

    template<typename Mesh, typename ValueType>
    static void init_container(const Mesh& _mesh, const OpenMesh::BasePropHandleT<ValueType>& _ph, std::vector<ValueType>& _container)
    {
        auto ph = static_cast<OpenMesh::FPropHandleT<ValueType>>(_ph);
        ISM_ASSERT(ph.is_valid());

        _container.clear();
        _container.resize(_mesh.n_faces());
        for (const auto fh : _mesh.facehandles())
            _container[fh.idx] = _mesh.property(ph, fh);
    }

    template <typename Mesh, typename Function>
    static void foreach(const Mesh& _mesh, Function _f)
    {
        for (auto f : _mesh.faces())
            _f(f);
    }
};
}

/**
 * An ExternalProperty allows to conveniently assign values to
 * Mesh elements, i.e. vertices, halfedges, edges and faces.
 *
 * In contrast to standard OpenMesh properties, no write operations
 * on the mesh are performed. This allows creating and deleting
 * properties from within multiple threads.
 *
 * This class performs no lifecycle management whatsoever.
 * The lifetime of the property is purely defined by the
 * instantiation and deletion of this class and may even
 * exceed the lifetime of the mesh. Moreover, it is possible
 * to use the same ExternalProperty with different meshes,
 * having the same numbers of elements.
 */
template <typename HandleType, typename ValueType>
class ExternalProperty
{
    private:
        using TW = TemplateWorkaround<HandleType>;

    public:
        /**
         * Uninitialized default constructor.
         * Call init() before usage.
         */
        ExternalProperty() { }

        /**
         * Copies all values.
         */
        ExternalProperty(const ExternalProperty& _other)
        {
            init(_other);
        }

        /**
         * Initializes the property with the correct number of elements,
         * calling the default constructor of ValueType.
         * Warning: Some types, e.g. OpenMesh::Vec3d, remain uninitialized.
         */
        template <typename Mesh>
        ExternalProperty(const Mesh& _mesh)
        {
            init(_mesh);
        }

        /**
         * Initializes the property with the given number of elements,
         * calling the default constructor of ValueType.
         * Warning: Some types, e.g. OpenMesh::Vec3d, remain uninitialized.
         */
        ExternalProperty(const size_t _size)
        {
            init(_size);
        }

        /**
         * Initializes all elements with the given value.
         */
        template <typename Mesh>
        ExternalProperty(const Mesh& _mesh, const ValueType& _value)
        {
            init(_mesh, _value);
        }

        /**
         * Initializes all elements with the given value.
         */
        ExternalProperty(const size_t _size, const ValueType& _value)
        {
            init(_size, _value);
        }

        template <typename Mesh>
        ExternalProperty(const Mesh& _mesh, const std::vector<ValueType>& _vector)
        {
            init(_mesh, _vector);
        }

        template <typename Mesh>
        ExternalProperty(const Mesh& _mesh, const Eigen::VectorX<ValueType>& _vector)
        {
            init(_mesh, _vector);
        }

        template <typename Mesh>
        ExternalProperty(const Mesh& _mesh, const OpenMesh::BasePropHandleT<ValueType> _ph)
        {
            init(_mesh, _ph);
        }

        template <typename Mesh>
        ExternalProperty(const Mesh& _mesh, std::function<ValueType(HandleType)> _f)
        {
            init(_mesh, _f);
        }

        virtual ~ExternalProperty() { }

        /**
         * Initializes the property with the correct number of elements,
         * calling the default constructor of ValueType.
         * Warning: Some types, e.g. OpenMesh::Vec3d, remain uninitialized.
         */
        template <typename Mesh>
        void init(const Mesh& _mesh)
        {
            container_ = std::vector<ValueType>(TW::nElements(_mesh));
        }

        /**
         * Initializes the property with the given number of elements,
         * calling the default constructor of ValueType.
         * Warning: Some types, e.g. OpenMesh::Vec3d, remain uninitialized.
         */
        void init(const size_t _size)
        {
            container_ = std::vector<ValueType>(_size);
        }

        /**
         * Initializes all elements with the given value.
         */
        template <typename Mesh>
        void init(const Mesh& _mesh, const ValueType& _value)
        {
            container_ = std::vector<ValueType>(TW::nElements(_mesh), _value);
        }

        /**
         * Initializes all elements with the given value if the property is empty.
         */
        template <typename Mesh>
        void init_if_empty(const Mesh& _mesh, const ValueType& _value)
        {
            ISM_ASSERT(empty() || size_okay(_mesh));
            if (empty())
                init(_mesh, _value);
        }

        /**
         * Initializes all elements with the given value.
         */
        void init(const size_t _size, const ValueType& _value)
        {
            container_ = std::vector<ValueType>(_size, _value);
        }

        /**
         * Initializes the property via copy.
         */
        void init(const ExternalProperty& _other)
        {
            container_ = _other.container_;
        }

        template <typename Mesh>
        void init(const Mesh& _mesh, const OpenMesh::BasePropHandleT<ValueType>& _ph)
        {
            TW::init_container(_mesh, _ph, container_);
        }

        template <typename Mesh>
        void init(const Mesh& _mesh, const std::vector<ValueType>& _vector)
        {
            ISM_ASSERT_EQ(TW::nElements(_mesh), _vector.size());
            container_ = _vector;
        }

        template <typename Mesh>
        void init(const Mesh& _mesh, const Eigen::VectorX<ValueType>& _vector)
        {
            ISM_ASSERT_EQ(TW::nElements(_mesh), _vector.size());
            container_.resize(_vector.size());
            Eigen::VectorX<ValueType>::Map(container_.data(), _vector.size()) = _vector;
        }

        template <typename Mesh>
        void init(const Mesh& _mesh, std::function<ValueType(HandleType)> _f)
        {
            container_ = std::vector<ValueType>(TW::nElements(_mesh));
            TW::foreach(_mesh, [&] (HandleType handle)
            {
                container_[handle.idx()] = _f(handle);
            });
        }

        void set_constant(const ValueType& _value)
        {
            for (int i = 0; i < container_.size(); ++i)
                container_[i] = _value;
        }

        /**
         * Calls reserve on the container.
         */
        void reserve(const int _n)
        {
            container_.reserve(_n);
        }

        /**
         * Resizes the container.
         */
        void resize(const int _n)
        {
            ISM_ASSERT_GEQ(_n, 0);
            container_.resize(_n);
        }

        /**
         * Resizes the container.
         */
        void resize(const int _n, const ValueType& _value)
        {
            ISM_ASSERT_GEQ(_n, 0);
            container_.resize(_n, _value);
        }

        /**
         * Resizes the container to match the mesh.
         */
        template <typename Mesh>
        void resize(const Mesh& _mesh)
        {
            container_.resize(TW::nElements(_mesh));
        }

        /**
         * Resizes the container to match the mesh.
         */
        template <typename Mesh>
        void resize(const Mesh& _mesh, const ValueType& _value)
        {
            container_.resize(TW::nElements(_mesh), _value);
        }

        size_t size() const
        {
            return container_.size();
        }

        bool empty() const
        {
            return container_.empty();
        }

        template <typename Mesh>
        bool size_okay(const Mesh& _mesh) const
        {
            return TW::nElements(_mesh) == container_.size();
        }

        void push_back(const ValueType& _value)
        {
            container_.push_back(_value);
        }

        void clear()
        {
            container_.clear();
        }

        typename std::vector<ValueType>::const_reference operator[] (const HandleType _h) const
        {
            EXTERNAL_PROPERTY_RANGE_CHECK(_h);
            return container_[_h.idx()];
        }

        typename std::vector<ValueType>::reference operator[] (const HandleType _h)
        {
            EXTERNAL_PROPERTY_RANGE_CHECK(_h);
            return container_[_h.idx()];
        }

        template <typename T> using VecX = Eigen::Matrix<T, Eigen::Dynamic, 1>;

        const Eigen::Map<const VecX<ValueType>> as_eigen() const
        {
            return Eigen::Map<const VecX<ValueType>>(container_.data(), container_.size());
        }

        Eigen::Map<VecX<ValueType>> as_eigen()
        {
            return Eigen::Map<VecX<ValueType>>(container_.data(), container_.size());
        }

        template <typename Mesh>
        OpenMesh::BasePropHandleT<ValueType> to_open_mesh(Mesh& _mesh)
        {
            OpenMesh::BasePropHandleT<ValueType> ph;
            TW::toProperty(_mesh, ph, container_);

            return ph;
        }

        std::vector<ValueType> container_;
};

}
