/*
  Copyright (C) 2017-2019 Sven Willner <sven.willner@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as published
  by the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef NETCDF_FILE_H
#define NETCDF_FILE_H

#include <cassert>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include "netcdftools.h"
#include "nvector.h"
#include "version.h"

namespace netCDF {

template<typename T>
struct NetCDFType {};
template<>
struct NetCDFType<double> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_DOUBLE;
};
template<>
struct NetCDFType<float> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_FLOAT;
};
template<>
struct NetCDFType<std::int16_t> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_USHORT;
};
template<>
struct NetCDFType<int> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_INT;
};
template<>
struct NetCDFType<const char*> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_STRING;
};
template<>
struct NetCDFType<std::string> {
    static const netCDF::NcType::ncType type = netCDF::NcType::nc_STRING;
};

template<typename T>
class DimVar : public std::vector<T> {
  protected:
    const std::string name_m;
    std::vector<std::tuple<std::string, NcType, std::size_t, std::vector<char>>> attributes;

  public:
    using std::vector<T>::size;

    DimVar(const NcDim& dim, const NcVar& var) : name_m(dim.getName()), std::vector<T>(dim.getSize()) {
        var.getVar(&(*this)[0]);
        for (const auto& att : var.getAtts()) {
            std::size_t size = att.second.getAttLength();
            std::vector<char> value(size);
            att.second.getValues(&value[0]);
            attributes.emplace_back(std::make_tuple(att.first, att.second.getType(), size, value));
        }
    }
    DimVar(const DimVar&) = delete;
    DimVar(DimVar&&) = default;  // NOLINT(performance-noexcept-move-constructor,hicpp-noexcept-move) [cannot use noexpect here]
    NcDim write_to(NcFile& file) const {
        netCDF::NcDim result = file.addDim(name_m, size());
        netCDF::NcVar var = file.addVar(name_m, NetCDFType<T>::type, {result});
        for (const auto& att : attributes) {
            var.putAtt(std::get<0>(att), std::get<1>(att), std::get<2>(att), &std::get<3>(att)[0]);
        }
        var.putVar(&(*this)[0]);
        return result;
    }
    const std::string& name() const { return name_m; }
};

class File : public netCDF::NcFile {
  protected:
    template<typename T, std::size_t c, std::size_t dim, std::size_t... Ns>
    struct get_helper {
        static inline nvector::Vector<T, dim> reserve(const netCDF::NcVar& var_p) { return get_helper<T, c + 1, dim, Ns..., c>::reserve(var_p); }
    };

    template<typename T, std::size_t dim, std::size_t... Ns>
    struct get_helper<T, dim, dim, Ns...> {
        static inline nvector::Vector<T, dim> reserve(const netCDF::NcVar& var_p) { return nvector::Vector<T, dim>(0, var_p.getDim(Ns).getSize()...); }
    };

    template<typename T, std::size_t c, std::size_t dim, typename Storage, std::size_t... Ns>
    struct set_helper {
        static inline void set(netCDF::NcVar& var, const nvector::Vector<T, dim, Storage>& data) {
            set_helper<T, c + 1, dim, Storage, Ns..., c>::set(var, data);
        }
    };

    template<typename T, std::size_t dim, typename Storage, std::size_t... Ns>
    struct set_helper<T, dim, dim, Storage, Ns...> {
        static inline void set(netCDF::NcVar& var, const nvector::Vector<T, dim, Storage>& data) {
            var.putVar({(0 * Ns)...}, {data.size(Ns)...}, &data.data()[0]);
        }
    };

    template<typename T, std::size_t c, std::size_t dim, typename Storage, std::size_t... Ns>
    struct var_set_helper {
        template<typename... Indices>
        static inline void set(netCDF::NcVar& var, const nvector::Vector<T, dim, Storage>& data, Indices&&... indices) {
            var_set_helper<T, c + 1, dim, Storage, Ns..., c>::set(var, data, std::forward<Indices>(indices)...);
        }
    };

    template<typename T, std::size_t dim, typename Storage, std::size_t... Ns>
    struct var_set_helper<T, dim, dim, Storage, Ns...> {
        template<typename... Indices>
        static inline void set(netCDF::NcVar& var, const nvector::Vector<T, dim, Storage>& data, Indices&&... indices) {
            var.putVar({std::forward<Indices>(indices)..., (0 * Ns)...}, {(0 * sizeof(Indices) + 1)..., data.size(Ns)...}, &data.data()[0]);
        }
    };

  public:
    File() = default;

    File(const std::string& filename, char mode) { open(filename, mode); }

    using NcFile::open;

    void open(const std::string& filename, char mode) {
        switch (mode) {
            case 'r':
                open(filename, netCDF::NcFile::read);
                break;
            case 'w':
                open(filename, netCDF::NcFile::replace, netCDF::NcFile::nc4);
                {
                    auto now = std::time(nullptr);
                    std::ostringstream ss;
                    ss << std::put_time(std::localtime(&now), "%Y-%m-%d %H:%M:%S");
                    putAtt("created_at", ss.str());
                }
                putAtt("created_with", "https://github.com/swillner/flood-processing by S. Willner");
                putAtt("citation", "https://doi.org/10.5281/zenodo.1241051");
                putAtt("contact", "sven.willner@pik-potsdam.de");
                putAtt("flood_processing_version", FLOOD_PROCESSING_VERSION);
                break;
            default:
                throw std::runtime_error("unknown file mode");
        }
    }

    netCDF::NcVar var(const std::string& varname) { return getVar(varname); }

    template<typename T>
    netCDF::NcVar var(const std::string& varname, const std::vector<netCDF::NcDim>& dims) {
        netCDF::NcVar res = addVar(varname, NetCDFType<T>::type, dims);
        res.setCompression(false, true, 7);
        res.setFill<T>(true, std::numeric_limits<T>::quiet_NaN());
        return res;
    }

    template<typename T>
    std::vector<T> get(const std::string& varname) {
        auto variable = getVar(varname);
        std::vector<T> result(variable.getDim(0).getSize());
        variable.getVar(&result[0]);
        return result;
    }

    template<typename T, std::size_t dim>
    nvector::Vector<T, dim> get(const std::string& varname) {
        auto variable = getVar(varname);
        auto result = get_helper<T, 0, dim>::reserve(variable);
        variable.getVar(&result.data()[0]);
        return result;
    }

    template<typename T, std::size_t dim>
    nvector::Vector<T, dim> get(const netCDF::NcVar& var) {
        auto result = get_helper<T, 0, dim>::reserve(var);
        var.getVar(&result.data()[0]);
        return result;
    }

    template<typename T>
    DimVar<T> dimvar(const netCDF::NcDim& dim) {
        return DimVar<T>(dim, getVar(dim.getName()));
    }

    template<typename T>
    netCDF::NcDim dimvar(const DimVar<T>& dimvar) {
        return dimvar.write_to(*this);
    }

    template<typename T, std::size_t dim>
    nvector::Vector<T, dim> reserve(const std::string& varname) {
        return get_helper<T, 0, dim>::reserve(getVar(varname));
    }

    template<typename T, std::size_t dim>
    nvector::Vector<T, dim> reserve(const netCDF::NcVar& var) {
        return get_helper<T, 0, dim>::reserve(var);
    }

    template<unsigned char dim>
    inline std::size_t size(const netCDF::NcVar& var) const {
        return var.getDim(dim).getSize();
    }

    netCDF::NcDim lat(std::size_t lat_count, double from_lat = -90.0, double to_lat = 90.0) {
        netCDF::NcDim lat_dim = addDim("lat", lat_count);
        netCDF::NcVar lat_var = addVar("lat", NetCDFType<double>::type, {lat_dim});
        lat_var.putAtt("standard_name", "latitude");
        lat_var.putAtt("long_name", "latitude");
        lat_var.putAtt("units", "degrees_north");
        lat_var.putAtt("axis", "Y");
        for (std::size_t ilat = 0; ilat < lat_count; ++ilat) {
            lat_var.putVar({ilat}, to_lat - (to_lat - from_lat) * (ilat + 0.5) / lat_count);
        }
        return lat_dim;
    }

    netCDF::NcDim lon(std::size_t lon_count, double from_lon = -180.0, double to_lon = 180.0) {
        netCDF::NcDim lon_dim = addDim("lon", lon_count);
        netCDF::NcVar lon_var = addVar("lon", NetCDFType<double>::type, {lon_dim});
        lon_var.putAtt("standard_name", "longitude");
        lon_var.putAtt("long_name", "longitude");
        lon_var.putAtt("units", "degrees_east");
        lon_var.putAtt("axis", "X");
        for (std::size_t ilon = 0; ilon < lon_count; ++ilon) {
            lon_var.putVar({ilon}, from_lon + (to_lon - from_lon) * (ilon + 0.5) / lon_count);
        }
        return lon_dim;
    }

    template<typename T>
    void set(netCDF::NcVar var, const std::vector<T>& data) {
        var.putVar(&data[0]);
    }

    // template<typename T, std::size_t dim, typename Storage>
    // void set(netCDF::NcVar var, const nvector::Vector<T, dim, Storage>& data) { // TODO Really needed?
    //     set_helper<T, 0, dim, Storage>::set(var, data); // TODO Fix issue with dimension missmatch
    // }

    template<typename T, std::size_t dim, typename Storage, typename... Indices>
    void set(netCDF::NcVar var, const nvector::Vector<T, dim, Storage>& data, Indices&&... indices) {
        var_set_helper<T, 0, dim, Storage>::set(var, data, std::forward<Indices>(indices)...);  // TODO Fix issue with dimension missmatch
    }

    template<typename T>
    T get_attribute(const std::string& name) {
        T res;
        getAtt(name).getValues(res);
        return res;
    }

    template<typename T>
    T get_var_attribute(const netCDF::NcVar& var_p, const std::string& name) {
        T res;
        var_p.getAtt(name).getValues(res);
        return res;
    }

    netCDF::NcDim dim(const netCDF::NcDim& dim) {
        if (dim.isUnlimited()) {
            return addDim(dim.getName());
        }
        return addDim(dim.getName(), dim.getSize());
    }

    netCDF::NcVar var(const netCDF::NcVar& var_p) {
        std::vector<netCDF::NcDim> dims;
        std::vector<std::size_t> sizes;
        std::vector<std::size_t> indices;
        auto s = var_p.getType().getSize();
        for (const auto& dim : var_p.getDims()) {
            dims.emplace_back(getDim(dim.getName()));
            sizes.emplace_back(dim.getSize());
            indices.emplace_back(0);
            s *= dim.getSize();
        }
        netCDF::NcVar res = addVar(var_p.getName(), var_p.getType().getTypeClass(), dims);
        res.setCompression(false, true, 7);
        for (const auto& att : var_p.getAtts()) {
            const auto len = att.second.getAttLength();
            std::vector<char> value(len);
            att.second.getValues(&value[0]);
            res.putAtt(att.first, att.second.getType(), len, &value[0]);
        }
        std::vector<char> fill_value(var_p.getType().getSize());
        std::vector<char> value(s);
        bool fill_mode;
        var_p.getFillModeParameters(fill_mode, &fill_value[0]);
        if (fill_mode) {
            res.setFill(fill_mode, &fill_value[0]);
        }
        var_p.getVar(&value[0]);
        res.putVar(indices, sizes, &value[0]);
        return res;
    }
};

template<>
inline netCDF::NcVar File::var<const char*>(const std::string& varname, const std::vector<netCDF::NcDim>& dims) {
    netCDF::NcVar res = addVar(varname, NetCDFType<const char*>::type, dims);
    res.setCompression(false, true, 7);
    return res;
}

}  // namespace netCDF

#endif
