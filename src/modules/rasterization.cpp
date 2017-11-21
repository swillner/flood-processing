/*
  Copyright (C) 2017 Sven Willner <sven.willner@pik-potsdam.de>

  This file is part of flood-processing.

  flood-processing is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of
  the License, or (at your option) any later version.

  flood-processing is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with flood-processing.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "modules/rasterization.h"
#ifdef FLOOD_PROCESSING_WITH_TQDM
#include "tqdm/tqdm.h"
#endif
#include <gdal_alg.h>
#include <gdal_priv.h>
#include <ogrsf_frmts.h>
#include <algorithm>
#include <vector>
#undef CPL_CVSID
#define CPL_CVSID(string)
#include "../src/gdal/gdalrasterize.cpp"
#include "../src/gdal/llrasterize.cpp"

namespace flood_processing {
namespace modules {

template<typename T>
Rasterization<T>::Rasterization(const settings::SettingsNode& settings) {
    resolution_mask_name = settings["resolution_from"].as<std::string>();
    shapefilename = settings["shapefile"].as<std::string>();
    layername = settings["layer"].as<std::string>();
    fieldname = settings["field"].as<std::string>();
    max_advance = settings["advance"].as<std::size_t>(0);
    adjust_scale = settings["adjust_scale"].as<std::size_t>(1);
    adjust_max = settings["adjust_max"].as<bool>(false);
    invalid_value = settings["invalid_value"].as<T>(std::numeric_limits<T>::quiet_NaN());
    GDALAllRegister();
}

template<typename T>
static T most_common(std::vector<T>& vec) {
    std::sort(std::begin(vec), std::end(vec), [](T a, T b) { return (a < b || std::isnan(b)) && !std::isnan(a); });
    T result = std::numeric_limits<T>::quiet_NaN();
    std::size_t result_count = 0;
    T current = std::numeric_limits<T>::quiet_NaN();
    std::size_t current_count = 0;
    for (const auto& v : vec) {
        if (std::isnan(v)) {
            break;
        }
        if (v == current) {
            ++current_count;
        } else {
            if (current_count > result_count) {
                result = current;
                result_count = current_count;
            }
            current = v;
            current_count = 1;
        }
    }
    if (current_count > result_count) {
        result = current;
        result_count = current_count;
    }
    return result;
}

template<typename T>
inline void Rasterization<T>::combine_fine(const nvector::View<T, 2>& fine_raster, nvector::View<T, 2>& raster) {
#ifdef FLOOD_PROCESSING_WITH_TQDM
    const auto size = raster.template size<0>() * raster.template size<1>();
    tqdm::Params p;
    p.desc = "Combining";
    p.ascii = "";
    p.f = stdout;
    tqdm::RangeTqdm<std::size_t> it{tqdm::RangeIterator<std::size_t>(size), tqdm::RangeIterator<std::size_t>(size, size), p};
#endif
    raster.foreach_parallel([&](std::size_t lat, std::size_t lon, T& n) {
        std::vector<T> inner;
        inner.reserve(adjust_scale * adjust_scale);
        for (std::size_t inner_lat = adjust_scale * lat; inner_lat < adjust_scale * (lat + 1); ++inner_lat) {
            for (std::size_t inner_lon = adjust_scale * lon; inner_lon < adjust_scale * (lon + 1); ++inner_lon) {
                inner.push_back(fine_raster(inner_lat, inner_lon));
            }
        }
        if (adjust_max) {
            const auto& m = std::max_element(std::begin(inner), std::end(inner), [](T a, T b) { return (a < b || std::isnan(b)) && !std::isnan(a); });
            if (m != std::end(inner)) {
                n = *m;
            }
        } else {
            n = most_common(inner);
        }
#ifdef FLOOD_PROCESSING_WITH_TQDM
#pragma omp critical(output)
        { ++it; }
#endif
    });
#ifdef FLOOD_PROCESSING_WITH_TQDM
    it.close();
#endif
}

template<typename T>
inline void Rasterization<T>::advance(nvector::View<T, 2>& result, std::size_t max_advance) {
#ifdef FLOOD_PROCESSING_WITH_TQDM
    tqdm::Params p;
    p.desc = "Advancing";
    p.ascii = "";
    p.f = stdout;
    tqdm::RangeTqdm<std::size_t> it{tqdm::RangeIterator<std::size_t>(max_advance), tqdm::RangeIterator<std::size_t>(max_advance, max_advance), p};
#endif
    nvector::Vector<T, 2> last(0, result.template size<0>(), result.template size<1>());
    for (std::size_t i = 0; i < max_advance; ++i) {
        std::copy(std::begin(result), std::end(result), std::begin(last));
        bool found = false;
        foreach_view_parallel(std::make_tuple(result, last), [&](std::size_t lat, std::size_t lon, T& n, const T& l) {
            if (std::isnan(l)) {
                std::vector<T> neighbours;
                if (lon > 0) {
                    if (lat > 0) {
                        neighbours.push_back(last(lat - 1, lon - 1));
                    }
                    neighbours.push_back(last(lat, lon - 1));
                    if (lat < result.template size<0>() - 1) {
                        neighbours.push_back(last(lat + 1, lon - 1));
                    }
                }
                if (lon < result.template size<1>() - 1) {
                    if (lat > 0) {
                        neighbours.push_back(last(lat - 1, lon + 1));
                    }
                    neighbours.push_back(last(lat, lon + 1));
                    if (lat < result.template size<0>() - 1) {
                        neighbours.push_back(last(lat + 1, lon + 1));
                    }
                }
                if (lat > 0) {
                    neighbours.push_back(last(lat - 1, lon));
                }
                neighbours.push_back(last(lat, lon));
                if (lat < result.template size<0>() - 1) {
                    neighbours.push_back(last(lat + 1, lon));
                }
                n = most_common(neighbours);
                if (!std::isnan(n)) {
                    found = true;
                }
            }
        });
        if (!found) {
            break;
        }
#ifdef FLOOD_PROCESSING_WITH_TQDM
        ++it;
#endif
    }
#ifdef FLOOD_PROCESSING_WITH_TQDM
    it.close();
#endif
}

template<typename T>
template<typename Function>
inline void Rasterization<T>::rasterize(nvector::View<T, 2>& result, Function&& func) {
    auto infile = static_cast<GDALDataset*>(GDALOpenEx(shapefilename.c_str(), GDAL_OF_VECTOR, nullptr, nullptr, nullptr));
    if (!infile) {
        throw std::runtime_error("could not open shape file");
    }
    OGRLayer* inlayer = infile->GetLayerByName(layername.c_str());
    if (!inlayer) {
        throw std::runtime_error("could not read layer from shape file");
    }

    const std::size_t lat_count = result.template size<0>();
    const std::size_t lon_count = result.template size<1>();
    double geotransform[] = {-180., 360. / lon_count, 0, 90., 0, -180. / lat_count};
    // this is how geotransform is used by gdal:
    //   *pdfGeoX = padfGeoTransform[0] + dfPixel * padfGeoTransform[1] - dfLine * padfGeoTransform[2];
    //   *pdfGeoY = padfGeoTransform[3] + dfPixel * padfGeoTransform[4] - dfLine * padfGeoTransform[5];

    void* transform;
    {
        char* projection = nullptr;
        OGRSpatialReference* spatial_reference = inlayer->GetSpatialRef();
        if (!spatial_reference) {
            throw std::runtime_error("could not get spatial reference");
        } else {
            spatial_reference->exportToWkt(&projection);
        }
        transform = GDALCreateGenImgProjTransformer3(projection, nullptr, nullptr, geotransform);
        CPLFree(projection);
    }

    std::vector<double> data(lat_count * lon_count, std::numeric_limits<double>::quiet_NaN());
    const std::size_t size = inlayer->GetFeatureCount();
#ifdef FLOOD_PROCESSING_WITH_TQDM
    tqdm::Params p;
    p.desc = "Rasterizing";
    p.ascii = "";
    p.f = stdout;
    tqdm::RangeTqdm<std::size_t> it{tqdm::RangeIterator<std::size_t>(size), tqdm::RangeIterator<std::size_t>(size, size), p};
#endif
    inlayer->ResetReading();
#pragma omp parallel for default(shared) schedule(dynamic)
    for (std::size_t i = 0; i < size; ++i) {
        OGRFeature* infeature;
#pragma omp critical(infeature)
        { infeature = inlayer->GetNextFeature(); }
        OGRGeometry* geometry = infeature->GetGeometryRef();  //->SimplifyPreserveTopology(pixel_size / 16);
        double value = func(infeature, i);
        gv_rasterize_one_shape(static_cast<unsigned char*>(static_cast<void*>(&data.data()[0])),  // unsigned char *pabyChunkBuf,
                               0,                                                                 // int nYOff,
                               lon_count,                                                         // int nXSize,
                               lat_count,                                                         // int nYSize,
                               1,                                                                 // int nBands,
                               GDT_Float64,                                                       // GDALDataType eType,
                               0,                                                                 // int bAllTouched,
                               geometry,                                                          // OGRGeometry *poShape,
                               &value,                                                            // double *padfBurnValue,
                               GBV_UserBurnValue,                                                 // GDALBurnValueSrc eBurnValueSrc,
                               GRMA_Replace,                                                      // GDALRasterMergeAlg eMergeAlg,
                               GDALGenImgProjTransform,                                           // GDALTransformerFunc pfnTransformer,
                               transform                                                          // void *pTransformArg
        );
        OGRFeature::DestroyFeature(infeature);
        // OGRGeometryFactory::destroyGeometry(geometry);
#ifdef FLOOD_PROCESSING_WITH_TQDM
#pragma omp critical(output)
        { ++it; }
#endif
    }
    GDALClose(infile);
#ifdef FLOOD_PROCESSING_WITH_TQDM
    it.close();
#endif

    GDALDestroyTransformer(transform);

    std::move(std::begin(data), std::end(data), std::begin(result));
}

template<typename T>
void Rasterization<T>::run(pipeline::Pipeline* p) {
    const auto resolution_mask = p->consume<nvector::View<T, 2>>(resolution_mask_name);
    auto raster =
        std::make_shared<nvector::Vector<T, 2>>(std::numeric_limits<T>::quiet_NaN(), resolution_mask->template size<0>(), resolution_mask->template size<1>());
    if (adjust_scale > 1) {
        nvector::Vector<T, 2> fine_raster(std::numeric_limits<T>::quiet_NaN(), adjust_scale * resolution_mask->template size<0>(),
                                          adjust_scale * resolution_mask->template size<1>());
        rasterize(fine_raster, [&](OGRFeature* feature, std::size_t index) {
            (void)index;
            return feature->GetFieldAsDouble(fieldname.c_str());
        });
        combine_fine(fine_raster, *raster);
    } else {
        rasterize(*raster, [&](OGRFeature* feature, std::size_t index) {
            (void)index;
            return feature->GetFieldAsDouble(fieldname.c_str());
        });
    }
    if (max_advance > 0) {
        advance(*raster, max_advance);
    }
    if (!std::isnan(invalid_value)) {
        raster->foreach_parallel([&](std::size_t lat, std::size_t lon, T& n) {
            if (n == invalid_value) {
                n = std::numeric_limits<T>::quiet_NaN();
            }
        });
    }
    p->provide<nvector::Vector<T, 2>>("raster", raster);
}

template<typename T>
void RegionIndexRasterization<T>::run(pipeline::Pipeline* p) {
    const auto regions = p->consume<std::vector<std::string>>("regions");
    const auto resolution_mask = p->consume<nvector::View<T, 2>>(resolution_mask_name);
    auto raster =
        std::make_shared<nvector::Vector<T, 2>>(std::numeric_limits<T>::quiet_NaN(), resolution_mask->template size<0>(), resolution_mask->template size<1>());
    rasterize(*raster, [&](OGRFeature* feature, std::size_t index) {
        (void)index;
        const auto& id = feature->GetFieldAsString(this->fieldname.c_str());
        const auto& r = std::find(regions->begin(), regions->end(), id);
        if (r == regions->end()) {
            if (!altfieldname.empty()) {
                const auto& id2 = feature->GetFieldAsString(altfieldname.c_str());
                const auto& r2 = std::find(regions->begin(), regions->end(), id2);
                if (r2 == regions->end()) {
                    return std::numeric_limits<double>::quiet_NaN();
                } else {
                    return static_cast<double>(std::distance(regions->begin(), r2));
                }
            } else {
                return std::numeric_limits<double>::quiet_NaN();
            }
        } else {
            return static_cast<double>(std::distance(regions->begin(), r));
        }
    });
    if (max_advance > 0) {
        advance(*raster, max_advance);
    }
    p->provide<nvector::Vector<T, 2>>("region_index_raster", raster);
}

template class Rasterization<float>;
template class Rasterization<double>;
template class RegionIndexRasterization<float>;
template class RegionIndexRasterization<double>;

}  // namespace modules
}  // namespace flood_processing
