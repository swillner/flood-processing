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

#ifndef FLOOD_PROCESSING_FORTRANGRID_H
#define FLOOD_PROCESSING_FORTRANGRID_H

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <iomanip>
#include <string>
#include <vector>

template<typename T>
class FortranGrid {
  protected:
    int fd;
    std::unique_ptr<T[]> data;
    std::size_t lon_count_;
    std::size_t lat_count_;

  public:
    FortranGrid(const std::string& filename, std::size_t lon_count_p, std::size_t lat_count_p, char mode) : lat_count_(lat_count_p), lon_count_(lon_count_p) {
        const std::size_t size = lon_count_ * lat_count_ * sizeof(T);

        switch (mode) {
            case 'r': {
                fd = open(filename.c_str(), O_RDONLY | O_CLOEXEC);  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg,hicpp-signed-bitwise)
                if (fd < 0) {
                    throw std::runtime_error("could not open file " + filename);
                }
                data.reset(static_cast<T*>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0)));  // NOLINT(hicpp-signed-bitwise)
            } break;

            case 'w': {
                fd = open(filename.c_str(), O_RDWR | O_CREAT | O_CLOEXEC,  // NOLINT(cppcoreguidelines-pro-type-vararg,hicpp-vararg,hicpp-signed-bitwise)
                          static_cast<mode_t>(0600));
                if (fd < 0) {
                    throw std::runtime_error("could not create file " + filename);
                }
                int rc = lseek(fd, size - 1, SEEK_SET);
                if (rc < 0) {
                    throw std::runtime_error("lseek failed");
                }
                rc = write(fd, "", 1);
                if (rc < 0) {
                    throw std::runtime_error("write failed");
                }
                data.reset(static_cast<T*>(mmap(nullptr, size, PROT_WRITE, MAP_SHARED, fd, 0)));
            } break;

            default:
                throw std::runtime_error("unknown file mode");
        }

        if (data.get() == MAP_FAILED) {  // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
            throw std::runtime_error("mmap failed");
        }
        int rc = madvise(data.get(), size, MADV_WILLNEED | MADV_SEQUENTIAL);  // NOLINT(hicpp-signed-bitwise)
        if (rc < 0) {
            throw std::runtime_error("madvice failed");
        }
    }

    ~FortranGrid() {
        munmap(data.get(), lon_count_ * lat_count_ * sizeof(T));
        close(fd);
        data.release();
    }

    operator T*() { return data.get(); }        // NOLINT(hicpp-explicit-conversions,google-explicit-constructor)
    operator T*() const { return data.get(); }  // NOLINT(hicpp-explicit-conversions,google-explicit-constructor)

    T& operator()(const std::size_t& lon, const std::size_t& lat) noexcept { return data[lat * lon_count_ + lon]; }

    const T& operator()(const std::size_t& lon, const std::size_t& lat) const noexcept { return data[lat * lon_count_ + lon]; }

    constexpr std::size_t lat_count() const { return lat_count_; }
    constexpr std::size_t lon_count() const { return lon_count_; }
};

#endif
