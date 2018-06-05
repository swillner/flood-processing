/*
  Copyright (C) 2017 Sven Willner <sven.willner@gmail.com>

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

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#ifdef FLOOD_PROCESSING_WITH_TQDM
#include <string>
#include "tqdm/tqdm.h"
#endif

namespace flood_processing {

class ProgressBar {
  protected:
#ifdef FLOOD_PROCESSING_WITH_TQDM
    std::unique_ptr<tqdm::RangeTqdm<std::size_t>> it;
#endif

  public:
    ProgressBar(std::string description, std::size_t length) {
#ifdef FLOOD_PROCESSING_WITH_TQDM
        tqdm::Params p;
        p.desc = description;
        p.ascii = "";
        p.f = stdout;
        it.reset(new tqdm::RangeTqdm<std::size_t>{tqdm::RangeIterator<std::size_t>(length), tqdm::RangeIterator<std::size_t>(length, length), p});
#endif
    }
    inline void tick() {
#ifdef FLOOD_PROCESSING_WITH_TQDM
#pragma omp critical(output)
        { ++(*it); }
#endif
    }
    ~ProgressBar() {
#ifdef FLOOD_PROCESSING_WITH_TQDM
        it->close();
#endif
    }
};

}  // namespace flood_processing

#endif
