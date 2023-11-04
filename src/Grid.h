#pragma once

#include<iostream>
#include <stdint.h>
#include<memory>
#include "Vector.h"
#include <thread>
#include <atomic>
#include "Task.h"
#include <atomic>
#include <vector>
#include"MedianCounter.h"
#include<functional>

using FieldCell = bool;

namespace fieldCell {
    static constexpr FieldCell cellAlive = true;
    static constexpr FieldCell cellDead = false;

    inline FieldCell nextGeneration(FieldCell const cell, uint8_t const aliveNeighboursCount) {
        if(!cell && aliveNeighboursCount == 3) return true;
        else if(cell && aliveNeighboursCount != 2 && aliveNeighboursCount != 3) return false;
        else return cell;
    }
    inline constexpr char const *asString(const FieldCell cell) {
        return cell ? "alive" : "dead";
    }
}

struct FieldModification {
    uint32_t startIndex_int, size_int;
    const uint32_t *data;
};

struct FieldOutput { //TODO: remove
    virtual void write(FieldModification fm) {
        this->batched()->write(fm);
    }

    virtual std::unique_ptr<FieldOutput> batched() const = 0;
    virtual ~FieldOutput() = default;
}; //must be used as output only in one thread


struct Cell {
    FieldCell cell;
    int32_t index;
};

class Field final {
public:
    struct FieldPimpl;
    struct GridData;
private:
    std::unique_ptr<FieldPimpl> gridPimpl;
    bool isStopped;
    std::unique_ptr<FieldOutput> const current_output;
    std::unique_ptr<FieldOutput> const buffer_output;

    const uint32_t numberOfTasks;
    std::unique_ptr<std::unique_ptr<Task<GridData>>[/*numberOfTasks*/]> gridTasks;
    std::atomic_bool interrupt_flag;
    std::vector<uint32_t> indecesToBrokenCells;
public:
    Field(
        const uint32_t gridWidth, const uint32_t gridHeight, const size_t numberOfTasks_, 
        std::function<std::unique_ptr<FieldOutput>()> current_outputs, 
        std::function<std::unique_ptr<FieldOutput>()> buffer_outputs
    );
    ~Field();

    Field(Field const&) = delete;
    Field& operator=(Field const&) = delete;
public:
    bool tryFinishGeneration();
    void startCurGeneration();
    void startNewGeneration();

    void fill(const FieldCell cell);

    FieldCell cellAtIndex(const uint32_t index) const;

    void setCellAtIndex(const uint32_t index, FieldCell cell);
    void setCellAtCoord(const vec2i& coord, FieldCell cell);
    void setCells(Cell const *const cells, size_t const count);

    FieldCell cellAtCoord(const vec2i& coord) const;
    FieldCell cellAtCoord(const int32_t column, const int32_t row) const;

    vec2i indexAsCoord(const int32_t index) const;

    int32_t coordAsIndex(const vec2i& coord) const;
    int32_t coordAsIndex(const int32_t column, const int32_t row) const;
    
    int32_t normalizeIndex(const int32_t index) const;
    vec2i normalizeCoord(const vec2i& coord) const;

    uint32_t width() const;
    uint32_t height() const;
    uint32_t size() const;

    uint32_t size_bytes() const;
    //uint32_t size_actual() const;
    uint32_t width_actual() const;
    uint32_t *rawData() const;
private:
    void waitForGridTasks();
    void deployGridTasks();
};

inline void Field::setCellAtCoord(const vec2i& coord, FieldCell cell) {
    setCellAtIndex(coordAsIndex(coord), cell);
}

inline FieldCell Field::cellAtCoord(const vec2i& coord) const {
    return cellAtIndex(static_cast<int32_t>(coordAsIndex(coord)));
}

inline FieldCell Field::cellAtCoord(const int32_t column, const int32_t row) const {
    return cellAtCoord(vec2i(column, row));
}

inline vec2i Field::indexAsCoord(const int32_t index) const {
    const auto index_n = normalizeIndex(index);
    const auto x = index_n % width();
    const auto y = index_n / width();
    return vec2i(x, y);
}

inline int32_t Field::coordAsIndex(const vec2i& coord) const {
    const auto coord_n = normalizeCoord(coord);
    return coord_n.x + coord_n.y * width();
}

inline int32_t Field::coordAsIndex(const int32_t column, const int32_t row) const {
    return coordAsIndex(vec2i(column, row));
}

inline int32_t Field::normalizeIndex(const int32_t index) const {
    return misc::mod(index, static_cast<int32_t>(size()));
}
inline vec2i Field::normalizeCoord(const vec2i& coord) const {
    return vec2i(misc::mod(coord.x, width()), misc::mod(coord.y, height()));
}

